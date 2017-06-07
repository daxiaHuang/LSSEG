//===========================================================================
// The Level-Set Segmentation Library (LSSEG)
//
//
// Copyright (C) 2000-2005 SINTEF ICT, Applied Mathematics, Norway.
//
// This program is free software; you can redistribute it and/or          
// modify it under the terms of the GNU General Public License            
// as published by the Free Software Foundation version 2 of the License. 
//
// This program is distributed in the hope that it will be useful,        
// but WITHOUT ANY WARRANTY; without even the implied warranty of         
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
// GNU General Public License for more details.                           
//
// You should have received a copy of the GNU General Public License      
// along with this program; if not, write to the Free Software            
// Foundation, Inc.,                                                      
// 59 Temple Place - Suite 330,                                           
// Boston, MA  02111-1307, USA.                                           
//
// Contact information: e-mail: tor.dokken@sintef.no                      
// SINTEF ICT, Department of Applied Mathematics,                         
// P.O. Box 124 Blindern,                                                 
// 0314 Oslo, Norway.                                                     
// 
//
// Other licenses are also available for this software, notably licenses
// for:
// - Building commercial software.                                        
// - Building software whose source code you wish to keep private.        
//
//===========================================================================
#include "LevelSetFunction.h"
#include "NormalDistributionForce.h"
#include "ParzenDistributionForce.h"
#include "SingleRegionAlgorithm.h"
#include "MultiRegionAlgorithm.h"
#include "Region.h"
#include "simple_tools.h"
#include "cimg_dependent.h"
#include "level_set.h"
#include "Filters.h"
#include "image_utils.h"
#include "colordefs.h"
#include <time.h>
#include <string>
#include <assert.h>
#include <stdexcept>

using namespace std;
using namespace lsseg;

/// lowest resolution of downsampled image
int LOWEST_RES            = 70; 
/// boundary smoothness factor
double NU_FACTOR          = 1.5; 
/// nonlinear smoothing factor
double P                  = 1.5; 
/// threshold for when to stop when carrying out anisotropic presmoothing (lower value means more smoothing.
double SMOOTHING_FACTOR   = 0.3; 
/// by which factor should we downsample the image (in x-resolution and y-resolution)
double DOWNSCALING_FACTOR = sqrt(double(2)); 
/// number of iterations for segmentation on lowest level
int FIRST_ITERATIONS      = 1800;
/// number of iterations on higher levels 
int LATER_ITERATIONS      = 500; 
/// 2->use color info, 1->use grayscale, 0->do not use original image
int USE_ORIG_IMAGE        = 2; 
/// use the \ref anchor_StructureTensor "structure tensor" in the segmentation process, 1->yes, 2->no
int USE_STRUCTURE_TENSOR  = 1; 
/// use the scale measure in the segmentation process, 1->yes, 2->no
int USE_SCALE_INFO        = 1; 
/// number of segmented regions.  Default is normal, two-region segmentation
int NUM_REG               = 2; 

///  centered circle for two regions; 0 -> centered circle for two-region segmentation, or random voronoi for multiregions.
///  n > 0 -> a given number of horizontal stripes, n < 0 -> (for multiregions) number of fragmented voronoi regions.
int INIT_SEG_TYPE         = 0;  
/// show intermediary results.  Turn off if you do not want to push a key to continue from time to time
bool SHOW_INTERMEDIARY    = true;  
/// save the resulting segmentation (currently works for two-region segmentation only)
bool SAVE_SEG             = false; 
/// filename to use when saving results
string SAVE_FILENAME; 
/// whether or not to use a \ref Mask
bool USE_MASK = false;
/// \ref Mask to use
Mask mask_object;   
/// define if the \ref Mask or its inverse should be used for masking.
int MFLIP                 = 0; 
/// timestep size for computing \ref anchor_ScaleMeasure "scale measure" 
const double SCALE_FAC_DT = 1; 
/// number of timesteps to use when computing the \ref anchor_ScaleMeasure "scale measure"
const double SCALE_FAC_T = 15; 
/// how often \ref anchor_Reinitialization "reinitialization" should be applied on the coarsest level.
const int REINIT_MODULO_LOWEST = 100;
/// how often \ref anchor_Reinitialization "reinitialization" should be applied on higher levels 
const int REINIT_MODULO_OTHER = 20;  
/// timestep size when carrying out anisotropic smoothing
const int ANISOTROPIC_DT = 3;
/// how high the image resolution should get before applying the Parzen probability estimate rather than using a normal distribution
const int PARZEN_LIMIT = 150 * 100; 
/// number of available colors (excluding the background color, black)
const int TOTAL_NUM_COLORS = 9; 
const int* COLORS[] = {YELLOW, CYAN, MAGENTA, WHITE, GREY, BROWN, GREEN, BLUE, RED, BLACK}; /// list of colors

void show_usage();
void parse_command_line(int varnum, char** vararg);
void initialize_levelsets(vector<LevelSetFunction>& phi);


//===========================================================================
int main(int varnum, char** vararg)
//===========================================================================
{
    if (varnum == 1) {
	show_usage();
	return 0;
    }
    parse_command_line(varnum, vararg);
    
    // loading image
    Image<double> img;
    const bool convert_to_greyscale = USE_ORIG_IMAGE != 2;
    load_image(vararg[1], img, convert_to_greyscale);

    // generating and initializing level-set-functions
    const bool multireg = NUM_REG > 2;
    vector<LevelSetFunction> phi(multireg ? NUM_REG : 1, LevelSetFunction(img.dimx(), img.dimy()));
    initialize_levelsets(phi);

    // making multiresolution version of image and level-set functions
    vector<Image<double> > img_stack;
    vector<Mask> mask_stack;
    vector<const Mask*> mask_stack_ptr;
    vector<vector<LevelSetFunction> > phi_stack(phi.size());
    const bool downsample_z = false;
    downsample_series(img, img_stack, LOWEST_RES, downsample_z, DOWNSCALING_FACTOR, false);

    // setting up mask
    mask_stack_ptr.resize(img_stack.size(), 0);
    if (USE_MASK) {
	if (MFLIP) {
	    // inverting mask
	    for (char* p = mask_object.begin(); p != mask_object.end(); ++p) {
		*p = (*p) ? 0 : 1;
	    }
	}
	downsample_series(mask_object, mask_stack, LOWEST_RES, downsample_z, DOWNSCALING_FACTOR, false);
	for (size_t i = 0; i != mask_stack.size(); ++i) {
	    mask_stack_ptr[i] = &(mask_stack[i]);
	}
    }  

    for (size_t i = 0; i != phi.size(); ++i) {
	downsample_series(phi[i], phi_stack[i], LOWEST_RES, downsample_z, DOWNSCALING_FACTOR, false);
    }
    cout << "Downsampling finished" << endl;

    const int num_levels = img_stack.size();
    cout << "Number of levels: " << num_levels << endl;
    

    vector<Region> regs(phi.size());
    Image<double> multichan;
    Image<double> G, S, T_acc;
    UpdatableImage status_visualisator(img, "viewer");

    vector<ParzenDistributionForce> parzen_force(phi.size(), ParzenDistributionForce(multireg)); 
    vector<NormalDistributionForce> ndist_force(phi.size(), NormalDistributionForce(multireg)); 

    // running segmentation on each level
    for (int level = num_levels - 1; level >= 0; --level) {
	
	// initialize regions
	for (int r = 0; r != regs.size(); ++r) {
	    regs[r].mu = compute_nu_Brox(phi_stack[r][level].size(), NU_FACTOR);
	    regs[r].phi.swap(phi_stack[r][level]);
	}
	
	// making composite channel image to use for generating driving force
	vector<const Image<double>*> ch_vec;
	ch_vec.push_back(&img_stack[level]);
	Image<double> grayimage(img_stack[level]);
	if (grayimage.numChannels() > 1) {
	    to_grayscale(grayimage); 
	}

	if (USE_STRUCTURE_TENSOR) {
	    compute_structure_tensor_2D(grayimage, G, true);
	    if (SHOW_INTERMEDIARY) {
		display_image(G);
	    }
	    ch_vec.push_back(&G);
	} 
	if (USE_SCALE_INFO) {
	    compute_scale_factor_2D(grayimage, S, T_acc, SCALE_FAC_DT, SCALE_FAC_T);
	    if (SHOW_INTERMEDIARY) {
		display_image(S);
	    }
	    ch_vec.push_back(&S);
	}

	combine_channel_images(&ch_vec[0], ch_vec.size(), multichan);

	// anisotropic blurring
	cout << "doing anisotropic blurring" << endl;
	Image<double> tmp, diff;
 	while (level != 0) { // no smoothing at the highest level
	    anisotropic_smoothing(multichan, tmp, ANISOTROPIC_DT, 0, P);
	    diff = tmp;
	    diff -= multichan;
	    diff.setAbsolute();
	    const double average = diff.getAverage();
	    cout << "average pixel variation: " << average << endl;
	    multichan.swap(tmp);
	    if (average < SMOOTHING_FACTOR) {
		break; // this is the only way to get out of here
	    }
	}
	rescale_channels(multichan, 0, 255);

	const bool use_parzen = ((multichan.dimx() * multichan.dimy()) >= PARZEN_LIMIT) || (level == 0);
	for (int r = 0; r < regs.size(); ++r) {
	    if (use_parzen) {
		cout << "Using PARZEN distribution on level: " << level << endl;
		regs[r].fgen = &parzen_force[r];
	    } else {
		cout << "Using NORMAL distribution on level: " << level << endl;
		regs[r].fgen = &ndist_force[r];
	    }
	    regs[r].fgen->init(&multichan, mask_stack_ptr[level]);
	}
	
	const int num_iter =      (level == num_levels - 1) ? FIRST_ITERATIONS : LATER_ITERATIONS;
	const int reinit_modulo = (level == num_levels - 1) ? REINIT_MODULO_LOWEST : REINIT_MODULO_OTHER;
	
	Image<double> visualization;
	visualization.resize(regs[0].phi);

	for (int i = 0; i < num_iter; i+= reinit_modulo) {
	    
	    // develop region(s)
	    if (NUM_REG > 2) { // multiregion
		develop_multiregion_2D(&regs[0], regs.size(), reinit_modulo, 0, mask_stack_ptr[level]);
	    } else {           // single region
		develop_single_region_2D(regs[0], reinit_modulo, 0, mask_stack_ptr[level]);
	    } 
	    cout << "We are currently at iteration: " << i << ", level " << level << endl;

	    if (SHOW_INTERMEDIARY) {
		// display temporary result
		Image<double> tmp(img_stack[level]);
		Image<double> visu(regs[0].phi, false);
		if (NUM_REG > 2) {
		    vector<const LevelSetFunction*> phipointers(regs.size());
		    for (size_t r = 0; r != regs.size(); ++r) {
			phipointers[r] = &(regs[r].phi);
		    }
		    const int offset = TOTAL_NUM_COLORS - regs.size();
		    visualize_multisets(&phipointers[0], phipointers.size(), visu, COLORS + offset);
		} else {
		    visualize_level_set(regs[0].phi, visu, 0);
		}
		if (visu.numChannels() > tmp.numChannels()) {
		    visu.swap(tmp);
		}
		tmp.superpose(visu);
		if (mask_stack_ptr[level]) {
		    tmp.applyMask(*mask_stack_ptr[level]);
		}
		status_visualisator.update(tmp);
	    }
	    
	    for (size_t r = 0; r != regs.size(); ++r) {
		regs[r].phi.reinitialize2D();
	    }
	}
	if (level > 0) {
	    for (size_t r = 0; r != regs.size(); ++r) {
		resample_into(regs[r].phi, phi_stack[r][level-1]);
	    }
	} else {
	    for (size_t r = 0; r != regs.size(); ++r) {
		regs[r].phi.swap(phi_stack[r][0]);
	    }
	}
    }

    Image<double> visu(phi_stack[0][0], false);
    Image<double> tmp(img_stack[0]);
    if (NUM_REG > 2) {

	vector<const LevelSetFunction*> phipointers(regs.size());
	for (size_t r = 0; r != regs.size(); ++r) {
	    phipointers[r] = &(phi_stack[r][0]);
	}
	const int offset = TOTAL_NUM_COLORS - regs.size();
	visualize_multisets(&phipointers[0], phipointers.size(), visu, COLORS + offset);

    } else {
	visualize_level_set(phi_stack[0][0], visu, 0);
    }
    tmp.superpose(visu);
    display_image(tmp);
    if (SAVE_SEG) {
	Mask m;
	mask_from_segmentation(phi_stack[0][0], m, SEG_NEGATIVE);
	for (char* p = m.begin(); p != m.end(); ++p) {
	    if (*p) {
		*p = 255; // makes it easier to see the mask with standard imaging software
	    }
	}
	save_image(SAVE_FILENAME.c_str(), m);
    }
    return 0;
};

//===========================================================================
void parse_command_line(int varnum, char** vararg)
//===========================================================================
{
    for (int i = 2; i < varnum-1; ++i) {
	string opt(vararg[i]);
	if (opt == string("-m")) { 
	    Image<int> tmp;
	    load_image(vararg[++i], tmp, true);
	    mask_object = Mask(tmp, true);
	    USE_MASK = true;
	} else if (opt == string("-mflip")) {
	    const int tmp = atoi(vararg[++i]);
	    MFLIP = (tmp ? 1 : 0);
	} else if (opt == string("-save")) {
	    SAVE_SEG = true;
	    SAVE_FILENAME = string(vararg[++i]);
	    SAVE_FILENAME += string(".png");
	} else if (opt == string("-lr")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp > 0) {
		LOWEST_RES = tmp;
	    } else {
		cerr << "Ignored invalid value for -lr option (must be positive)" << endl;
	    }
	} else if (opt == string("-n")) {
	    const float tmp  = atof(vararg[++i]);
	    if (tmp > 0) {
		NU_FACTOR = tmp; 
	    } else {
		cerr << "Ignored invalid value for -n option (must be positive)" << endl;
	    }
	} else if (opt == string("-p")) {
	    const float tmp = atof(vararg[++i]);
	    if (tmp > 0) {
		P = tmp;
	    } else {
		cerr << "Ignored invalid value for -p option (must be positive)" << endl;
	    }
	} else if (opt == string("-s")) {
	    const float tmp = atof(vararg[++i]);
	    if (tmp > 0) {
		SMOOTHING_FACTOR = tmp;
	    } else {
		cerr << "Ignored invalid value for -s option (must be positive)" << endl;
	    }
	} else if (opt == string("-d")) {
	    const float tmp = atof(vararg[++i]);
	    if (tmp > 1) {
		DOWNSCALING_FACTOR = tmp;
	    } else {
		cerr << "Ignored invalid value for -d option (must be > 1)" << endl;
	    }
	} else if (opt == string("-i1")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp > 0) {
		FIRST_ITERATIONS = tmp;
	    } else {
		cerr << "Ignored invalid value for -i1 option (must be positive)" << endl;
	    }
	} else if (opt == string("-i2")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp > 0) {
		LATER_ITERATIONS = tmp;
	    } else {
		cerr << "Ignored invalid value for -i2 option (must be positive)" << endl;
	    } 
	} else if (opt == string("-c")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp == 0 || tmp == 1 || tmp == 2) {
		USE_ORIG_IMAGE = tmp;		
	    } else {
		cerr <<" Ignored invalid value for -c option (must be 0, 1 or 2)" << endl;
	    }
	} else if (opt == string("-g")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp == 0 || tmp == 1) {
		USE_STRUCTURE_TENSOR = tmp;
	    } else {
		cerr << "Ignored invalid value for -g oiption (must be 0 or 1)" << endl;
	    }
	} else if (opt == string("-scale")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp == 0 || tmp == 1){
		USE_SCALE_INFO = tmp;
	    } else {
		cerr << "Ignored invalid value for -scale option (must be 0 or 1)" << endl;
	    }
	} else if (opt == string("-show")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp == 0 || tmp == 1) {
		SHOW_INTERMEDIARY = tmp;
	    } else { 
		cerr << "Ignored invalid value for -show option (must be 0 or 1)" << endl;
	    }
	} else if (opt == string("-nr")) {
	    const int tmp = atoi(vararg[++i]);
	    if (tmp >= 2 && tmp <= TOTAL_NUM_COLORS) {
		NUM_REG = tmp;
	    } else {
		cout << "Ignored invalid value for -nr option (must be >= 2 and not more than) " << TOTAL_NUM_COLORS << endl;
	    }
	} else if (opt == string("-init")) {
	    const int tmp = atoi(vararg[++i]);
	    INIT_SEG_TYPE = tmp;
	} else {
	    cerr << "Unknown option: " << opt << endl;
	}
    }
}


//===========================================================================
void initialize_levelsets(vector<LevelSetFunction>& phi)
//===========================================================================
{
    if (phi.size() > 1) {
	// multiregion
	if (INIT_SEG_TYPE == 0) {
	    // initialize using voronoi regions
	    random_scattered_voronoi(&phi[0], phi.size(), 1);
	} else if (INIT_SEG_TYPE < 0) {
	    // initialize using voronoi regions
	    random_scattered_voronoi(&phi[0], phi.size(), -INIT_SEG_TYPE);
	} else {
	    // initialize using horizontal bands
	    const int total_num_bands = phi.size() * INIT_SEG_TYPE;
	    const int bandwidth = phi[0].dimy() / total_num_bands;
	    multiregion_bands(&phi[0], phi.size(), bandwidth);
	}
    } else {
	// single region
	assert(phi.size() == 1);
	if (INIT_SEG_TYPE == 0) {
	    sphere(phi[0], 0.3, 0.5, 0.5);
	} else {
	    const int num_bands = INIT_SEG_TYPE > 0 ? INIT_SEG_TYPE : -INIT_SEG_TYPE;
	    horizontal_sinusoidal_bands(phi[0], num_bands);
	}
    }
}

//===========================================================================
void show_usage()
//===========================================================================
{
    cerr << "This program carries out 2-region or multiregion segmentation of arbitrary images." << endl;
    cerr << endl;
    cerr << "Usage: segmentation2D <filename> [option list]" << endl;
    cerr << endl;
    cerr << "Where options are: " << endl;
    cerr << endl;
    cerr << "-m <mask name> mask to describe active part of image to segment (must have\n";
    cerr << "               same resolution as image)  \n";
    cerr << "-mflip <bool> decide if the zero or nonzero part of a mask should be the active part\n";
    cerr << "              of segmentation.  Only matters if a mask is actually used.  \n";
    cerr << "    -- default: " << (MFLIP ? " true " : " false " ) << endl;
    cerr << "-lr <positive integer> lowest level resolution\n";
    cerr << "    -- default: " << LOWEST_RES << endl;
    cerr << "-n <pos. floating-point number> boundary smoothness factor 'nu'\n";
    cerr << "    -- default: " << NU_FACTOR << endl;
    cerr << "-p <pos. floating-point number> anisotropic factor for pre-blurring\n";
    cerr << "    -- default: " << P << endl;
    cerr << "-s <pos. floating-point number> level of presmoothing (smaller value meens \n";
    cerr << "                                more smoothing)\n";
    cerr << "    -- default: " << SMOOTHING_FACTOR << endl;
    cerr << "-d <floating-point number greater than 1> downscaling factor\n";
    cerr << "    -- default: " << DOWNSCALING_FACTOR << endl;
    cerr << "-i1 <pos. integer> number of iterations at lowest level\n";
    cerr << "    -- default: " << FIRST_ITERATIONS << endl;
    cerr << "-i2 <pos. integer> number of iterations at higher levels\n";
    cerr << "    -- default: " << LATER_ITERATIONS << endl;
    cerr << "-c <number> 2->use color info, 1->use only grayscale, 0-> use neither\n";
    cerr << "    -- default: " << USE_ORIG_IMAGE << endl;
    cerr << "-g <bool> use structure tensor information, 1->yes, 0->no \n";
    cerr << "    -- default: " << (USE_STRUCTURE_TENSOR ? "use " : "do not use ") << "structure tensor\n";
    cerr << "-scale <bool> use scale information, 1->yes, 0->no \n";
    cerr << "    -- default: " << (USE_SCALE_INFO ? "use " : "do not use ") << "scale info\n";
    cerr << "-show <bool> whether intermediate segmentations should be shown\n";			
    cerr << "    -- default: " << (SHOW_INTERMEDIARY ? "show " : "do not show ") << "intermediary results.\n";
    cerr << "-nr <number> number of distinct regions in segmentation (positive integer)\n";
    cerr << "    -- default: " << NUM_REG << endl;
    cerr << "-init <number> shape of initial segmentation.  The meaning of the number \n";
    cerr << "               depends on whether\n";
    cerr << "               we are segmenting using TWO regions or SEVERAL regions.\n";
    cerr << "               * for TWO regions, we have: \n";
    cerr << "                  0       -> centered circle\n";
    cerr << "                  n, n>0  -> number of horizontal bands\n";
    cerr << "               * for MULTI-regions, we have \n";
    cerr << "                  0       -> randomly placed voronoi-type regions\n";
    cerr << "                  n, n>0  -> number of horizontal bands\n";
    cerr << "                  n, n<0  -> randomly placed voronoi-type fragmented regions,\n";
    cerr << "                             where |n| is the number of fragments of each region\n";
    cerr << "    -- default: " << INIT_SEG_TYPE << endl;
    cerr << "-save <mask filename> Save the segmented region as a mask (image containing only \n";
    cerr << "                      white or black pixels).  Only works for 2-region segmentation.\n";
    cerr << "                      The mask will be saved as a PNG file, and the extension '.png' will\n";
    cerr << "                      be added to the filename given.\n";
    cerr << "    -- default: do not save" << endl;
    cerr << endl;
};
