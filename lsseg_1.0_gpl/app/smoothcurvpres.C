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
#include <fstream>
#include <limits>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <boost/tokenizer.hpp>

#include "Image.h"
#include "Filters.h"
#include "LIC.h"
#include "simple_tools.h"
#include"cimg_dependent.h"

using namespace std;
using namespace lsseg;

const double PI = 3.1415926535897932384;

// argument that may be overridden by command line
/// timestep to use for \ref anchor_LIC "line integral convolution".
double DT = 30; 
/// spatial step to use for \ref anchor_LIC "line integral convolution".
double DL = 0.8; 
/// amount of isotropic blurring of image \em before computing structure tensor
double ALPHA = 1.2; 
/// amount of gaussian blur of structure tensor before computing the smoothing geometry T.
double SIGMA = 1; 
/// number of timesteps to use
int NUM_TIMESTEPS = 1;
/// factor to use for determining amount of smoothing along major smoothing direction when computing the 
/// smoothing geometry T.
double P1 = 0.3; 
/// factor to use for determining amount of smoothing along second smoothing direction when computing the 
/// smoothing geometry T.
double P2 = 0.7;
/// factor to use for determining amount of smoothing along third smoothing direction when computing the
/// smoothing geometry T (only used for 3D). 
double P3 = 0.9; 
/// length of window for \ref anchor_LIC "line integral convolution"
double L = 2 * sqrt(DT); 
/// how many longitude angular steps to use when decomposing the smoothing geometry tensor field T.
const int DIV_THETA = 6; 
/// size of longitude step when decomposing the smoothing geometry T.
const double DTHETA = PI / DIV_THETA;
/// how many latitude angular steps to use when decomposing the smoothing geometry tensor field T (only used in 3D)
const int DIV_PHI = 6; 
/// size of latitude step when decomposing the smoothing geometry T (3D only)
const double DPHI = PI / DIV_PHI;
/// where to save the result
const char SAVED_3D_FILE[] = "greycstoration_result.stack";  

double C = double(1) / sqrt(4 * PI * DT); // factor used in gauss
double D = double(1) / (4 * DT);

//===========================================================================
void make_vector_field_2D(double angle, const Image<double>& tensor, Image<double>& vec);
//===========================================================================

//===========================================================================
void make_vector_field_3D(double phi, double theta, const Image<double>& tensor, Image<double>& vec);
//===========================================================================

//===========================================================================
void show_command_info();
//===========================================================================

//===========================================================================
void show_current_vars();
//===========================================================================

//===========================================================================
void read_command_info(int varnum, char** vararg);
//===========================================================================

//===========================================================================
double gauss(double t) 
//===========================================================================
{
    //return 1; // @@
    const double res = C * exp(-(t*t) * D);
    return res;
}

//===========================================================================
bool refers_to_stack(string str)
//===========================================================================
{
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(".");
    tokenizer tokens(str, sep);
    string saved;
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it) {
	saved = *it;
    }
    if (saved == string("stack") || saved == string("Stack") || saved == string("STACK")) {
	return true;
    }
    return false;
}

//===========================================================================
void to_field(const Image<double>& T, Image<double>& W)
//===========================================================================
{
    const int X = T.dimx();
    const int Y = T.dimy();
    const int Z = T.dimz();
    W.resize(X, Y, Z, 4);
    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		double wx = T(x, y, z, 0);
		double wy = T(x, y, z, 1);
		double wz = T(x, y, z, 2);
		double norm = sqrt(wx * wx + wy * wy + wz*wz);
		double norm_inv = (norm == 0) ? 1 : double(1) / norm;
		W(x, y, z, 0) = wx * norm_inv;
		W(x, y, z, 1) = wy * norm_inv;
		W(x, y, z, 2) = wz * norm_inv;
		W(x, y, z, 3) = norm;
	    }
	}
    }
}

void make_debug_field(const Image<double>& T, Image<double>& W)
{
    //to_field(T, W); return; //@

    const int X = T.dimx();
    const int Y = T.dimy();
    const int Z = T.dimz();
    W.resize(X, Y, Z, 4);
    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		double vx = (x - X/2+0.5);
		double vy = (y - Y/2 +0.5);
		double vz =(z - Z/2 + 0.5);
		double norm = sqrt(vx * vx + vy * vy + vz * vz);
		vx /= norm;
		vy /= norm;
		vz /= norm;
		
		W(x, y, z, 0) = vx;
		W(x, y, z, 1) = vy;
		W(x, y, z, 2) = vz;
		W(x, y, z, 3) = 1;
	    }
	}
    }
    
}

//===========================================================================
int main(int varnum, char** vararg)  
//===========================================================================
{
    if (varnum == 1) {
	show_command_info();
	return -1;
    } 
    read_command_info(varnum, vararg);
    cout << "Now running with the following options: " << endl;
    show_current_vars();
    
    // Setting global values necessary for the gauss function or LIC
    L = 2 * sqrt(DT); // // parametric length for LIC 
    C = double(1) / sqrt(4 * PI * DT);
    D = double(1) / (4 * DT);

    clock_t start_time = clock();
    // loading source image
    Image<double> img;
    if (refers_to_stack(vararg[1])) {
	cout << "Reading 3D image stack...";
	ifstream is(vararg[1]);
	if (!is) {
	    cerr << "Could not open file " << vararg[1] << endl;
	    return 1;
	}
	img.read(is);
	cout << "Finished!" << endl;
    } else {
	cout << "Image is 2D " << endl;
	load_image(vararg[1], img, false);
    }
    permanent_display(img);

    const int X = img.dimx();
    const int Y = img.dimy();
    const int Z = img.dimz();
    const bool three_dimensional = Z > 1; 
    const int tensor_components = three_dimensional ? 6 : 3;
    const int vector_components = three_dimensional ? 4 : 3;

    Image<double> target(img, false); // no copying of contents
    target = 0;

    cout << "Allocating memory.." << endl;
    Image<double> G(X, Y, Z, tensor_components); // structure tensor
    Image<double> T(X, Y, Z, tensor_components); // smoothing tensor
    Image<double> W(X, Y, Z, vector_components); // smoothing vector field
    Image<double> img_blurred(img, false);
    cout << "Finished allocating memory" << endl;
    const double norm_fac = 
	three_dimensional ? double(1) / (DIV_THETA * DIV_PHI) : double(1) / DIV_THETA;

    for (int tstep = 0; tstep < NUM_TIMESTEPS; ++tstep) {
	img_blurred = img;
	cout << "Bluring image" << endl;
	blur_image(img_blurred, ALPHA);
	cout << "Rescaling image" << endl;
	rescale(img_blurred, 0, 255);
	//permanent_display(img_blurred); //@@

	if (three_dimensional) {
	    cout << "Computing structure tensor" << endl;
	    compute_structure_tensor_3D(img_blurred, G); 
	    cout << "Bluring structure tensor" << endl;
	    blur_image(G, SIGMA); 
	    
	    //permanent_display(G); //@@
	    cout << "Computing smoothing geometry " << endl;
	    compute_smoothing_geometry_3D(G, P1, P2, P3, T, true);
	    //permanent_display(T);
	    // decomposing T into a sum of vectors and carry out smoothing

	    
	    for (int phistep = 0; phistep < DIV_PHI; ++phistep) {
		for (int thetastep = 0; thetastep < DIV_THETA; ++thetastep) {
		    const double phi = phistep * DPHI - PI/2;
		    const double theta = thetastep * DTHETA;
		    make_vector_field_3D(phi, theta, T, W);
		    
		    // line integral convolution
		    cout << "Convoluting with PHI = " << phi << " and THETA = " << theta << endl;
		    LIC_3D_FS(img, W, target, &gauss, DL, L);
		}
	    }
	    // normalize target
	    target *= norm_fac;
	} else {

	    // 2D case
	    compute_structure_tensor_2D(img_blurred, G);       
	    blur_image(G, SIGMA);
	    compute_smoothing_geometry_2D(G, P1, P2, T, true); 
	    // decomposing T into a sum of vectors and carry out smoothing
	    for (int a = 0; a < DIV_THETA; ++a) {
		const double angle = a * DTHETA;
		make_vector_field_2D(angle, T, W);
		
		// line integral convolution
		LIC_2D_FS(img, W, target, &gauss, DL, L);
	    }
	    target *= norm_fac;
	}
	// swap images
	target.swap(img);
    }
    clock_t end_time = clock();
    //permanent_display(img);
    display_image(img);
    if (three_dimensional) {
	ofstream os(SAVED_3D_FILE);
	if (!os) {
	    cerr << "Warning: not able to save 3D result because unable to open file." << endl;
	} else {
	    img.write(os);
	    os.close();
	}
    }
    cout << "Total CPU time was: " << double(end_time - start_time) / 1000 << " millisecs." << endl;
    return 1;
};

//===========================================================================
void make_vector_field_3D(double phi, double theta, const Image<double>& tensor, Image<double>& vec)
//===========================================================================
{
    const double POLE_CORRECTION_FACTOR = sqrt(0.5);
    const double TINY = numeric_limits<double>::epsilon();
    assert(tensor.numChannels() == 6);
    const int X = tensor.dimx();
    const int Y = tensor.dimy();
    const int Z = tensor.dimz();
    const double v_x = cos(phi) * cos(theta);
    const double v_y = cos(phi) * sin(theta);
    const double v_z = POLE_CORRECTION_FACTOR * sin(phi);
    vec.resize(X, Y, Z, 4);
    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		const double txx = tensor(x, y, z, 0);
		const double txy = tensor(x, y, z, 1);
		const double tyy = tensor(x, y, z, 2);
		const double txz = tensor(x, y, z, 3);
		const double tyz = tensor(x, y, z, 4);
		const double tzz = tensor(x, y, z, 5);

		const double vx = txx * v_x + txy * v_y + txz * v_z;
		const double vy = txy * v_x + tyy * v_y + tyz * v_z;
		const double vz = txz * v_x + tyz * v_y + tzz * v_z; 
		const double norm = sqrt(vx * vx + vy * vy + vz * vz);
		const double norm_inv = (norm > TINY) ? double(1) / norm : 0;
		vec(x, y, z, 0) = vx * norm_inv;
		vec(x, y, z, 1) = vy * norm_inv;
		vec(x, y, z, 2) = vz * norm_inv;
		vec(x, y, z, 3) = norm;

		assert(vec(x, y, z, 0) <= 1 && vec(x, y, z, 0) >= -1);
		assert(vec(x, y, z, 1) <= 1 && vec(x, y, z, 1) >= -1);
		assert(vec(x, y, z, 2) <= 1 && vec(x, y, z, 2) >= -1);
	    }
	}
    }
}

//===========================================================================
void make_vector_field_2D(double theta, const Image<double>& tensor, Image<double>& vec)
//===========================================================================
{
    assert(tensor.dimz() == 1 && tensor.numChannels() == 3);
    const int X = tensor.dimx();
    const int Y = tensor.dimy();
    const double sin_theta = sin(theta);
    const double cos_theta = cos(theta);

    vec.resize(X, Y, 1, 3);

    for (int y = 0; y < Y; ++y) {
	for (int x = 0; x < X; ++x) {
	    const double a = tensor(x, y, 0, 0);
	    const double b = tensor(x, y, 0, 1);
	    const double c = tensor(x, y, 0, 2);
	    const double vx = a * cos_theta + b * sin_theta;
	    const double vy = b * cos_theta + c * sin_theta;
	    const double norm = sqrt(vx*vx + vy * vy);
	    const double norm_inv = double(1) / norm;
	    vec(x, y, 0, 0) = vx * norm_inv;
	    vec(x, y, 0, 1) = vy * norm_inv;
	    vec(x, y, 0, 2) = norm;
	}
    }
}


//===========================================================================
void show_command_info()
//===========================================================================
{
    cerr << "Usage: greycstoration2D <image> <option(s)>" << endl;
    cerr << "Where options are: " << endl;
    cerr << "-DT <timestep value (pos. double)>  -- default: " << DT << endl;
    cerr << "-ALPHA <smooth value(pos. double)>" << endl;
    cerr << "       -> how much to smooth the image prior to computing" << endl;
    cerr << "          the structure tensor G.   -- default: " << ALPHA << endl;
    cerr << "-SIGMA <smooth value (pos.double)>" << endl;
    cerr << "       -> how much to smooth the structure tensor G prior" << endl;
    cerr << "          to line integral convolution." << endl;
    cerr << "                                    -- default: " << SIGMA << endl;
    cerr << "-TIMESTEPS <positive integer>" << endl;
    cerr << "       -> how many timesteps (of length DT) that should be " << endl;
    cerr << "          computed.                 -- default: " << NUM_TIMESTEPS << endl;
    cerr << "-P1 <fac 1 (pos. double)> " << endl;
    cerr << "       -> specify how to smooth along the main smoothing direction." << endl;
    cerr << "          (approx. along isophote).  Higher values yield less" << endl;
    cerr << "          smoothing.                -- default: " << P1 << endl;
    cerr << "-P2 <fac 2 (pos. double)> " << endl;
    cerr << "       -> specify how to smooth along the minor smoothing direction. " << endl;
    cerr << "          (approx. along image gradient).  Higher values yield less" << endl;
    cerr << "          smoothing.  Should be higher than P1." << endl;
    cerr << "                                    -- default: " << P2 << endl;
    cerr << "-P3 <fac 3 (pos. double)> " << endl;
    cerr << "       -> specify how to smooth along the minor smoothing direction (3D). " << endl;
    cerr << "          Higher values yield less smoothing.  Should be higher than P2." << endl;
    cerr << "                                    -- default: " << P3 << endl;
    cerr << endl;
    cerr << "The <image> argument can be the name of a picture in any of the most common ";
    cerr << "formats (jpeg, tiff, png, etc.).  In that case, the 2D image will be read ";
    cerr << "and treated by the algorithm in a 2D way.  You can also give the algorithm ";
    cerr << "a stack of images, and it will process this stack as one 3D image.  A ";
    cerr << "stack file has 'stack' as its suffix, and can be made using the ";
    cerr << "'generate_stack' utility program (run it to get an explanation of " << endl;
    cerr << "its use.  In case a 3D stack is to be processed, the result is ALWAYS saved ";
    cerr << "to the file \"" << SAVED_3D_FILE << "\"." << endl;
}

//===========================================================================
void read_command_info(int varnum, char** vararg)
//===========================================================================
{
    for (int i = 2; i != varnum; i+=2) {
	string arg(vararg[i]);
	if (arg=="-DT") {
	    DT = atof(vararg[i+1]);
	} else if (arg=="-ALPHA") {
	    ALPHA = atof(vararg[i+1]);
	} else if (arg=="-SIGMA") { 
	    SIGMA = atof(vararg[i+1]);
	} else if (arg=="-TIMESTEPS") {
	    NUM_TIMESTEPS = atoi(vararg[i+1]);
	} else if (arg=="-P1") {
	    P1 = atof(vararg[i+1]);
	} else if (arg=="-P2") {
	    P2 = atof(vararg[i+1]);
	} else if (arg=="-P3") {
	    P3 = atof(vararg[i+1]);
	} else {
	    cerr << "Unrecognized option: " << arg << endl;
	    exit(-1);
	}
    }
}

//===========================================================================
void show_current_vars()
//===========================================================================
{
    cout << "--------------------------" << endl;
    cout << "DT:             " << DT << endl;
    cout << "ALPHA:          " << ALPHA << endl;
    cout << "SIGMA:          " << SIGMA << endl;
    cout << "NUM_TIMESTEPS:  " << NUM_TIMESTEPS << endl;
    cout << "P1:             " << P1 << endl;
    cout << "P2:             " << P2 << endl;
    cout << "P3:             " << P3 << endl;
    cout << "--------------------------" << endl;
}
