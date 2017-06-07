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
//===========================================================================
//                                                                           
// File: simple_tools.C                                                      
//                                                                           
// Created: Tue Oct 25 13:56:16 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: simple_tools.C,v 1.32 2006/11/25 20:08:28 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements simple_tools.h.
//                                                                           
//===========================================================================

// recent changes:
// 2006/03/30:	added function read_image_sequence_multichn() (kgre)

#include <iostream>
#include <vector>
#include "simple_tools.h"
#include "cimg_dependent.h"
#include "errormacros.h"
#include "Image.h"
//#include "GoBorrowedMVGrid.h"
//#include "GoTensorProductSpline.h"
//#include "randomnoise.h"

using namespace std;
//using namespace Go;

namespace {
    inline bool in_zero_one_interval(double d) { return (d >= 0) && (d <= 1);}

    // function taken from the Sintef Applied Mathematics GoTools library ('utils' module)
    inline void uniformNoise(double* res, double lval, double uval, int num_samples)
    {
	if (uval <= lval) {
	    throw runtime_error("uniformNoise(...) : erroneous range.");
	}
	double range = uval - lval;
	double scale_factor = range / double(RAND_MAX);
	for (int i = 0; i < num_samples; ++i) {
	    res[i] = double(rand()) * scale_factor + lval;
	}
    }

}; // end anonymous namespace

namespace lsseg {

//===========================================================================
void negate(Image<double>& img, double min, double max)
//===========================================================================
{
    double* it;
    for (it = img.begin(); it != img.end(); ++it) {
	*it = max - (*it-min);
    }
}

//===========================================================================
void rescale(Image<double>& img, 
	     double cur_min, double cur_max, double to_min, double to_max)
//===========================================================================
{
    double cur_int_inv = double(1) / (cur_max - cur_min);
    double* data = img.begin();
    for (size_t i = 0; i < img.size(); ++i) {
	double t = (data[i] - cur_min) * cur_int_inv;   
	assert(t >= 0);
	data[i] = t * to_max + (1-t) * to_min;
    }
}

//===========================================================================
void rescale_channels(Image<double>& img, double to_min, double to_max)
//===========================================================================
{
    int chan_size = img.channelSize();

    // looping over channels
    for (int c = 0; c < img.numChannels(); ++c) {
	double cur_min = *min_element(img.begin() + c * chan_size,
				      img.begin() + (c+1) * chan_size);
	double cur_max = *max_element(img.begin() + c * chan_size,
				      img.begin() + (c+1) * chan_size);

	double cur_int_inv = double(1) / (cur_max - cur_min);

	double* data = img.begin();
	for (int i = 0; i < chan_size; ++i) {
	    double& cur = data[c * chan_size + i];
	    double t = (cur - cur_min) * cur_int_inv;
	    cur = t * to_max + (1-t) * to_min;
	    data[c * chan_size + i] = cur;
	}
    }
}

//===========================================================================
void clip(Image<double>& img, double min, double max)
//===========================================================================
{
    double* data = img.begin();
    for (int i = 0; i < int(img.size()); ++i) {
	double tmp = data[i];
	if (tmp < min) {
	    data[i] = min;
	} else if (tmp > max) {
	    data[i] = max;
	}
    }
}

//===========================================================================
void transpose(Image<double>& img)
//===========================================================================
{
    Image<double> tmp;
    transpose(img, tmp);
    img.swap(tmp);
}

//===========================================================================
void transpose(const Image<double>& img, Image<double>& target)
//===========================================================================
{
    target.resize(img);
    for (int c = 0; c < img.numChannels(); ++c) {
	for (int z = 0; z < img.dimz(); ++z) {
	    for (int y = 0; y < img.dimy(); ++y) {
		for (int x = 0; x < img.dimx(); ++x) {
		    target(y, x, z, c) = img(x, y, z, c);
		}
	    }
	}
    }
}

//===========================================================================
void horizontal_sinusoidal_bands(LevelSetFunction& img, int num_bands, double phase)
//===========================================================================
{
    const double PI = 3.1415926535897932384;
    for (int h = 0; h < img.dimy(); ++h) {
	double angle = (double(h) / double(img.dimy()) + phase) * 2 * PI;
	double value = sin(angle * num_bands);
	for (int z = 0; z < img.dimz(); ++z) {
	    for (int x = 0; x < img.dimx(); ++x) {
		img(x, h, z, 0) = value; // value > 0 ? 1 : -1
	    }
	}
    }
}

//===========================================================================
void rectangle(LevelSetFunction& img, 
	       double xmin_ratio, 
	       double xmax_ratio,
	       double ymin_ratio,
	       double ymax_ratio,
	       double zmin_ratio,
	       double zmax_ratio)
//===========================================================================
{
    assert(img.numChannels() == 1);

    // setting image white
    std::fill(img.begin(), img.end(), 1);//127);

    // determining square corners
    int xmin = int(img.dimx() * xmin_ratio);
    int xmax = int(img.dimx() * xmax_ratio);
    int ymin = int(img.dimy() * ymin_ratio);
    int ymax = int(img.dimy() * ymax_ratio);
    int zmin = int(img.dimz() * zmin_ratio);
    int zmax = int(img.dimz() * zmax_ratio);

    for (int z = zmin; z < zmax; ++z) {
	for (int y = ymin; y < ymax; ++y) {
	    for (int x = xmin; x < xmax; ++x) {
		img(x, y, z) = -1;
	    }
	}
    }
}

//===========================================================================
void init_voronoi_regions(LevelSetFunction* regs,
			  const double* center_coords,
			  int num_regions,
			  bool three_d)
//===========================================================================
{
    assert(num_regions > 0);
    const int X = regs[0].dimx();
    const int Y = regs[0].dimy();
    const int Z = regs[0].dimz();
    assert(three_d || Z == 1);
    for (int i = 0; i < num_regions; ++i) {
 	assert(regs[i].dimx() == X && regs[i].dimy() == Y && regs[i].dimz() == Z);
    }
    vector<double> dists(num_regions);
    const int stride = three_d ? 3 : 2;
    vector<double> abs_coords(num_regions * stride);
    for (int i = 0; i < num_regions; ++i) {
	abs_coords[i * stride] = center_coords[i * stride] * X;
	abs_coords[i * stride + 1] = center_coords[i * stride + 1] * Y;
	if (stride == 3) {
	    abs_coords[i * stride + 2] = center_coords[i * stride + 2] * Z;
	}
    }

    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		int closest_reg = -1; // uninitialized
		double closest_dist2 = (X+Y+Z) * (X+Y+Z); // surely a number exceeding the image sizes
		for (int r = 0; r < num_regions; ++r) {
		    double dx = x - abs_coords[r * stride + 0];
		    double dy = y - abs_coords[r * stride + 1];
		    double dz = (three_d) ? z - abs_coords[r * stride + 2] : 0;
		    double cur_dist2 = dx * dx + dy * dy + dz * dz;
		    
		    if (cur_dist2 < closest_dist2) {
			closest_dist2 = cur_dist2;
			closest_reg = r;
		    }
		}
		assert(closest_reg != -1);
		// filling this pixel in all functions
		for (int r = 0; r < num_regions; ++r) {
		    regs[r](x, y) = (r == closest_reg) ? -1 : 1;
		}
	    }
	}
    }
    for (int r = 0; r < num_regions; ++r) {
	if (three_d) {
	    regs[r].reinitialize3D();
	} else {
	    regs[r].reinitialize2D();
	}
    }
}

//===========================================================================
void multiregion_bands(LevelSetFunction* regs,
		       int num_regs,
		       int pixel_bandwith)
//===========================================================================
{
    assert(num_regs > 1);
    const int X = regs[0].dimx();
    const int Y = regs[0].dimy();
    const int Z = regs[0].dimz();
    for (int i = 0; i < num_regs; ++i) {
	assert(regs[i].dimx() == X && regs[i].dimy() == Y && regs[i].dimz() == Z);
    }
    int current_band = 0;
    int b_count = 0;
    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int r = 0; r < num_regs; ++r) {
		const double val = (r == current_band) ? -1 : 1;
		for (int x = 0; x < X; ++x) {
		    regs[r](x, y, z) = val;
		}
	    }
	    ++b_count;
	    if (b_count == pixel_bandwith) {
		b_count = 0;
		++current_band;
		if (current_band == num_regs) {
		    current_band = 0;
		}
	    }
	}
    }
    for (int r = 0; r < num_regs; ++r) {
	if (Z == 1) {
	    regs[r].reinitialize2D();
	} else {
	    regs[r].reinitialize3D();
	}
    }
}

//===========================================================================
void random_scattered_voronoi(LevelSetFunction* regs,
			      int num_regs,
			      int num_fragments) // per region
//===========================================================================
{
    assert(num_regs > 1);
    const int X = regs[0].dimx();
    const int Y = regs[0].dimy();
    const int Z = regs[0].dimz();
    for (int i = 0; i < num_regs; ++i) {
	assert(regs[i].dimx() == X && regs[i].dimy() == Y && regs[i].dimz() == Z);
    }
    const bool three_d = (Z == 1);
    const int stride = three_d ? 3 : 2;
    
    vector<vector<double> > seed_abs_coords(num_regs);
    for (int i = 0; i < num_regs; ++i) {
	for (int fr = 0; fr < num_fragments; ++fr) {
	    double seed_x, seed_y, seed_z;
	    uniformNoise(&seed_x, 0, X, 1);
	    uniformNoise(&seed_y, 0, Y, 1);
	    seed_abs_coords[i].push_back(seed_x);
	    seed_abs_coords[i].push_back(seed_y);
	    if (three_d) {
		uniformNoise(&seed_z, 0, Z, 1);
		seed_abs_coords[i].push_back(seed_z);
	    }
	}
    }

    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		int closest_reg = -1; // uninitialized;
		double closest_dist2 = (X+Y+Z) * (X+Y+Z); // surely a number exceeding the image sizes
		for (int r = 0; r < num_regs; ++r) {
		    for (int fr = 0; fr < num_fragments; ++fr) {
			const double dx = x - seed_abs_coords[r][stride * fr + 0];
			const double dy = y - seed_abs_coords[r][stride * fr + 1];
			const double dz = three_d ? z - seed_abs_coords[r][stride * fr + 2] : 0;
			double cur_dist2 = dx * dx + dy * dy + dz * dz;
			if (cur_dist2 < closest_dist2) {
			    closest_dist2 = cur_dist2;
			    closest_reg = r;
			}
		    }
		}
		assert(closest_reg != -1);
		// filling this pixel in all functions
		for (int r = 0; r < num_regs; ++r) {
		    regs[r](x, y, z) = (r == closest_reg) ? -1 : 1;
		}
	    }
	}
    }
    for (int r = 0; r < num_regs; ++r) {
	if (three_d) {
	    regs[r].reinitialize3D();
	} else {
	    regs[r].reinitialize2D();
	}
    }
}

//===========================================================================
void set_from_parzen(LevelSetFunction& phi, 
		      const ParzenDistributionForce pf,
		      const Mask* m)
//===========================================================================
{
    assert(pf.baseImage()->spatial_compatible(phi));

    if (!m) {
	for (size_t ix = 0; ix < phi.size(); ++ix) {
	    if (pf.force(ix) < 0) {
		phi[ix] = 1;
	    } else {
		phi[ix] = -1;
	    }
	}
    } else {
	phi = 1;
	for (size_t ix = 0; ix < phi.size(); ++ix) {
	    if ((*m)[ix] && pf.force(ix) > 0) {
		phi[ix] = -1;
	    }
	}
    }
    phi.reinitialize3D(m);
}

//===========================================================================
void sphere(LevelSetFunction& img,
	    double relrad,  // between 0 and 1
	    double xrelpos, // between 0 and 1
	    double yrelpos, // between 0 and 1
	    double zrelpos)
//===========================================================================
{
    assert(img.numChannels() == 1);
    
    const int xmax = img.dimx();
    const int ymax = img.dimy();
    const int zmax = img.dimz();

    const int xpos = int (xmax * xrelpos);
    const int ypos = int (ymax * yrelpos);
    const int zpos = int (zmax * zrelpos);

    int maxim = xmax > ymax ? xmax : ymax;
    maxim = maxim > zmax ? maxim : zmax;
    int radius2 = int(relrad * maxim);
    radius2 *= radius2; // working with square of radius
    
    for (int z = 0; z < zmax; ++z) {
	for (int y = 0; y < ymax; ++y) {
	    for (int x = 0; x < xmax; ++x) {
		double dx = (xpos - x);
		double dy = (ypos - y);
		double dz = (zpos - z);
		double dist2 = dx * dx + dy * dy + dz * dz;
		img(x, y, z) = (dist2 > radius2) ? sqrt(dist2 - radius2) : -sqrt(radius2 - dist2); //(dist2 < radius2) ? -1 : 1;
	    }
	}
    }
}



// //===========================================================================
// // (kgre) changed type of first parameter from Image<double>& to
// // LevelSetFunction& as in the fuction declaration
// void distance_from_point(LevelSetFunction& img,
// 			 double zero_level, 
// 			 double point_x_relpos,
// 			 double point_y_relpos,
// 			 double point_z_relpos)
// //===========================================================================
// {
//     assert(in_zero_one_interval(zero_level));
//     assert(in_zero_one_interval(point_x_relpos));
//     assert(in_zero_one_interval(point_y_relpos));
//     assert(in_zero_one_interval(point_z_relpos));

//     // setting image white
//     std::fill(img.begin(), img.end(), 1);

//     // determining longest direction, central point position and zero level distance
//     int X = img.dimx();
//     int Y = img.dimy();
//     int Z = img.dimz();
//     int longest_dir = (X > Y) ? X : Y;
//     longest_dir = (longest_dir < Z) ? Z : longest_dir;
//     int px = int(point_x_relpos * X);
//     int py = int(point_y_relpos * Y);
//     int pz = int(point_z_relpos * Z);
//     (px == X) ? px-- : px; // adjust case where px is too big...
//     (py == Y) ? py-- : py;
//     (pz == Z) ? pz-- : pz;
//     zero_level *= 0.5 *  longest_dir;
//     for (int z = 0; z < Z; ++z) {
// 	for (int y = 0; y < Y; ++y) {
// 	    for (int x = 0; x < X; ++x) {
// 		double x2 = (x - px); x2 *= x2;
// 		double y2 = (y - py); y2 *= y2;
// 		double z2 = (z - pz); z2 *= z2;
// 		img(x, y, z) = sqrt(x2 + y2 + z2) - zero_level;
// 	    }
// 	}
//     }
// }

//===========================================================================
void make_border_mask(const LevelSetFunction& phi, Mask& target, int bandwidth, const Mask* geom_mask)
//===========================================================================
{
    ALWAYS_ERROR_IF(phi.numChannels() != 1, "wrong number of channels");
    // setting complete mask to 0 before starting
    target.resize(phi.dimx(), phi.dimy(), phi.dimz(), 1);
    fill(target.begin(), target.end(), 0);

    ALWAYS_ERROR_IF(bandwidth < 1, "zero or negative bandwidth");
	   
    // detecting zero-level
    const int X = phi.dimx();
    const int Y = phi.dimy();
    const int Z = phi.dimz();
    const int bm1 = bandwidth - 1;
    const int Xmb = X - bandwidth;
    const int Ymb = Y - bandwidth;
    const int Zmb = Z - bandwidth;

    const Mask::value_type* m_it = geom_mask ? geom_mask->begin() : 0;

    for (int z = 0; z < Z; ++z) {
	const int zn =  (z < Z-1) ? z+1 : z;
	for (int y = 0; y < Y; ++y) {
	    const int yn = (y < Y-1) ? y+1 : y;
	    for (int x = 0; x < X; ++x) {
		if (geom_mask && !(*m_it++)) {
		    continue; // geometry mask is zero for this pixel
		}
		const int xn = (x < X-1) ? x+1 : x;

		const double Iccc = phi( x,  y,  z, 0);
		const double Incc = phi(xn,  y,  z, 0);
		const double Icnc = phi( x, yn,  z, 0);
		const double Iccn = phi( x,  y, zn, 0);

		if (Iccc * Incc <= 0 || Iccc * Icnc <= 0 || Iccc * Iccn <= 0) {

		    const int xstart = (x >= bm1) ? x - bm1 : 0;
		    const int xend = (x < Xmb) ? x + bandwidth : X - 1;

		    const int ystart = (y >= bm1) ? y - bm1 : 0;
		    const int yend = (y < Ymb) ? y + bandwidth : Y-1;

		    const int zstart = (z >= bm1) ? z - bm1 : 0;
		    const int zend = (z < Zmb) ? z + bandwidth : Z-1;

		    for (int zz = zstart; zz <= zend; ++zz) {
			for (int yy = ystart; yy <= yend; ++yy) {
			    for (int xx = xstart; xx <= xend; ++xx) {
				target(xx, yy, zz) = 1;
			    }
			}
		    }
		}
	    }
	}
    }
}

//===========================================================================
void mask_from_segmentation(const LevelSetFunction& phi, 
			    Mask& target, 
			    SEG_REGION reg)
//===========================================================================
{
    target.resize(phi);
    const double* in = phi.begin();
    char* out = target.begin();
    
    switch (reg) {
    case SEG_NEGATIVE:
	for (;in != phi.end(); ++in, ++out) {
	    *out = *in < 0;
	}
	break;
    case SEG_POSITIVE:
	for (;in != phi.end(); ++in, ++out) {
	    *out = *in > 0;
	}
	break;
    default:
        ALWAYS_ERROR_IF(true, "logical error in mask_from_segmentation");
    }
}


//===========================================================================
void read_image_sequence(istream& image_list, 
			 Image<double>& result,
			 bool convert_to_grayscale)
//===========================================================================
{
    vector<string> filenames;
    
// (kgre) changed type of cur_name to char* for Windows
#if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
	char cur_name[_MAX_PATH];
#else
    string cur_name;
#endif

    while (image_list >> cur_name) {
	filenames.push_back(cur_name);
    }
    if (filenames.size() == 0) {
	return;
    }
    cout << "Total number of images to read: " << filenames.size() << endl;

    // opening first image in sequence
    Image<double> im;
    load_image(filenames[0].c_str(), im, convert_to_grayscale);
    assert(im.dimz() == 1);
    
    result.resize(im.dimx(), im.dimy(), int(filenames.size()), im.numChannels());

    int total_size = im.size();
    double* cur_pos = result.begin();
    copy(im.begin(), im.end(), cur_pos);
    cur_pos += total_size;
    for (size_t fnum = 1; fnum < filenames.size(); ++fnum) {
	cout << "Now reading image: " << fnum << endl;
	load_image(filenames[fnum].c_str(), im, convert_to_grayscale);
	assert(im.size() == total_size);
	copy(im.begin(), im.end(), cur_pos);
	cur_pos += total_size;
    }
}


// //===========================================================================
// void read_image_sequence_multichn(istream& image_list, 
// 				  Image<double>& result,
// 				  bool convert_to_grayscale)
// //===========================================================================
// {
// #if defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
//     char cur_name[_MAX_PATH];
// #else
//     string cur_name;
// #endif
//     vector<string> filenames;
    
//     while (image_list >> cur_name) {
// 	filenames.push_back(cur_name);
//     }
//     cout << "Total number of images to read: " << filenames.size() << endl;

//     // opening first image in sequence
//     Image<double> im;
//     load_image(filenames[0].c_str(), im, convert_to_grayscale);
//     assert(im.dimz() == 1);
    
//     result.resize(im.dimx(), im.dimy(), 1, int(filenames.size()));

//     int total_size = im.size();
//     copy(im.begin(), im.end(), result.channelBegin(0));

//     for (int fnum = 1; fnum < filenames.size(); ++fnum) {
// 	cout << "Now reading image: " << fnum << endl;
// 	load_image(filenames[fnum].c_str(), im, convert_to_grayscale);
// 	assert(im.size() == total_size);
// 	copy(im.begin(), im.end(), result.channelBegin(fnum));
//     }
// }


//===========================================================================
unsigned long int nonzeroes(const Image<int>& img)
//===========================================================================
{
    unsigned long int res = 0;
    for (const int* dp = img.begin(); dp != img.end(); ++dp) {
	if (*dp != 0) {
	    ++res;
	}
    }
    return res;
}

//===========================================================================
double nonzero_ratio(const Image<int>& img)
//===========================================================================
{
    unsigned long int total_size = img.size();
    
    return double(nonzeroes(img)) / double(total_size);
}

//===========================================================================
unsigned long int positives(const Image<double>& img)
//===========================================================================
{
    unsigned long int res = 0;
    for (unsigned long int i =0 ; i < img.size(); ++i) {
	if (img[i] >= 0) {
	    ++res;
	}
    }
    return res;
}

//===========================================================================
unsigned long int negatives(const Image<double>& img)
//===========================================================================
{
    return img.size() - positives(img);
}

//===========================================================================
double positive_ratio(const Image<double>& img)
//===========================================================================
{
    return double(positives(img)) / double(img.size());
}

//===========================================================================
double negative_ratio(const Image<double>& img)
//===========================================================================
{
    return double(negatives(img)) / double(img.size());
}

}; // end namespace lsseg
