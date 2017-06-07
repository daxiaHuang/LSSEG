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
// File: simple_tools_templates.h                                            
//                                                                           
// Created: Wed Dec 14 09:16:15 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: simple_tools_templates.h,v 1.6 2006/11/25 20:08:24 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements the inline functions specified in \ref simple_tools.h.
//                                                                           
//===========================================================================

#ifndef _SIMPLE_TOOLS_TEMPLATES_H
#define _SIMPLE_TOOLS_TEMPLATES_H

//#include "GoTensorProductSpline.h"
#include <stdexcept>

namespace {
//===========================================================================
template<typename T>
T interpolate_value(const lsseg::Image<T>& img, 
		    int channel,
		    int x_pixel,
		    int y_pixel,
		    int z_pixel,
		    double x_off,
		    double y_off,
		    double z_off,
		    bool linear);
//===========================================================================
};

namespace lsseg {

//==============================================================================
template<typename ImgType>
void downsample_series(const ImgType& input, 
		       std::vector<ImgType >& result,
		       int min_num_pixels,
		       bool downscale_z,
		       double factor,
		       bool only_grayscale,
		       bool linear)
//==============================================================================
{
    // downsample in x and y (and possibly z)
    assert(input.dimx() >= min_num_pixels &&  input.dimy() >= min_num_pixels);
    if (downscale_z) {
	assert(input.dimz() >= min_num_pixels);
    }
    if (only_grayscale) {
	assert(input.numChannels() == 1 || input.numChannels() == 3);
    }
    assert(factor > 1);

    result.clear();
    result.push_back(input);
    if (only_grayscale && result[0].numChannels() != 1) {
	assert(result[0].numChannels() == 3);
	to_grayscale(result[0]);
    }
    if (factor <= 1) {
	return; // nothing more to do
    }
    
    int img_size[4];
    img_size[0] = result.front().dimx();
    img_size[1] = result.front().dimy();
    img_size[2] = result.front().dimz();
    img_size[3] = result.front().numChannels();
    // downscaling
    img_size[0] = int(img_size[0] / factor);
    img_size[1] = int(img_size[1] / factor);
    if (downscale_z) {
	img_size[2] = int(img_size[2] / factor);
    }
    while (img_size[0] >= min_num_pixels && 
	   img_size[1] >= min_num_pixels && 
	   (!downscale_z || img_size[2] >= min_num_pixels)) {
	result.push_back(ImgType(img_size[0], img_size[1], img_size[2], img_size[3]));
	resample_into(result.front(), result.back(), linear);
	img_size[0] = int(img_size[0] / factor);
	img_size[1] = int(img_size[1] / factor);
	if (downscale_z) {
	    img_size[2] = int(img_size[2] / factor);
	}
    }
}

//===========================================================================
template<typename ImgType>
void resample_into(const ImgType& input, ImgType& target, bool linear)
//===========================================================================
{
    assert(target.numChannels() == input.numChannels());
    double input_len[3];
    double target_len_inv[3];
    input_len[0] = input.dimx() - 1;
    input_len[1] = input.dimy() - 1;
    input_len[2] = input.dimz() - 1;
    target_len_inv[0] = (target.dimx() > 1) ? double(1) / (target.dimx() - 1) : 0;
    target_len_inv[1] = (target.dimy() > 1) ? double(1) / (target.dimy() - 1) : 0;
    target_len_inv[2] = (target.dimz() > 1) ? double(1) / (target.dimz() - 1) : 0;

    typename ImgType::value_type *dp = target.begin();

    for (int c = 0; c != target.numChannels(); ++c) {
	for (int z = 0; z != target.dimz(); ++z) {
	    const double z_rel_pos = z * target_len_inv[2];
	    const int z_pix = int(input_len[2] * z_rel_pos);
	    const double z_off = (input_len[2] * z_rel_pos) - z_pix;
	    for (int y = 0; y != target.dimy(); ++y) {
		const double y_rel_pos = y * target_len_inv[1];
		const int y_pix = int(input_len[1] * y_rel_pos);
		const double y_off = (input_len[1] * y_rel_pos) - y_pix;
		for (int x = 0; x != target.dimx(); ++x) {
		    const double x_rel_pos = x * target_len_inv[0]; // in interval [0, 1]
		    const int x_pix = int(input_len[0] * x_rel_pos);
		    const double x_off = (input_len[0] * x_rel_pos) - x_pix;
		    *dp++ = interpolate_value(input, c, x_pix, y_pix, z_pix, x_off, y_off, z_off, linear);
		}
	    }
	}
    }
}

//===========================================================================
template<typename T>
void to_grayscale(Image<T>& img)
//===========================================================================
{
    ALWAYS_ERROR_IF(img.numChannels() != 1 && 
		    img.numChannels() != 3, "Image dimension must be 1 or 3");
    if (img.numChannels() == 1) {
	return;
    }
    // img.dim == 3
    Image<T> res(img.dimx(), img.dimy(), img.dimz(), 1);
    for (int z = 0; z < img.dimz(); ++z) {
	for (int y = 0; y < img.dimy(); ++y) {
	    for (int x = 0; x < img.dimx(); ++x) {
		double r = img(x, y, z, 0);
		double g = img(x, y, z, 1);
		double b = img(x, y, z, 2);
		res(x, y, z, 0) = T(0.3 * r + 0.59 * g + 0.11 * b);
	    }
	}
    }
    img.swap(res);
}

}; // end namespace lsseg

namespace {

//===========================================================================
template<typename T>
T interpolate_value(const lsseg::Image<T>& img, 
		    int ch,
		    int x_px,
		    int y_px,
		    int z_px,
		    double x_off,
		    double y_off,
		    double z_off,
		    bool linear)
//===========================================================================
{
    if (!linear) {
	// piecewise constant
	x_px = (x_off >= 0.5 && x_px < (img.dimx() - 1)) ? x_px + 1 : x_px;
	y_px = (y_off >= 0.5 && y_px < (img.dimy() - 1)) ? y_px + 1 : y_px;
	z_px = (z_off >= 0.5 && z_px < (img.dimz() - 1)) ? z_px + 1 : z_px;
	return img(x_px, y_px, z_px, ch);
    }
    // we will apply linear interpolation
    const int x_top = (x_px < (img.dimx() - 1)) ? x_px + 1 : x_px;
    const int y_top = (y_px < (img.dimy() - 1)) ? y_px + 1 : y_px;
    const int z_top = (z_px < (img.dimz() - 1)) ? z_px + 1 : z_px;

    // m/M ~ min/Max in y and z
    const T x_mm = x_off * img(x_top, y_px, z_px, ch) + (1 - x_off) * img(x_px, y_px, z_px, ch) ; 
    const T x_mM = x_off * img(x_top, y_px, z_top, ch) + (1 - x_off) * img(x_px, y_px, z_top, ch) ; 
    const T x_Mm = x_off * img(x_top, y_top, z_px, ch) + (1 - x_off) * img(x_px, y_top, z_px, ch) ; 
    const T x_MM = x_off * img(x_top, y_top, z_top, ch) + (1 - x_off) * img(x_px, y_top, z_top, ch) ; 

    const T y_m = y_off * x_Mm + (1 - y_off) * x_mm; // m ~ min, M ~ max in z coordinate
    const T y_M = y_off * x_MM + (1 - y_off) * x_mM; 
    
    const T z_res = z_off * y_M + (1 - z_off) * y_m;

    return z_res;
}

};

#endif // _SIMPLE_TOOLS_TEMPLATES_H

