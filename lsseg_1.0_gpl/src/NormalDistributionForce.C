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
// File: NormalDistributionForce.C                                           
//                                                                           
// Created: Fri Mar  3 10:56:16 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: NormalDistributionForce.C,v 1.4 2006/11/25 20:08:27 oan Exp $
//                                                                           
// Description:
/// \file 
/// \brief Implements NormalDistributionForce.h.
//                                                                           
//===========================================================================

#include "NormalDistributionForce.h"

namespace lsseg {

//===========================================================================
NormalDistributionForce::NormalDistributionForce(const Image<double>* img,
						 Mask* m,
						 bool multireg_mode)
//===========================================================================
    : multi_region_mode_(multireg_mode)
{
    init(img, m);
}

//===========================================================================
void NormalDistributionForce::init(const Image<double>* img, const Mask* mask)
//===========================================================================
{
    // verification of contract
    assert(img);
    assert(!mask || mask->numChannels() == 1);
    assert(!mask || img->spatial_compatible(*mask));

    img_ = img;
    mask_ = mask;
    const int nchan = img->numChannels();
    mu_in_.resize(nchan);
    mu_out_.resize(nchan);
    inv_sigma_in_.resize(nchan);
    inv_sigma_out_.resize(nchan);
    precalc_log_.resize(nchan);
}

//===========================================================================
void NormalDistributionForce::update(const LevelSetFunction& phi) 
//===========================================================================
{
    // verification of contract
    assert(phi.spatial_compatible(*img_));
    
    compute_averages(phi);
    compute_deviations(phi);
}

//===========================================================================
double NormalDistributionForce::force2D(int x, int y) const
//===========================================================================
{
    return force(img_->indexOf(x, y));
}

//===========================================================================
double NormalDistributionForce::force3D(int x, int y, int z) const
//===========================================================================
{
    return force(img_->indexOf(x, y, z));
}

//===========================================================================
double NormalDistributionForce::force(size_t ix) const
//===========================================================================
{
    // force only defined inside mask
    assert(!mask_ || (*mask_)[ix]); 

    double res = 0;
    const int CNUM = img_->numChannels();
    const int CSIZE = img_->channelSize();

    if (multi_region_mode_) {
	for (int c = 0; c < CNUM; ++c) {
	    const size_t c_offset = c * CSIZE;
	    const double s1_inv = inv_sigma_in_[c];
	    const double img_val = (*img_)[ix + c_offset];
	    const double tmp = (img_val - mu_in_[c]) * s1_inv;
	    res += precalc_log_[c] - 0.5 * tmp * tmp;
	}
    } else {
	for (int c = 0; c < CNUM; ++c) {
	    const size_t c_offset = c * CSIZE;
	    const double s1_inv = inv_sigma_in_[c];
	    const double s2_inv = inv_sigma_out_[c];
	    const double img_val = (*img_)[ix + c_offset];
	    const double tmp1 = (img_val - mu_in_[c]) * s1_inv;
	    const double tmp2 = (img_val - mu_out_[c]) * s2_inv;
	    res += precalc_log_[c] + 0.5 * (tmp2 * tmp2 - tmp1 * tmp1);
	}
    }
    return res;
}    

//===========================================================================
void NormalDistributionForce::force(LevelSetFunction& res, const Mask* m) const
//===========================================================================
{
    assert(!m || img_->spatial_compatible(*m));
    assert(img_->spatial_compatible(res));
    const size_t SIZE = res.size();

    // explicit checking of mask presence for optimization reasons
    if (m) {
	if (mask_) {
	    for (size_t it = 0; it < SIZE; ++it) {
		res[it] = ((*m)[it] && (*mask_)[it]) ? force(it) : 0;
	    }
	} else {
	    for (size_t it = 0; it < SIZE; ++it) {
		res[it] = (*m)[it] ? force(it) : 0;
	    }
	}
    } else {
	if (mask_) {
	    for (size_t it = 0; it < SIZE; ++it) {
		res[it] = (*mask_)[it] ? force(it) : 0;
	    }
	} else {
	    for (size_t it = 0; it < SIZE; ++it) {
		res[it] = force(it);
	    }
	}
    }
}

//===========================================================================
void NormalDistributionForce::compute_averages(const LevelSetFunction& phi)
//===========================================================================
{
    assert(img_->spatial_compatible(phi));
    const size_t CSIZE = img_->channelSize();
    
    for (int c = 0; c < img_->numChannels(); ++c) {
	size_t num_in = 0; // count pixels inside region
	size_t num_out = 0; // count pixels outside region
	
	// compute statistics for each channel
	mu_in_[c] = mu_out_[c] = 0;
	const size_t c_offset = c * CSIZE;
	
	if (mask_) {
	    for (size_t i = 0; i < CSIZE; ++i) {
		if ((*mask_)[i]) {
		    const double img_val = (*img_)[i + c_offset];
		    if (phi[i] <= 0) {
			++num_in;
			mu_in_[c] += img_val;
		    } else {
			++num_out;
			mu_out_[c] += img_val;
		    }
		}
	    }
	} else {
	    // no mask to worry about
	    for (size_t i = 0; i < CSIZE; ++i) {
		const double img_val = (*img_)[i + c_offset];
		if (phi[i] <= 0) {
		    ++num_in;
		    mu_in_[c] +=img_val;
		} else {
		    ++num_out;
		    mu_out_[c] += img_val;
		}
	    }
	}
	if (num_in > 0) {
	    mu_in_[c] /= double(num_in);
	}
	if (num_out > 0) {
	    mu_out_[c] /= double(num_out);
	}
    }
}

//===========================================================================
void NormalDistributionForce::compute_deviations(const LevelSetFunction& phi)
//===========================================================================
{
    // this function supposes that the averages (mu_in_ and mu_out_)
    // have already been calculated for all channels
    const double PI = 3.1415926535897932384;
    const double log_inv_root_2pi = log(double(1)/sqrt(2 * PI));

    assert(img_->spatial_compatible(phi));
    const double EPS = 1.0e-2;
    const size_t CSIZE = img_->channelSize();

    for (int c = 0; c < img_->numChannels(); ++c) {
	size_t num_in = 0;
	size_t num_out = 0;

	const size_t c_offset = c * CSIZE;

	if (mask_) {
	    for (size_t i = 0; i < CSIZE; ++i) {
		if((*mask_)[i]) {
		    const double img_val = (*img_)[i + c_offset];
		    if (phi[i] <= double(0)) {
			++num_in;
			const double diff = mu_in_[c] - img_val;
			inv_sigma_in_[c] += diff * diff;
		    } else {
			++num_out;
			const double diff = mu_out_[c] - img_val;
			inv_sigma_out_[c] += diff * diff;
		    }
		}
	    }
	} else {
	    // no mask to worry about
	    for (size_t i = 0; i < CSIZE; ++i) {
		const double img_val = (*img_)[i + c_offset];
		if (phi[i] <= double(0)) {
		    ++num_in;
		    const double diff = mu_in_[c] - img_val;
		    inv_sigma_in_[c] += diff * diff;
		} else {
		    ++num_out;
		    const double diff = mu_out_[c] - img_val;
		    inv_sigma_out_[c] += diff * diff;
		}
	    }
	}
	if (num_in > 0) {
	    inv_sigma_in_[c] /= double(num_in);
	    inv_sigma_in_[c] = sqrt(inv_sigma_in_[c]);
	    inv_sigma_in_[c] = double(1) / (inv_sigma_in_[c] > EPS ? inv_sigma_in_[c] : EPS);
	} else {
	    inv_sigma_in_[c] = double(1) / EPS;
	}
	if (num_out > 0) {
	    inv_sigma_out_[c] /= double(num_out);
	    inv_sigma_out_[c] = sqrt(inv_sigma_out_[c]);
	    inv_sigma_out_[c] = double(1) / (inv_sigma_out_[c] > EPS ? inv_sigma_out_[c] : EPS);
	} else {
	    inv_sigma_out_[c] = double(1) / EPS;
	}
	precalc_log_[c] = multi_region_mode_ ?
	    log(inv_sigma_in_[c]) + log_inv_root_2pi :
	    log(inv_sigma_in_[c] / inv_sigma_out_[c]);
    }
}


}; // end namespace lsseg
