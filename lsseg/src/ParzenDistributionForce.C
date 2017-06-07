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
// File: ParzenDistributionForce.C                                          
//                                                                           
// Created: Thu Mar  2 15:10:49 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: ParzenDistributionForce.C,v 1.2 2006/11/13 02:29:28 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements ParzenDistributionForce.h.
//                                                                           
//===========================================================================

// recent changes:
// 2006/03/30:	changed functions update(), fixChannelDistributionTo(), saveChannelDistribution(),
//				init() and get_histogram() so that they can handle histograms for inner and outer 
//				regions separately.
//				removed function dumpHistogram() -> functionality is now covert by saveChannelDistribution()

#include "ParzenDistributionForce.h"
#include "cimg_dependent.h" // for blurring function
#include "errormacros.h"
#include <assert.h>
#include <vector>
#include <fstream> //@@ debug

using namespace std;

namespace lsseg {

//===========================================================================
ParzenDistributionForce::ParzenDistributionForce(const Image<double>* img, 
						   Mask* m,
						   bool multireg_mode)
//===========================================================================
    : multi_region_mode_(multireg_mode), num_force_bins_(256)
{
    init(img, m);
}

//===========================================================================
void ParzenDistributionForce::init(const Image<double>* img, const Mask* m)
//===========================================================================
{
    // verification of contract
    assert(img); // there should be an image to segment
    assert(!m || m->numChannels() == 1); // if mask is present, it shall have 1 channel
    assert(!m || img->spatial_compatible(*m)); // if mask is present, it shall have same shape
                                               // as the image to segment

    img_ = img;
    mask_ = m;

    const int nchan = img->numChannels();
    hist_.resize(nchan);
    is_fixed_.resize(nchan);
    precalc_force_.resize(nchan);
    force_bin_factor_.resize(nchan);
    channel_ranges_.resize(nchan);

    for (int c = 0; c != nchan; ++c) {
	double& rmin = channel_ranges_[c].first;
	double& rmax = channel_ranges_[c].second;
	rmin = *std::min_element(img_->channelBegin(c), img_->channelEnd(c));
	rmax = *std::max_element(img_->channelBegin(c), img_->channelEnd(c)); 
	
	force_bin_factor_[c] = num_force_bins_ / (rmax - rmin);

	hist_[c].first = hist_[c].second = makeDefaultHistogram(c);
	precompute_force(c, hist_[c].first, hist_[c].second, precalc_force_[c]);
	is_fixed_[c].first = is_fixed_[c].second = false;
    }
}

//===========================================================================
void ParzenDistributionForce::saveChannelDistribution(int channel, 
						       std::ostream& os, 
						       bool inside) const
//===========================================================================
{
    assert(channel >= 0 && channel < img_->numChannels());
    if (inside) {
	hist_[channel].first.write(os);
    } else {
	hist_[channel].second.write(os);
    }
}

//===========================================================================
void ParzenDistributionForce::fixDistributionToUniform(int channel, bool inside)
//===========================================================================
{
    const int num_bins = 1; // use only one bin
    Histogram h(channel_ranges_[channel].first, channel_ranges_[channel].second, num_bins);
    h.binValues()[0] = 1;
    if (inside) {
	hist_[channel].first = h;
	is_fixed_[channel].first = true;
    } else {
	hist_[channel].second = h;
	is_fixed_[channel].second = true;
    }
    precompute_force(channel, hist_[channel].first, hist_[channel].second, precalc_force_[channel]);
}

//===========================================================================
void ParzenDistributionForce::fixChannelDistributionTo(int channel, std::istream& is, bool inside)
//===========================================================================
{
    // reading distribution
    Histogram h;
    h.read(is);
    
    if (inside) {
	hist_[channel].first = h;
	is_fixed_[channel].first = true;
    } else {
	hist_[channel].second = h;
	is_fixed_[channel].second = true;
    }

    // recompute forces
    precompute_force(channel, hist_[channel].first, hist_[channel].second, precalc_force_[channel]);
}

//===========================================================================
void ParzenDistributionForce::unfixChannelDistribution(int channel, bool inside)
//===========================================================================
{
    if (inside) {
	is_fixed_[channel].first = false;
	hist_[channel].first = makeDefaultHistogram(channel);
    } else {
	is_fixed_[channel].second = false;
	hist_[channel].second = makeDefaultHistogram(channel);
    }
}

//===========================================================================
Histogram ParzenDistributionForce::makeDefaultHistogram(int channel) const 
//===========================================================================
{
    const int DEFAULT_NUM_BINS = 256;

    return Histogram(channel_ranges_[channel].first, 
		     channel_ranges_[channel].second,
		     DEFAULT_NUM_BINS);
}

////===========================================================================
//void ParzenDistributionForce::dumpHistogram(char* filename, 
//					     int channel, 
//					     bool inside) const
////===========================================================================
//{
//    ofstream os(filename);
//    for (int i = 0; i < num_bins_; ++i) {
//	if (inside) {
//	    os << idist_[channel][i] << " ";
//	} else {
//	    os << odist_[channel][i] << " ";
//	}
//    }
//    os.close();
//}

//===========================================================================
void ParzenDistributionForce::update(const LevelSetFunction& phi) 
//===========================================================================
{
    assert(phi.spatial_compatible(*img_));
    const double sigma_h = 8; // blur factor
    const int num_channels = img_->numChannels();

    for (int ch = 0; ch < num_channels; ++ch) {
	if (is_fixed_[ch].first && is_fixed_[ch].second) {
	    // nothing to do for this channel
	    continue;
	}
	Histogram dummy; // even though either of 'idist' and 'odist' below
	                 // might be assigned to this vector, we know that they
                         // will not both be assigned to it at the same time,
	                 // due to the test above.
	Histogram& ihist = is_fixed_[ch].first ? dummy : hist_[ch].first;
	Histogram& ohist = is_fixed_[ch].second ? dummy : hist_[ch].second;

	get_histogram(ch, phi, ihist, ohist);
	    
	// blurring histogram (a little hacked)
	ihist.blur(sigma_h);
	ohist.blur(sigma_h);

	// precalculating the induced force for this value in this channel
	precompute_force(ch, hist_[ch].first, hist_[ch].second, precalc_force_[ch]);
    }
}

//===========================================================================
void ParzenDistributionForce::get_histogram(int channel, 
					     const LevelSetFunction& phi,
					     Histogram& ihist,
					     Histogram& ohist) const
//===========================================================================
{
    ihist.init(channel_ranges_[channel].first, channel_ranges_[channel].second, num_force_bins_);
    ohist.init(channel_ranges_[channel].first, channel_ranges_[channel].second, num_force_bins_);

    vector<double>& idist = ihist.binValues();
    vector<double>& odist = ohist.binValues();
    std::fill(idist.begin(), idist.end(), 1); // by starting with one 'free pixel' in each bin,
    std::fill(odist.begin(), odist.end(), 1); // we assure that no bin will end up with a zero value,
    double area_inside = 0;                   // which would give rise to possibly infinite forces later.
    double area_outside = 0;
    
    const double* data = img_->channelBegin(channel);
    const double* end = phi.end();

    // computing distributions
    if (mask_) {
	const char* m_ptr = mask_ ? mask_->begin() : 0;
	for (const double* it_phi = phi.begin(); it_phi != end; ++it_phi, ++data) {
	    if (*m_ptr++) {
		const int bin_ix = getForceBin(*data, channel);
		assert(bin_ix >= 0 && bin_ix < num_force_bins_);
		if (*it_phi <= 0) {
		    idist[bin_ix] += 1;
		    area_inside += 1;
		} else {
		    odist[bin_ix] += 1;
		    area_outside += 1;
		}
	    }
	}
    } else {
	for (const double* it_phi = phi.begin(); it_phi != end; ++it_phi, ++data) {
	    const int bin_ix = getForceBin(*data, channel);
	    assert(bin_ix >= 0 && bin_ix < num_force_bins_);
	    if (*it_phi <= 0) {
		idist[bin_ix] += 1;
		area_inside += 1;
	    } else {
		odist[bin_ix] += 1;
		area_outside += 1;
	    }
	}
    }
    // normalizing
    double area_inside_inv = area_inside != 0 ? double(1) / area_inside : 1;
    double area_outside_inv = area_outside != 0 ? double(1) / area_outside : 1;
    for (int i = 0; i < num_force_bins_; ++i) {
	odist[i] *= area_outside_inv;
	idist[i] *= area_inside_inv;
    }
}

//===========================================================================
int ParzenDistributionForce::getForceBin(double val, int channel) const
//===========================================================================
{
    int res = int((val - channel_ranges_[channel].first) * force_bin_factor_[channel]);
    if (res >= num_force_bins_) {
	--res;
    }
    return res;
}

//===========================================================================
double ParzenDistributionForce::force2D(int x, int y) const
//===========================================================================
{
    return force(img_->indexOf(x, y));    
}

//===========================================================================
double ParzenDistributionForce::force3D(int x, int y, int z) const 
//===========================================================================
{
    return force(img_->indexOf(x, y, z));
}

//===========================================================================
double ParzenDistributionForce::force(size_t ix) const
//===========================================================================
{
    // force only defined inside mask
    assert(!mask_ || (*mask_)[ix]);

    const size_t chan_size = img_->channelSize();
    const int num_channels = img_->numChannels();

    double res = 0;

    for (int c = 0; c < num_channels; ++c) {
	
	int bin_ix = getForceBin((*img_)[c * chan_size + ix], c);
	res += precalc_force_[c][bin_ix];
    }
    return res;
}

//===========================================================================
void ParzenDistributionForce::force(LevelSetFunction& res, const Mask* m) const
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
double ParzenDistributionForce::getCenterBinValue(int channel, int ix)
//===========================================================================
{
    double min = channel_ranges_[channel].first;
    double max = channel_ranges_[channel].second;
    return min + (double(ix) + 0.5) / force_bin_factor_[channel];
}

//===========================================================================
void ParzenDistributionForce::precompute_force(int channel,
					       const Histogram& ihist,
					       const Histogram& ohist,
					       std::vector<double>& forcevec)
//===========================================================================
{
    forcevec.resize(num_force_bins_);
    
    // precalculating the induced force for this value
    std::fill(forcevec.begin(), forcevec.end(), 0);
    const double min_possible = double(1) / double(img_->size());
    if (multi_region_mode_) {
	for (int i = 0; i < num_force_bins_; ++i) {
	    double val = ihist.valueFor(getCenterBinValue(channel, i));
	    val = (val > min_possible) ? val : min_possible;
	    forcevec[i] = log(val);
	}
    } else {
	for (int i = 0; i < num_force_bins_; ++i) {

	    double val = getCenterBinValue(channel, i);
	    double ival = ihist.valueFor(val);
	    double oval = ohist.valueFor(val);
	    ival = (ival > min_possible) ? ival : min_possible;
	    oval = (oval > min_possible) ? oval : min_possible;
	    forcevec[i] = log(ival / oval);
	}
    }
}


}; // end namespace lsseg
