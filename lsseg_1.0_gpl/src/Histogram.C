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
// File: Histogram.C                                                         
//                                                                           
// Created: Fri Mar 31 17:54:22 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: Histogram.C,v 1.5 2006/11/13 02:29:28 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements Histogram.h.
//                                                                           
//===========================================================================

#include "Histogram.h"
#include "cimg_dependent.h"

using namespace std;

namespace lsseg {

void Histogram::write(std::ostream& os) const
{
    // writing number of bins
    os << numBins() << endl;

    // writing range
    os << rangeMin() << " " << rangeMax() << endl;

    // writing distribution
    for (int i = 0; i < numBins(); ++i) {
	os << dist_[i] << " ";
    }
    os << endl;
}

void Histogram::read(std::istream& is)
{
    // reading number of bins
    int num_bins ;
    is >> num_bins;
    dist_.resize(num_bins);

    // reading range
    is >> range_min_;
    is >> range_max_;

    // reading distribution
    for (int i = 0; i < numBins(); ++i) {
	is >> dist_[i];
    }

    bin_factor_ = double(numBins()) / (rangeMax() - rangeMin());
}

int Histogram::getBin(double val) const
{
    int res = int((val - range_min_) * bin_factor_);
    if (res == numBins()) {
	--res;
    }
    return res;
}

void Histogram::blur(double sigma)
{
    if (numBins() > 1) {
	blur_1D(&dist_[0], numBins(), sigma);
    }
}

double Histogram::valueFor(double val) const 
{
    int ix = getBin(val);
    if (ix < 0 || ix >= numBins()) {
	return 0; // outside histogram range
    }
    return dist_[ix] * bin_factor_;
}

}; // end namespace lsseg
