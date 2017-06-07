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
// File: FG_ChanVese.h                                                       
//                                                                           
// Created: Wed Jan 11 09:51:32 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: ChanVeseForce.h,v 1.6 2006/09/20 22:55:47 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Header containing the definition of the
///  \ref lsseg::ForceGenerator "ForceGenerator" called \ref lsseg::ChanVeseForce "ChanVeseForce".
//                                                                           
//===========================================================================

#ifndef _FG_CHANVESE_H
#define _FG_CHANVESE_H

#include "errormacros.h"
#include "ForceGenerator.h"

namespace lsseg{


//===========================================================================
/// \brief This ForceGenerator creates a \ref anchor_NormalForceField "normal force field",
/// that leads to a segmentation similar to the model proposed by Chan and Vese.
///
/// This model seeks to minimize the 
/// \ref anchor_MumfordShahFunctional "Mumford-Shah functional" by splitting the image into
/// two regions where the pixel distribution is modeled by \em constant \em values for each
/// region.  It can be seen as the simplest model (in our setting, at least) that use a
/// statistical model to describe each region.  (The statistical model being, in our case, the
/// average pixel value in each region).  The relevant papers are \ref anchor_Chan99 "[Chan99]" 
/// and \ref anchor_Chan01 "[Chan01]".
/// Thomas Brox also gives a nice overview of this and other methods in 
/// \ref anchor_Brox05 "[Brox05]".
class ChanVeseForce : public ForceGenerator
//===========================================================================
{
public:
    /// Constructor, taking a pointer to the image to be segmented.
     ChanVeseForce(const Image<double>* img) : img_(img) 
    {
	assert(img_->numChannels() == 1); // @@ temporary restriction...
    }

    virtual ~ChanVeseForce() { }

    virtual void init(const Image<double>* img, const Mask* mask) 
    {
	MESSAGE("Warning!  Mask ignored in ChanVeseForce::init()");
	img_ = img;
    }

    virtual void update(const LevelSetFunction& phi) 
    {
	assert(phi.spatial_compatible(*img_));
	const int SIZE = phi.size();
	unsigned int nb_in = 0;
	unsigned int nb_out = 0;
	mu_inside_ = mu_outside_ = 0;
	for (int i = 0; i < SIZE; ++i) {
	    if (phi[i] > 0) {
		++nb_out;
		mu_outside_ += (*img_)[i];
	    } else {
		++nb_in;
		mu_inside_ += (*img_)[i];
	    }
	}
	if (nb_in) {
	    mu_inside_ /= nb_in;
	}
	if (nb_out) {
	    mu_outside_ /= nb_out;
	}
    }

    virtual void force(LevelSetFunction& phi, const Mask* mask = 0) const
    {
	assert(!m || img_->spatial_compatible(*m));
	assert(img_->spatial_compatible(phi));
	const size_t SIZE = phi.size();
	
	// explicit checking of mask presence for optimization reasons
	if (m) {
	    if (mask_) {
		for (size_t it = 0; it < SIZE; ++it) {
		    phi[it] = ((*m)[it] && (*mask_)[it]) ? force(it) : 0;
		}
	    } else {
		for (size_t it = 0; it < SIZE; ++it) {
		    phi[it] = (*m)[it] ? force(it) : 0;
		}
	    }
	} else {
	    if (mask_) {
		for (size_t it = 0; it < SIZE; ++it) {
		    phi[it] = (*mask_)[it] ? force(it) : 0;
		}
	    } else {
		for (size_t it = 0; it < SIZE; ++it) {
		    phi[it] = force(it);
		}
	    }
	}
    }
    

    virtual double force3D(int x, int y, int z) const 
    {
	int ix = img_->indexOf(x, y, z);
	return force(ix);
    }

    virtual double force2D(int x, int y) const
    {
	int ix = img_->indexOf(x,y);
	return force(ix);
    }

    virtual double force(size_t ix) const 
    {
	double tmp = (*img_)[ix];
	double d_out = mu_outside_ - tmp;
	double d_in = mu_inside_ - tmp;
	return d_out * d_out - d_in * d_in;
    }

private:
    /// \brief average pixel value inside the Image's closed region, for the last segmentation
    double mu_inside_;
    /// \brief average pixel value inside the Image's open region, for the last segmentation
    double mu_outside_;
    /// \brief pointer to the image to be segmented; specified in the constructor
    const Image<double>* const img_;
};

}; // end namespace lsseg

#endif // _FG_CHANVESE_H

