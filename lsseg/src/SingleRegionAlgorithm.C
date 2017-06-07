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
// File: SingleRegionAlgorithm.C                                             
//                                                                           
// Created: Tue Feb 21 13:14:07 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: SingleRegionAlgorithm.C,v 1.9 2006/11/13 02:29:28 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements SingleRegionAlgorithm.h.
//                                                                           
//===========================================================================

#include "SingleRegionAlgorithm.h"

#include "BorderCondition.h"
#include <assert.h>
#include <limits>
#include "Mask.h"
#include "simple_tools.h"
#include "level_set.h"
#include "cimg_dependent.h" // for debug purposes

using namespace std;

namespace lsseg {

//===========================================================================
double develop_single_region_2D(Region& reg, int num_iter, int reinit_modulo, const Mask* geom_mask)
//===========================================================================
{
    // constants that are candidates to be user-defined
    const int MASK_WIDTH = 2;
    const BorderCondition bcond = NEUMANN;

    assert(int(reg.phi.dimz()) == 1); // this function is written for the 2D case
    assert(!geom_mask || geom_mask->numChannels() == 1);
    assert(!geom_mask || reg.phi.size_compatible(*geom_mask));

    LevelSetFunction curv_term(reg.phi, false);
    LevelSetFunction norm_term(reg.phi, false);
    Mask old_mask, mask(reg.phi, false);
    if (geom_mask) {
	mask = *geom_mask;
    } else {
	mask = 1;
    }
    const double mu = reg.mu;
    double time_accum = 0;
    for (int i = 0; i < num_iter; ++i) {
	
	// compute border mask
	mask.swap(old_mask);
	make_border_mask(reg.phi, mask, MASK_WIDTH, &old_mask);
	if (geom_mask) {
	    mask.intersect(*geom_mask);
	}

	// compute curvature-driven force
	reg.phi.curvatureTimesGrad2D(curv_term, &mask);

	// compute normal force
	reg.fgen->update(reg.phi);
	reg.fgen->force(norm_term, &mask);
	
	double H1, H2;
	//display_image(norm_term);
	normal_direction_flow_2D(reg.phi, 
				 norm_term, 
				 bcond, 
				 H1, 
				 H2, 
				 &mask, 
				 geom_mask); // change to 2D

	// computing maximum allowed time step
	const double dt = double(1) / (4 * reg.mu + H1 + H2);
	
	const double* nt = norm_term.begin();
	const double* ct = curv_term.begin();
	const char* msk = mask.begin();
	const double* end = reg.phi.end();
	for (double* cur = reg.phi.begin(); cur != end; ++cur, ++nt, ++msk, ++ct) {
	    if (*msk) {
		*cur += dt * (mu * (*ct) + (*nt));
	    }
	}
	time_accum += dt;
	if (reinit_modulo && ((i+1) % reinit_modulo == 0)) {
	    reg.phi.reinitialize2D(&mask);
	}
    }
    return time_accum;
}

//===========================================================================
double develop_single_region_3D(Region& reg, int num_iter, int reinit_modulo, const Mask* geom_mask)
//===========================================================================
{
    // constants that are candidates to be user-defined
    const int MASK_WIDTH = 2;
    const BorderCondition bcond = NEUMANN;

    assert(!geom_mask || geom_mask->numChannels() == 1);
    assert(!geom_mask || reg.phi.size_compatible(*geom_mask));

    LevelSetFunction curv_term(reg.phi, false);
    LevelSetFunction norm_term(reg.phi, false);
    Mask old_mask, mask(reg.phi, false);
    if (geom_mask) {
	mask = *geom_mask;
    } else {
	mask = 1;
    }

    const double mu = reg.mu;
    double time_accum = 0;
    for (int i = 0; i < num_iter; ++i) {
	
	// compute border mask
	//make_border_mask(reg.phi, mask, MASK_WIDTH, geom_mask);
	mask.swap(old_mask);
	make_border_mask(reg.phi, mask, MASK_WIDTH, &old_mask);
	if (geom_mask) {
	    mask.intersect(*geom_mask);
	}

	// compute curvature-driven force
	reg.phi.curvatureTimesGrad3D(curv_term, &mask);

	// compute normal force
	reg.fgen->update(reg.phi);
	reg.fgen->force(norm_term, &mask);
	
	double H1, H2, H3;
	//display_image(norm_term);
	normal_direction_flow_3D(reg.phi, 
				 norm_term, 
				 bcond, 
				 H1, 
				 H2, 
				 H3,
				 &mask, 
				 geom_mask); // change to 3D

	// computing maximum allowed time step
	const double dt = double(1) / (4 * reg.mu + H1 + H2 + H3);
	
	const double* nt = norm_term.begin();
	const double* ct = curv_term.begin();
	const char* msk = mask.begin();
	const double* end = reg.phi.end();
	for (double* cur = reg.phi.begin(); cur != end; ++cur, ++nt, ++msk, ++ct) {
	    if (*msk) {
		*cur += dt * (mu * (*ct) + (*nt));
	    }
	}
	time_accum += dt;
	if (reinit_modulo && ((i+1) % reinit_modulo == 0)) {
	    reg.phi.reinitialize3D(&mask);
	}
    }
    return time_accum;
}



}; // end namespace lsseg


