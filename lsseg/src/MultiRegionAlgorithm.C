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
// File: MultiRegionAlgorithm.C                                              
//                                                                           
// Created: Fri Mar  3 11:13:18 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: MultiRegionAlgorithm.C,v 1.15 2006/11/25 20:08:26 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements MultiRegionAlgorithm.h.
//                                                                           
//===========================================================================

#include "BorderCondition.h"
#include <assert.h>
#include <limits>
#include <time.h>
#include "level_set.h"
#include "Mask.h"
#include "simple_tools.h"
#include "Region.h"
#include "cimg_dependent.h" // for debug purposes

namespace { // anonymous namespace

void region_compete(const std::vector<lsseg::LevelSetFunction>& cterms,
		    std::vector<lsseg::LevelSetFunction>& fterms, // will be updated
		    const std::vector<lsseg::Mask>& masks,
		    const lsseg::Region* regs);

void resolve_competing(double* vals, int num_vals);

inline double smooth_dirac(double val)
{
    const double FAC = 3;
    val *= FAC;
    return double(1) / (1 + val * val);
}

inline bool regional_overlap(size_t k, const lsseg::Region* regs, int num_regs)
{
    bool one_region_found = false;
    for (int r = 0; r < num_regs; ++r) {
	if (regs[r].phi[k] <= 0) {
	    if (one_region_found) {
		return true;
	    } else {
		one_region_found = true;
	    }
	}
    }
    return false;
}

}; // end anonymous namespace

using namespace std;

namespace lsseg {

//===========================================================================
double develop_multiregion_2D(Region* regs, 
			      int num_regs,
			      int num_iter, 
			      int reinit_modulo,
			      const Mask* geom_mask)
//===========================================================================
{
    // initial assertions
    assert(num_regs > 0 && num_iter >= 0);
    assert(regs[0].phi.dimz() == 1); // function written for the 2D case
    assert(!geom_mask || regs[0].phi.size_compatible(*geom_mask));
    for (int i = 1; i < num_regs; ++i) {
	assert(regs[i].phi.size_compatible(regs[0].phi)); // all regions should be of same size
    }

    // constants that are candidates to be user-defined
    const int MASK_WIDTH = 2;
    const BorderCondition bcond = NEUMANN;
    const int Y = regs[0].phi.dimy();
    const int X = regs[0].phi.dimx();

    std::vector<LevelSetFunction> curv_terms(num_regs);
    std::vector<LevelSetFunction> force_terms(num_regs); 
    std::vector<Mask> region_masks(num_regs);
    std::vector<Mask> old_region_masks(num_regs);
    std::vector<double> max_allowed_timestep(num_regs);
    double time_accum = 0;
    for (int r = 0; r < num_regs; ++r) {
	curv_terms[r].resize(X, Y);
	force_terms[r].resize(X, Y);
	old_region_masks[r].resize(X, Y);
	if (geom_mask) {
	    region_masks[r] = *geom_mask;
	} else {
	    region_masks[r].resize(X, Y);
	    region_masks[r] = 1;
	}
    }
    

    for (int i = 0; i < num_iter; ++i) {
	for (int r = 0; r < num_regs; ++r) {
	    region_masks[r].swap(old_region_masks[r]);
	    make_border_mask(regs[r].phi, region_masks[r], MASK_WIDTH, &old_region_masks[r]);
	    if (geom_mask) {
		region_masks[r].intersect(*geom_mask);
	    }
	    //region_masks[r] = 1; //@@ KRULL
	    //display_image(region_masks[r]);

	    // compute curvature-driven force
	    regs[r].phi.curvature2D(curv_terms[r], &(region_masks[r]));
	    curv_terms[r] *= regs[r].mu;
	    //display_image(curv_terms[r]);

	    // compute individual normal force term
	    regs[r].fgen->update(regs[r].phi);
	    regs[r].fgen->force(force_terms[r], &(region_masks[r]));
	}    

	// let regions compete for influencence
	region_compete(curv_terms, force_terms, region_masks, regs); 

 	//display_image(force_terms[0]);
 	//display_image(force_terms[1]);
// 	display_image(force_terms[2]);
// 	Image<double> diff;
// 	diff = force_terms[0];
// 	diff -= force_terms[1];
// 	display_image(diff);

	for (int r = 0; r < num_regs; ++r) {
	    double H1, H2;
	    //display_image(force_terms[r]);
	    normal_direction_flow_2D(regs[r].phi,
				     force_terms[r],
				     bcond,
				     H1,
				     H2,
				     &(region_masks[r]),
				     geom_mask); 
	    // computing maximum allowed time step
	    //cout << "H1 and H2: " << H1 << " " << H2 << endl;
	    //cout << "Mu: " << regs[r].mu << endl;
	    max_allowed_timestep[r] = double(1) / (4 * regs[r].mu + H1 + H2);
	}

	double dt = *min_element(max_allowed_timestep.begin(), max_allowed_timestep.end());
	time_accum += dt;
	LevelSetFunction grad;
	for (int r = 0; r < num_regs; ++r) {
	    regs[r].phi.gradientNorm2D(grad, &(region_masks[r]));
	    for (int y = 0; y < Y; ++y) {
 		for (int x = 0; x < X; ++x) {
		    if (region_masks[r](x, y) ) {
			const double gradphi = grad(x, y);
			const double uval = 
			    (gradphi * curv_terms[r](x, y)) + force_terms[r](x, y); 
			regs[r].phi(x, y) += dt * uval;
		    }
		}
	    }
	    if (reinit_modulo && ((i+1) % reinit_modulo == 0)) {
		regs[r].phi.reinitialize2D(&region_masks[r]);
		//regs[r].phi.reinitialize2D();
	    }
	}
	cout << "dt: " << dt << endl;
    }
    return time_accum;
}

//===========================================================================
double develop_multiregion_3D(Region* regs, 
			      int num_regs,
			      int num_iter, 
			      int reinit_modulo,
			      const Mask* geom_mask)
//===========================================================================
{
    // initial assertions
    assert(num_regs > 0 && num_iter >= 0);
    assert(!geom_mask || regs[0].phi.size_compatible(*geom_mask));
    for (int i = 1; i < num_regs; ++i) {
	assert(regs[i].phi.size_compatible(regs[0].phi)); // all regions should be of same size
    }

    // constants that are candidates to be user-defined
    const int MASK_WIDTH = 2;
    const BorderCondition bcond = NEUMANN;
    const int Z = regs[0].phi.dimz();
    const int Y = regs[0].phi.dimy();
    const int X = regs[0].phi.dimx();

    std::vector<LevelSetFunction> curv_terms(num_regs);
    std::vector<LevelSetFunction> force_terms(num_regs); 
    std::vector<Mask> region_masks(num_regs);
    std::vector<Mask> old_region_masks(num_regs);
    std::vector<double> max_allowed_timestep(num_regs);
    double time_accum = 0;
    for (int r = 0; r < num_regs; ++r) {
	curv_terms[r].resize(X, Y, Z);
	force_terms[r].resize(X, Y, Z);
	region_masks[r].resize(X, Y, Z);
	old_region_masks[r].resize(X, Y, Z);
	if (geom_mask) {
	    region_masks[r] = *geom_mask;
	} else {
	    region_masks[r].resize(X, Y, Z);
	    region_masks[r] = 1;
	}

    }

    for (int i = 0; i < num_iter; ++i) {
	for (int r = 0; r < num_regs; ++r) {
	    region_masks[r].swap(old_region_masks[r]);
	    make_border_mask(regs[r].phi, region_masks[r], MASK_WIDTH, &old_region_masks[r]);
	    if (geom_mask) {
		region_masks[r].intersect(*geom_mask);
	    }

	    // compute curvature-driven force
	    regs[r].phi.curvature3D(curv_terms[r], &(region_masks[r]));
	    curv_terms[r] *= regs[r].mu;

	    // compute individual normal force term
	    regs[r].fgen->update(regs[r].phi);
	    regs[r].fgen->force(force_terms[r], &(region_masks[r]));
	}    

	// let regions compete for influencence
	region_compete(curv_terms, force_terms, region_masks, regs); 

	for (int r = 0; r < num_regs; ++r) {
	    double H1, H2, H3;
	    //display_image(force_terms[r]);
	    normal_direction_flow_3D(regs[r].phi,
				     force_terms[r],
				     bcond,
				     H1,
				     H2,
				     H3,
				     &(region_masks[r]),
				     geom_mask); 
	    // computing maximum allowed time step
	    max_allowed_timestep[r] = double(1) / (4 * regs[r].mu + H1 + H2 + H3);
	}

	double dt = *min_element(max_allowed_timestep.begin(), max_allowed_timestep.end());
	time_accum += dt;
	LevelSetFunction grad;
	for (int r = 0; r < num_regs; ++r) {
	    regs[r].phi.gradientNorm3D(grad, &(region_masks[r]));
	    for (int z= 0; z < Z; ++z) {
		for (int y = 0; y < Y; ++y) {
		    for (int x = 0; x < X; ++x) {
			if (region_masks[r](x, y, z) ) {
			    const double gradphi = grad(x, y, z);
			    const double uval = 
				(gradphi * curv_terms[r](x, y, z)) + force_terms[r](x, y, z); 
			    regs[r].phi(x, y, z) += dt * uval;
			}
		    }
		}
	    }
	    if (reinit_modulo && ((i+1) % reinit_modulo == 0)) {
		regs[r].phi.reinitialize3D(&region_masks[r]);
	    }
	}
	cout << "dt: " << dt << endl;
    }
    return time_accum;
}


}; // end namespace lsseg

namespace {

//===========================================================================
void region_compete(const std::vector<lsseg::LevelSetFunction>& cterms,
		    std::vector<lsseg::LevelSetFunction>& fterms, // will be updated
		    const std::vector<lsseg::Mask>& masks,
		    const lsseg::Region* regs)
//===========================================================================
{
    const double MIN_INFTY = -std::numeric_limits<double>::infinity();
    const size_t num_regs = cterms.size();
    assert(fterms.size() == num_regs && masks.size() == num_regs);
    const size_t REG_SIZE = fterms[0].size();
    const double OUT_CT = 1;
    std::vector<double> competing_vals(num_regs);

    for (size_t k = 0; k < REG_SIZE; ++k) {

	int overlaps = 0;
	for (size_t r = 0; r < num_regs; ++r) {
	    if (masks[r][k]) {
		++overlaps;
	    }
	}

	if (overlaps > 1) {
	    for (size_t r = 0; r < num_regs; ++r) {
		competing_vals[r] =
		    (masks[r][k]) ? fterms[r][k] - cterms[r][k] : MIN_INFTY;
	    }
	    resolve_competing(&competing_vals[0], num_regs);
	    
	    for (size_t r = 0; r < num_regs; ++r) {
		if (masks[r][k]) {
		    fterms[r][k] = competing_vals[r] + cterms[r][k]; // @@ ?? is this correct?
		}
	    }
	} else {
	    for (size_t r = 0; r < num_regs; ++r) {
		if (masks[r][k] && !regional_overlap(k, regs, num_regs)) {
		    if (fterms[r][k] <= 0) {
			// empty area will "suck" the function into it.
			fterms[r][k] = OUT_CT; 
		    }
		}
	    }
	}
    }
}

//===========================================================================
void resolve_competing(double* vals, int num_vals)
//===========================================================================
{
    const double MIN_VAL = -std::numeric_limits<double>::infinity();
    assert(num_vals >= 2);

    // finding highest and next highest value, as well as their indexes;
    int ix_h = -1; // index of highest value
    double hsf = MIN_VAL;
    double nhsf = hsf;
    for (int i = 0; i < num_vals; ++i) {
	if (vals[i] > hsf) {
	    nhsf = hsf;
	    hsf = vals[i];
	    ix_h = i;
	} else if (vals[i] > nhsf) {
	    nhsf = vals[i];
	}
    }
    assert(ix_h >= 0);
    assert(nhsf > MIN_VAL);
    assert(hsf > MIN_VAL);
    for (int i = 0; i < num_vals; ++i) {
	if (i != ix_h) {
	    vals[i] -= hsf;
	    assert(vals[i] <= 0);
	} else {
	    // highest value
	    vals[i] -= nhsf;
	    assert(vals[i] >= 0);
	}
    }
}



}; // end anonymous namespace
