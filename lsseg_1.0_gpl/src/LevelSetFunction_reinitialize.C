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
// File: LevelSetFunction_reinitialize.C                                     
//                                                                           
// Created: Mon Feb 20 11:50:38 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: LevelSetFunction_reinitialize.C,v 1.9 2006/11/22 15:23:00 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements the \ref lsseg::LevelSetFunction::reinitialize2D "reinitialize2D()" and 
///        \ref lsseg::LevelSetFunction::reinitialize3D "reinitialize3D()" functions of the
///        \ref lsseg::LevelSetFunction "LevelSetFunction class, as defined in LevelSetFunction.h.
//                                                                           
//===========================================================================

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <stdexcept>
#include "Image.h"
#include "LevelSetFunction.h"
#include "Mask.h"
#include "cimg_dependent.h" // @@ for debug purposes
#include <iostream> // @@ for debug purposes

// UNDEF_VAL should be LARGER than any number that we can expect from the
// computations.  It should under NO CIRCUMSTANCES be negative.

const double UNDEF_VAL = 100000; 

using namespace std;
using namespace lsseg;
using namespace boost;

namespace { // anonymous namespace


//typedef tuple<int, int> Coord2D;
//typedef tuple<int, int, int> Coord3D;
//===========================================================================
struct Coord2D
//===========================================================================
{
    Coord2D(int x, int y) : x_(x), y_(y) {}
    template<int I> const int& get() const ;

    int x_, y_;
};
inline bool operator<(const Coord2D& c1, const Coord2D& c2) { 
    if (c1.y_ != c2.y_) 
	return (c1.y_ < c2.y_); 
    return (c1.x_ < c2.x_); 
}
inline bool operator==(const Coord2D& c1, const Coord2D& c2) {
    return (c1.x_ == c2.x_) && (c1.y_ == c2.y_);
}
template<> const int& Coord2D::get<0>() const {return x_;}
template<> const int& Coord2D::get<1>() const {return y_;}

template<size_t MUL>
class HashCoord2D
{
public:
    size_t operator()(const Coord2D& c) const {return c.y_ * MUL + c.x_;}
};

//===========================================================================
struct Coord3D
//===========================================================================
{
    Coord3D(int x, int y, int z) : x_(x), y_(y), z_(z) {}
    template<int I> const int& get() const ;

    int x_, y_, z_;
};
inline bool operator<(const Coord3D& c1, const Coord3D& c2) {
    if (c1.z_ != c2.z_)	return (c1.z_ < c2.z_); 
    if (c1.y_ != c2.y_) return (c1.y_ < c2.y_);
    return (c1.x_ < c2.x_);
}

inline bool operator==(const Coord3D& c1, const Coord3D& c2) {
    return (c1.x_ == c2.x_) && (c1.y_ == c2.y_) && (c1.z_ == c2.z_);
}

template<> const int& Coord3D::get<0>() const {return x_;}
template<> const int& Coord3D::get<1>() const {return y_;}
template<> const int& Coord3D::get<2>() const {return z_;}

template<size_t MUL>
class HashCoord3D
{
public:
    size_t operator()(const Coord3D& c) const {return c.x_ + MUL * (c.y_ + MUL * c.z_);}
};

const size_t MAX_ANTICIPATED_RES=10000;
typedef multimap<double, Coord2D > dist_map_2D; // sparse distance map in 2D
typedef multimap<double, Coord3D > dist_map_3D; // sparse distance map in 3D
typedef map<Coord2D, dist_map_2D::iterator> coord_map_2D; // inverse of dist_map_2D
typedef map<Coord3D, dist_map_3D::iterator> coord_map_3D; // inverse of dist_map_3D


void init_zero_level_interface_2D(const LevelSetFunction&, LevelSetFunction&, const Mask* m);
void init_zero_level_interface_3D(const LevelSetFunction&, LevelSetFunction&, const Mask* m);
    
void init_tentative_values_2D(const LevelSetFunction&, dist_map_2D&, dist_map_2D&, const Mask* m);
void init_tentative_values_3D(const LevelSetFunction&, dist_map_3D&, dist_map_3D&, const Mask* m);

double compute_dist_estimate_2D(const double x1, const double x2, 
				const double y1, const double y2);
double compute_dist_estimate_3D(const double x1, const double x2, 
				const double y1, const double y2, 
				const double z1, const double z2);

void march_out_2D(dist_map_2D& tentatives, 
		  coord_map_2D& inv_tent, 
		  LevelSetFunction& target,
		  const Mask* m,
		  bool negative);
void march_out_3D(dist_map_3D& tentatives, 
		  coord_map_3D& inv_tent, 
		  LevelSetFunction& target,
		  const Mask* m,
		  bool negative);
inline double distance_poly(double d1, double d2);
inline double distance_poly(double d1, double d2, double d3);
}; // end anonymous namespace

namespace lsseg {

//===========================================================================
void LevelSetFunction::reinitialize2D(const Mask* m)
//===========================================================================
{
    // for the moment, this routine uses a first-order accurate fast marching
    // method.  (A different approach could be the reinitialization equation
    // using a higher-order HJ WENO scheme). 

    MESSAGE_IF(dimz() > 1, "Warning: applying 2D method on a grid with dim(z) > 1");
    assert(!m || size_compatible(*m));

    // detecting zero-level interface and computing neighbour distance values
    // (1-band)

    static LevelSetFunction res; // made static because memory allocation is expensive,
                                 // and this function is typically called (very) frequently!
    res.resize(*this);
    res = UNDEF_VAL;

    init_zero_level_interface_2D(*this, res, m);

    // initializing tentative values (those that have a defined neighbour while themselves 
    // being undefined)
    dist_map_2D pos_tentatives, neg_tentatives;
    coord_map_2D pos_tent_inv, neg_tent_inv;
    init_tentative_values_2D(res, pos_tentatives, neg_tentatives, m);
    
    // computing inverse maps  (there might be a more efficient way of carrying
    // out the task solved by these maps, possible optimization @@)
    for (dist_map_2D::iterator it = pos_tentatives.begin(); it != pos_tentatives.end(); ++it) {
	pos_tent_inv[it->second] = it;
    }
    for (dist_map_2D::iterator it = neg_tentatives.begin(); it != neg_tentatives.end(); ++it) {
	neg_tent_inv[it->second] = it;
    }

    // growing outwards
    march_out_2D(pos_tentatives, pos_tent_inv, res, m, false);
    // growing inwards
    march_out_2D(neg_tentatives, neg_tent_inv, res, m, true);

    assert(pos_tentatives.empty());
    assert(neg_tentatives.empty());
    assert(pos_tent_inv.empty());
    assert(neg_tent_inv.empty());
    
    if (!m) { // res has been updated everywhere, and we can just swap it
	swap(res);
    } else {
	// copying the defined values, leaving the rest as it is
	const double* res_pt = res.begin();
	const char* mask_pt = m->begin();
	for (double* phi_pt = begin(); phi_pt != end(); ++phi_pt, ++res_pt, ++mask_pt) {
	    if (*mask_pt) {
		if (*res_pt != UNDEF_VAL) {
		    *phi_pt = *res_pt;
		} else {
		    // this part has not been reached by the marching alrogithm.  It must
		    // be separated from the interface by the mask (the whole area it 
		    // belongs to has the same sign
		    *phi_pt = (*phi_pt > 0) ? 1: -1;
		}
	    }
	}
    }
}

//===========================================================================
void LevelSetFunction::reinitialize3D(const Mask* m)
//===========================================================================
{
    // for the moment, this routine uses a first-order accurate fast marching
    // method.  (A different approach could be the reinitialization equation
    // using a higher-order HJ WENO scheme). 
    assert(!m || size_compatible(*m));

    if (dimz() == 1) {
	return reinitialize2D(m); // not necessary to call the 3D routine
    }

    // detecting zero-level interface and computing neighbour distance values
    // (1-band)

    static LevelSetFunction res; // made static because memory allocation is expensive,
                                 // and this function is typically called (very) frequently!
    res.resize(*this);
    res = UNDEF_VAL;
    init_zero_level_interface_3D(*this, res, m);

    // initializing tentative values (those that have a defined neighbour while themselves 
    // being undefined)
    dist_map_3D pos_tentatives, neg_tentatives;
    coord_map_3D pos_tent_inv, neg_tent_inv;
    init_tentative_values_3D(res, pos_tentatives, neg_tentatives, m);
    
    // computing inverse maps  (there might be a more efficient way of carrying
    // out the task solved by these maps, possible optimization @@)
    for (dist_map_3D::iterator it = pos_tentatives.begin(); it != pos_tentatives.end(); ++it) {
	pos_tent_inv[it->second] = it;
    }
    for (dist_map_3D::iterator it = neg_tentatives.begin(); it != neg_tentatives.end(); ++it) {
	neg_tent_inv[it->second] = it;
    }

    // growing outwards
    march_out_3D(pos_tentatives, pos_tent_inv, res, m, false);
    // growing inwards
    march_out_3D(neg_tentatives, neg_tent_inv, res, m, true);

    assert(pos_tentatives.empty());
    assert(neg_tentatives.empty());
    assert(pos_tent_inv.empty());
    assert(neg_tent_inv.empty());
    
    if (!m) { // res has been updated everywhere, and we can just swap it
	swap(res);
    } else {
	// copying the defined values, leaving the rest as it is
	double* res_pt = res.begin();
	for (double* phi_pt = begin(); phi_pt != end(); ++phi_pt, ++res_pt) {
	    if (*res_pt != UNDEF_VAL) {
		*phi_pt = *res_pt;
	    }
	}
    }
}


}; // end anonymous namespace

namespace { // anonymous namespace

//===========================================================================
void init_zero_level_interface_2D(const lsseg::LevelSetFunction& src, 
				  lsseg::LevelSetFunction& trg,
				  const Mask* m)
//===========================================================================
{
    // @@ NB: Mask not taken into consideration when checking for a sign change.
    //    It seems to work OK though, and it saves some checks, but this might
    //    have to be rewritten.
    const double EPS = 1.0e-10;
    const int X = src.dimx();
    const int Y = src.dimy();
    vector<double> scratch(4);
    trg.resize(src);
    size_t mask_it = 0;
    for (int y = 0; y < Y; ++y) {
	for (int x = 0; x < X; ++x) {
	    if (m && !(*m)[mask_it++]) {
		continue; // outside mask
	    }
	    scratch.clear();
	    const double cur = src(x, y);
	    if (cur == 0.0) { // necessary when a pixel is found to be EXACTLY on the interface
		trg(x, y) = 0;
		continue;
	    }
	    // we know now that cur is not exactly zero and thus MUST have a sign
	    if (x > 0 && cur * src(x - 1, y) <= 0) {
		// estimating position of interface by linear interpolation
		const double neigh = src(x - 1, y);
		const double val = cur / (cur - neigh);
		assert (val >= 0 && val <= 1);
		scratch.push_back((val > EPS) ? val : EPS);
	    }
	    if (x < X - 1 && cur * src(x + 1, y) <= 0) {
		// estimating position of interface by linear interpolation
		const double neigh = src(x + 1, y);
		const double val = cur / (cur - neigh);
		assert (val >= 0 && val <= 1);
		scratch.push_back((val > EPS) ? val : EPS);
	    }
	    if (y > 0 && cur * src(x, y - 1) <= 0) {
		// estimating position of interface by linear interpolation
		const double neigh = src(x, y - 1);
		const double val = cur / (cur - neigh);
		assert (val >= 0 && val <= 1);
		scratch.push_back((val > EPS) ? val : EPS);
	    }
	    if (y < Y - 1 && cur * src(x, y + 1) <= 0 ) {
		// estimating position of interface by linear interpolation
		const double neigh = src(x, y + 1);
		const double val = cur / (cur - neigh);
		assert (val >= 0 && val <= 1);
		scratch.push_back((val > EPS) ? val : EPS);
	    }
	    if (!scratch.empty()) {
		trg(x, y) = *min_element(scratch.begin(), scratch.end());
		if (cur < 0) {
		    trg(x, y) *= -1;
		}
	    }
	}
    }
}

//===========================================================================
void init_zero_level_interface_3D(const lsseg::LevelSetFunction& src, 
				  lsseg::LevelSetFunction& trg,
				  const Mask* m)
//===========================================================================
{
    const double EPS = 1.0e-10;
    const int X = src.dimx();
    const int Y = src.dimy();
    const int Z = src.dimz();
    vector<double> scratch(6);
    trg.resize(src);
    size_t mask_it = 0;
    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		if (m && !(*m)[mask_it++]) {
		    continue; // outside mask
		}
		scratch.clear();
		const double cur = src(x, y, z);
		if (cur == 0.0) { // necessary when a pixel is found to be EXACTLY on the interface
		    trg(x, y, z) = 0;
		    continue;
		}
		// we know now that cur is not exactly zero and thus MUST have a sign
		if (x > 0 && cur * src(x - 1, y, z) <= 0) {
		    // estimating position of interface by linear interpolation
		    const double neigh = src(x - 1, y, z);
		    const double val = cur / (cur - neigh);
		    assert (val >= 0 && val <= 1);
		    scratch.push_back((val > EPS) ? val : EPS);
		}
		if (x < X - 1 && cur * src(x + 1, y, z) <= 0) {
		    // estimating position of interface by linear interpolation
		    const double neigh = src(x + 1, y, z);
		    const double val = cur / (cur - neigh);
		    assert (val >= 0 && val <= 1);
		    scratch.push_back((val > EPS) ? val : EPS);
		}
		if (y > 0 && cur * src(x, y - 1, z) <= 0) {
		    // estimating position of interface by linear interpolation
		    const double neigh = src(x, y - 1, z);
		    const double val = cur / (cur - neigh);
		    assert (val >= 0 && val <= 1);
		    scratch.push_back((val > EPS) ? val : EPS);
		}
		if (y < Y - 1 && cur * src(x, y + 1, z) <= 0 ) {
		    // estimating position of interface by linear interpolation
		    const double neigh = src(x, y + 1, z);
		    const double val = cur / (cur - neigh);
		    assert (val >= 0 && val <= 1);
		    scratch.push_back((val > EPS) ? val : EPS);
		}
		if (z > 0 && cur * src(x, y ,z - 1) <= 0) {
		    // estimating position of interface by linear interpolation
		    const double neigh = src(x, y, z - 1);
		    const double val = cur / (cur - neigh);
		    assert(val >= 0 && val <= 1);
		    scratch.push_back((val > EPS) ? val : EPS);
		}
		if (z < Z - 1 && cur * src(x, y, z + 1) <= 0) {
		    // estimating position of interface by linear interpolation
		    const double neigh = src(x, y, z+1);
		    const double val = cur / (cur - neigh);
		    assert(val >= 0 && val <= 1);
		    scratch.push_back((val > EPS) ? val : EPS);
		}
		if (!scratch.empty()) {
		    trg(x, y, z) = *min_element(scratch.begin(), scratch.end());
		    if (cur < 0) {
			trg(x, y, z) *= -1;
		    }
		}
	    }
	}
    }
}

//===========================================================================
void init_tentative_values_2D(const lsseg::LevelSetFunction& phi, 
			      dist_map_2D& positive_map,
			      dist_map_2D& negative_map, 
			      const Mask* m)
//===========================================================================
{
    positive_map.clear();
    negative_map.clear();
    const int X = phi.dimx();
    const int Y = phi.dimy();

    // mapping all coordinate pair that are undefined but having a defined neighbour
    set<Coord2D > candidates;
    if (!m) { // no mask
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		if (phi(x, y) != UNDEF_VAL) {
		    // checking neighbours
		    if (x > 0 && phi(x-1, y) == UNDEF_VAL) {
			candidates.insert(Coord2D(x-1, y));
		    }
		    if (x < X-1 && phi(x+1, y) == UNDEF_VAL) {
			candidates.insert(Coord2D(x+1, y));
		    }
		    if (y > 0 && phi(x, y-1) == UNDEF_VAL) {
			candidates.insert(Coord2D(x, y-1));
		    }
		    if (y < Y-1 && phi(x, y+1) == UNDEF_VAL) {
			candidates.insert(Coord2D(x, y+1));
		    }
		}
	    }
	}
    } else { // mask is present
	size_t mask_it = 0;
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		if ((*m)[mask_it++]) {
		    if (phi(x, y) != UNDEF_VAL) {
			// checking neighbours
			if (x > 0 && phi(x-1, y) == UNDEF_VAL && (*m)(x-1, y)) {
			    candidates.insert(Coord2D(x-1, y));
			}
			if (x < X-1 && phi(x+1, y) == UNDEF_VAL && (*m)(x+1, y)) {
			    candidates.insert(Coord2D(x+1, y));
			}
			if (y > 0 && phi(x, y-1) == UNDEF_VAL && (*m)(x, y-1)) {
			    candidates.insert(Coord2D(x, y-1));
			}
			if (y < Y-1 && phi(x, y+1) == UNDEF_VAL && (*m)(x, y+1)) {
			    candidates.insert(Coord2D(x, y+1));
			}
		    }
		}
	    }
	}
    }

    // looping through candidates and computing their distance to the interface
    set<Coord2D >::const_iterator it;
    double l, r, b, t; // left, right, bottom, top...
    for (it = candidates.begin(); it != candidates.end(); ++it) {
	const int x = it->get<0>();
	const int y = it->get<1>();
	l = (x > 0)   ? phi(x-1, y) : UNDEF_VAL;
	r = (x < X-1) ? phi(x+1, y) : UNDEF_VAL;
	b = (y > 0)   ? phi(x, y-1) : UNDEF_VAL;
	t = (y < Y-1) ? phi(x, y+1) : UNDEF_VAL;

	double estimate = compute_dist_estimate_2D(l, r, b, t);

	if (estimate > 0) {
	    positive_map.insert(pair<double, Coord2D >(estimate, *it)); 
	} else {
	    negative_map.insert(pair<double, Coord2D >(estimate, *it)); 
	}	
    }
}

//===========================================================================
void init_tentative_values_3D(const lsseg::LevelSetFunction& phi, 
			      dist_map_3D& positive_map,
			      dist_map_3D& negative_map, 
			      const Mask* m)
//===========================================================================
{
    positive_map.clear();
    negative_map.clear();
    const int X = phi.dimx();
    const int Y = phi.dimy();
    const int Z = phi.dimz();

    // mapping all coordinate pair that are undefined but having a defined neighbour
    set<Coord3D > candidates;
    if (!m) { // no mask
	for (int z = 0; z < Z; ++z) {
	    for (int y = 0; y < Y; ++y) {
		for (int x = 0; x < X; ++x) {
		    if (phi(x, y, z) != UNDEF_VAL) {
			// checking neighbours
			if (x > 0 && phi(x-1, y, z) == UNDEF_VAL) {
			    candidates.insert(Coord3D(x-1, y, z));
			}
			if (x < X-1 && phi(x+1, y, z) == UNDEF_VAL) {
			    candidates.insert(Coord3D(x+1, y, z));
			}
			if (y > 0 && phi(x, y-1, z) == UNDEF_VAL) {
			    candidates.insert(Coord3D(x, y-1, z));
			}
			if (y < Y-1 && phi(x, y+1, z) == UNDEF_VAL) {
			    candidates.insert(Coord3D(x, y+1, z));
			}
			if (z > 0 && phi(x, y, z-1) == UNDEF_VAL) {
			    candidates.insert(Coord3D(x, y, z-1));
			}
			if (z < Z-1 && phi(x, y, z+1) == UNDEF_VAL) {
			    candidates.insert(Coord3D(x, y, z+1));
			}
		    }
		}
	    }
	}
    } else { // mask is present
	size_t mask_it = 0;
	for (int z = 0; z < Z; ++z) {
	    for (int y = 0; y < Y; ++y) {
		for (int x = 0; x < X; ++x) {
		    if ((*m)[mask_it++]) {
			if (phi(x, y, z) != UNDEF_VAL) {
			    // checking neighbours
			    if (x > 0 && phi(x-1, y, z) == UNDEF_VAL && (*m)(x-1, y, z)) {
				candidates.insert(Coord3D(x-1, y, z));
			    }
			    if (x < X-1 && phi(x+1, y, z) == UNDEF_VAL && (*m)(x+1, y, z)) {
				candidates.insert(Coord3D(x+1, y, z));
			    }
			    if (y > 0 && phi(x, y-1, z) == UNDEF_VAL && (*m)(x, y-1, z)) {
				candidates.insert(Coord3D(x, y-1, z));
			    }
			    if (y < Y-1 && phi(x, y+1, z) == UNDEF_VAL && (*m)(x, y+1, z)) {
				candidates.insert(Coord3D(x, y+1, z));
			    }
			    if (z > 0 && phi(x, y, z-1) == UNDEF_VAL && (*m)(x, y, z-1)) {
				candidates.insert(Coord3D(x, y, z-1));
			    }
			    if (z < Z-1 && phi(x, y, z+1) == UNDEF_VAL && (*m)(x, y, z+1)) {
				candidates.insert(Coord3D(x, y, z+1));
			    }
			}
		    }
		}
	    }
	}
    }

    // looping through candidates and computing their distance to the interface
    set<Coord3D >::const_iterator it;
    double l, r, f, b, u, d; // left, right, front, back, up, down
    for (it = candidates.begin(); it != candidates.end(); ++it) {
	const int x = it->get<0>();
	const int y = it->get<1>();
	const int z = it->get<2>();
	l = (x > 0) ? phi(x-1, y, z) : UNDEF_VAL;
	r = (x < X-1) ? phi(x+1, y, z) : UNDEF_VAL;
	f = (y > 0) ? phi(x, y-1, z) : UNDEF_VAL;
	b = (y < Y-1) ? phi(x, y+1, z) : UNDEF_VAL;
	u = (z > 0) ? phi(x, y, z-1) : UNDEF_VAL;
	d = (z < Z-1) ? phi(x, y, z+1) : UNDEF_VAL;

	double estimate = compute_dist_estimate_3D(l, r, f, b, u, d);

	if (estimate > 0) {
	    positive_map.insert(pair<double, Coord3D >(estimate, *it)); 
	} else {
	    negative_map.insert(pair<double, Coord3D >(estimate, *it)); 
	}	
    }
}

//===========================================================================
double compute_dist_estimate_2D(const double x1, const double x2,
				const double y1, const double y2)
//===========================================================================
{
    assert(x1 * x2 >= 0 || (x1 == UNDEF_VAL || x2 == UNDEF_VAL));
    assert(y1 * y2 >= 0 || (y1 == UNDEF_VAL || y2 == UNDEF_VAL)); // should have same sign
    assert(x1 * y1 >= 0 || (x1 == UNDEF_VAL || y1 == UNDEF_VAL)); // should have same sign

    const double fac = (x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0) ? -1 : 1; // sign
    const double x_min = fabs(x1) < fabs(x2) ? fabs(x1) : fabs(x2);
    const double y_min = fabs(y1) < fabs(y2) ? fabs(y1) : fabs(y2);

    if (x_min == UNDEF_VAL) {
	assert (y_min != UNDEF_VAL);
	return fac * (y_min + 1);
    } 
    if (y_min == UNDEF_VAL) {
	assert (x_min != UNDEF_VAL);
	return fac * (x_min + 1);
    }

    // both x_min and y_min have defined values, we have to solve a quadratic
    // equation with a non-zero discriminant  (ref. book "Level Set Methods and 
    // Dynamic Implicit Surfaces" (Osher/Sethian), page 73.)

    // no numerical trouble, proceed with resolution of quadratic equation
    return fac * distance_poly(x_min, y_min);
}

//===========================================================================
double compute_dist_estimate_3D(const double x1, const double x2,
				const double y1, const double y2,
				const double z1, const double z2)
//===========================================================================
{
    assert(x1 * x2 >= 0 || (x1 == UNDEF_VAL || x2 == UNDEF_VAL));
    assert(y1 * y2 >= 0 || (y1 == UNDEF_VAL || y2 == UNDEF_VAL)); // should have same sign
    assert(z1 * z2 >= 0 || (z1 == UNDEF_VAL || z2 == UNDEF_VAL));
    assert(x1 * y1 >= 0 || (x1 == UNDEF_VAL || y1 == UNDEF_VAL)); // should have same sign
    assert(x1 * z1 >= 0 || (x1 == UNDEF_VAL || z1 == UNDEF_VAL));

    const double fac = 
	(x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0 || z1 < 0 || z2 < 0) ? -1 : 1; // sign

    const double x_min = fabs(x1) < fabs(x2) ? fabs(x1) : fabs(x2);
    const double y_min = fabs(y1) < fabs(y2) ? fabs(y1) : fabs(y2);
    const double z_min = fabs(z1) < fabs(z2) ? fabs(z1) : fabs(z2); 

    // for an explanation of the below, read the text in the book
    // "Level Set Methods and Dynamic Implicit Surfaces" by Osher/Sethian, page 73.

    if (y_min == UNDEF_VAL && z_min == UNDEF_VAL) {
	assert(x_min != UNDEF_VAL);
	return fac * (x_min + 1);
    } else if (x_min == UNDEF_VAL && z_min == UNDEF_VAL) {
	assert(y_min != UNDEF_VAL);
	return fac * (y_min + 1);
    } else if (x_min == UNDEF_VAL && y_min == UNDEF_VAL) {
	assert(z_min != UNDEF_VAL);
	return fac * (z_min + 1);
    } else if (x_min == UNDEF_VAL) {
	assert(y_min != UNDEF_VAL && z_min != UNDEF_VAL);
	return fac * distance_poly(y_min, z_min); 
    } else if (y_min == UNDEF_VAL) {
	assert(x_min != UNDEF_VAL && z_min != UNDEF_VAL);
	return fac * distance_poly(x_min, z_min);
    } else if (z_min == UNDEF_VAL) {
	assert(x_min != UNDEF_VAL && y_min != UNDEF_VAL);
	return fac * distance_poly(x_min, y_min);
    } 

    // if we got here, then x_min, y_min and z_min should all be defined.
    return fac * distance_poly(x_min, y_min, z_min);
}

//===========================================================================
double distance_poly(double d1, double d2)
//===========================================================================
{
    // consistency check
    double ccheck = fabs(d1 - d2);
    if (ccheck > 1) {
	// numerical trouble, resort to "backup" solution, which is
	// the smallest distance + 1 (which is considered pixel distance)
	return d1 < d2 ? (d1 + 1) : (d2 + 1);
    }

    const double b = -2 * (d1 + d2);
    const double c = (d1 * d1) + (d2 * d2) - 1;
    const double discriminant = b * b - 8 * c;
    assert(discriminant >= 0); // should be assured before calling this function
    return 0.25 * (-b + sqrt(discriminant));
}


//===========================================================================
double distance_poly(double d1, double d2, double d3)
//===========================================================================
{
    // bubble-sorting d1, d2 and d3 according to increasing value
    if (d1 > d2) swap(d1, d2);
    if (d2 > d3) swap(d2, d3);
    if (d1 > d2) swap(d1, d2);
    assert(d1 <= d2 && d2 <= d3);

    // consistency check
    const double tmp1 = d3 - d1;
    const double tmp2 = d3 - d2;
    const double ccheck = (tmp1 * tmp1) + (tmp2 * tmp2);
    if (ccheck > 1) {
	// try solving poly with two variables
	return distance_poly(d1, d2);
    } 
    
    const double six_inv = double(1) / double(6);
    // (a is equal to 3, but this is hardcoded into the solution)
    const double b = -2 * (d1 + d2 + d3);
    const double c = (d1 * d1) + (d2 * d2) + (d3 * d3) - 1;
    const double discriminant = (b * b - 12 * c);
    assert(discriminant >= 0); // should be assured before calling this function
    return six_inv * (-b + sqrt(discriminant));
}


//===========================================================================
void march_out_2D(dist_map_2D& tentatives, 
		  coord_map_2D& inv_tent, 
		  lsseg::LevelSetFunction& target, 
		  const Mask* m,
		  bool negative)
//===========================================================================
{
    const int X = target.dimx();
    const int Y = target.dimy();
    
    dist_map_2D::iterator tmp_it, cur_it;
    coord_map_2D::iterator c_it;
    
    double l, r, b, t; // left, right, bottom, top, ...

    while (!tentatives.empty()) {
	if (negative) {
	    cur_it = tentatives.end();
	    --cur_it;
	} else {
	    cur_it = tentatives.begin();
	}

	const int x = cur_it->second.get<0>();
	const int y = cur_it->second.get<1>();
	target(x, y) = cur_it->first;
	
	// updating neighbours
	if (x > 0 && target(x-1, y) == UNDEF_VAL && (!m || (*m)(x-1, y))) {  // value to the left of current
	    r = cur_it->first;
	    l = (x > 1) ? target(x-2, y) : UNDEF_VAL;
	    b = (y > 0) ? target(x-1, y-1) : UNDEF_VAL;
	    t = (y < Y-1) ? target(x-1, y+1) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_2D(l, r, b, t);
	    const Coord2D xy(x-1, y);
	    c_it = inv_tent.find(xy);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord2D >( dist, xy));
	    inv_tent[xy]= tmp_it;
	}
	if (x < X-1 && target(x+1, y) == UNDEF_VAL && (!m || (*m)(x+1, y))) { // value to the right of current
	    r = (x < X-2) ? target(x+2, y) : UNDEF_VAL;
	    l = cur_it->first;
	    b = (y > 0) ? target(x+1,y-1) : UNDEF_VAL;
	    t = (y < Y-1) ? target(x+1, y+1) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_2D(l, r, b, t);
	    const Coord2D xy(x+1, y);
	    c_it = inv_tent.find(xy);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord2D >(dist,xy));
	    inv_tent[xy] = tmp_it;
	}
	if (y > 0 && target(x, y-1) == UNDEF_VAL && (!m || (*m)(x, y-1))) { // value on the bottom of current
	    r = (x < X-1) ? target(x+1, y-1) : UNDEF_VAL;
	    l = (x > 0) ? target(x-1, y-1) : UNDEF_VAL;
	    b = (y > 1) ? target(x, y-2) : UNDEF_VAL;
	    t = cur_it->first;
	    const double dist = compute_dist_estimate_2D(l, r, b, t);
	    const Coord2D xy(x, y-1);
	    c_it = inv_tent.find(xy);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord2D >(dist, xy));
	    inv_tent[xy] = tmp_it;
	}
	if (y < Y-1 && target(x, y+1) == UNDEF_VAL && (!m || (*m)(x, y+1))) { // value on top of current
	    r = (x < X-1) ? target(x+1, y+1) : UNDEF_VAL;
	    l = (x > 0) ? target(x-1, y+1): UNDEF_VAL;
	    b = cur_it->first;
	    t = (Y < Y-2) ? target(x, y+2) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_2D(l, r, b, t);
	    const Coord2D xy(x, y+1);
	    c_it= inv_tent.find(xy);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord2D >(dist, xy));
	    inv_tent[xy] = tmp_it;
	}
	assert(inv_tent.size() == tentatives.size());
	
	// removing element
	assert(inv_tent[cur_it->second] == cur_it);
	inv_tent.erase(cur_it->second);
	tentatives.erase(cur_it);
    }
}

//===========================================================================
void march_out_3D(dist_map_3D& tentatives, 
		  coord_map_3D& inv_tent, 
		  LevelSetFunction& target,
		  const Mask* m,
		  bool negative)
//===========================================================================
{
    const int X = target.dimx();
    const int Y = target.dimy();
    const int Z = target.dimz();
    
    dist_map_3D::iterator tmp_it, cur_it;
    coord_map_3D::iterator c_it;

    double l, r, f, b, u, d; // left, right, front, back, up, down

    while (!tentatives.empty()) {
	if (negative) {
	    cur_it = tentatives.end();
	    --cur_it;
	} else {
	    cur_it = tentatives.begin();
	}

	const int x = cur_it->second.get<0>();
	const int y = cur_it->second.get<1>();
	const int z = cur_it->second.get<2>();
	target(x, y, z) = cur_it->first;
	
	// updating neighbours

	// value to the left of current (x-axis)
	if (x > 0 && target(x-1, y, z) == UNDEF_VAL && (!m || (*m)(x-1, y, z))) {  
	    r = cur_it->first;
	    l = (x > 1) ? target(x-2, y, z) : UNDEF_VAL;
	    b = (y > 0) ? target(x-1, y-1, z) : UNDEF_VAL;
	    f = (y < Y-1) ? target(x-1, y+1, z) : UNDEF_VAL;
	    d = (z > 0) ? target(x-1, y, z-1) : UNDEF_VAL;
	    u = (z < Z-1) ? target(x-1, y, z+1) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_3D(l, r, b, f, d, u);
	    const Coord3D xyz(x-1, y, z);
	    c_it = inv_tent.find(xyz);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord3D >( dist, xyz));
	    inv_tent[xyz]= tmp_it;
	}

	// value to the right of current (x-axis)
	if (x < X-1 && target(x+1, y, z) == UNDEF_VAL && (!m || (*m)(x+1, y, z))) { 
	    r = (x < X-2) ? target(x+2, y, z) : UNDEF_VAL;
	    l = cur_it->first;
	    b = (y > 0) ? target(x+1,y-1, z) : UNDEF_VAL;
	    f = (y < Y-1) ? target(x+1, y+1, z) : UNDEF_VAL;
	    d = (z > 0) ? target(x+1, y, z-1) : UNDEF_VAL;
	    u = (z < Z-1) ? target(x+1, y, z+1) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_3D(l, r, b, f, d, u);
	    const Coord3D xyz(x+1, y, z);
	    c_it = inv_tent.find(xyz);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord3D >(dist,xyz));
	    inv_tent[xyz] = tmp_it;
	}

	// value on the back of current (y-axis)
	if (y > 0 && target(x, y-1, z) == UNDEF_VAL && (!m || (*m)(x, y-1, z))) { 
	    r = (x < X-1) ? target(x+1, y-1, z) : UNDEF_VAL;
	    l = (x > 0) ? target(x-1, y-1, z) : UNDEF_VAL;
	    b = (y > 1) ? target(x, y-2, z) : UNDEF_VAL;
	    f = cur_it->first;
	    d = (z > 0) ? target(x, y-1, z-1) : UNDEF_VAL;
	    u = (z < Z-1) ? target(x, y-1, z+1) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_3D(l, r, b, f, d, u);
	    const Coord3D xyz(x, y-1, z);
	    c_it = inv_tent.find(xyz);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord3D >(dist, xyz));
	    inv_tent[xyz] = tmp_it;
	}

	// value on front of current (y-axis)
	if (y < Y-1 && target(x, y+1, z) == UNDEF_VAL && (!m || (*m)(x, y+1, z))) { 
	    r = (x < X-1) ? target(x+1, y+1, z) : UNDEF_VAL;
	    l = (x > 0) ? target(x-1, y+1, z): UNDEF_VAL;
	    b = cur_it->first;
	    f = (Y < Y-2) ? target(x, y+2, z) : UNDEF_VAL;
	    d = (z > 0) ? target(x, y+1, z-1) : UNDEF_VAL;
	    u = (z < Z-1) ? target(x, y+1, z+1) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_3D(l, r, b, f, d, u);
	    const Coord3D xyz(x, y+1, z);
	    c_it= inv_tent.find(xyz);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord3D >(dist, xyz));
	    inv_tent[xyz] = tmp_it;
	}
	// value below current (z-axis)
	if (z > 0 && target(x, y, z-1) == UNDEF_VAL && (!m || (*m)(x, y, z-1))) {
	    r = (x < X-1) ? target(x+1, y, z-1) : UNDEF_VAL;
	    l = (x > 0) ? target(x-1, y, z-1) : UNDEF_VAL;
	    b = (y > 0) ? target(x, y-1, z-1) : UNDEF_VAL;
	    f = (y < Y-1) ? target(x, y+1, z-1) : UNDEF_VAL;
	    d = (z > 1) ? target(x, y, z-2) : UNDEF_VAL;
	    u = cur_it->first;
	    const double dist = compute_dist_estimate_3D(l, r, b, f, d, u);
	    const Coord3D xyz(x, y, z-1);
	    c_it = inv_tent.find(xyz);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord3D >(dist, xyz));
	    inv_tent[xyz] = tmp_it;
	}

	// value above current (z-axis)
	if (z < Z-1 && target(x, y, z+1) == UNDEF_VAL && (!m || (*m)(x, y, z+1))) {
	    r = (x < X-1) ? target(x+1, y, z+1) : UNDEF_VAL;
	    l = (x > 0) ? target(x-1, y, z+1) : UNDEF_VAL;
	    b = (y > 0) ? target(x, y-1, z+1) : UNDEF_VAL;
	    f = (y < Y-1) ? target(x, y+1, z+1) : UNDEF_VAL;
	    d = cur_it->first;
	    u = (z < Z-2) ? target(x, y, z+2) : UNDEF_VAL;
	    const double dist = compute_dist_estimate_3D(l, r, b, f, d, u);
	    const Coord3D xyz(x, y, z+1);
	    c_it = inv_tent.find(xyz);
	    if (c_it != inv_tent.end()) {
		tentatives.erase(c_it->second);
	    }
	    tmp_it = tentatives.insert(pair<double, Coord3D>(dist, xyz));
	    inv_tent[xyz] = tmp_it;
	}
	assert(inv_tent.size() == tentatives.size());
	
	// removing element
	tentatives.erase(cur_it);
	assert(inv_tent[cur_it->second] == cur_it);
	inv_tent.erase(cur_it->second);
    }
}

}; // end anonymous namespace

