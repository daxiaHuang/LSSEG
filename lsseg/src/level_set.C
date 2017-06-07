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
// File: level_set.C                                                         
//                                                                           
// Created: Tue Oct 25 14:02:40 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: level_set.C,v 1.35 2006/11/25 20:08:28 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements level_set.h.
//                                                                           
//===========================================================================

#include <functional>
#include <stdexcept>
#include <time.h>
#include <vector>
#include <boost/lambda/lambda.hpp>
#include "Region.h"
#include "errormacros.h"
#include "Image.h"
#include "LevelSetFunction.h"
#include "level_set.h"
#include "simple_tools.h"
#include "cimg_dependent.h"
//#include "GoTensorProductSpline.h"
//#include "GLviewUtils.h"
#include "Mask.h"
#include "Filters.h"
#include "DiscreteApproximations.h"

using namespace boost::lambda;
using namespace std;
//using namespace Go;
using namespace lsseg;

namespace  // anonymous namespace 
{
    void visualize_multisets_impl(vector<const LevelSetFunction*> phi_pts,
				  Image<double>& target,
				  const int** const rgb_color);

    inline void godunov_phi_2D(double a_term, 
			       double phi_xl, double phi_xr,
			       double phi_yl, double phi_yr,
			       double& phi_x,
			       double& phi_y);

    inline void godunov_phi_3D(double a_term, 
			       double phi_xl, double phi_xr,
			       double phi_yl, double phi_yr,
			       double phi_zl, double phi_zr,
			       double& phi_x,
			       double& phi_y,
			       double& phi_z);

//===========================================================================
// differential approximation schemes
//===========================================================================

inline bool defined(int x, int y, int z, int X, int Y, int Z, const Mask* mask) 
{
    return (x >= 0 && x < X &&
	    y >= 0 && y < Y &&
	    z >= 0 && z < Z &&
	    (!mask || (*mask)(x, y, z)));
}

double neumann_curvature_at(int x, 
			    int y, 
			    int z, 
			    const Image<double>& phi, 
			    const Mask* mask);

double dirichlet_curvature_at(int x, 
			      int y,
			      int z, 
			      const Image<double>& phi,
			      const Mask* mask);

void neumann_d1_leftright_2D(int x, 
			     int y, 
			     const Image<double>& phi, 
			     const Mask* mask,
			     double& dxl, 
			     double& dxr, 
			     double& dyl, 
			     double& dyr);

void dirichlet_d1_leftright_2D(int x, 
			       int y, 
			       const Image<double>& phi, 
			       const Mask* mask,
			       double& dxl, 
			       double& dxr, 
			       double& dyl, 
			       double& dyr);

void neumann_d1_leftright_3D(int x, 
			     int y, 
			     int z, 
			     const Image<double>& phi, 
			     const Mask* mask,
			     double& dxl, 
			     double& dxr, 
			     double& dyl, 
			     double& dyr, 
			     double& dzl, 
			     double& dzr);

void dirichlet_d1_leftright_3D(int x, 
			       int y, 
			       int z, 
			       const Image<double>& phi, 
			       const Mask* mask,
			       double& dxl, 
			       double& dxr, 
			       double& dyl, 
			       double& dyr, 
			       double& dzl, 
			       double& dzr);
    
}; // end anonymous namespace 

namespace lsseg {


//==============================================================================
void normal_direction_flow_2D(const LevelSetFunction& phi, 
			      LevelSetFunction& advect,
			      BorderCondition bcond,
			      double& H1,
			      double& H2,
			      const Mask* const changes_allowed_mask,
			      const Mask* const defined_region_mask)
//==============================================================================
{
    assert(phi.dimz() == 1); // this should be the 2D case
    assert(advect.size_compatible(phi));
    assert(phi.numChannels() == 1);
    const int X = phi.dimx();
    const int Y = phi.dimy();
    H1 = H2 = 0;

    void (*der_estimate)(int, int, const Image<double>&, const Mask*,
			 double&, double&, double&, double&);
    if (bcond == NEUMANN) {
	der_estimate = neumann_d1_leftright_2D;
    } else {
	assert(bcond == DIRICHLET);
	der_estimate = dirichlet_d1_leftright_2D;
    }

    const Mask::value_type* m_it = changes_allowed_mask ? changes_allowed_mask->begin() : 0;

    for (int yc = 0; yc < Y; ++yc) {
	for (int xc = 0; xc < X; ++xc) {
	    if (changes_allowed_mask && !(*m_it++)) {
		// the changes_allowed_mask is defined and is false at this pixel, => no advection
		advect(xc, yc) = 0;
		continue;
	    }
	    const double a_term = advect(xc, yc);
	    double phi_xl, phi_xr, phi_yl, phi_yr;
	    
	    der_estimate(xc, yc, phi, 
			 defined_region_mask, 
			 phi_xl, phi_xr,
			 phi_yl, phi_yr);
	    
	    double phi_x, phi_y;
	    godunov_phi_2D(a_term, 
			   phi_xl, phi_xr, 
			   phi_yl, phi_yr, 
			   phi_x, 
			   phi_y);
	    
	    const double norm_phi2 = phi_x * phi_x + phi_y * phi_y;
	    const double norm_phi = sqrt(norm_phi2);
	    advect(xc, yc) = -a_term * norm_phi; 
	    
	    if (norm_phi > 0) {
		const double norm_inv = double(1) / norm_phi;
		const double fac = a_term * norm_inv;
		double tmp = fabs(phi_x * fac);
		H1 = H1 > tmp ? H1 : tmp;
		tmp = fabs(phi_y * fac);
		H2 = H2 > tmp ? H2 : tmp;
	    }
	}
    }
}

//==============================================================================
void normal_direction_flow_3D(const LevelSetFunction& phi, 
			      LevelSetFunction& advect,
			      BorderCondition bcond,
			      double& H1,
			      double& H2,
			      double& H3,
			      const Mask* const changes_allowed_mask,
			      const Mask* const defined_region_mask)
//==============================================================================
{
    assert(advect.size_compatible(phi));
    assert(phi.numChannels() == 1);
    const int X = phi.dimx();
    const int Y = phi.dimy();
    const int Z = phi.dimz();
    H1 = H2 = H3 = 0;

    void (*der_estimate)(int, int, int, const Image<double>&, const Mask*,
			 double&, double&, double&, double&, double&, double&);
    if (bcond == NEUMANN) {
	der_estimate = neumann_d1_leftright_3D;
    } else {
	assert(bcond == DIRICHLET);
	der_estimate = dirichlet_d1_leftright_3D;
    }

    const Mask::value_type* m_it = changes_allowed_mask ? changes_allowed_mask->begin() : 0;

    for (int zc = 0; zc < Z; ++zc) {
	for (int yc = 0; yc < Y; ++yc) {
	    for (int xc = 0; xc < X; ++xc) {
		if (changes_allowed_mask && !(*m_it++)) {
		    // the changes_allowed_mask is defined and is false at this pixel, => no advection
		    advect(xc, yc, zc) = 0;
		    continue;
		}
		const double a_term = advect(xc, yc, zc);
		
		double phi_xl, phi_xr, phi_yl, phi_yr, phi_zl, phi_zr;
		
		der_estimate(xc, yc, zc, phi, 
			     defined_region_mask, 
			     phi_xl, phi_xr,
			     phi_yl, phi_yr,
			     phi_zl, phi_zr);
		double phi_x, phi_y, phi_z;
		godunov_phi_3D(a_term, 
			       phi_xl, phi_xr, 
			       phi_yl, phi_yr, 
			       phi_zl, phi_zr,
			       phi_x, 
			       phi_y,
			       phi_z);
		const double norm_phi2 = phi_x * phi_x + phi_y * phi_y + phi_z * phi_z;
		const double norm_phi = sqrt(norm_phi2);
		advect(xc, yc, zc) = -a_term * norm_phi;
		
		if (norm_phi > 0) {
		    const double norm_inv = double(1) / norm_phi;
		    const double fac = a_term * norm_inv;
		    double tmp = fabs(phi_x * fac);
		    
		    H1 = H1 > tmp ? H1 : tmp;
		    tmp = fabs(phi_y * fac);
		    H2 = H2 > tmp ? H2 : tmp;
		    tmp = fabs(phi_z * fac);
		    H3 = H3 > tmp ? H3 : tmp;
		}
	    }
	}
    }
}

//===========================================================================
void visualize_multisets(const LevelSetFunction** const images,
			 int num_images,
			 Image<double>& target,
			 const int** const rgb_color) 
//===========================================================================
{
    vector<const LevelSetFunction* > phi_pointers(images, images + num_images);
//     for (int i = 0; i < num_images; ++i) {
// 	phi_pointers[i] = const_cast<Image<double>* >(images + i);
//     }
    visualize_multisets_impl(phi_pointers, target, rgb_color);
}

//===========================================================================
int visualize_level_set(const LevelSetFunction& img, 
			 Image<double>& target, 
			 double threshold,
			 double outside_intensity,
			 double curve_intensity,
			 double interior_intensity,
			 const Mask* const mask,
			 double mask_intensity)
//===========================================================================
{
    int num_border_pixels = 0;
    ALWAYS_ERROR_IF(img.numChannels() != 1, "dim should be 1 here");
    assert(img.size_compatible(target));
    assert(!mask || img.size_compatible(*mask)); // check for size if not null pointer
    bool no_mask = (mask == 0);
    for (int z = 0; z < img.dimz(); ++z) {
	for (int y = 0; y < img.dimy(); ++y) {
	    for (int x = 0; x < img.dimx(); ++x) {
		if (no_mask || (*mask)(x, y, z) > 0) {
		    double Iccc = img(x, y, z);
		    double Icnc = img(x, (y == img.dimy() - 1) ? y : y+1, z);
		    double Incc = img((x == img.dimx() - 1) ? x : x+1, y, z);
		    double Iccn = img(x, y, (z == img.dimz() - 1) ? z : z+1);
		    double Icpc = img(x, (y==0) ? y : y-1, z);
		    double Ipcc = img((x==0) ? 0 : x-1, y, z);
		    double Iccp = img(x, y, (z==0) ? 0 : z-1);

		    Iccc -= threshold;
		    Icnc -= threshold;
		    Incc -= threshold;
		    Iccn -= threshold;
		    Icpc -= threshold;
		    Ipcc -= threshold;
		    Iccp -= threshold;
		    
		    if (Iccc * Incc < 0 || (Iccc == 0 && Incc != 0) ||
			Iccc * Icnc < 0 || (Iccc == 0 && Icnc != 0) ||
			Iccc * Iccn < 0 || (Iccc == 0 && Iccn != 0) ||
			Iccc * Ipcc < 0 || (Iccc == 0 && Ipcc != 0) ||
			Iccc * Icpc < 0 || (Iccc == 0 && Icpc != 0) ||
			Iccc * Iccp < 0 || (Iccc == 0 && Iccp != 0)) {
			target(x, y, z) = curve_intensity;
			++num_border_pixels;
		    } else {
			target(x, y, z) = Iccc < 0 ? interior_intensity : outside_intensity;
		    }
		} else {
		    target(x, y, z) = mask_intensity;
		}
	    }
	}
    }
    return num_border_pixels;
}

}; // end namespace lsseg

namespace  // anonymous namespace 
{

//===========================================================================
void neumann_d1_leftright_2D(int x, 
			     int y, 
			     const Image<double>& phi, 
			     const Mask* mask,
			     double& dxl, 
			     double& dxr, 
			     double& dyl, 
			     double& dyr)
//===========================================================================
{
    assert(phi.dimz() == 1);
    const int X = phi.dimx();
    const int Y = phi.dimy();

    // negative value means that the index has not been set
    int xp(-1), yp(-1), xn(-1), yn(-1);

    if (!mask) {
	if (x >   0) { xp = x-1; }
	if (x < X-1) { xn = x+1; }
	if (y >   0) { yp = y-1; }
	if (y < Y-1) { yn = y+1; }
    } else {
	if (x >   0 && (*mask)(x-1, y)) { xp = x-1; }
	if (x < X-1 && (*mask)(x+1, y)) { xn = x+1; }
	if (y >   0 && (*mask)(x, y-1)) { yp = y-1; }
	if (y < Y-1 && (*mask)(x, y+1)) { yn = y+1; }
    }
    const double Icc = phi(x, y);
    const double Ipc = (xp >= 0) ? phi(xp, y) : (xn >= 0) ? 2 * Icc - phi(xn, y) : Icc;
    const double Inc = (xn >= 0) ? phi(xn, y) : (xp >= 0) ? 2 * Icc - phi(xp, y) : Icc;
    const double Icp = (yp >= 0) ? phi(x, yp) : (yn >= 0) ? 2 * Icc - phi(x, yn) : Icc;
    const double Icn = (yn >= 0) ? phi(x, yn) : (yp >= 0) ? 2 * Icc - phi(x, yp) : Icc;

    dxl = Icc - Ipc;
    dxr = Inc - Icc;
    dyl = Icc - Icp;
    dyr = Icn - Icc;

}

//===========================================================================
void dirichlet_d1_leftright_2D(int x, 
			       int y, 
			       const Image<double>& phi, 
			       const Mask* mask,
			       double& dxl, 
			       double& dxr, 
			       double& dyl, 
			       double& dyr)
//===========================================================================
{
    const int X = phi.dimx();
    const int Y = phi.dimy();

    // negative value means that the index has not been set
    int xp, yp, xn, yn;
    if (!mask) {
	xp = (x > 0) ? x-1 : x;
	yp = (y > 0) ? y-1 : x;
	xn = (x < X-1) ? x+1 : x;
	yn = (y < Y-1) ? y+1 : y;
    } else {
	xp = (x > 0 && (*mask)(x-1,   y)) ? x-1 : x;
	yp = (y > 0 && (*mask)(  x, y-1)) ? y-1 : y;
	xn = (x < X-1 && (*mask)(x+1,   y)) ? x+1 : x;
	yn = (y < Y-1 && (*mask)(  x, y+1)) ? y+1 : y;
    }
    
    double Icc = phi(x, y);
    double Ipc = phi(xp, y);
    double Inc = phi(xn, y);
    double Icp = phi(x, yp);
    double Icn = phi(x, yn);

    dxl = Icc - Ipc;
    dxr = Inc - Icc;
    dyl = Icc - Icp;
    dyr = Icn - Icc;
}


//===========================================================================
void dirichlet_d1_leftright_3D(int x, 
			       int y, 
			       int z, 
			       const Image<double>& phi, 
			       const Mask* mask,
			       double& dxl, 
			       double& dxr, 
			       double& dyl, 
			       double& dyr, 
			       double& dzl, 
			       double& dzr)
//===========================================================================
{
    const int X = phi.dimx();
    const int Y = phi.dimy();
    const int Z = phi.dimz();

    // negative value means that the index has not been set
    int xp, yp, zp, xn, yn, zn;
    if (!mask) {
	xp = (x > 0) ? x-1 : x;
	yp = (y > 0) ? y-1 : x;
	zp = (z > 0) ? z-1 : x;
	xn = (x < X-1) ? x+1 : x;
	yn = (y < Y-1) ? y+1 : y;
	zn = (z < Z-1) ? z+1 : z;
    } else {
	xp = (x > 0 && (*mask)(x-1,   y,   z)) ? x-1 : x;
	yp = (y > 0 && (*mask)(  x, y-1,   z)) ? y-1 : y;
	zp = (z > 0 && (*mask)(  x,   y, z-1)) ? z-1 : z;
	xn = (x < X-1 && (*mask)(x+1,   y,   z)) ? x+1 : x;
	yn = (y < Y-1 && (*mask)(  x, y+1,   z)) ? y+1 : y;
	zn = (z < Z-1 && (*mask)(  x,   y, z+1)) ? z+1 : z;
    }
    
    double Iccc = phi(x, y, z);
    double Ipcc = phi(xp, y, z);
    double Incc = phi(xn, y, z);
    double Icpc = phi(x, yp, z);
    double Icnc = phi(x, yn, z);
    double Iccp = phi(x, y, zp);
    double Iccn = phi(x, y, zn);

    dxl = Iccc - Ipcc;
    dxr = Incc - Iccc;
    dyl = Iccc - Icpc;
    dyr = Icnc - Iccc;
    dzl = Iccc - Iccp;
    dzr = Iccn - Iccc;
}


//===========================================================================
void neumann_d1_leftright_3D(int x, 
			     int y, 
			     int z, 
			     const Image<double>& phi, 
			     const Mask* mask,
			     double& dxl, 
			     double& dxr, 
			     double& dyl, 
			     double& dyr, 
			     double& dzl, 
			     double& dzr)
//===========================================================================
{
    const int X = phi.dimx();
    const int Y = phi.dimy();
    const int Z = phi.dimz();

    // negative value means that the index has not been set
    int xp(-1), yp(-1), zp(-1), xn(-1), yn(-1), zn(-1);

    if (!mask) {
	if (x >   0) { xp = x-1; }
	if (x < X-1) { xn = x+1; }
	if (y >   0) { yp = y-1; }
	if (y < Y-1) { yn = y+1; }
	if (z >   0) { zp = z-1; }
	if (z < Z-1) { zn = z+1; }
    } else {
	if (x >   0 && (*mask)(x-1, y, z)) { xp = x-1; }
	if (x < X-1 && (*mask)(x+1, y, z)) { xn = x+1; }
	if (y >   0 && (*mask)(x, y-1, z)) { yp = y-1; }
	if (y < Y-1 && (*mask)(x, y+1, z)) { yn = y+1; }
	if (z >   0 && (*mask)(x, y, z-1)) { zp = z-1; }
	if (z < Z-1 && (*mask)(x, y, z+1)) { zn = z+1; }
    }
    const double Iccc = phi(x, y, z);
    const double Ipcc = (xp >= 0) ? phi(xp, y, z) : (xn >= 0) ? 2 * Iccc - phi(xn, y, z) : Iccc;
    const double Incc = (xn >= 0) ? phi(xn, y, z) : (xp >= 0) ? 2 * Iccc - phi(xp, y, z) : Iccc;
    const double Icpc = (yp >= 0) ? phi(x, yp, z) : (yn >= 0) ? 2 * Iccc - phi(x, yn, z) : Iccc;
    const double Icnc = (yn >= 0) ? phi(x, yn, z) : (yp >= 0) ? 2 * Iccc - phi(x, yp, z) : Iccc;
    const double Iccp = (zp >= 0) ? phi(x, y, zp) : (zn >= 0) ? 2 * Iccc - phi(x, y, zn) : Iccc;
    const double Iccn = (zn >= 0) ? phi(x, y, zn) : (zp >= 0) ? 2 * Iccc - phi(x, y, zp) : Iccc;


    dxl = Iccc - Ipcc;
    dxr = Incc - Iccc;
    dyl = Iccc - Icpc;
    dyr = Icnc - Iccc;
    dzl = Iccc - Iccp;
    dzr = Iccn - Iccc;
}

//===========================================================================
double dirichlet_curvature_at(int x, 
			      int y,
			      int z, 
			      const Image<double>& phi,
			      const Mask* mask)
//===========================================================================
{
    MESSAGE("dirichlet_curvature_at() unimplemented.");
    exit(0);
}

//===========================================================================
double neumann_curvature_at(int x, 
			    int y, 
			    int z, 
			    const Image<double>& phi, 
			    const Mask* mask)
//===========================================================================
{
    const double EPS=1.0e-10;
    const int X = phi.dimx();
    const int Y = phi.dimy();
    const int Z = phi.dimz();
    
    // negative value means that the index has not been set
   int xp(-1), yp(-1), zp(-1), xn(-1), yn(-1), zn(-1);

    if (!mask) {
	if (x >   0) { xp = x-1; }
	if (x < X-1) { xn = x+1; }
	if (y >   0) { yp = y-1; }
	if (y < Y-1) { yn = y+1; }
	if (z >   0) { zp = z-1; }
	if (z < Z-1) { zn = z+1; }
    } else {
	if (x >   0 && (*mask)(x-1, y, z)) { xp = x-1; }
	if (x < X-1 && (*mask)(x+1, y, z)) { xn = x+1; }
	if (y >   0 && (*mask)(x, y-1, z)) { yp = y-1; }
	if (y < Y-1 && (*mask)(x, y+1, z)) { yn = y+1; }
	if (z >   0 && (*mask)(x, y, z-1)) { zp = z-1; }
	if (z < Z-1 && (*mask)(x, y, x+1)) { zn = z+1; }
    }
    double Iccc = phi(x, y, z);
    double Ipcc = (xp >= 0) ? phi(xp, y, z) : (xn >= 0) ? 2 * Iccc - phi(xn, y, z) : Iccc;
    double Incc = (xn >= 0) ? phi(xn, y, z) : (xp >= 0) ? 2 * Iccc - phi(xp, y, z) : Iccc;
    double Icpc = (yp >= 0) ? phi(x, yp, z) : (yn >= 0) ? 2 * Iccc - phi(x, yn, z) : Iccc;
    double Icnc = (yn >= 0) ? phi(x, yn, z) : (yp >= 0) ? 2 * Iccc - phi(x, yp, z) : Iccc;
    double Iccp = (zp >= 0) ? phi(x, y, zp) : (zn >= 0) ? 2 * Iccc - phi(x, y, zn) : Iccc;
    double Iccn = (zn >= 0) ? phi(x, y, zn) : (zp >= 0) ? 2 * Iccc - phi(x, y, zp) : Iccc;

    double dx = 0.5 * (Incc - Ipcc);
    double dy = 0.5 * (Icnc - Icpc);
    double dz = 0.5 * (Iccn - Iccp);
    
    double norm2 = (dx * dx + dy * dy + dz * dz);
    
    if (norm2 < EPS) return 0; // close to zero gradient norm - no significant curvature

    //double Ippp = defined(x-1,y-1,z-1,X,Y,Z,mask) ? phi(x-1,y-1,z-1) : Ipcc+Icpc+Iccp-2*Iccc;
    double Ippc = defined(x-1,y-1,z  ,X,Y,Z,mask) ? phi(x-1,y-1,z )  : Ipcc+Icpc     -  Iccc;
    //double Ippn = defined(x-1,y-1,z+1,X,Y,Z,mask) ? phi(x-1,y-1,z+1) : Ipcc+Icpc+Iccn-2*Iccc;
    double Ipcp = defined(x-1,y  ,z-1,X,Y,Z,mask) ? phi(x-1,y  ,z-1) : Ipcc+     Iccp-  Iccc;
    double Ipcn = defined(x-1,y  ,z+1,X,Y,Z,mask) ? phi(x-1,y  ,z+1) : Ipcc+     Iccn-  Iccc;
    //double Ipnp = defined(x-1,y+1,z-1,X,Y,Z,mask) ? phi(x-1,y+1,z-1) : Ipcc+Icnc+Iccp-2*Iccc;
    double Ipnc = defined(x-1,y+1,z  ,X,Y,Z,mask) ? phi(x-1,y+1,z  ) : Ipcc+Icnc     -  Iccc;
    //double Ipnn = defined(x-1,y+1,z+1,X,Y,Z,mask) ? phi(x-1,y+1,z+1) : Ipcc+Icnc+Iccn-2*Iccc;
    double Icpp = defined(x  ,y-1,z-1,X,Y,Z,mask) ? phi(x  ,y-1,z-1) :      Icpc+Iccp-  Iccc;
    double Icpn = defined(x  ,y-1,z+1,X,Y,Z,mask) ? phi(x  ,y-1,z+1) :      Icpc+Iccn-  Iccc;
    double Icnp = defined(x  ,y+1,z-1,X,Y,Z,mask) ? phi(x  ,y+1,z-1) :      Icnc+Iccp-  Iccc;
    double Icnn = defined(x  ,y+1,z+1,X,Y,Z,mask) ? phi(x  ,y+1,z+1) :      Icnc+Iccn-  Iccc;
    //double Inpp = defined(x+1,y-1,z-1,X,Y,Z,mask) ? phi(x+1,y-1,z-1) : Incc+Icpc+Iccp-2*Iccc;
    double Inpc = defined(x+1,y-1,z  ,X,Y,Z,mask) ? phi(x+1,y-1,z  ) : Incc+Icpc     -  Iccc;
    //double Inpn = defined(x+1,y-1,z+1,X,Y,Z,mask) ? phi(x+1,y-1,z+1) : Incc+Icpc+Iccn-2*Iccc;
    double Incp = defined(x+1,y  ,z-1,X,Y,Z,mask) ? phi(x+1,y  ,z-1) : Incc+     Iccp-  Iccc;
    double Incn = defined(x+1,y  ,z+1,X,Y,Z,mask) ? phi(x+1,y  ,z+1) : Incc+     Iccn-  Iccc;
    //double Innp = defined(x+1,y+1,z-1,X,Y,Z,mask) ? phi(x+1,y+1,z-1) : Incc+Icnc+Iccp-2*Iccc;
    double Innc = defined(x+1,y+1,z  ,X,Y,Z,mask) ? phi(x+1,y+1,z  ) : Incc+Icnc     -  Iccc;
    //double Innn = defined(x+1,y+1,z+1,X,Y,Z,mask) ? phi(x+1,y+1,z+1) : Incc+Icnc+Iccn-2*Iccc;

    double dxx = Incc + Ipcc - 2 * Iccc;
    double dyy = Icnc + Icpc - 2 * Iccc;
    double dzz = Iccn + Iccp - 2 * Iccc;

    double dxy = (dx * dy < 0) ? 
	0.5 * (2 * Iccc + Ippc + Innc - Ipcc - Incc - Icpc - Icnc) :
	0.5 * (Ipcc + Incc + Icpc + Icnc - 2 * Iccc - Ipnc - Inpc);

    double dxz =  (dx * dz < 0) ? 
	0.5 * (2 * Iccc + Ipcp + Incn - Ipcc - Incc - Iccp - Iccn) :
	0.5 * (Ipcc + Incc + Iccp + Iccn - 2 * Iccc - Ipcn - Incp);

    double dyz = (dy * dz < 0) ?
	0.5 * (2 * Iccc + Icpp + Icnn - Icpc - Icnc - Iccp - Iccn) :
	0.5 * (Icpc + Icnc + Iccp + Iccn - 2 * Iccc - Icpn - Icnp);

    double norm = sqrt(norm2);
    
    return (dx * dx * (dyy + dzz) + 
	    dy * dy * (dxx + dzz) +
	    dz * dz * (dxx + dyy) - 
	    2 * (dx * dy * dxy + 
		 dx * dz * dxz +
		 dy * dz * dyz)
	    ) / (norm * norm2);
}

//===========================================================================
void visualize_multisets_impl(vector<const LevelSetFunction* > images,
			      Image<double>& target,
			      const int** const rgb_color)
//===========================================================================
{
    const int num_images = images.size();
    assert(num_images > 0);
    assert(images[0]->dimz() == 1); // currently only implemented in 2D
    const int X = images[0]->dimx();
    const int Y = images[0]->dimy();
    target.resize(X, Y, 1, 3); // 3 color channels

    int blend[3]; 
    int num_covering = 0;
    for (int y = 0; y < Y; ++y) {
	for (int x = 0; x < X; ++x) {
	    blend[0] = blend[1] = blend[2] = 0;
	    num_covering = 0;
	    for (int i = 0; i < num_images; ++i) {
		if ((*images[i])(x, y) <= 0) {
		    // this function is in the overlap
		    blend[0] += rgb_color[i][0];
		    blend[1] += rgb_color[i][1];
		    blend[2] += rgb_color[i][2];
		    ++num_covering;
		}
	    }
	    if (num_covering > 0) {
		target(x, y, 0, 0) = blend[0] / num_covering;
		target(x, y, 0, 1) = blend[1] / num_covering;
		target(x, y, 0, 2) = blend[2] / num_covering;
	    } else {
 		target(x, y, 0, 0) = rgb_color[num_images][0];
 		target(x, y, 0, 1) = rgb_color[num_images][1];
 		target(x, y, 0, 2) = rgb_color[num_images][2];
	    }
	}
    }
}

//===========================================================================
void godunov_phi_2D(double a_term, 
		    double phi_xl, double phi_xr,
		    double phi_yl, double phi_yr,
		    double& phi_x,
		    double& phi_y)
//===========================================================================
{
    // using Godunov to choose correct phi_x and phi_y 
    if (a_term * phi_xr > 0) {
	phi_x = (a_term * phi_xl > 0) ? phi_xl : 0;
    } else { // a_term * phi_xr < 0
	if (a_term * phi_xl > 0) {
	    phi_x = (fabs(phi_xl) > fabs(phi_xr)) ? phi_xl : phi_xr;
	} else {
	    phi_x = phi_xr;
	}
    }
    if (a_term * phi_yr > 0) {
	phi_y = (a_term * phi_yl > 0) ? phi_yl : 0;
    } else { // a_term * phi_yr < 0
	if (a_term * phi_yl > 0) {
	    phi_y = (fabs(phi_yl) > fabs(phi_yr)) ? phi_yl : phi_yr;
	} else {
	    phi_y = phi_yr;
	}
    }
}

//===========================================================================
void godunov_phi_3D(double a_term, 
		    double phi_xl, double phi_xr,
		    double phi_yl, double phi_yr,
		    double phi_zl, double phi_zr,
		    double& phi_x,
		    double& phi_y,
		    double& phi_z)
//===========================================================================
{
    // using Godunov to choose correct phi_x, phi_y and phi_z
    if (a_term * phi_xr > 0) {
	phi_x = (a_term * phi_xl > 0) ? phi_xl : 0;
    } else { // a_term * phi_xr < 0
	if (a_term * phi_xl > 0) {
	    phi_x = (fabs(phi_xl) > fabs(phi_xr)) ? phi_xl : phi_xr;
	} else {
	    phi_x = phi_xr;
	}
    }
    if (a_term * phi_yr > 0) {
	phi_y = (a_term * phi_yl > 0) ? phi_yl : 0;
    } else { // a_term * phi_yr < 0
	if (a_term * phi_yl > 0) {
	    phi_y = (fabs(phi_yl) > fabs(phi_yr)) ? phi_yl : phi_yr;
	} else {
	    phi_y = phi_yr;
	}
    }
    if (a_term * phi_zr > 0) {
	phi_z = (a_term * phi_zl > 0) ? phi_zl : 0;
    } else { // a_term * phi_zr < 0
	if (a_term * phi_zl > 0) {
	    phi_z = (fabs(phi_zl) > fabs(phi_zr)) ? phi_zl : phi_zr;
	} else {
	    phi_z = phi_zr;
	}
    }
}

}; // end anonymous namespace 

