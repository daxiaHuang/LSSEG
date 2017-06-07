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
// File: Filters.C                                                           
//                                                                           
// Created: Wed Feb 22 13:34:51 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: Filters.C,v 1.12 2006/11/25 20:08:26 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements Filters.h.
//                                                                           
//===========================================================================

#include <cmath>
#include <stdexcept>
#include <assert.h>
#include "EigValComp3x3.h"
#include "cimg_dependent.h"
#include "Filters.h"
#include "DiscreteApproximations.h"

using namespace std;

namespace {

void diffusion_1D(const double* const u_old, 
		  const double* const g,
		  double * u_new, 
		  double dt,
		  int dim);
void tridiag(const double* const a, 
	     const double* const b, 
	     const double* const c,
	     const double* const r, 
	     double* const sol, int n);

void accumulate_scale(const lsseg::Image<double>& img1,
		      const lsseg::Image<double>& img2,
		      lsseg::Image<double>& inv_scale, 
		      lsseg::Image<double>& time_span,
		      double dt);

inline double p_func(double u_squared, double p);

}; // end anonymous namespace 


namespace lsseg {

//===========================================================================
void compute_structure_tensor_2D(const Image<double>& img, 
				 Image<double>& G, 
				 bool square_root)
//===========================================================================
{
    assert(img.dimz() == 1);
    G.resize(img.dimx(), img.dimy(), 1, 3);
    G = 0;
    const int X = img.dimx();
    const int Y = img.dimy();
    for (int ch = 0; ch != img.numChannels(); ++ch) {
	for (int y = 0; y < Y; ++y) {
	    const int yp = (y > 0) ? y-1 : y;
	    const int yn = (y < Y-1) ? y+1 : y;
	    const double dy_inv = (yn-yp == 2) ? 0.5 : 1;

	    for (int x = 0; x < X; ++x) {
		const int xp = (x > 0) ? x-1 : x;
		const int xn = (x < X-1) ? x+1 : x;
		const double dx_inv = (xn-xp == 2) ? 0.5 : 1;

		const double Ipc = img(xp, y, 0, ch);
		const double Inc = img(xn, y, 0, ch);
		const double Icp = img(x, yp, 0, ch);
		const double Icn = img(x, yn, 0, ch);
		const double ix = dx_inv * (Inc - Ipc);
		const double iy = dy_inv * (Icn - Icp);
		
		double norm_inv = 1;
		if (square_root) {
		    const double norm = sqrt(ix * ix + iy * iy);
		    norm_inv = (norm > 0) ? double(1) / norm : 1;
		}
		G(x, y, 0, 0) += ix * ix * norm_inv;
		G(x, y, 0, 1) += ix * iy * norm_inv;
		G(x, y, 0, 2) += iy * iy * norm_inv;
	    }
	}
    }
}

//===========================================================================
void compute_structure_tensor_3D(const Image<double>& img,
				 Image<double>& G,
				 bool square_root)
//===========================================================================
{
    G.resize(img.dimx(), img.dimy(), img.dimz(), 6);
    G = 0;
    const int X = img.dimx();
    const int Y = img.dimy();
    const int Z = img.dimz();
    for (int ch = 0; ch != img.numChannels(); ++ch) {
	for (int z = 0; z < Z; ++z) {
	    const int zp = (z > 0) ? z-1 : z;
	    const int zn = (z < Z-1) ? z+1 : z;
	    const double dz_inv = (zn - zp == 2) ? 0.5 : 1;

	    for (int y = 0; y < Y; ++y) {
		const int yp = (y > 0) ? y-1 : y;
		const int yn = (y < Y-1) ? y+1 : y;
		const double dy_inv = (yn - yp == 2) ? 0.5 : 1;
		
		for (int x = 0; x < X; ++x) {
		    const int xp = (x > 0) ? x-1 : x;
		    const int xn = (x < X-1) ? x+1 : x;
		    const double dx_inv = (xn-xp == 2) ? 0.5 : 1;
		    
		    const double Ipcc = img(xp, y, z, ch);
		    const double Incc = img(xn, y, z, ch);
		    const double Icpc = img(x, yp, z, ch);
		    const double Icnc = img(x, yn, z, ch);
		    const double Iccp = img(x, y, zp, ch);
		    const double Iccn = img(x, y, zn, ch);
		    const double ix = dx_inv * (Incc - Ipcc);
		    const double iy = dy_inv * (Icnc - Icpc);
		    const double iz = dz_inv * (Iccn - Iccp);

		    double norm_inv = 1;
		    if(square_root) {
			const double norm = sqrt(ix * ix + iy * iy + iz * iz);
			norm_inv = (norm > 0) ? double(1)/ norm : 1;
		    }
		    G(x, y, z, 0) = ix * ix * norm_inv;  // x x
		    G(x, y, z, 1) = ix * iy * norm_inv;  // x y
		    G(x, y, z, 2) = iy * iy * norm_inv;  // y y
		    G(x, y, z, 3) = ix * iz * norm_inv;  // x z
		    G(x, y, z, 4) = iy * iz * norm_inv;  // y z
		    G(x, y, z, 5) = iz * iz * norm_inv;  // z z
		}
	    }
	}
    }
}


//===========================================================================
void compute_scale_factor_2D(const Image<double>& img,
			     Image<double>& scale_accum,
			     Image<double>& time_accum,
			     double dt,
			     double T)
//===========================================================================
{
    assert(img.numChannels() == 1);
    assert(img.dimz() == 1);
    // applying isotropic, nonlinear smoothing on image, and measuring the 
    // rate of change
    
    Image<double> tmp1(img); // copy content
    Image<double> tmp2(img, false); // do not copy content
    
    scale_accum.resize(img);
    scale_accum = 0;
    time_accum.resize(img);
    time_accum = 0;
    for (double t = 0; t < T; t+= dt) {
	nonlinear_gauss_filter_2D(tmp1, tmp2, dt, 0, 1);
	accumulate_scale(tmp1, tmp2, scale_accum, time_accum, dt);
	tmp1.swap(tmp2);
    }
}

//===========================================================================
void compute_scale_factor_3D(const Image<double>& img,
			     Image<double>& scale_accum,
			     Image<double>& time_accum,
			     double dt,
			     double T)
//===========================================================================
{
    assert(img.numChannels() == 1);
    // applying isotropic, nonlinear smoothing on image, and measuring the 
    // rate of change
    
    Image<double> tmp1(img); // copy content
    Image<double> tmp2(img, false); // do not copy content
    
    scale_accum.resize(img);
    scale_accum = 0;
    time_accum.resize(img);
    time_accum = 0;
    for (double t = 0; t < T; t+= dt) {
	nonlinear_gauss_filter_3D(tmp1, tmp2, dt, 0, 1);
	accumulate_scale(tmp1, tmp2, scale_accum, time_accum, dt);
	tmp1.swap(tmp2);
    }
}

//===========================================================================
void nonlinear_gauss_filter_2D(const Image<double>& img_old,  
			       Image<double>& img_new,
			       double dt,
			       double sigma,
			       double p)
//===========================================================================
{
    assert(img_old.size_compatible(img_new));
    
    Image<double> img_cpy(img_old);

    if (sigma > 0) {
	blur_image(img_cpy, sigma);
    }

    // computing the squared gradient matrix
    Image<double> g; 
    compute_squared_gradient_sum_2D(img_cpy, g);
    for (double* d = g.begin(); d != g.end(); ++d) {
	*d = p_func(*d, p);
    }

    // temporary result image
    Image<double> tmp_res;

    int PERMXY[] = {1, 0, 2, 3};
    for (int run = 0; run < 2; ++run) {
	int len = img_cpy.dimx();
	double* gp = g.begin();
	double* trg = (run == 0) ? img_new.begin() : tmp_res.begin();
	for (double* cp = img_cpy.begin(); 
	     cp != img_cpy.end(); cp += len, gp += len, trg += len) {
	    if (gp == g.end()) { // reset pointer to g data
		gp = g.begin();
	    }
	    diffusion_1D(cp, gp, trg, dt, len);
	}
	img_cpy.permute(PERMXY); // permuting x, y and z component
	g.permute(PERMXY);
	img_new.permute(PERMXY);
	if (run == 0) {
	    tmp_res.resize(img_new.dimx(), 
			   img_new.dimy(), 
			   img_new.dimz(), 
			   img_new.numChannels());
	} else {
	    tmp_res.permute(PERMXY);
	    img_new += tmp_res;
	}
    }
    img_new /= 2;
}


//===========================================================================
void nonlinear_gauss_filter_3D(const Image<double>& img_old,
			       Image<double>& img_new,
			       double dt,
			       double sigma,
			       double p)
//===========================================================================
{
    const int PERM[] = {1, 2, 0, 3};
    assert(img_old.size_compatible(img_new));
    
    Image<double> img_cpy(img_old);

    if (sigma > 0) {
	blur_image(img_cpy, sigma);
    }

    // computing the squared gradient matrix
    Image<double> g; 
    compute_squared_gradient_sum_3D(img_cpy, g);
    for (double* d = g.begin(); d != g.end(); ++d) {
	*d = p_func(*d, p);
    }

    // temporary result image
    Image<double> tmp_res;

    // treating image lines in first run, columns in the second loop and pillars in
    // third run
    for (int run = 0; run < 3; ++run) {
	int len = img_cpy.dimx();
	double* gp = g.begin();
	double* trg = (run == 0) ? img_new.begin() : tmp_res.begin();
	for (double* cp = img_cpy.begin(); 
	     cp != img_cpy.end(); cp += len, gp += len, trg += len) {
	    if (gp == g.end()) { // reset pointer to g data
		gp = g.begin();
	    }
	    diffusion_1D(cp, gp, trg, dt, len);
	}
	img_cpy.permute(PERM); // permuting x, y and z component
	g.permute(PERM);
	img_new.permute(PERM);
	if (run == 0) {
	    tmp_res.resize(img_new.dimx(), 
			   img_new.dimy(), 
			   img_new.dimz(), 
			   img_new.numChannels());
	} else {
	    tmp_res.permute(PERM);
	    img_new += tmp_res;
	}
    }

    img_new /= 3;
}

//===========================================================================
void anisotropic_smoothing(const Image<double>& img_old, 
			   Image<double>& img_new, 
			   double dt, 
			   double sigma,
			   double p)
//===========================================================================
{
    assert(p >= 0 && sigma >= 0 && dt >= 0);
    img_new.resize(img_old);
    if (img_old.dimz() > 1) {
	nonlinear_gauss_filter_3D(img_old, img_new, dt, sigma, p);
    } else {
	nonlinear_gauss_filter_2D(img_old, img_new, dt, sigma, p);
    }
}

//===========================================================================
void compute_smoothing_geometry_3D(const Image<double>& G,
				   double p1,
				   double p2,
				   double p3,
				   Image<double>& T,
				   bool take_square_root)
//===========================================================================
{
    const double EPS = 1.0e-8; // @@ only used in assert below
    assert(G.numChannels() == 6);
    T.resize(G);
    const int X = G.dimx(), Y = G.dimy(), Z = G.dimz();
    if (take_square_root) {
	p1 *= 0.5;
	p2 *= 0.5;
	p3 *= 0.5;
    }
    vector<double> v1(3), v2(3),v3(3);
    double lambda1, lambda2, lambda3;
    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		// finding eigenvalues and vectors of current structure tensor
		const double alpha1 = G(x, y, z, 0);  // dx*dx
		const double alpha2 = G(x, y, z, 2);  // dy*dy
		const double alpha3 = G(x, y, z, 5);  // dz*dz
		const double beta = G(x, y, z, 1);    // dx*dy
		const double gamma = G(x, y, z, 3);   // dx*dz
		const double delta = G(x, y, z, 4);   // dy*dz

		numeric_eigsys(alpha1, alpha2, alpha3, 
			       beta, gamma, delta, 
			       lambda1, lambda2, lambda3, 
			       &v1[0], &v2[0], &v3[0]);

		// sorting eigen-elements

		if (fabs(lambda3) > fabs(lambda2)) {
		    swap(lambda3, lambda2);
		    swap(v3, v2);
		}
		if (fabs(lambda2) > fabs(lambda1)) {
		    swap(lambda2, lambda1);
		    swap(v2, v1);
		    if (fabs(lambda3) > fabs(lambda2)) {
			swap(lambda3, lambda2);
			swap(v3, v2);
		    }
		}
 		assert(fabs(lambda1) >= fabs(lambda2) && fabs(lambda2) >= fabs(lambda3));

 		// constructing T
		assert(lambda1 >= -EPS && lambda2 >= -EPS && lambda3 >= -EPS);
		const double tmp = 1 + lambda1 + lambda2 + lambda3;
		const double fac_3 = double(1) / pow(tmp, p1); // biggest @@
		const double fac_2 = double(1) / pow(tmp, p2); 
		const double fac_1 = double(1) / pow(tmp, p3); // smallepst@@

		assert(fac_3 >= fac_2 && fac_2 >= fac_1);

 		T(x, y, z, 0) = 
		    (fac_3 * v3[0] * v3[0]) + 
		    (fac_2 * v2[0] * v2[0]) + 
		    (fac_1 * v1[0] * v1[0]);
		T(x, y, z, 1) = 
		    (fac_3 * v3[0] * v3[1]) + 
		    (fac_2 * v2[0] * v2[1]) + 
		    (fac_1 * v1[0] * v1[1]);
		T(x, y, z, 2) = 
		    (fac_3 * v3[1] * v3[1]) + 
		    (fac_2 * v2[1] * v2[1]) + 
		    (fac_1 * v1[1] * v1[1]);
		T(x, y, z, 3) = 
		    (fac_3 * v3[0] * v3[2]) + 
		    (fac_2 * v2[0] * v2[2]) + 
		    (fac_1 * v1[0] * v1[2]);
		T(x, y, z, 4) = 
		    (fac_3 * v3[1] * v3[2]) + 
		    (fac_2 * v2[1] * v2[2]) + 
		    (fac_1 * v1[1] * v1[2]);
		T(x, y, z, 5) = 
		    (fac_3 * v3[2] * v3[2]) + 
		    (fac_2 * v2[2] * v2[2]) + 
		    (fac_1 * v1[2] * v1[2]);
	    }
	}
    }
}


//===========================================================================
void compute_smoothing_geometry_2D(const Image<double>& G,
				   double p1,
				   double p2,
				   Image<double>& T,
				   bool take_square_root)
//===========================================================================
{
    const double EPS = 1.0e-8;  // used in computing eigenvectors
    assert(G.numChannels() == 3 && G.dimz() == 1);
    T.resize(G);
    const int X = G.dimx();
    const int Y = G.dimy();

    if (take_square_root) {
	p1 *= 0.5;
	p2 *= 0.5;
    }
    double lm, lp, xm, xp, ym, yp; // lambda plus/minus, x and y components of theta plus/minus
    for (int y = 0; y < Y; ++y) {
	for (int x = 0; x < X; ++x) {
	    const double a = G(x, y, 0, 0);
	    const double b = G(x, y, 0, 1);
	    const double c = G(x, y, 0, 2);

	    // computing eigenvalues/vecs
	    if (fabs(b) < EPS) {
		// we consider the off-diagonal values to be zero
		if (a > c) { 
		    lp = a;  lm = c;
		    xp = 1;  xm = 0;
		    yp = 0;  ym = 1;
		} else {
		    lp = c;  lm = a;
		    xp = 0;  xm = 1;
		    yp = 1;  ym = 0;
		}
	    } else {
		const double amc = a - c;
		const double apc = a + c;
		const double b2 = b * b;
		const double discr = sqrt(amc*amc + 4 * b2); // root of discriminant
		assert(discr > 0); // should only happen if b=0, which is taken care of above

		lp = 0.5 * (apc + discr);  
		const double amlp = a - lp;
		double norm_fac_p = sqrt(b2 + amlp * amlp);
		assert(norm_fac_p > 0); // same reason as above
		norm_fac_p = double(1) / norm_fac_p;
		xp = -b * norm_fac_p;
		yp = amlp * norm_fac_p;

		lm = 0.5 * (apc - discr);
		const double amlm = a - lm;
		double norm_fac_m = sqrt(b2 + amlm * amlm);
		assert(norm_fac_m > 0); // same reason as above
		norm_fac_m = double(1) / norm_fac_m;
		xm = -b * norm_fac_m;
		ym = amlm * norm_fac_m;
	    }

	    // constructing T
	    const double tmp = 1 + lp + lm;
	    const double fac_m = double(1) / pow(tmp, p1);
	    const double fac_p = double(1) / pow(tmp, p2); 
	    T(x, y, 0, 0) = (fac_p * xp * xp) + (fac_m * xm * xm); // upper left value
	    T(x, y, 0, 1) = (fac_p * xp * yp) + (fac_m * xm * ym); // symmetric off-diagonal value
	    T(x, y, 0, 2) = (fac_p * yp * yp) + (fac_m * ym * ym); // lower right value	
 
	}
    }
}



}; // end namespace lsseg

namespace {

//===========================================================================
void accumulate_scale(const lsseg::Image<double>& img1, 
		      const lsseg::Image<double>& img2, 
 		      lsseg::Image<double>& inv_scale, 
 		      lsseg::Image<double>& time_span,
 		      double dt)
//===========================================================================
{
    const double EPS = 1.0e-5;
    assert(img1.size_compatible(img2));
    assert(img1.size_compatible(inv_scale));
    assert(img1.size_compatible(time_span));
    assert(img1.numChannels() == 1);

    const double* img1_it = img1.begin();
    const double* img2_it = img2.begin();
    double* scale_it = inv_scale.begin();
    double* time_it = time_span.begin();

    for ( ; img1_it != img1.end(); ++img1_it, ++img2_it, ++scale_it, ++time_it) {
	double norm = fabs(*img1_it - *img2_it);

	if (norm > EPS) {
	    *scale_it += norm * dt;
	    *time_it += 4 * dt;
	}
    }
}

//==============================================================================
double p_func(double u_squared, double p)
//==============================================================================
{
    static const double EPS_SQUARED = 1.0e-4; // regularising constant
    
    if (fabs(p-1) < EPS_SQUARED) {
	return double(1)/sqrt(u_squared + EPS_SQUARED);
    } else {
	return double(1) / pow(u_squared + EPS_SQUARED, 0.5 * p);
    }
}

//==============================================================================
void tridiag(const double* const a, 
	     const double* const b, 
	     const double* const c,
	     const double* const r, 
	     double* const sol, int n)
//==============================================================================
{
    double bet; // temporary variable
    vector<double> gam(n); // workspace
    if (b[0] == 0.0) throw runtime_error("Error 1 in tridiag.");
    sol[0] = r[0] / (bet = b[0]);
    for (int j = 1; j < n; ++j) {
	gam[j] = c[j-1] / bet;
	bet = b[j] - a[j] * gam[j];
	if (bet == 0.0) throw runtime_error("Error 2 in tridiag.");
	sol[j] = (r[j] - a[j] * sol[j-1]) / bet;
    }
    for (int j = (n-2); j >= 0; --j) {
	sol[j] -= gam[j+1] * sol[j+1];
    }
}

//===========================================================================
void diffusion_1D(const double* const u_old, 
		  const double* const g,
		  double * u_new, 
		  double dt,
		  int dim)
//===========================================================================
{
    if (dim == 1) {
	// no diffusion possible, just copy input data
	*u_new = *u_old;
	return;
    }
    static vector<double> a, b, c;
    a.resize(dim);
    b.resize(dim);
    c.resize(dim);
    
    // border cases
    a[0] = 0;
    c[0] = -dt * (g[0] + g[1]);
    b[0] = 1 - c[0];

    a[dim-1] = -dt * (g[dim-2] + g[dim-1]);
    c[dim-1] = 0;
    b[dim-1] = 1 - a[dim-1];

    // interior
    for (int i = 1; i < dim-1; ++i) {
	a[i] = c[i-1];
	c[i] = -dt * (g[i] + g[i+1]);// * 0.5;
	b[i] = 1 - (a[i] + c[i]);
    }

    // solving linear system
    tridiag(&a[0], &b[0], &c[0], u_old, u_new, dim);

}


}; // end anonymous namespace 




// //===========================================================================
// void compute_scale_factor(const Image& img,
// 			  Image& scale_accum,
// 			  Image& time_accum,
// 			  double dt,
// 			  double T)
// //===========================================================================
// {
//     assert(img.dim == 1);
//     // applying isotropic, nonlinear smoothing on image, and measuring the
//     // rate of change

//     Image tmp1(img); // copy content 
//     Image tmp2(img, false); // do not copy content
//     scale_accum.resize(img);
//     scale_accum.fill(0);
//     time_accum.resize(img);
//     time_accum.fill(0);
//     for (double t = 0; t < T; t+= dt) {
// 	cout << "Time: " << t << endl;
// 	nonlinear_gauss_filter(tmp1, tmp2, dt, 0, 1);
// 	accumulate_scale(tmp1, tmp2, scale_accum, time_accum, dt);
// 	tmp1.swap(tmp2);
//     }
// }



// //===========================================================================
// void compute_structure_tensor(const Image& img, Image& G, bool square_root)
// //===========================================================================
// {
//     assert(img.dim == 1); // if not, offset will not work
//     G.resize(img.width, img.height, 3);
//     fill(G.data.begin(), G.data.end(), 0);
    
//     double Ipc, Inc, Icp, Icn;
//     if (square_root) {
// 	int N = img.width;
// 	int offset = 0;
// 	for (int y = 0; y < img.height; ++y) {
// 	    for (int x = 0; x < img.width; ++x, ++offset) {
// 		Ipc = (x>0) ? img.data[offset-1] : 2 * img.data[offset] - img.data[offset+1];
// 		Inc = (x<img.width-1) ? img.data[offset+1] : 2*img.data[offset]-img.data[offset-1];
// 		Icp = (y>0) ? img.data[offset - N] : 2 * img.data[offset] - img.data[offset+N];
// 		Icn = (y<img.height-1) ? img.data[offset+N] : 2*img.data[offset]-img.data[offset-N];
		
// 		double ix = 0.5 * (Inc - Ipc);
// 		double iy = 0.5 * (Icn - Icp);
// 		double norm = sqrt(ix * ix + iy * iy);
// 		if (norm > 0) {
// 		    double norm_inv = double(1) / norm;
// 		    G.data[offset + 0 * img.size()] = ix * ix * norm_inv;
// 		    G.data[offset + 1 * img.size()] = ix * iy * norm_inv;
// 		    G.data[offset + 2 * img.size()] = iy * iy * norm_inv;
// 		} else {
// 		    G.data[offset + 0 * img.size()] = 0;
// 		    G.data[offset + 1 * img.size()] = 0;
// 		    G.data[offset + 2 * img.size()] = 0;
// 		}
// 	    }
// 	}
//     } else {
// 	int N = img.width;
// 	int offset = 0;
// 	for (int y = 0; y < img.height; ++y) {
// 	    for (int x = 0; x < img.width; ++x, ++offset) {
// 		Ipc = (x>0) ? img.data[offset-1] : 2 * img.data[offset] - img.data[offset+1];
// 		Inc = (x<img.width-1) ? img.data[offset+1] : 2*img.data[offset]-img.data[offset-1];
// 		Icp = (y>0) ? img.data[offset - N] : 2 * img.data[offset] - img.data[offset+N];
// 		Icn = (y<img.height-1) ? img.data[offset+N] : 2*img.data[offset]-img.data[offset-N];
		
// 		double ix = 0.5 * (Inc - Ipc);
// 		double iy = 0.5 * (Icn - Icp);
// 		G.data[offset + 0 * img.size()] = ix * ix;
// 		G.data[offset + 1 * img.size()] = ix * iy;
// 		G.data[offset + 2 * img.size()] = iy * iy;
// 	    }
// 	}
//     }
// }

