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
// File: EigValComp3x3.C                                                     
//                                                                           
// Created: Wed May  3 14:55:25 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: EigValComp3x3.C,v 1.7 2006/11/22 15:23:00 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements EigValComp3x3.h.
//                                                                           
//===========================================================================

#include <assert.h>
#include <math.h>
#include "EigValComp3x3.h"
#include <iostream>
#include <limits>
using namespace std;

namespace { // anonymous namespace
    const double EPS = numeric_limits<double>::epsilon();
    const double PI = 3.1415926535897932384;
    const double THIRD = double(1)/3;
    const double TW7TH = double(1)/27;

    inline void xprod(const double* const u, const double* const v, double* res)
    {
	res[0] = u[1] * v[2] - u[2] * v[1];
	res[1] = u[2] * v[0] - u[0] * v[2];
	res[2] = u[0] * v[1] - u[1] * v[0];
    }

    inline double norm2(const double* v) {return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];}

    // do not call for a zero vector
    inline void normalize(double* v) 
    {
	double invnorm = double(1) / sqrt(norm2(v));
	for (int i = 0; i < 3; ++i) {
	    v[i] *= invnorm;
	}
    }

    inline void find_non_collinear(const double* const v, double* res)
    {
	// heuristic for finding a unit vector that we know is not parallel to v
	int min_ix = fabs(v[0]) < fabs(v[1]) ? 0 : 1;
	min_ix = fabs(v[min_ix]) < fabs(v[2]) ? min_ix : 2;

	for (int i = 0; i < 3; ++i) {
	    res[i] = (i == min_ix) ? 1 : 0;
	}

    }

    // only call if dimension of kernel is known to be 1 (not 0, 2 or 3)
    void find_kernel(const double* const col1, 
		     const double* const col2, 
		     const double* const col3, 
		     double* ker);

inline void givens(const double a, const double b, double& c, double& s) 
{
    if (b == 0) {
	c = 1;
	s = 0;
    } else {
	if (fabs(b) > fabs(a)) {
	    const double tau = -a / b;
	    s = 1 / sqrt(1 + tau * tau);
	    c = tau * s;
	} else {
	    const double tau = -b / a;
	    c = 1 / sqrt(1 + tau * tau);
	    s = tau * c;
	}
    }
}

void apply_givens(const double c, const double s, 
		  double& alpha1, double& alpha2, double& alpha3,
		  double& beta, double& delta, double& corner)
{
    const double tmp1 = c*alpha1 - s*beta;
    const double tmp2 = c*beta - s*alpha2;
    const double tmp3 = s*alpha1 + c*beta;
    const double tmp4 = s*beta + c*alpha2;
    const double ctemp = corner;
    
    corner = -s * delta + c * ctemp;
    delta = c * delta + s * ctemp;
    
    alpha1 = c * tmp1 - s * tmp2;
    alpha2 = s * tmp3 + c * tmp4;
    beta = s * tmp1 + c * tmp2;
}

void apply_givens(const double c, const double s,
		  double& alpha1, double& alpha2, double& beta)
{
    const double tmp1 = c*alpha1 - s*beta;
    const double tmp2 = c*beta - s*alpha2;
    const double tmp3 = s*alpha1 + c*beta;
    const double tmp4 = s*beta + c*alpha2;

    alpha1 = c * tmp1 - s * tmp2;
    alpha2 = s * tmp3 + c * tmp4;
    beta = s * tmp1 + c * tmp2;
}

void post_mul_givens(const double c, const double s, double* u, double* v)
{
    const double u0 = u[0];
    const double u1 = u[1];
    const double u2 = u[2];
    const double v0 = v[0];
    const double v1 = v[1];
    const double v2 = v[2];

    u[0] = c * u0 - s * v0;
    u[1] = c * u1 - s * v1;
    u[2] = c * u2 - s * v2;
    v[0] = s * u0 + c * v0;
    v[1] = s * u1 + c * v1;
    v[2] = s * u2 + c * v2;
}

void construct_eigvecs(double alpha1, double alpha2, double alpha3,
		       double beta, double gamma, double delta,
		       double& lambda1, 
		       double& lambda2, 
		       double& lambda3,
		       double* v1, double* v2, double* v3);

}; // end anonymous namespace 

namespace lsseg {

//===========================================================================
void analytic_eigvals(double alpha1, double alpha2, double alpha3,
		      double beta, double gamma, double delta,
		      double& lambda1, double& lambda2, double& lambda3)
//===========================================================================
{
    // computing coefficients of characteristic polynomial, a_2, a_1 and a_0
    const double bgd = beta * gamma * delta;
    const double beta2 = beta * beta;
    const double gamma2 = gamma * gamma;
    const double delta2 = delta * delta;
   
    const double a_2 = -(alpha1 + alpha2 + alpha3);  // coef. of x^2
    const double a_1 = 
	-1*(beta2 + gamma2 + delta2 - alpha1*alpha2 - alpha2*alpha3 - alpha1*alpha3);  // coef. of x
    const double a_0 = 
	-1*(alpha1*alpha2*alpha3 - alpha1*delta2 - alpha2*gamma2 - alpha3*beta2 + 2*bgd);  // coef. of 1

    
    //cout << a_0 << " " << a_1 << " " <<a_2 << endl;
    // computing the three real roots.  Formula is taken from 
    // http://mathworld.wolfram.com/CubicFormula.html

    const double a_2_sq = a_2 * a_2;
    const double a_2_cub = a_2_sq * a_2;

    const double p = a_1 - (a_2_sq) * THIRD;

    //cout << "p: " << p << endl;
    const double q = -a_0 + (9 * a_1 * a_2 - 2 * a_2_cub) * TW7TH;

    const double fac = fabs(p) * THIRD;
    
    double C = (fac > 0) ? 0.5 * q * pow(double(1)/fac, 1.5) : 1;
    //cout << "C: " << C << endl;
    
    // according to theory, |C| should be less or equal to one (since our symmetric matrix
    // guarantee us no complex roots).  To make sure numerical
    // rounding errors do not cause pain in the acos function below, we assure this here:

    const double acos_C = (fabs(C) <= 1) ? acos(C) : C > 0 ? 0 : PI;
    const double y_1 = cos (THIRD * acos_C);
    const double y_2 = cos (THIRD * (acos_C + 2 * PI));
    const double y_3 = cos (THIRD * (acos_C + 4 * PI));

    //cout << "y1, y2, y3: " << y_1 << " " << y_2 << " " << y_3 << endl;

    const double tmp1 = 2 * sqrt(fac);
    const double tmp2 = THIRD * a_2;

    lambda1 = tmp1 * y_1 - tmp2;
    lambda2 = tmp1 * y_2 - tmp2;
    lambda3 = tmp1 * y_3 - tmp2;
}

//===========================================================================
void analytic_eigsys(double alpha1, double alpha2, double alpha3,
		     double beta, double gamma, double delta,
		     double& lambda1, double& lambda2, double& lambda3,
		     double* v1, double* v2, double* v3)
//===========================================================================
{
    // computing eigenvalues (analytically)
    analytic_eigvals(alpha1, alpha2, alpha3, beta, gamma, delta, lambda1, lambda2, lambda3);
    construct_eigvecs(alpha1, alpha2, alpha3, beta, gamma, delta, lambda1, lambda2, lambda3,
		      v1, v2, v3);
}

//===========================================================================
void numeric_eigsys(double alpha1, double alpha2, double alpha3,
		    double beta, double gamma, double delta,
		    double& lambda1, double& lambda2, double& lambda3,
		    double* u, double* v, double* w)
//===========================================================================
{
    const bool compute_eigvecs = u && v && w;

    // reduce to tridiagonal form ("get rid of gamma") by using Householder
    bool apply_tridiagonalization = fabs(gamma) > EPS * (fabs(alpha1) + fabs(alpha2) + fabs(alpha3)) * THIRD;
    double root_sigma, u_1, u_2, H_inv,p_1, p_2, p_3, K, q_1, q_2, q_3;
    if (apply_tridiagonalization) {
	root_sigma = sqrt(gamma*gamma + delta*delta);
	u_1 = gamma; // @@
	u_2 = (delta > 0) ? delta + root_sigma : delta - root_sigma;
	// u_3 is 0 and we do not take it into consideration
	H_inv = double(2) / (u_1 * u_1 + u_2 * u_2);
	p_1 = H_inv * (alpha1 * u_1 + beta * u_2);
	p_2 = H_inv * (beta * u_1 + alpha2 * u_2);
	p_3 = H_inv * (gamma * u_1 + delta * u_2);
	K = 0.5 * H_inv * (u_1 * p_1 + u_2 * p_2);
	q_1 = p_1 - K * u_1;
	q_2 = p_2 - K * u_2;
	q_3 = p_3;
	
	// A_tridiag = A - q*u' - u*q'
	alpha1 -= 2 * (q_1 * u_1);
	alpha2 -= 2 * (q_2 * u_2);
	// alpha3 not modified since u_3 = 0
	beta -= (q_1 * u_2 + u_1 * q_2);
	gamma -= q_3 * u_1;  // (u_3 is 0, so omit term q_1 * u_3...
	delta -= q_3 * u_2;  // idem
	
	//cout << alpha1 << " " << alpha2 << " " << alpha3 << endl;
	//cout << beta << " " << gamma << " " << delta << endl;
    } else {
	gamma = 0; // from now, we consider gamma identically zero
    }

    double corner = 0;
    bool beta_zero = fabs(beta) <= EPS * (fabs(alpha1) + fabs(alpha2));
    bool delta_zero = fabs(delta) <= EPS * fabs(alpha3);
    double d, sign_d, mu, c1, s1, c2, s2;
    if (compute_eigvecs) {
	u[0] = 1; u[1] = 0; u[2] = 0;
	v[0] = 0; v[1] = 1; v[2] = 0;
	w[0] = 0; w[1] = 0; w[2] = 1;
    }

    while(true) {
	if (!beta_zero && !delta_zero) { // neither beta nor delta are zero

	    d = (alpha2 - alpha3) * 0.5;
	    sign_d = d > 0 ? 1 : -1;
	    mu = alpha3 - delta * delta / (d + sign_d * sqrt(d*d + delta*delta));
	    
	    givens(alpha1-mu, beta, c1, s1);
	    apply_givens(c1,s1,alpha1, alpha2, alpha3, beta, delta, corner);
	    beta_zero = fabs(beta) <= EPS * (fabs(alpha1) + fabs(alpha2));
	    
	    givens(beta, corner, c2, s2);
	    apply_givens(c2, s2, alpha2, alpha3, alpha1, delta, corner, beta);
	    delta_zero = fabs(delta) <= EPS * (fabs(alpha3));

	    if (compute_eigvecs) {
		post_mul_givens(c1, s1, u, v);
		post_mul_givens(c2, s2, v, w);
	    }

	} else if (!beta_zero) { // delta is zero, but not beta

	    d = (alpha1 - alpha2) * 0.5;
	    sign_d = d > 0 ? 1 : -1;
	    mu = alpha2 - beta * beta / (d + sign_d * sqrt(d*d + beta*beta));
	    givens(alpha1 - mu, beta, c1, s1);
	    apply_givens(c1, s1, alpha1, alpha2, beta);
	    beta_zero = fabs(beta) <= EPS * (fabs(alpha1) + fabs(alpha2));
	    if (compute_eigvecs) {
		post_mul_givens(c1, s1, u, v);
	    }
	    
	} else if (!delta_zero) { // beta is zero, but not delta

	    d = (alpha2 - alpha3) * 0.5;
	    sign_d = d > 0 ? 1 : -1;
	    mu = alpha3 - delta * delta / (d + sign_d * sqrt(d*d + delta*delta));
	    
	    givens(alpha2 - mu, delta, c2, s2);
	    apply_givens(c2, s2, alpha2, alpha3, delta);
	    delta_zero = fabs(delta) <= EPS * (fabs(alpha3));
	    if (compute_eigvecs) {
		post_mul_givens(c2, s2, v, w);
	    }
	} else { // both delta and beta are zero
	    break;
	}
    }
    lambda1 = alpha1;
    lambda2 = alpha2;
    lambda3 = alpha3;

    // transforming eigenvectors for the tridiagonal matrix to the eigenvectors
    // for the original matrix
    if (compute_eigvecs && apply_tridiagonalization) {
	const double fac_u = H_inv * (u[0] * u_1 + u[1] * u_2);
	u[0] -= fac_u * u_1;
	u[1] -= fac_u * u_2;

	const double fac_v = H_inv * (v[0] * u_1 + v[1] * u_2);
	v[0] -= fac_v * u_1;
	v[1] -= fac_v * u_2;

	const double fac_w = H_inv * (w[0] * u_1 + w[1] * u_2);
	w[0] -= fac_w * u_1;
	w[1] -= fac_w * u_2;
    }
}

}; // end namespace lsseg

namespace {

double tempvec1[3];
double tempvec2[3];
double tempvec3[3];

// only call if dimension of kernel is known to be 1 (not 0, 2 or 3)
void find_kernel(const double* const col1, 
		 const double* const col2, 
		 const double* const col3, 
		 double* ker)
{
    xprod(col1, col2, tempvec1);
    xprod(col2, col3, tempvec2);
    xprod(col1, col3, tempvec3);

    double n1 = norm2(tempvec1);
    double n2 = norm2(tempvec2);
    double n3 = norm2(tempvec3);

    if (n1 > n2) {
	if (n1 > n3) {
	    // n1 is biggest
	    n1 = double(1) / sqrt(n1);
	    for (int i = 0; i < 3; ker[i] = tempvec1[i++] * n1);
	} else {
	    // n3 is biggest
	    n3 = double(1) / sqrt(n3);
	    for (int i = 0; i < 3; ker[i] = tempvec3[i++] * n3);
	}
    } else {
	if (n2 > n3) {
	    // n2 is biggest
	    n2 = double(1) / sqrt(n2);
	    for (int i = 0; i < 3; ker[i] = tempvec2[i++] * n2);
	} else {
	    // n3 is biggest
	    n3 = double(1) / sqrt(n3);
	    for (int i = 0; i < 3; ker[i] = tempvec3[i++] * n3);
	}
    }
}

//===========================================================================
void construct_eigvecs(double alpha1, double alpha2, double alpha3,
		       double beta, double gamma, double delta,
		       double& lambda1, 
		       double& lambda2, 
		       double& lambda3,
		       double* v1, double* v2, double* v3)
//===========================================================================
{
    // sorting eigenvals
    if (fabs(lambda3) > fabs(lambda2)) swap(lambda3, lambda2);
    if (fabs(lambda2) > fabs(lambda1)) {
	swap(lambda2, lambda1);
	if (fabs(lambda3) > fabs(lambda2)) swap(lambda3, lambda2);
    }
    // we have now fabs(lambda1) >= fabs(lambda2) >= fabs(lambda3)

    // constructing eigenvectors
    double col1[3], col2[3], col3[3];
    double tmp[3];
    const double numeric_threshold =  EPS * (fabs(alpha1) + fabs(alpha2) + fabs(alpha3));
    
    //if (fabs(lambda1 - lambda2) > EPS * (fabs(lambda1) + fabs(lambda2)) * 0.5) {
    if (fabs(lambda1 - lambda2) > numeric_threshold) {

	col1[0] = alpha1-lambda1; col1[1] = beta;            col1[2] = gamma;
	col2[0] = beta;           col2[1] = alpha2-lambda1 ; col2[2] = delta;
	col3[0] = gamma;          col3[1] = delta;           col3[2] = alpha3-lambda1;
	find_kernel(col1, col2, col3, v1);

	if (fabs(lambda2 - lambda3) > numeric_threshold) {

	    // lambda1 > lambda2 > lambda3
	    col1[0] = alpha1-lambda2; col2[1] = alpha2-lambda2 ; col3[2] = alpha3-lambda2;
	    find_kernel(col1, col2, col3, v2);

	    col1[0] = alpha1-lambda3; col2[1] = alpha2-lambda3 ; col3[2] = alpha3-lambda3;
	    find_kernel(col1, col2, col3, v3);

	} else {
	    //lambda1 > lambda2 == lambda3

	    find_non_collinear(v1, tmp);
	    
	    // finding vector of space perpendicular to v1
	    xprod(v1, tmp, v2);
	    normalize(v2);
	    xprod(v1, v2, v3);
	    normalize(v3);
	}
    } else if (fabs(lambda1 - lambda3) > numeric_threshold) {
	// lambda1 == lambda2 > lambda3

	col1[0] = alpha1-lambda3; col1[1] = beta;            col1[2] = gamma;
	col2[0] = beta;           col2[1] = alpha2-lambda3 ; col2[2] = delta;
	col3[0] = gamma;          col3[1] = delta;           col3[2] = alpha3-lambda3;
	find_kernel(col1, col2, col3, v3);

	find_non_collinear(v3, tmp);
	
	// finding vector of space perpendicular to v3
	xprod(v3, tmp, v1);
	normalize(v1);
	xprod(v3, v1, v2);
	normalize(v2);

    } else {
	// lambda1 == lambda2 == lambda3
	v1[0] = 1; v1[1] = 0; v1[2] = 0;
	v2[0] = 0; v2[1] = 1; v2[2] = 0;
	v3[0] = 0; v3[1] = 0; v3[2] = 1;
    }
}

};
