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
// File: EigValComp3x3.h                                                     
//                                                                           
// Created: Wed May  3 14:51:38 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: EigValComp3x3.h,v 1.6 2006/11/13 02:29:24 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief This file contains functions for computing the eigensystem of a 
/// 3x3 symmetric matrix.  
//                                                                           
//===========================================================================

#ifndef _EIGVALCOMP3X3_H
#define _EIGVALCOMP3X3_H

namespace lsseg {

/// \brief Compute the eigenvalues of a 3x3 symmetric matrix by analytically
/// solving the corresponding characteristic polynomial.
///
/// The terms in the symmetric matrix are given like this:
/// \f[ \left(\begin{array}{ccc}  \alpha_1 & \beta & \gamma \\  \beta & \alpha_2 & \delta \\ \gamma & \delta & \alpha_3 \end{array} \right)\f]
/// \note This function might not be sufficiently robust in nearly-degenerate cases.
/// It is better to use the numerical algorithm numeric_eigsys().  This function
/// is mostly here for theoretical (and historical) reasons.
/// \param[in] alpha1 \f$\alpha_1\f$ in the matrix above (first diagonal term)
/// \param[in] alpha2 \f$\alpha_2\f$ in the matrix above (second diagonal term)
/// \param[in] alpha3 \f$\alpha_3\f$ in the matrix above (third diagonal term)
/// \param[in] beta   \f$\beta\f$ in the matrix above
/// \param[in] gamma \f$\gamma\f$ in the matrix above
/// \param[in] delta \f$\delta\f$ in the matrix above
/// \param[out] lambda1 the first of the computed eigenvalues (in arbitrary order)
/// \param[out] lambda2 the second of the computed eigenvalues (in arbitrary order)
/// \param[out] lambda3 the third of the computed eigenvalues (in arbitrary order)
void analytic_eigvals(double alpha1, double alpha2, double alpha3,
		      double beta, double gamma, double delta,
		      double& lambda1, double& lambda2, double& lambda3);

/// \brief Compute the eigenvalues and eigenvectors of a 3x3 symmetric matrix by analytically
/// solving the corresponding characteristic polynomial to find the eigenvalues, 
/// and then explicitly constructing the corresponding eigenvectors.
///
/// The terms in the symmetric matrix are given like this:
/// \f[ \left(\begin{array}{ccc}  \alpha_1 & \beta & \gamma \\  \beta & \alpha_2 & \delta \\ \gamma & \delta & \alpha_3 \end{array} \right)\f]
/// \note This function might not be sufficiently robust in nearly-degenerate cases.
/// It is better to use the numerical algorithm numeric_eigsys().  This function
/// is mostly here for theoretical (and historical) reasons.
/// \param[in] alpha1 \f$\alpha_1\f$ in the matrix above (first diagonal term)
/// \param[in] alpha2 \f$\alpha_2\f$ in the matrix above (second diagonal term)
/// \param[in] alpha3 \f$\alpha_3\f$ in the matrix above (third diagonal term)
/// \param[in] beta   \f$\beta\f$ in the matrix above
/// \param[in] gamma \f$\gamma\f$ in the matrix above
/// \param[in] delta \f$\delta\f$ in the matrix above
/// \param[out] lambda1 the first of the computed eigenvalues (the largest one)
/// \param[out] lambda2 the second of the computed eigenvalues (the middle one)
/// \param[out] lambda3 the third of the computed eigenvalues (the smallest one)
/// \param[out] v1 pointer to an array where the eigenvector corresponding to \c lambda1
///                has been stored.  (The array contains 3 elements).
/// \param[out] v2 pointer to an array where the eigenvector corresponding to \c lambda2
///                has been stored.  (The array contains 3 elements).
/// \param[out] v3 pointer to an array where the eigenvector corresponding to \c lambda3
///                has been stored.  (The array contains 3 elements).
void analytic_eigsys(double alpha1, double alpha2, double alpha3,
		     double beta, double gamma, double delta,
		     double& lambda1, double& lambda2, double& lambda3,
		     double* v1, double* v2, double* v3);

/// \brief Compute the eigenvalues and eigenvectors of a 3x3 symmetric matrix by 
/// means of an iterative numerical algorithm using Householder transformations and
/// Givens rotations.
///
/// The terms in the symmetric matrix are given like this:
/// \f[ \left(\begin{array}{ccc}  \alpha_1 & \beta & \gamma \\  \beta & \alpha_2 & \delta \\ \gamma & \delta & \alpha_3 \end{array} \right)\f]
/// \note This function is more numerically robust than analytic_eigsys() and should be prefered to 
/// that one.
/// \param[in] alpha1 \f$\alpha_1\f$ in the matrix above (first diagonal term)
/// \param[in] alpha2 \f$\alpha_2\f$ in the matrix above (second diagonal term)
/// \param[in] alpha3 \f$\alpha_3\f$ in the matrix above (third diagonal term)
/// \param[in] beta   \f$\beta\f$ in the matrix above
/// \param[in] gamma \f$\gamma\f$ in the matrix above
/// \param[in] delta \f$\delta\f$ in the matrix above
/// \param[out] lambda1 the first of the computed eigenvalues (arbitrary order)
/// \param[out] lambda2 the second of the computed eigenvalues (arbitrary order)
/// \param[out] lambda3 the third of the computed eigenvalues (arbitrary order)
/// \param[out] v1 pointer to an array where the eigenvector corresponding to \c lambda1
///                has been stored.  (The array contains 3 elements).
/// \param[out] v2 pointer to an array where the eigenvector corresponding to \c lambda2
///                has been stored.  (The array contains 3 elements).
/// \param[out] v3 pointer to an array where the eigenvector corresponding to \c lambda3
///                has been stored.  (The array contains 3 elements).
void numeric_eigsys(double alpha1, double alpha2, double alpha3,
		    double beta, double gamma, double delta,
		    double& lambda1, double& lambda2, double& lambda3,
		    double* v1, double* v2, double* v3);


}; // end namespace

#endif // _EIGVALCOMP3X3_H

