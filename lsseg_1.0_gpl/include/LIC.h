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
// File: LIC.h                                                               
//                                                                           
// Created: Mon Apr 24 17:21:25 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: LIC.h,v 1.7 2006/11/20 03:56:25 oan Exp $
//                                                                           
// Description: line integral convolution
/// \file
/// \brief Contains functionality for \ref anchor_LIC "line integral convolution".
///                                                                      
//===========================================================================

#ifndef _LIC_H
#define _LIC_H

#include "Image.h"

namespace lsseg {

//===========================================================================
/// Carry out \ref anchor_LIC "line integral convolution" of a 2D-image with 
/// a kernel function, along streamlines given by a user-specified vector field.
/// The result is an image smoothed along the streamlines, with a smoothing
/// defined by the kernel function.
/// For more on line integral convolution, refer to  \ref anchor_Cabral93 "this paper [Cabral93]",
/// or (more adapted for the \ref section_GreycStoration "curvature-preserving smoothing algorithm"): 
/// \ref anchor_Tschumperle06 "this paper [Tschumperle06]".
/// \param[in] src the Image that will participate in the
///                \ref anchor_LIC "line integral convolution"
/// \param[in] vec a vector field defining the streamlines along which convolution
///                will take place.  The vector field is specified as an Image with
///                three channels; the two first containing the x- and y-components
///                of a normalized vector field, the third containing the \em magnitudes
///                of the vectors in the field.
/// \param[out] target the Image obtained when convoluting 'src' with 'kernel_func' along the streamlines of 'vec' (the 
///                    result of the \ref anchor_LIC "line integral convolution".
/// \param[in] kernel_func pointer to the kernel function (a function from \f$\mathcal{R}\f$ to \f$\mathcal{R}\f$.  
///                        This is typically a symmetrical function with compact support, or with infinite support but
///                        with a value that goes off to zero relatively quickly.
/// \param[in] dl steplength used for integration purposes
/// \param[in] L total length along which to integrate (in each direction from zero).  The important part of 
///              the support of 'kernel_func' should fall within this range.
/// \note This function implements line integral convolution with \em fixed steplength.  
/// A \ref LIC_2D "version using \em adaptive steplengths" is also provided in this library, but it can 
/// \ref anchor_LIC_2D_TROUBLE "run into problems".
/// For that reason, it is recommended to use this function instead, until those problems are fixed.
void LIC_2D_FS(const Image<double>& src,
	       const Image<double>& vec,
	       Image<double>& target,
	       double (*kernel_func)(double),
	       const double dl = 0.8, // steplength
	       const double L = 10); // total length along which to integrate
//===========================================================================


//===========================================================================
/// Carry out \ref anchor_LIC "line integral convolution" of a 3D-image with 
/// a kernel function, along streamlines given by a user-specified vector field.
/// The result is an image smoothed along the streamlines, with a smoothing
/// defined by the kernel function.
/// For more on line integral convolution, refer to  \ref anchor_Cabral93 "this paper [Cabral93]",
/// or (more adapted for the \ref section_GreycStoration "curvature-preserving smoothing algorithm"): 
/// \ref anchor_Tschumperle06 "this paper [Tschumperle06]".
/// \param[in] src the Image that will participate in the
///                \ref anchor_LIC "line integral convolution"
/// \param[in] vec a vector field defining the streamlines along which convolution
///                will take place.  The vector field is specified as an Image with
///                four channels; the three first containing the x- , y- and z-components
///                of a normalized vector field, the third containing the \em magnitudes
///                of the vectors in the field.
/// \param[out] target the Image obtained when convoluting 'src' with 'kernel_func' along the streamlines of 'vec' (the 
///                    result of the \ref anchor_LIC "line integral convolution".
/// \param[in] kernel_func pointer to the kernel function (a function from \f$\mathcal{R}\f$ to \f$\mathcal{R}\f$.  
///                        This is typically a symmetrical function with compact support, or with infinite support but
///                        with a value that goes off to zero relatively quickly.
/// \param[in] dl steplength used for integration purposes
/// \param[in] L total length along which to integrate (in each direction from zero).  The important part of 
///              the support of 'kernel_func' should fall within this range.
/// \note This function implements line integral convolution with \em fixed steplength.  
/// A \ref LIC_3D "version using \em adaptive steplengths" is also provided in this library, but it can 
/// \ref anchor_LIC_3D_TROUBLE "run into problems".
/// For that reason, it is recommended to use this function instead, until those problems are fixed.
void LIC_3D_FS(const Image<double>& src,
	       const Image<double>& vec,
	       Image<double>& target,
	       double (*kernel_func)(double),
	       const double dl = 0.8, // steplength
	       const double L = 10); // total length along which to integrate
//===========================================================================


// NB: The below two functions were written along the lines of the paper
// "Imaging Vector Fields Using Line Integral Convolution" (Cabral).  However,
// even though the adaptive steplength is likely to allow for a more accurate
// tracing of stream lines, there have been cases where the tracing "gets stuck"
// in a cycle of 4 cells.  A quick fix for this has not been found.  Moreover,
// the accuracy of the steplength might not be needed for our purposes, which is
// smoothing of images.  Therefore, the above two functions should be used instead.
// The two functions below are however kept for future reference.
// NB!  The vector field is here specified directly by its x, y, and z coordinate ('vec')
//===========================================================================



//===========================================================================
/// Carry out \ref anchor_LIC "line integral convolution" of a 2D-image with 
/// a kernel function, along streamlines given by a user-specified vector field.
/// The result is an image smoothed along the streamlines, with a smoothing
/// defined by the kernel function.  During integration, this function uses an \em adaptive steplength.
/// For more on line integral convolution, refer to  \ref anchor_Cabral93 "this paper [Cabral93]",
/// or (more adapted for the \ref section_GreycStoration "curvature-preserving smoothing algorithm"): 
/// \ref anchor_Tschumperle06 "this paper [Tschumperle06]".
/// \param[in] src the Image that will participate in the
///                \ref anchor_LIC "line integral convolution"
/// \param[in] vec a vector field defining the streamlines along which convolution
///                will take place.  The vector field is specified as an Image with
///                two channels, expressing the x- and y-components of the vector field.
/// \param[out] target the Image obtained when convoluting 'src' with 'kernel_func' along the streamlines of 'vec' (the 
///                    result of the \ref anchor_LIC "line integral convolution".
/// \param[in] kernel_func pointer to the kernel function (a function from \f$\mathcal{R}\f$ to \f$\mathcal{R}\f$.  
///                        This is typically a symmetrical function with compact support, or with infinite support but
///                        with a value that goes off to zero relatively quickly.
/// \param[in] L total length along which to integrate (in each direction from zero).  The important part of 
///              the support of 'kernel_func' should fall within this range.
/// 
/// \anchor anchor_LIC_2D_TROUBLE
/// \note This function was written along the lines of the paper
/// \ref anchor_Cabral93 "Imaging Vector Fields Using Line Integral Convolution".  However
/// even though the adaptive steplength is likely to allow for a more accurate tracing of stream
/// lines, there have been testcases where the tracing "gets stuck" in a cycle of 4 cells.
/// A quick fix for this has not been found at the time of writing of this documentation.
/// Moreover, the accuracy of the steplength might not be needed for the main purposes of the
/// lsseg library (namely, the \ref section_GreycStoration "curvature-preserving smoothing algorithm").
/// Therefore, <span style="color: red ">use of this function should be avoided</span>, and the 
/// \ref LIC_2D_FS "fixed-steplength version" should be used instead.  However, this function
/// is kept in the library for future reference and perhaps fix.
void LIC_2D(const Image<double>& src,
	    const Image<double>& vec,
	    Image<double>& target,
	    double (*kernel_func)(double),
	    const double L = 10); // parametric length in each direction
//===========================================================================


//===========================================================================
/// Carry out \ref anchor_LIC "line integral convolution" of a 3D-image with 
/// a kernel function, along streamlines given by a user-specified vector field.
/// The result is an image smoothed along the streamlines, with a smoothing
/// defined by the kernel function. During integration, this function uses an \em adaptive steplength.
/// For more on line integral convolution, refer to  \ref anchor_Cabral93 "this paper [Cabral93]",
/// or (more adapted for the \ref section_GreycStoration "curvature-preserving smoothing algorithm"): 
/// \ref anchor_Tschumperle06 "this paper [Tschumperle06]".
/// \param[in] src the Image that will participate in the
///                \ref anchor_LIC "line integral convolution"
/// \param[in] vec a vector field defining the streamlines along which convolution
///                will take place.  The vector field is specified as an Image with
///                three channels, expressing the x-, y- and z-components of the vector field.
/// \param[out] target the Image obtained when convoluting 'src' with 'kernel_func' along the streamlines of 'vec' (the 
///                    result of the \ref anchor_LIC "line integral convolution".
/// \param[in] kernel_func pointer to the kernel function (a function from \f$\mathcal{R}\f$ to \f$\mathcal{R}\f$.  
///                        This is typically a symmetrical function with compact support, or with infinite support but
///                        with a value that goes off to zero relatively quickly.
/// \param[in] L total length along which to integrate (in each direction from zero).  The important part of 
///              the support of 'kernel_func' should fall within this range.
/// 
/// \anchor anchor_LIC_3D_TROUBLE
/// \note This function was written along the lines of the paper
/// \ref anchor_Cabral93 "Imaging Vector Fields Using Line Integral Convolution".  However
/// even though the adaptive steplength is likely to allow for a more accurate tracing of stream
/// lines, there have been testcases where the tracing "gets stuck" in a cycle of 4 cells.
/// A quick fix for this has not been found at the time of writing of this documentation.
/// Moreover, the accuracy of the steplength might not be needed for the main purposes of the
/// lsseg library (namely, the \ref section_GreycStoration "curvature-preserving smoothing algorithm").
/// Therefore, <span style="color: red">use of this function should be avoided</span>, and the 
/// \ref LIC_2D_FS "fixed-steplength version" should be used instead.  However, this function
/// is kept in the library for future reference and perhaps fix.
void LIC_3D(const Image<double>& src,
	    const Image<double>& vec,
	    Image<double>& target,
	    double (*kernel_func)(double),
	    const double L = 10); // parametric length in each direction
//===========================================================================

} // end namespace lsseg	    

#endif // _LIC_H
