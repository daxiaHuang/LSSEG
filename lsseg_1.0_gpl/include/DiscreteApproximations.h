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
// File: DiscreteApproximations.h                                            
//                                                                           
// Created: Wed Feb 22 17:49:15 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: DiscreteApproximations.h,v 1.2 2006/09/20 22:55:48 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief This file contains functions for computing the squared gradient norm
/// over an \ref lsseg::Image "Image<double>".
///
//                                                                           
//===========================================================================

#ifndef _DISCRETEAPPROXIMATIONS_H
#define _DISCRETEAPPROXIMATIONS_H

#include "Image.h"

namespace lsseg {

//===========================================================================
/// \brief Compute an estimate of the squared gradient norm, for each pixel in 
/// a two-dimensional image.  
///
/// If the image contains more than one channel, the \em sum of the squared gradient
/// norms of each channel will be computed for each image pixel.
///
/// \param[in] image The input image.  It is supposed to be 2D (only 1 pixel along
///                  <tt>z</tt>-direction).  It can have an arbitrary number of channels.
/// \param[out] result Upon completion of the function, this image will have been
///                    resized to the same resolution as \c image, but will contain 
///                    only \em one channel.  Each pixel in this channel corresponds
///                    to the estimated squared gradient norm of the corresponding 
///                    pixel in \c image (or, in the case that \c image is multichanneled,
///                    it will correspond to the \em sum of the squared gradient norms of 
///                    this pixel for each image channel in \c image.
void compute_squared_gradient_sum_2D(const Image<double>& image, 
				     Image<double>& result);
//===========================================================================

//===========================================================================
/// \brief Compute an estimate of the squared gradient norm, for each pixel in 
/// a three-dimensional image.  
///
/// If the image contains more than one channel, the \em sum of the squared gradient
/// norms of each channel will be computed for each image pixel.
///
/// \param[in] image The input image.  It is supposed to be 3D (more than 1 pixel along
///                  <tt>z</tt>-direction).  \c image can in theory also be a 2D image,
///                  but in that case it is better to call the slightly faster function
///                  compute_squared_gradient_sum_2D(). 
///                  \c image can have an arbitrary number of channels.
/// \param[out] result Upon completion of the function, this image will have been
///                    resized to the same resolution as \c image, but will contain 
///                    only \em one channel.  Each pixel in this channel corresponds
///                    to the estimated squared gradient norm of the corresponding 
///                    pixel in \c image (or, in the case that \c image is multichanneled,
///                    it will correspond to the \em sum of the squared gradient norms of 
///                    this pixel for each image channel in \c image.
void compute_squared_gradient_sum_3D(const Image<double>& image, 
				     Image<double>& result);
//===========================================================================


};

#endif // _DISCRETEAPPROXIMATIONS_H

