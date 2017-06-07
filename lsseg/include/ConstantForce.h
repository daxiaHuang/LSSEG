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
// File: ConstantForce.h                                                     
//                                                                           
// Created: Fri Mar  3 14:21:01 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: ConstantForce.h,v 1.4 2006/09/20 22:55:47 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Header containing the definition of the
///  \ref lsseg::ForceGenerator "ForceGenerator" called \ref lsseg::ConstantForce "ConstantForce".
//                                                                           
//===========================================================================

#ifndef _CONSTANTFORCE_H
#define _CONSTANTFORCE_H

#include "ForceGenerator.h"

namespace lsseg {

//===========================================================================
/// \brief This ForceGenerator creates a constant \ref anchor_NormalForceField "normal force", 
/// which does not depend on the underlying image or its current segmentation
/// 
/// This ForceGenerator might not be very useful for segmentation purposes,
/// but can be used to make predictable \ref anchor_NormalForceField "normal force field"
/// for testing or debug purposes.
class ConstantForce : public ForceGenerator
//===========================================================================
{
public:
    /// \brief constructor.  The value of the constant normal force is specified with the
    /// argument \c c.
    ///
    /// A positive value for \c c represents a force pointing \em outwards of the domain
    /// enclosed by the curve/surface.  A negative value represents a force pointing \em inwards.
    ConstantForce(double c) : constant_(c) {}

    /// virtual destructor
    virtual ~ConstantForce() {}

    /// \brief Initializes the ConstantForce ForceGenerator.
    ///
    /// inherited from ForceGenerator.  For this class, however, the image information
    /// will not contribute to the determination of the force field, which is constant over
    /// the whole image domain..
    /// 
    /// \param img (pointer to) the image that the ForceGenerator should use to 
    ///            derive the force field (i.e. the Image to be segmented).  For this class,
    ///            the image information will \em not participate in the determination of the
    ///            (constant) force field.  However, the \em domain of the defined force field 
    ///            will be that ot the image.
    /// \param mask (pointer to) an optional \ref Mask that specifies which part(s) of the
    ///             image are active.  (No forces will be computed for the inactive
    ///             part(s) of an image).  If this pointer is left at zero, the whole
    ///             of \c img is considered to be active.
    virtual void init(const Image<double>* img, const Mask* mask = 0) {}
    virtual void update(const LevelSetFunction& phi) {}
    virtual double force2D(int x, int y) const {return constant_;}
    virtual double force3D(int x, int y, int z) const {return constant_;}
    virtual double force(size_t ix) const {return constant_;}
    virtual void force(LevelSetFunction& phi, const Mask* mask) const { phi = constant_;}

private:
    
    /// storage of the value for the constant force
    double constant_;

};

}; // end namespace lsseg

#endif // _CONSTANTFORCE_H

