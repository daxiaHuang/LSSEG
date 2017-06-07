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
// File: DiscreteApproximations.C                                            
//                                                                           
// Created: Wed Feb 22 17:50:23 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: DiscreteApproximations.C,v 1.2 2006/11/13 02:29:27 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements DiscreteApproximations.h
//                                                                           
//===========================================================================

#include "DiscreteApproximations.h"

namespace lsseg {

//===========================================================================
void compute_squared_gradient_sum_2D(const Image<double>& image, Image<double>& result)
//===========================================================================
{
    assert(image.dimz() == 1);
    result.resize(image.dimx(), image.dimy(), 1, 1);
    result = 0;
    double Inc, Ipc, Icn, Icp;
    const int X = image.dimx();
    const int Y = image.dimy();
    int xp, xn, yp, yn;

    for (int c = 0; c < image.numChannels(); ++c) {
	for (int y = 0; y < Y; ++y) {
	    yp = (y > 0) ? y - 1 : y;
	    yn = (y+1 < Y) ? y + 1 : y;
	    double y_div = (yp + 2 == yn) ? 0.5 : (yp + 1 == yn) ? 1 : 0;
	    for (int x = 0; x < X; ++x) {
		xp = (x > 0) ? x - 1 : x;
		xn = (x+1 < X) ? x + 1 : x;
		double x_div = (xp + 2 == xn) ? 0.5 : (xp + 1 == xn) ? 1 : 0;
		
		Ipc = image(xp,  y,  0, c);
		Inc = image(xn,  y,  0, c);
		Icp = image( x, yp,  0, c);
		Icn = image( x, yn,   0, c);
		
		double dx = (Inc - Ipc) * x_div;
		double dy = (Icn - Icp) * y_div;
		
		result(x, y) += (dx * dx) + (dy * dy);
	    }
	}
    }
}

//===========================================================================
void compute_squared_gradient_sum_3D(const Image<double>& image, Image<double>& result)
//===========================================================================
{
    result.resize(image.dimx(), image.dimy(), image.dimz(), 1);
    result = 0;
    double Incc, Ipcc, Icnc, Icpc, Iccn, Iccp;
    const int X = image.dimx();
    const int Y = image.dimy();
    const int Z = image.dimz();
    int xp, xn, yp, yn, zp, zn;

    for (int c = 0; c < image.numChannels(); ++c) {
	for (int z = 0; z < Z; ++z) {
	    zp = (z > 0) ? z - 1 : z;
	    zn = (z+1 < Z) ? z + 1 : z;
	    double z_div = (zp + 2 == zn) ? 0.5 : (zp + 1 == zn) ? 1 : 0;
	    for (int y = 0; y < Y; ++y) {
		yp = (y > 0) ? y - 1 : y;
		yn = (y+1 < Y) ? y + 1 : y;
		double y_div = (yp + 2 == yn) ? 0.5 : (yp + 1 == yn) ? 1 : 0;
		for (int x = 0; x < X; ++x) {
		    xp = (x > 0) ? x - 1 : x;
		    xn = (x+1 < X) ? x + 1 : x;
		    double x_div = (xp + 2 == xn) ? 0.5 : (xp + 1 == xn) ? 1 : 0;

		    Ipcc = image(xp,  y,  z, c);
		    Incc = image(xn,  y,  z, c);
		    Icpc = image( x, yp,  z, c);
		    Icnc = image( x, yn,  z, c);
		    Iccp = image( x,  y, zp, c);
		    Iccn = image( x,  y, zn, c);
		    
		    double dx = (Incc - Ipcc) * x_div;
		    double dy = (Icnc - Icpc) * y_div;
		    double dz = (Iccn - Iccp) * z_div;

		    result(x, y, z) += (dx * dx) + (dy * dy) + (dz * dz);
		}
	    }
	}
    }
}





}; // end namespace lsseg
