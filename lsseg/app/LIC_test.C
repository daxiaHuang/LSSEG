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
#include "LIC.h"
#include "cimg_dependent.h"
#include <stdexcept>

using namespace std;
using namespace lsseg;

// function taken from the Sintef Applied Mathematics GoTools library ('utils' module)
inline void uniformNoise(double* res, double lval, double uval, int num_samples)
{
    if (uval <= lval) {
	throw runtime_error("uniformNoise(...) : erroneous range.");
    }
    double range = uval - lval;
    double scale_factor = range / double(RAND_MAX);
    for (int i = 0; i < num_samples; ++i) {
	res[i] = double(rand()) * scale_factor + lval;
    }
}


double gauss(double t) 
{
    //return 1;
    return exp(-(t*t) / 1000);
}

int main()
{
    const int X = 200;
    const int Y = 200;

    Image<double> noise(X, Y);
    uniformNoise(noise.begin(), -128, 128, noise.size());
    blur_image(noise, 0.5);
    permanent_display(noise);
    Image<double> vecfield(X, Y, 1, 2);
    for (int y = 0; y < Y; ++y) {
	for (int x = 0; x < X; ++x) {

	    //vecfield(x, y, 0, 0) = 0.1 * 1;
	    //vecfield(x, y, 0, 1) = 0.1 * sin(double(x)/90);

// 	    vecfield(x, y, 0, 0) = double(x)/100;
// 	    vecfield(x, y, 0, 1) = 0;

// 	    double xx = double(x-50)/50;
// 	    double yy = double(y-50)/50;
// 	    double krull = 1;//sqrt(xx * xx + yy * yy);
// 	    vecfield(x, y, 0, 0) = xx / krull;
// 	    vecfield(x, y, 0, 1) = yy / krull;

// 	    vecfield(x, y, 0, 0) = x < 50 ? 1 : -1;;
// 	    vecfield(x, y, 0, 1) = 0;
	    

//   	    vecfield(x, y, 0, 0) = cos(double(y)/10);
//  	    vecfield(x, y, 0, 1) = sin(double(x)/10);

	    double dx = x - X/2;
	    double dy = y - Y/2;
	    double tmp = double(1) / sqrt(dx * dx + dy * dy);
	    vecfield(x, y, 0, 0) = -dy * tmp;
	    vecfield(x, y, 0, 1) = dx * tmp;
	}
    }
//    display_image(krull);
    Image<double> target(noise,false);
    LIC_2D(noise, vecfield, target, &gauss, 10);
    display_image(target);

    return 1;
}
