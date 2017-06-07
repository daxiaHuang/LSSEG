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
#include <math.h>
#include <time.h>
#include <stdexcept>
#include <iostream>
#include "Image.h"
#include "LIC.h"
#include "cimg_dependent.h"

using namespace lsseg;
using namespace std;

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

void monopole_field(double xpos, double ypos, double zpos, double q, Image<double>& I);

double gauss(double t) 
{
    //return 1;
    return exp(-(t*t) / 10);
}

int main()
{
    const int X = 200;
    const int Y = 200;
    const int Z = 200;
    const double FIELD_STRENGTH = 1;
    
    cout << "Now allocating memory for field..." << endl;
    Image<double> field(X, Y, Z, 3);
    cout << "Finished!" << endl;

    // deciding pole positions
    const double p1x = double(X-1)/2;   
    const double p1y = double(Y-1)/2; 
    const double p1z = double(Z-1)/3;
    const double p2x = double(X-1)/2; 
    const double p2y = double(Y-1)/2;
    const double p2z = double(2*(Z-1))/3;

    cout << "Now generating field..." << endl;
    monopole_field(p1x, p1y, p1z, FIELD_STRENGTH, field);
    cout << "halfway..." << endl;
    monopole_field(p2x, p2y, p2z, -FIELD_STRENGTH, field);
    cout << "Finished!" << endl;

    cout << "Now allocating memory for noise image..." << endl;
    Image<double> noise(X, Y, Z);
    cout << "Finished!" << endl;
    cout << "Now making noise image..." << endl;
    uniformNoise(noise.begin(), -128, 128, noise.size());
    blur_image(noise, 0.5);
    cout << "Finished!" << endl;

    cout << "Now allocating memory for target image" << endl;
    Image<double> target(X, Y, Z);
    cout << "Finished!" << endl;
    target = 0;
    clock_t tstart = clock();
    cout << "Now starting processing" << endl;
    LIC_3D(noise, field, target, &gauss, 2);
    cout << "Finished!" << endl;
    clock_t tend = clock();

    int perm[] = {1, 2, 0, 3}; // permutation of spatial coordinates 
    target.permute(perm);
    
    for (int z = 0; z < Z; z += (Z/10 < 1) ? 1 : Z/10) {
	cout << "now showing z slice: " << z << " out of " << Z << endl;
	display_image(target, z);
    }
    cout << "Total processing time: " << double(tend-tstart)/1000 << " ms." << endl;

    return 1;
};

double maxof(double a, double b, double c) 
{
    double res = a > b ? a : b;
    res = res > c ? res : c;
    return res;
}

void monopole_field(double xpos, double ypos, double zpos, double q, Image<double>& I)
{
    const int X = I.dimx();
    const int Y = I.dimy();
    const int Z = I.dimz();
    const double K = 0.001; // regularising constant
    const double M = maxof(X, Y, Z);

    for (int z = 0; z < Z; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		const double dx = (x - xpos)/M;
		const double dy = (y - ypos)/M;
		const double dz = (z - zpos)/M;
		const double n = sqrt(dx*dx + dy*dy + dz*dz);
		
		const double fac = double(q) / (K + (n * n * n));
		
		I(x,y,z,0) += fac * dx;
		I(x,y,z,1) += fac * dy;
		I(x,y,z,2) += fac * dz;
	    }
	}
    }
}
