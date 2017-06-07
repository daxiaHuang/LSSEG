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
#include <time.h>
#include <iostream>
#include "Image.h"
#include "Filters.h"
#include "cimg_dependent.h"


using namespace lsseg;
using namespace std;

//const char filename[] = "data/zebra.jpg";

int main(int varnum, char** vararg)
{
    if (varnum < 2) {
	cerr << endl;
	cerr << "This program takes an image as input, converts it to greycsale,  and computes " << endl;
	cerr << "and displays the following derived images : " << endl;
	cerr << "1) The structure tensor, with its 3 different components " << endl;
	cerr << "2) The image smoothed with the total variation diminishing scheme " << endl;
	cerr << "3) The scale measure (gives an idea about the magnitude of the structures " << endl;
	cerr << "   to which each pixel belongs" << endl;
	cerr << endl;
	cerr << "Push spacebar to go from one visualization to the next. " << endl;
	cerr << endl;
	cerr << "Usage: filter2DTests <filename> " << endl;
	cerr << endl;
	return -1;
    }


    const char* filename = vararg[1];

    Image<double> img;
    load_image(filename, img, true); // to grey
    display_image(img);

    // ---- TEST OF STRUCTURE TENSOR ---- 

    Image<double> G;
    cout << "Entering compute_structure_tensor_2D" << endl;
    clock_t start = clock();
    compute_structure_tensor_2D(img, G);
    clock_t end = clock();
    cout << "Exited from compute_structure_tensor_2D" << endl;
    cout << "Computation took: " << (end - start) / double(1000) << " ms." << endl;

    cout << "Displaying all three channels: " << endl;
    display_image(G);

    Image<double> G_ch;
    for (int ch = 0; ch < G.numChannels(); ++ch) {
	G.makeChannelImage(G_ch, ch);
	cout << "Displaying separate channel: " << ch << endl;
	display_image(G_ch);
    }

    // ---- TEST OF ANISOTROPIC SMOOTHING ---- 
    const double dt = 1;
    const double T = 20;

    Image<double> NGF(img);
    Image<double> NGF2(img, false);
    for (double t = 0; t < T; t+= dt) {
	cout << t << endl;
	nonlinear_gauss_filter_2D(NGF, NGF2, dt, 0, 1);
	NGF2.swap(NGF);
    }
    cout << "Displaying total variation image: " << endl;
    display_image(NGF);
		       
    // ---- TEST OF SCALE MEASURE ---- 
    
    Image<double> S, Tacc;
    compute_scale_factor_2D(img, S, Tacc, dt, T);
    cout << "Displaying scale image: " << endl;
    display_image(S);



};

