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

//mconst char filename[] = "data/zebra.jpg";

int main(int varnum, char** vararg)
{
    if (varnum < 2) {
	cerr << endl;
	cerr << "This program takes as input an image, and computes and display its structure tensor. " << endl;
	cerr << "First, the original image will be displayed.  Then, a color image will be shown, where " << endl;
	cerr << "the RGB values represent the three different components of the structure tensor.  Then, " << endl;
	cerr << "each channel will be shown separately." << endl;
	cerr << endl;
	cerr << "Usage: structureTensorComputation <image filename> " << endl;
	cerr << endl;
	return -1;
    }

    const char* filename = vararg[1];

    Image<double> img;
    load_image(filename, img);
    display_image(img);
    Image<double> G;
    cout << "Entering compute_structure_tensor_2D" << endl;
    clock_t start = clock();
    compute_structure_tensor_2D(img, G);
    clock_t end = clock();
    cout << "Exited from compute_structure_tensor_2D" << endl;
    cout << "Computation took: " << (end - start) / double(1000) << " ms." << endl;
    display_image(G);
    Image<double> G_ch;
    for (int ch = 0; ch < G.numChannels(); ++ch) {
	G.makeChannelImage(G_ch, ch);
	display_image(G_ch);
    }
};

