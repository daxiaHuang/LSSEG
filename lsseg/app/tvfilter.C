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
#include <iostream>
#include "cimg_dependent.h"
#include "Filters.h"

using namespace std;
using namespace lsseg;

int main(int varnum, char** vararg)
{
    if (varnum < 4) {
	cerr << "This program will display the process of TV-flow (total-varaion-diminishing flow)" << endl;
	cerr << "on an image." << endl;
	cerr << endl;
	cerr << "Usage: tvfilter <filename> <dt> <num timesteps> [savefile]" << endl;
	cerr << "Suggestion:try dt = 1 and num_timesteps = 20" << endl;
	cerr << endl;
	return -1;
    }

    Image<double> img;
    load_image(vararg[1], img, false);

    permanent_display(img);

    Image<double> res;
    
    double dt = atof(vararg[2]);
    int iter = atoi(vararg[3]);
    
    const char* save = (varnum > 4) ? vararg[4] : 0;

    UpdatableImage uimg(img, "development");

    for (int i = 0; i < iter; ++i) {
	anisotropic_smoothing(img, res, dt, 0, 1);
	uimg.update(res);
	img = res;
	cout << i << endl;
    }

    display_image(res);
    if (save) {
	cout << "Now saving result." << endl;
	save_image(save, res);
    }
    return 1;
};

