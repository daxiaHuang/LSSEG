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
#include <vector>
#include "simple_tools.h"
#include "colordefs.h"
#include "LevelSetFunction.h"
#include "level_set.h"
#include "cimg_dependent.h"

using namespace std;
using namespace lsseg;

int main()
{
    const int num_regs = 4;        // number of regions
    const int num_fragments = 80;  // number of disjoint fragments
    const int X = 200; // resolution of image
    const int Y = 200; //     -  '' -
    vector<LevelSetFunction> ls_vec(num_regs, LevelSetFunction(X,Y));
    
    random_scattered_voronoi(&ls_vec[0], num_regs, num_fragments);

    for (int i = 0; i != num_regs; ++i) {
	display_image(ls_vec[i]);
    }

    vector<const int*> colors;
    colors.push_back(RED);
    colors.push_back(GREEN);
    colors.push_back(BLUE);
    colors.push_back(YELLOW);
    colors.push_back(CYAN);
    colors.push_back(MAGENTA);
    colors.push_back(GREY);
    
    Image<double> visu(X, Y, 1, 3);
    vector<const LevelSetFunction*> lsptr(num_regs);
    for (int i = 0; i < num_regs; lsptr[i] = &ls_vec[i++]);
    visualize_multisets(&lsptr[0], num_regs, visu, &colors[0]);
    display_image(visu);

    return 1;
};

