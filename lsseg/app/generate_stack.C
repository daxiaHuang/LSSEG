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
#include "Image.h"
#include "level_set.h"
#include "simple_tools.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lsseg;

int main(int varnum, char** vararg) 
{
    if (varnum < 2) {
	cerr << endl << "usage: volumesmoothing_test <filelist_name> [savefile]" << endl << endl;
	cerr << "The first argument (filelist_name) should be an ASCII file containing the name " << endl;
	cerr << "of one image file per line.  All the images mentioned are supposed to have the " << endl;
	cerr << "same x and y resolution and number of channels.  The fileformats are directly " << endl;
	cerr << "specified by the file suffixes (jpeg, png, etc.)" << endl << endl;
	cerr << "The saved file will contain a 3D image where the original images are stacked along" << endl;
	cerr << "the z coordinate.  This file can be read into an lsseg::Image object by using its" << endl;
	cerr << "read() member function.  The saved file should conventionally have the suffix .stack ." << endl;
	cerr << endl;
	return -1;
    }
    string savefile;
    if (varnum > 2) {
	savefile = string(vararg[2]);
    }

    string filelist_name(vararg[1]);
    ifstream is(filelist_name.c_str());
    if (!is) {
	cerr << "Unable to open file: " << filelist_name << endl;
	return -1;
    }

    Image<double> I;
    read_image_sequence(is, I, true);

    cout << "Image dimensions are: ";
    cout << I.dimx() << " " << I.dimy() << " " << I.dimz() << endl;
    cout << "Number of channels: " << I.numChannels() << endl;

    if (savefile.size() > 0) {
	ofstream os(savefile.c_str());
	I.write(os, true);
	os.close();
    }

    return 0;
}
