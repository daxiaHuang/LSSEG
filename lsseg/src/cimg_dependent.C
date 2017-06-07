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
// File: cimg_dependent.C                                                    
//                                                                           
// Created: Tue Oct 25 13:29:10 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: cimg_dependent.C,v 1.22 2006/11/25 20:08:28 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements cimg_dependent.h
//                                                                           
//===========================================================================

#include <set>
#include <string>
#include "errormacros.h"
#include "CImg.h"
#include "Image.h"
#include "cimg_dependent.h"
#include "simple_tools.h"
//#include "GoBorrowedMVGrid.h"
//#include "GoSelfcontainedGrid.h"

//using namespace Go;
using namespace std;
using namespace cimg_library;
using namespace lsseg;

//===========================================================================
namespace 
//===========================================================================
{


    void image2cimg(const Image<double>& img, CImg<double>& target, int z = -1)
    {
	if (z < 0) { // copy everything
	    CImg<double> tmpImg(img.dimx(), img.dimy(), img.dimz(), img.numChannels());
	    copy(img.begin(), img.end(), tmpImg.data);
	    tmpImg.swap(target);
	} else { // copy only one image
	    CImg<double> tmpImg(img.dimx(), img.dimy(), 1, img.numChannels());
	    for (int c = 0; c < img.numChannels(); ++c) {
		//cout << z << endl;
		unsigned int offset = c * img.channelSize() + z * img.graySize2D();
		//cout << offset << " " << offset + img.graySize2D() << endl;
		copy(img.begin() + offset, 
		     img.begin() + offset + img.graySize2D(),
		     tmpImg.data + c * img.graySize2D());
	    }
	    tmpImg.swap(target);
	}
    }

    template<typename T>
    void cimg2image(const CImg<T>& img, Image<T>& target) 
    {
	target.resize(img.dimx(), img.dimy(), img.dimz(), img.dimv());
	copy(img.data, img.data + img.size(), target.begin());
    }

}; // end anonymous namespace


namespace lsseg {

//===========================================================================
void load_image(const char* name, Image<double>& target, bool convert_to_grayscale)
//===========================================================================
{
    CImg<double> tmp(name);
    cimg2image(tmp, target);
    if (convert_to_grayscale) {
	to_grayscale(target);
    }
}

//===========================================================================
void load_image(const char* name, Image<int>& target, bool convert_to_grayscale)
//===========================================================================
{
    CImg<int> tmp(name);
    cimg2image(tmp, target);
    if (convert_to_grayscale) {
	to_grayscale(target);
    }
}

//===========================================================================
void save_image(const char* name, const Image<double>& img)
//===========================================================================
{
    CImg<double> tmp(1, 1, 1, 1);
    image2cimg(img, tmp);
    tmp.save(name);
}

//===========================================================================
void display_image(const Image<double>& img, int z)
//===========================================================================
{
    //   return;
//#ifdef NDEBUG
    CImg<double> tmp(1, 1, 1, 1);
    image2cimg(img, tmp, z);
    tmp.display();
//#endif
}



//===========================================================================
void permanent_display(const Image<double>& img)
//===========================================================================
{
    //  the pointer to the display window will be lost, so there WILL be 
    // a memory leak.  However, this function is only intended for debug
    // purposes, and the display window is expected to live throughout the
    // execution of whatever test program uses this function.
    CImg<double> tmp(1, 1, 1, 1);
    image2cimg(img, tmp);
    new CImgDisplay(tmp);
}

//===========================================================================
void blur_image(Image<double>& img, double rho)
//===========================================================================
{
    CImg<double> tmp(1,1,1,1);
    image2cimg(img,tmp);
    tmp.blur(rho);
    cimg2image(tmp, img);
}

//===========================================================================
void blur_1D(double* data, unsigned int data_size, double rho)
//===========================================================================
{
    CImg<double> tmp(data_size, 1, 1, 1);
    std::copy(data, data + data_size, tmp.data);
    tmp.blur(rho);
    std::copy(tmp.data, tmp.data + data_size, data);
}

//===========================================================================
void gaussian_noise(Image<double>& img, double sigma)
//===========================================================================
{
    CImg<double> tmp(1,1,1,1);
    image2cimg(img,tmp);
    tmp.noise(sigma, 0);
    cimg2image(tmp, img);
}
//
////===========================================================================
//void display_distribution(const Image<double>& img)
////===========================================================================
//{
//    // only applied on the first image in the stack represented by 'img'
//    MESSAGE_IF(img.dimz() != 1,
//	       "Warning: displaying distribution only on first image of stack.");
//    
//    CImg<double> tmp(1, 1, 1, 1);
//    const unsigned int level_size = img.dimx() * img.dimy();
//    for (int d = 0; d < img.numChannels(); ++d) {
//	Image<double> tmp_img(img.dimx(), img.dimy());
//	copy(img.begin() + d * level_size, img.begin() + (d+1) * level_size, tmp_img.begin());
//	image2cimg(tmp_img, tmp);
//	CImg<double> hist = tmp.get_histogram();
//	display_graph(hist.dimx(), hist.data);
//    }
//}

//===========================================================================
/// \cond OMIT_FROM_DOXYGEN
struct Pimpl  // implementation of UpdatableImage
//===========================================================================
{
    Pimpl(const Image<double>& img, const char* name) : 
	disp_(img.dimx(), img.dimy(), name), name_(name)
    {
	update(img, false);
    }

    void update(const Image<double>& img, bool reshape) {
	image2cimg(img, i_);
	if (disp_.is_resized) {
	    resize(); // this will also display the new i_
	} else {
	    disp_.display(i_);
	}
	if (reshape) {
	    resize(img.dimx(), img.dimy());
	}
	disp_.show();	
    }

    void callbackLoop() 
    {
	while(true) {
	    int mouse_state = disp_.button;
	    if (mouse_state & 1 ) { // left button or space
		int x, y;
		getMouseXY(x, y);
		cout << "(" << x << ", " << y << ") : " << endl;
		printColorInfo(x, y);
		set<Pimpl*>::iterator it;
		for (it = connecteds_.begin(); it != connecteds_.end(); ++it) {
		    (*it)->printColorInfo(x, y);
		}
		cout << endl;
		disp_.wait();
	    }
	    if (mouse_state & 2 || disp_.key == cimg::keySPACE) { // right button or space
		disp_.key =0;
		break; // handle control to program
	    }
	    if (mouse_state & 4) { // middle button
		cout << "no function for the middle button yet" << endl; // print friendly message
	    }
	    if (disp_.is_resized) {
		resize();
		set<Pimpl*>::iterator it;
		for (it = connecteds_.begin(); it != connecteds_.end(); ++it) {
		    (*it)->resize(disp_.window_width, disp_.window_height);
		}
	    }
	} 
    }

    void connect(Pimpl* rhs) {
	connecteds_.insert(rhs);
	rhs->connecteds_.insert(this);
    }

    void resize(int x = -1, int y = -1) 
    {
	if (x > 0 && y > 0) {
	    disp_.resize(x, y);
	} else {
	    disp_.resize();
	}
	disp_.display(i_);
    }

private:
    // data members 
    CImg<double> i_;
    CImgDisplay disp_;
    string name_;
    set<Pimpl*> connecteds_;

    void getMouseXY(int& x, int& y) 
    {
	x = disp_.mouse_x;
	y = disp_.mouse_y;
	
	x *= i_.dimx();
	y *= i_.dimy();
	
	x /= disp_.window_width;
	y /= disp_.window_height;
    }

    void printColorInfo(int x, int y) {
	cout << name_ << ": ";
	cout << "[ ";
	for (int i = 0; i < i_.dimv(); ++i) {
	    cout << i_(x, y, 0, i) << " ";
	}
	cout << "]" << endl;
    }
};
/// \endcond

//===========================================================================
UpdatableImage::UpdatableImage(const Image<double>& img, 
			       const char* name) : p_(new Pimpl(img, name)) {}
UpdatableImage::~UpdatableImage() { delete p_;}
void UpdatableImage::update(const Image<double>& img, bool reshape) { p_->update(img, reshape);}
void UpdatableImage::interact() { p_->callbackLoop();}
void UpdatableImage::connect(UpdatableImage& rhs) { p_->connect(rhs.p_);}
//===========================================================================
 
}; // end namespace lsseg
