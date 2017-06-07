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
// File: LIC.C                                                               
//                                                                           
// Created: Mon Apr 24 18:09:44 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: LIC.C,v 1.11 2006/11/25 20:08:26 oan Exp $
//                                                                           
// Description:
/// \file
/// \brief Implements LIC.h.
//                                                                           
//===========================================================================

#include <iostream> // debug
#include <cmath>
#include <limits>
#include "errormacros.h"
#include "LIC.h"

#include "simplethreads.h" // test this

const int NUM_PROC=4; // how many threads should run at the same time (>= nb. of processors)

using namespace lsseg;
using namespace std;
//using namespace Go;

namespace { // anonymous namespace

const double TAN_LIM_ANGLE = tan(0.05); // tan-function of 3 degrees, approx.
const double SIN2_LIM_ANGLE = sin(0.05) * sin(0.05); // sin2 function of 3 degrees, approx.
const double ROUNDOFF_TERM = 1.0e-5; // added to stepsize to ensure entry into neighbour cell
const double INFINITE = std::numeric_limits<double>::max();

 double steplength(const double posx, 
			 const double posy, 
			 const double vecx,
			 const double vecy);


inline double steplength(const double posx,
			 const double posy,
			 const double posz,
			 const double vecx,
			 const double vecy,
			 const double vecz);

struct LIC_idata {
    int start; // refers to y or z coord, depending on 2D or 3D
    int end;   // idem
    const Image<double>* src;
    const Image<double>* vec;
    double* target_ptr;
    double (*kernel_func)(double);
    double L;
};

struct LIC_FS_idata {
    int start; // refers to y or z coord, depending on 2D or 3D 
    int end; // idem
    const Image<double>* src;
    const Image<double>* vec;
    double* target_ptr;
    double (*kernel_func)(double);
    double L;
    double dl;
};

void LIC_2D_subpart(const int job,
		    const int thread,
		    void * const idata,
		    void * const odata);

void LIC_3D_subpart(const int job,
		    const int thread,
		    void * const idata,
		    void * const odata);

void LIC_2D_FS_subpart(const int job,
		       const int thread,
		       void * const idata,
		       void * const odata);
void LIC_3D_FS_subpart(const int job,
		       const int thread,
		       void * const idata,
		       void * const odata);

}; // end anonymous namespace

namespace lsseg {

void LIC_2D_FS(const Image<double>& src,
	       const Image<double>& vec,
	       Image<double>& target,
	       double (*kernel_func)(double),
	       const double dl, // steplength
	       const double L) // total length along which to integrate
{
    ALWAYS_ERROR_IF(!src.spatial_compatible(vec), 
		    "source image and vector field size mismatch");
    ALWAYS_ERROR_IF(vec.numChannels() != 3, "given vector field was not 2D + magnitude");
    ALWAYS_ERROR_IF(!target.size_compatible(src), "target not same size as input image");

    // preparing indata for the threads
    const int Y = src.dimy();
    const int USED_PROC = Y < NUM_PROC ? Y : NUM_PROC;
    vector<LIC_FS_idata> job_idata(USED_PROC); 
    for (int i = 0; i < USED_PROC; ++i) {
	int y1 = (i * Y) / USED_PROC;
	int y2 = ((i+1) * Y) / USED_PROC;
	job_idata[i].start = y1;
	job_idata[i].end = y2;

	job_idata[i].target_ptr = &(target(0, y1, 0, 0));
	job_idata[i].src = &src;
	job_idata[i].vec = &vec;
	job_idata[i].kernel_func = kernel_func;
	job_idata[i].L = L;
	job_idata[i].dl = dl;
    }
    // carry out the work
    do_threaded_jobs(LIC_2D_FS_subpart, &job_idata[0], USED_PROC, USED_PROC, false, NULL);  
}

void LIC_3D_FS(const Image<double>& src,
	       const Image<double>& vec,
	       Image<double>& target,
	       double (*kernel_func)(double),
	       const double dl, // steplength
	       const double L) // total length along which to integrate
{
    ALWAYS_ERROR_IF(!src.spatial_compatible(vec), 
		    "source image and vector field size mismatch");
    ALWAYS_ERROR_IF(vec.numChannels() != 4, "given vector field was not 3D + magnitude");
    ALWAYS_ERROR_IF(!target.size_compatible(src), "target not same size as input image");

    // preparing indata for the threads
    const int Z = src.dimz();
    const int USED_PROC = Z < NUM_PROC ? Z : NUM_PROC;
    vector<LIC_FS_idata> job_idata(USED_PROC);
    for (int i = 0; i < USED_PROC; ++i) {
	int z1 = (i * Z) / USED_PROC;;
	int z2 = ((i+1) * Z) / USED_PROC;
	job_idata[i].start = z1;
	job_idata[i].end = z2;

	job_idata[i].target_ptr = &(target(0, 0, z1, 0));
	job_idata[i].src = &src;
	job_idata[i].vec = &vec;
	job_idata[i].kernel_func = kernel_func;
	job_idata[i].L = L;
	job_idata[i].dl = dl;
    }
    // carry out the work
    do_threaded_jobs(LIC_3D_FS_subpart, &job_idata[0], USED_PROC, USED_PROC, false, NULL);
}




void LIC_3D(const Image<double>& src,
	    const Image<double>& vec,
	    Image<double>& target,
	    double (*kernel_func)(double),
	    const double L)
{
    ALWAYS_ERROR_IF(!src.spatial_compatible(vec), 
		    "source image and vector field size mismatch");
    ALWAYS_ERROR_IF(vec.numChannels() != 3, "given vector field was not 3D");
    ALWAYS_ERROR_IF(!target.size_compatible(src), "target not same size as input image");

    // preparing indata for the threads
    const int Z = src.dimz();
    const int USED_PROC = Z < NUM_PROC ? Z : NUM_PROC;
    vector<LIC_idata> job_idata(USED_PROC);
    for (int i = 0; i < USED_PROC; ++i) {
	int z1 = (i * Z) / USED_PROC;;
	int z2 = ((i+1) * Z) / USED_PROC;
	job_idata[i].start = z1;
	job_idata[i].end = z2;

	job_idata[i].target_ptr = &(target(0, 0, z1, 0));
	job_idata[i].src = &src;
	job_idata[i].vec = &vec;
	job_idata[i].kernel_func = kernel_func;
	job_idata[i].L = L;
    }
    // carry out the work
    do_threaded_jobs(LIC_3D_subpart, &job_idata[0], USED_PROC, USED_PROC, false, NULL);
}


void LIC_2D(const Image<double>& src,
	    const Image<double>& vec,
	    Image<double>& target,
	    double (*kernel_func)(double),
	    const double L)
{
    ALWAYS_ERROR_IF(!src.spatial_compatible(vec), 
		    "source image and vector field size mismatch");
    ALWAYS_ERROR_IF(vec.numChannels() != 2, "given vector field was not 2D");
    ALWAYS_ERROR_IF(!target.size_compatible(src), "target not same size as input image");

    // preparing indata for the threads
    const int Y = src.dimy();
    const int USED_PROC = Y < NUM_PROC ? Y : NUM_PROC;
    vector<LIC_idata> job_idata(USED_PROC); 
    for (int i = 0; i < USED_PROC; ++i) {
	int y1 = (i * Y) / USED_PROC;
	int y2 = ((i+1) * Y) / USED_PROC;
	job_idata[i].start = y1;
	job_idata[i].end = y2;

	job_idata[i].target_ptr = &(target(0, y1, 0, 0));
	job_idata[i].src = &src;
	job_idata[i].vec = &vec;
	job_idata[i].kernel_func = kernel_func;
	job_idata[i].L = L;
    }
    // carry out the work
    do_threaded_jobs(LIC_2D_subpart, &job_idata[0], USED_PROC, USED_PROC, false, NULL);  
}


// void LIC_2D(const Image<double>& src,
// 	    const Image<double>& vec,
// 	    Image<double>& target,
// 	    double (*kernel_func)(double),
// 	    const double L)
// {
//     ALWAYS_ERROR_IF(!src.spatial_compatible(vec), 
// 		    "source image and vector field size mismatch");
//     ALWAYS_ERROR_IF(vec.numChannels() != 2, "given vector field was not 2D");
//     ALWAYS_ERROR_IF(!target.size_compatible(src), "target not same size as input image");
//     //target.resize(src);
//     //target = 0;
//     const int X = src.dimx();
//     const int Y = src.dimy();
//     const int C = src.numChannels();
//     std::vector<double> px_val(C);
//     double posx, posy; // position of current point on the traced line
//     for (int y = 0; y < Y; ++y) {
// 	for (int x = 0; x < X; ++x) {
// 	    fill(px_val.begin(), px_val.end(), 0);
// 	    double accum = 0;
// 	    int sign = 1;
// 	    for (int pass = 0; pass != 2; ++pass) { // first forwards, then backwards
		
// 		double posx = x + 0.5; // position of current point on the traced line
// 		double posy = y + 0.5; //                  - " -
// 		double l = 0;
// 		bool end_reached = false;
// 		int prev_x = int(posx);
// 		int prev_y = int(posy);
// 		do {
// 		    const int iposx = int(posx); // integer position
// 		    const int iposy = int(posy); // integer position
// 		    const double w_x = sign * vec(iposx, iposy, 0, 0); // roundoff to integer
// 		    const double w_y = sign * vec(iposx, iposy, 0, 1); // roundoff to integer
// 		    double dl = steplength(posx, posy, w_x, w_y); 
// 		    if (dl < 0) break; // degeneracy detected
// 		    l += dl; // l measures parametric length traversed
// 		    if (l > L) {
// 			dl = (l - L);
// 			l = L;
// 			// this will be the last cycle of updates.
// 			end_reached = true;
// 		    }
		    
// 		    // contribution to integral
// 		    // integrate over param.
// 		    const double h = kernel_func(l) * dl;
// 		    for (int c = 0; c < C; ++c) {
// 			px_val[c] += h * src(iposx, iposy, 0, c); 
// 		    }
// 		    accum += h;
		    
// 		    // updating position
// 		    posx += dl * w_x;
// 		    posy += dl * w_y;
// 		    if (posx < 0 || posx > X) break; // outside image domain
// 		    if (posy < 0 || posy > Y) break; // outside image domain
// 		    if (int(posx) == prev_x && int(posy) == prev_y) break; // degen. point in field
// 		    prev_x = iposx; // NB: iposx and iposy have not yet been updated
// 		    prev_y = iposy; // from last iteration, so we here get the previous value
// 		} while (!end_reached);
// 		sign *= -1;
// 	    }
// 	    if (accum > 0) {
// 		for (int c = 0; c < C; ++c) 
// 		    target(x, y, 0, c) += px_val[c] / accum;
// 	    } else {
// 		for (int c = 0; c < C; ++c) 
// 		    target(x, y, 0, c) += src(x, y, 0, c);
// 	    }
// 	}
//     }
// }

// void LIC_3D(const Image<double>& src,
// 	    const Image<double>& vec,
// 	    Image<double>& target,
// 	    double (*kernel_func)(double),
// 	    const double L)
// {
//     ALWAYS_ERROR_IF(!src.spatial_compatible(vec), 
// 		    "source image and vector field size mismatch");
//     ALWAYS_ERROR_IF(vec.numChannels() != 3, "given vector field was not 2D");
//     ALWAYS_ERROR_IF(!target.size_compatible(src), "target not same size as input image");

//     const int X = src.dimx();
//     const int Y = src.dimy();
//     const int Z = src.dimz();
//     const int C = src.numChannels();
//     std::vector<double> px_val(C);
//     double posx, posy, posz; // position of current point on the traced line
//     for (int z = 0; z < Z; ++z) {
// 	for (int y = 0; y < Y; ++y) {
// 	    for (int x = 0; x < X; ++x) {
// 		fill(px_val.begin(), px_val.end(), 0);
// 		double accum = 0;
// 		int sign = 1;
// 		for (int pass = 0; pass != 2; ++pass) { // first forwards, then backwards
		    
// 		    double posx = x + 0.5; // position of current point on the traced line
// 		    double posy = y + 0.5; //                  - " -
// 		    double posz = z + 0.5; //                  - " -
// 		    double l = 0;
// 		    bool end_reached = false;
// 		    int prev_x = int(posx);
// 		    int prev_y = int(posy);
// 		    int prev_z = int(posz);
// 		    do {
// 			const int iposx = int(posx); // integer position
// 			const int iposy = int(posy); // integer position
// 			const int iposz = int(posz);
// 			const double w_x = sign * vec(iposx, iposy, iposz, 0); // roundoff to int
// 			const double w_y = sign * vec(iposx, iposy, iposz, 1); // roundoff to int
// 			const double w_z = sign * vec(iposx, iposy, iposz, 2); // roundoff to int
// 			double dl = steplength(posx, posy, posz, w_x, w_y, w_z); 
// 			if (dl < 0) break; // degeneracy detected
// 			l += dl; // l measures parametric length traversed
// 			if (l > L) {
// 			    dl = (l - L);
// 			    l = L;
// 			    // this will be the last cycle of updates.
// 			    end_reached = true;
// 			}
		    
// 			// contribution to integral
// 			// integrate over param.
// 			const double h = kernel_func(l) * dl;
// 			for (int c = 0; c < C; ++c) {
// 			    px_val[c] += h * src(iposx, iposy, iposz, c); 
// 			}
// 			accum += h;
			
// 			// updating position
// 			posx += dl * w_x;
// 			posy += dl * w_y;
// 			posz += dl * w_z;
// 			if (posx < 0 || posx > X) break; // outside image domain
// 			if (posy < 0 || posy > Y) break; // outside image domain
// 			if (posz < 0 || posz > Z) break; // outside image domain
// 			if (int(posx) == prev_x && 
// 			    int(posy) == prev_y &&
// 			    int(posz) == prev_z) break; // degen. point in field
// 		    prev_x = iposx; // NB: iposx and iposy have not yet been updated
// 		    prev_y = iposy; // from last iteration, so we here get the previous value
// 		    prev_z = iposz; //
// 		    } while (!end_reached);
// 		    sign *= -1;
// 		}
// 		if (accum > 0) {
// 		    for (int c = 0; c < C; ++c) 
// 			target(x, y, z, c) += px_val[c] / accum;
// 		} else {
// 		    for (int c = 0; c < C; ++c) 
// 			target(x, y, z, c) += src(x, y, z, c);
// 		}
// 	    }
// 	}
//     }
// }

};// End namespace lsseg


namespace { // anonymous namespace

void LIC_3D_subpart(const int job,
		    const int thread,
		    void * const idata,
		    void * const odata)
{
    // grabbing parameters
    LIC_idata& input = ((LIC_idata*)idata)[job];

    const int zstart = input.start;
    const int zend = input.end;
    const Image<double>& src = *(input.src);
    const Image<double>& vec = *(input.vec);
    double* t_ptr = input.target_ptr;
    const double L = input.L;
    double (*kernel_func)(double) = input.kernel_func;

    std::cout <<"Now started job " << job << " in thread " << thread << endl;

    const int X = src.dimx();
    const int Y = src.dimy();
    const int Z = src.dimz();
    const int C = src.numChannels();
    const unsigned long int channelsize = src.channelSize();
    std::vector<double> px_val(C);
    for (int z = zstart; z != zend; ++z) {
	for (int y = 0; y < Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		fill(px_val.begin(), px_val.end(), 0);
		double accum = 0;
		int sign = 1;
		for (int pass = 0; pass != 2; ++pass) { // first forwards, then backwards
		    
		    double posx = x + 0.5; // position of current point on the traced line
		    double posy = y + 0.5; //                  - " -
		    double posz = z + 0.5; //                  - " -
		    double l = 0;
		    bool end_reached = false;
		    int prev_x = int(posx);
		    int prev_y = int(posy);
		    int prev_z = int(posz);
		    do {
			const int iposx = int(posx); // integer position
			const int iposy = int(posy); // integer position
			const int iposz = int(posz);
			const double w_x = sign * vec(iposx, iposy, iposz, 0); // roundoff to int
			const double w_y = sign * vec(iposx, iposy, iposz, 1); // roundoff to int
			const double w_z = sign * vec(iposx, iposy, iposz, 2); // roundoff to int
			double dl = steplength(posx, posy, posz, w_x, w_y, w_z); 
			if (dl < 0) break; // degeneracy detected
			l += dl; // l measures parametric length traversed
			if (l > L) {
			    dl = (l - L);
			    l = L;
			    // this will be the last cycle of updates.
			    end_reached = true;
			}
		    
			// contribution to integral
			// integrate over param.
			const double h = kernel_func(l) * dl;
			for (int c = 0; c < C; ++c) {
			    px_val[c] += h * src(iposx, iposy, iposz, c); 
			}
			accum += h;
			
			// updating position
			posx += dl * w_x;
			posy += dl * w_y;
			posz += dl * w_z;
			if (posx < 0 || posx > X) break; // outside image domain
			if (posy < 0 || posy > Y) break; // outside image domain
			if (posz < 0 || posz > Z) break; // outside image domain
			if (int(posx) == prev_x && 
			    int(posy) == prev_y &&
			    int(posz) == prev_z) break; // degen. point in field
		    prev_x = iposx; // NB: iposx and iposy have not yet been updated
		    prev_y = iposy; // from last iteration, so we here get the previous value
		    prev_z = iposz; //
		    } while (!end_reached);
		    sign *= -1;
		}
		if (accum > 0) {
		    for (int c = 0; c < C; ++c) 
			*(t_ptr + c * channelsize) += px_val[c] / accum;
		} else {
		    for (int c = 0; c < C; ++c) 
			*(t_ptr + c * channelsize) += src(x, y, z, c);
		}
		++t_ptr;
	    }
	}
    }
}

void LIC_2D_FS_subpart(const int job,
		       const int thread,
		       void * const idata,
		       void * const odata)
{
    // grabbing parameters
    LIC_FS_idata& input = ((LIC_FS_idata*)idata)[job];
    const int ystart = input.start;
    const int yend = input.end;
    const Image<double>& src = *(input.src);
    const Image<double>& vec = *(input.vec);
    double* t_ptr = input.target_ptr;
    const double L = input.L;
    const double dl = input.dl;
    double (*kernel_func)(double) = input.kernel_func;

    std::cout << "Now started job " << job << " in thread: " << thread << endl;
    
    const int X = src.dimx();
    const int Y = src.dimy();
    const int C = src.numChannels();
    const unsigned long int channelsize = src.channelSize();
    std::vector<double> px_val(C);
    for (int y = ystart; y != yend; ++y) {
	for (int x = 0; x < X; ++x) {
	    fill(px_val.begin(), px_val.end(), 0);
	    double accum = 0;
	    for (int pass = 0; pass != 2; ++pass) { // first time forwards, second time backwards
		double posx = x + 0.5; // position of current point on the traced line
		double posy = y + 0.5; //                 - " -
		//double l = 0;
		//bool end_reached = false;
		double prev_wx = vec(int(posx), int(posy), 0, 0);
		double prev_wy = vec(int(posx), int(posy), 0, 1);
		if (pass != 0) {
		    prev_wx *= -1;
		    prev_wy *= -1;
		}
		double magnitude = vec(int(posx), int(posy), 0, 2);
		const double traverse_length = L * magnitude;
		for (double l = 0; l < traverse_length; l+=dl) {
		    const int iposx = int(posx);
		    const int iposy = int(posy);
		    double wx = vec(iposx, iposy, 0, 0);
		    double wy = vec(iposx, iposy, 0, 1);
		    if (wx * prev_wx + wy * prev_wy < 0) {
			wx *= -1;
			wy *= -1;
		    }
		    // constribution to integral
		    const double h= kernel_func(l) * dl;
		    for (int c = 0; c < C; ++c) {
			px_val[c] += h * src(iposx, iposy, 0, c);
		    }
		    accum += h;
		    
		    // updating position
		    posx += dl * wx;
		    posy += dl * wy;
		    if (posx < 0 || posx >= X) break; // outside image domain
		    if (posy < 0 || posy >= Y) break; // outside image domain
		    prev_wx = wx;
		    prev_wy = wy;
		}
	    }
	    if (accum > 0) {
		for (int c = 0; c < C; ++c) {
		    *(t_ptr + c * channelsize) += px_val[c] / accum;
		}
	    } else {
		for (int c = 0; c < C; ++c) {
		    *(t_ptr + c * channelsize) = src(x, y, 0, c);
		}
	    }
	    ++t_ptr;
	}
    }
}

void LIC_3D_FS_subpart(const int job,
		       const int thread,
		       void * const idata,
		       void * const odata)
{
    // grabbing parameters
    LIC_FS_idata& input = ((LIC_FS_idata*)idata)[job];
    const int zstart = input.start;
    const int zend = input.end;
    const Image<double>& src = *(input.src);
    const Image<double>& vec = *(input.vec);
    double* t_ptr = input.target_ptr;
    const double L = input.L;
    const double dl = input.dl;
    double (*kernel_func)(double) = input.kernel_func;

    std::cout << "Now started job " << job << " in thread: " << thread << endl;
    
    const int X = src.dimx();
    const int Y = src.dimy();
    const int Z = src.dimz();
    const int C = src.numChannels();
    const unsigned long int channelsize = src.channelSize();
    std::vector<double> px_val(C);
    for (int z = zstart; z != zend; ++z) {
	for (int y = 0; y != Y; ++y) {
	    for (int x = 0; x < X; ++x) {
		fill(px_val.begin(), px_val.end(), 0);
		double accum = 0;
		for (int pass = 0; pass != 2; ++pass) { // first time forwards, second time backwards
		    double posx = x + 0.5; // position of current point on the traced line
		    double posy = y + 0.5; //                 - " -
		    double posz = z + 0.5; //                 - " -
		    //double l = 0;
		    //bool end_reached = false;
		    double prev_wx = vec(int(posx), int(posy), int(posz), 0);
		    double prev_wy = vec(int(posx), int(posy), int(posz), 1);
		    double prev_wz = vec(int(posx), int(posy), int(posz), 2);
		    if (pass != 0) {
			prev_wx *= -1;
			prev_wy *= -1;
			prev_wz *= -1;
		    }
		    double magnitude = vec(int(posx), int(posy), int(posz), 3); 
		    const double traverse_length = L * magnitude;
		    for (double l = 0; l < traverse_length; l+=dl) {
			const int iposx = int(posx);
			const int iposy = int(posy);
			const int iposz = int(posz);
			double wx = vec(iposx, iposy, iposz, 0);
			double wy = vec(iposx, iposy, iposz, 1);
			double wz = vec(iposx, iposy, iposz, 2);
			assert(wx >= -1 && wx <= 1);// @@
			assert(wy >= -1 && wy <= 1);// @@ remove for optimizing!
			assert(wz >= -1 && wz <= 1);// @@
			if (wx * prev_wx + wy * prev_wy + wz * prev_wz < 0) {
			    wx *= -1;
			    wy *= -1;
			    wz *= -1;
			}
			// constribution to integral
			const double h= kernel_func(l) * dl;
			for (int c = 0; c < C; ++c) {
			    px_val[c] += h * src(iposx, iposy, iposz, c);
			}
			accum += h;
		    
			// updating position
			posx += dl * wx;
			posy += dl * wy;
			posz += dl * wz;
			if (posx < 0 || posx >= X) break; // outside image domain
			if (posy < 0 || posy >= Y) break; // outside image domain
			if (posz < 0 || posz >= Z) break; // outside image domain
			prev_wx = wx;
			prev_wy = wy;
			prev_wz = wz;
		    }
		}
		if (accum > 0) {
		    for (int c = 0; c < C; ++c) {
			*(t_ptr + c * channelsize) += px_val[c] / accum;
		    }
		} else {
		    for (int c = 0; c < C; ++c) {
			*(t_ptr + c * channelsize) = src(x, y, z, c);
		    }
		}
		++t_ptr;
	    }
	}
    }
}


void LIC_2D_subpart(const int job, const int thread,  void * const idata,  void * const odata)
{
    // grabbing parameters
    LIC_idata& input = ((LIC_idata*)idata)[job];

    const int ystart = input.start;
    const int yend = input.end;
    const Image<double>& src = *(input.src);
    const Image<double>& vec = *(input.vec);
    double* t_ptr = input.target_ptr;
    const double L = input.L;
    double (*kernel_func)(double) = input.kernel_func;

    std::cout << "Now started job " << job << " in thread " << thread << endl;

    const int X = src.dimx();
    const int Y = src.dimy();
    const int C = src.numChannels();
    const unsigned long int channelsize = src.channelSize();
    std::vector<double> px_val(C);
    for (int y = ystart; y != yend; ++y) {
	for (int x = 0; x < X; ++x) {
	    fill(px_val.begin(), px_val.end(), 0);
	    double accum = 0;
	    int sign = 1;
	    for (int pass = 0; pass != 2; ++pass) {
		double posx = x + 0.5; // position of current point on the traced line
		double posy = y + 0.5; //                  - " -
		double l = 0;
		bool end_reached = false;
		int prev_x = int(posx);
		int prev_y = int(posy);
		do {
		    const int iposx = int(posx); // integer position
		    const int iposy = int(posy); // integer position
		    const double w_x = sign * vec(iposx, iposy, 0, 0); // roundoff to integer
		    const double w_y = sign * vec(iposx, iposy, 0, 1); // roundoff to integer
		    double dl = steplength(posx, posy, w_x, w_y); 
		    if (dl < 0) break; // degeneracy detected
		    l += dl; // l measures parametric length traversed
		    if (l > L) {
			dl = (l - L);
			l = L;
			// this will be the last cycle of updates.
			end_reached = true;
		    }
		    
		    // contribution to integral
		    // integrate over param.
		    const double h = kernel_func(l) * dl;
		    for (int c = 0; c < C; ++c) {
			px_val[c] += h * src(iposx, iposy, 0, c); 
		    }
		    accum += h;
		    
		    // updating position
		    posx += dl * w_x;
		    posy += dl * w_y;
		    if (posx < 0 || posx > X) break; // outside image domain
		    if (posy < 0 || posy > Y) break; // outside image domain
		    if (int(posx) == prev_x && int(posy) == prev_y) break; // degen. point in field
		    prev_x = iposx; // NB: iposx and iposy have not yet been updated
		    prev_y = iposy; // from last iteration, so we here get the previous value
		} while (!end_reached);
		sign *= -1;
	    }
	    if (accum > 0) {
		for (int c = 0; c < C; ++c) {
		    //cout << px_val[c] << " " << accum << endl;
		    *(t_ptr + c * channelsize) += px_val[c] / accum;
		}
	    } else {
		for (int c = 0; c < C; ++c) 
		    *(t_ptr + c * channelsize) = src(x, y, 0, c);
	    }
	    ++t_ptr;
	}
    }
}

// negative return value indicates degenerated point
double steplength(const double posx, 
			 const double posy,
			 const double vecx,
			 const double vecy)
{
    const double EPS = 1.0e-10;
    if (fabs(vecx) < EPS && fabs(vecy) < EPS) {
	// quasi-zero vector.  This cell will be considered a degenerated point in
	// the vector field
	return -1; // 
    }
    const double rem_x = int(posx) - posx;
    const double rem_y = int(posy) - posy;
    
    // vecx must be negative, angle 
    const double x_limit = TAN_LIM_ANGLE * fabs(vecy);
    const double y_limit = TAN_LIM_ANGLE * fabs(vecx);

    double sx_left, sx_right, sy_left, sy_right;
    if (vecx < -x_limit) {
	sx_left = rem_x / vecx;
	sx_right = INFINITE;
    } else if (vecx > x_limit) {
	sx_left = INFINITE;
	sx_right = (rem_x+1)/vecx;
    } else {
	sx_left = sx_right = INFINITE;
    }
    if (vecy < -y_limit) {
	sy_left = rem_y / vecy;
	sy_right = INFINITE;
    } else if (vecy > y_limit) {
	sy_left = INFINITE;
	sy_right = (rem_y+1)/vecy;
    } else {
	sy_left = sy_right = INFINITE;
    }

    const double min_x = (sx_left < sx_right) ? sx_left : sx_right;
    const double min_y = (sy_left < sy_right) ? sy_left : sy_right;
    assert(min_x >=0);
    assert(min_y >=0);

    const double res = min_x < min_y ? min_x : min_y;
    if (res == INFINITE) {
	// problem, consider this point a degeneracy
	return -1;
    }
    return res + ROUNDOFF_TERM;
}

inline double steplength(const double posx, 
			 const double posy,
			 const double posz,
			 const double vecx,
			 const double vecy,
			 const double vecz)
{
    const double EPS = 1.0e-10;
    if ((fabs(vecx) < EPS) && (fabs(vecy) < EPS) && (fabs(vecz) < EPS)) {
	// quasi-zero vector.  This cell will be considered a degenerated point in
	// the vector field
	return -1; // 
    }
    const double rem_x = int(posx) - posx;
    const double rem_y = int(posy) - posy;
    const double rem_z = int(posz) - posz;
    
    const double vecx2 = vecx*vecx;
    const double vecy2 = vecy*vecy;
    const double vecz2 = vecz*vecz;
    const double vnorm2 = vecx2 + vecy2 + vecz2;
    const double limit2 = vnorm2 * SIN2_LIM_ANGLE;

    double sx_left, sx_right, sy_left, sy_right, sz_left, sz_right;
    sx_left = sx_right = sy_left = sy_right = sz_left = sz_right = INFINITE;
    if (vecx2 > limit2) {
	if (vecx < 0) {
	    sx_left = rem_x / vecx;
	} else {
	    sx_right = (rem_x + 1)/vecx;
	}
    } 
    if (vecy2 > limit2) {
	if (vecy < 0) {
	    sy_left = rem_y / vecy;
	} else {
	    sy_right = (rem_y + 1) / vecy;
	}
    }
    if (vecz2 > limit2) {
	if (vecz < 0) {
	    sz_left = rem_z / vecz;
	} else {
	    sz_right = (rem_z + 1) / vecz;
	}
    }
    const double min_x = (sx_left < sx_right) ? sx_left : sx_right;
    const double min_y = (sy_left < sy_right) ? sy_left : sy_right;
    const double min_z = (sz_left < sz_right) ? sz_left : sz_right;
    assert(min_x >=0);
    assert(min_y >=0);
    assert(min_z >=0);

    double res = min_x < min_y ? min_x : min_y;
    res = (min_z < res) ? min_z : res;
    if (res == INFINITE) {
	// problem, consider this point a degeneracy
	return -1;
    }
    return res + ROUNDOFF_TERM;
}

}; // end anonymous namespace
