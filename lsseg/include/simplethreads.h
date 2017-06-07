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
#ifndef THREADSTUFF_H_INCLUDED

#include <stdio.h>
//#ifdef JON_PTHR
#  include <pthread.h>
//#endif

#ifndef CRIT_ERR
#  define CRIT_ERR(stmnt) \
    printf("\nIn file %s, line %d:\n  ", __FILE__, __LINE__), (stmnt), exit(0)
#endif

/// \file
/// \brief File containing functionality for multiple threads.  Can be used for 
/// parallellization of computationally-heavy routine.  It is written by Jens-Olav
/// Nygaard, and is not part of \ref lsseg "the lsseg library" as such.


//
// The function passed through the first parameter is the function
// actually doing the work, that is, the work of a given 'job'. This
// function is to be supplied by the user.
//
// The function must process the data in the structure pointed to by
// the second parameter, and also return its results as members in
// this structure, or placed wherever this structure's data members
// tells it to.
//
// The function is passed an integer, the id of the job, and a pointer
// to the datastructure belonging to this job.
//
// Typically, that structure contains all the information specific to
// a particular job, and the function will place results in a location
// specified in this structure, which will not be used by any other
// thread.
//
// This means that enough memory for storing the results of all jobs
// must be available when 'do_all_jobs' is called.
//

// Overview: A number of threads will be started. Each of them will
// look for jobs not done yet, do it, report it done, and go on with
// another job, until all jobs are done. Then the threads are merged,
// and control is returned to the caller of 'do_threaded_jobs'.

// Note that a seemingly (perhaps) superfluous level of data
// structures is made use of. This is simply to be able to present to
// the user of this library an API in which the user makes a
// 'do_a_job' routine, instead of actually supplying the
// 'do_the_jobs_of_a_thread' routine, thus hiding almost all
// job-server/client, or thread-, specific stuff...

// A tip: When using this library: Make sure to maintain a
// non-threaded version, it will come in very handy when you need a
// fallback solution for the time when the threaded stuff gets broken,
// and you have no clue as to why!

//
// 020618: Missing information: Will 'do_a_job' get a pointer to
//         the base of 'job_data_struct_array', or directly to its
//         own member of that list?
//         Same question applies to the output.
//
//         Answer: Seems like both get the base of the array(s).
//         Is this really a good idea?
//




void do_threaded_jobs(void (*do_a_job)(const int,	// job
				       const int,	// thread
				       void * const,	// input data
				       void * const),	// output data
		      void * job_data_struct_array,
		      const int threads,
		      const int jobs,
		      const bool show_progress,
		      void * const results);






#define THREADSTUFF_H_INCLUDED
#endif
