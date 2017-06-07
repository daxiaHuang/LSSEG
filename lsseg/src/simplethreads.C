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
//#ifdef tor
//#  include <vector>
//#else
//#  include <vector.h>
//#endif

/// \file
/// \brief Implements simplethreads.h.
#include <vector>

#include "simplethreads.h"






//
// This module is more or less stable, so let's force this off.
//
#undef DEBUG  
#undef DEBUG_HEAVY





/// \cond
struct thread_data_struct
/// \endcond
{
  int thread_id;
  pthread_mutex_t *mutex;
  char *job_flag;		// Pointer to array of flags, where
				// 0=not started on, 1=done, or being done.
  int *jobs_done;
  int jobs;
  void (*do_a_job)(const int, const int, void * const, void * const);
  void *individual_job_data;	// Note that this is an array of structs,
				// which 'do_a_job' must be able to correctly
				// handle! If this is NULL, it means that there
				// is no individual job data, which might be
				// the case for simple jobs.
  void *results;		// The same applies here. This is an array
  				// which do_a_job must know how to use.
                                // Note that neither of 'individual_job_data'
				// and 'results' is locked/unlocked, so jobs
				// should not share their elemnts of these
  				// arrays with each other!
  bool show_progress; // 051101
};






//======================================================================
//
// This is not meant to be called directly by the user, instead the
// user supplies the 'do_a_job' function to 'do_threaded_jobs'.
//
// In this function, a thread tries to do as many "jobs" as possible,
// until all jobs are done (by itself or other threads.)
//
// 050620: So, if I understand my own comment correctly, 
//         'do_threaded_jobs' is rather a wrapper around 'do_jobs'...
//
//======================================================================

void *do_jobs(void *thread_data_ptr)
{
  //
  // Some convenient aliases:
  //
  const int &thread_id=((thread_data_struct *)thread_data_ptr)->thread_id;
  int *jobs_done=((thread_data_struct *)thread_data_ptr)->jobs_done;
  int jobs=((thread_data_struct *)thread_data_ptr)->jobs;
  void * const results=((thread_data_struct *)thread_data_ptr)->results;
  pthread_mutex_t *mutex=((thread_data_struct *)thread_data_ptr)->mutex;
  char *job_flag=((thread_data_struct *)thread_data_ptr)->job_flag;
  void (*do_a_job)(const int, const int, void * const, void * const)=
    ((thread_data_struct *)thread_data_ptr)->do_a_job;
  void *individual_job_data=
    ((thread_data_struct *)thread_data_ptr)->individual_job_data;
  bool show_progress=((thread_data_struct *)thread_data_ptr)->show_progress;

#ifdef DEBUG
  printf("simplethreads: do_jobs called for thread %d.\n", thread_id);
#endif
  
  //
  // Lock all shared data, so that it is not changed by another thread
  // while we look for a new job.
  //
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) waiting for lock.\n", thread_id);
#endif
  int err;
  if ((err=pthread_mutex_lock(mutex))!=0)
    CRIT_ERR(printf("mutex_lock error: %d (thread=%d)\n", err, thread_id));
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) got lock.\n", thread_id);
#endif

  int stride=std::max(jobs/20, 1);

  //
  // Look for a job, if not all are done. This we will continue to do,
  // as long as there are jobs to be done. Note the locking of the
  // shared data across the scope boundaries both at the start and end
  // of this while block!
  //
  while (*jobs_done<jobs)
    {
      int i;
      for (i=0; (i<jobs) && (job_flag[i]!=0); )
	i++;
      if (i==jobs)
	CRIT_ERR(puts("Shouldn't happen!"));

      //
      // Ok, we have a job to do, let's do it, after we have updated
      // the job_flag array, and released it, in effect letting the
      // other threads know that *we* are the ones doing this
      // particular job.
      //
      job_flag[i]=1; // Meaning "done, or being done".
      (*jobs_done)++;
      if ((show_progress) && ((*jobs_done)%stride==0))
	printf("%d ", (jobs-(*jobs_done))/stride), fflush(stdout);
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) unlocking.\n", thread_id);
#endif
      if ((err=pthread_mutex_unlock(mutex))!=0)
	CRIT_ERR(printf("mutex_unlock err: %d (thread %d)\n", err, thread_id));
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) unlocked.\n", thread_id);
#endif
      
      //
      // Do the work!
      //
      do_a_job(i, thread_id, individual_job_data, results);

      //
      // Now, since we entered the scope with shared data locked, and
      // we have to do that again, for the next "round" of this loop,
      // we must now lock the data.
      //
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) waiting for lock 2.\n", thread_id);
#endif
      int err;
      if ((err=pthread_mutex_lock(mutex))!=0)
	CRIT_ERR(printf("mutex_lock error: %d (thread=%d)\n", err, thread_id));
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) got lock 2.\n", thread_id);
#endif
    }
  if (show_progress)
    printf("\n");
  
  //
  // No more work to do. Not much to do than return. But we should
  // remember to unlock the data, since we left the loop with the data
  // locked!
  //
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) unlocking 2.\n", thread_id);
#endif
  if ((err=pthread_mutex_unlock(mutex))!=0)
    CRIT_ERR(printf("mutex_unlock err: %d (thread %d)\n", err, thread_id));
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_jobs(%d) unlocked 2.\n", thread_id);
#endif

  return NULL;

  //
  // Why shouldn't this return an int instead, since the value shall
  // be used as an "exit code"??? Is this some kind of bug, or relic,
  // in the pthreads library?
  //
}






void do_threaded_jobs(void (*do_a_job)(const int,
				       const int,
				       void * const,
				       void * const),
		      void * job_data_struct_array,
		      const int threads,
		      const int jobs,
		      const bool show_progress,
		      void * const results)
{
  // printf("Initializing %d threads for %d jobs.\n", threads, jobs);

  // 0=job not started, 1=done or being processed.
  std::vector<char> job_flag(jobs, 0);
  int jobs_done=0;

  // Used to lock the job-flag array by threads accessing it.
  pthread_mutex_t trace_admin_mutex;
  int err;
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_threaded_jobs initializing lock.\n");
#endif
  if ((err=pthread_mutex_init(&trace_admin_mutex, NULL))!=0)
    CRIT_ERR(printf("pthread_mutex_init error: %d\n", err));
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_threaded_jobs initialized lock.\n");
#endif

  //
  // Setting up a data-package for each thread.  Ok, there is some
  // unneccessary duplication of data here. Will fix some other time.
  //
  std::vector<thread_data_struct> thread_data(threads);
  int i;
  for (i=0; i<threads; i++)
    {
      thread_data[i].thread_id=i;
      thread_data[i].mutex=&trace_admin_mutex;
      thread_data[i].job_flag=&job_flag[0];
      thread_data[i].jobs_done=&jobs_done;
      thread_data[i].jobs=jobs;
      thread_data[i].do_a_job=do_a_job;
      thread_data[i].individual_job_data=job_data_struct_array;
      thread_data[i].results=results;
      thread_data[i].show_progress=show_progress;
    }

  //
  // Forking off the threads. They will immediately start working on
  // the individual jobs, and they will all return when there are no
  // more jobs to do.
  //
  std::vector<pthread_t> thread(threads);
  for (i=0; i<threads; i++)
    {
#ifdef DEBUG_HEAVY
      printf("simplethreads: do_threaded_jobs creating thread %d.\n", i);
#endif
      if ((err=pthread_create(&thread[i],
			      NULL,
			      do_jobs,
			      (void *)&(thread_data[i])))!=0)
	CRIT_ERR(printf("pthread_create error: %d\n", err));
#ifdef DEBUG_HEAVY
      printf("simplethreads: do_threaded_jobs created thread %d.\n", i);
#endif
    }

  //
  // Now we simply wait for all the threads to finish.
  //

//    {
//      int i;
//      double x=0.0;
//      for (i=0; i<1000000000; i++)
//        {
//  	x=sin(x+0.01);
//  	if (i%1000000==0) printf("z=%f\n", x);
//        }
//      puts("z ferdig");
//      for (i=0; i<1000000000; i++)
//        {
//  	x=sin(x+0.01);
//  	if (i%1000000==0) printf("zz=%f\n", x);
//        }
//      puts("zz ferdig");
//    }



  for (i=0; i<threads; i++)
    {
      void *retval;
#ifdef DEBUG
      printf("simplethreads: do_threaded_jobs joining thread %d.\n", i);
#endif
      if ((err=pthread_join(thread[i], &retval))!=0)
	CRIT_ERR(printf("pthread_join error on thread %d: %d\n", i, err));
#ifdef DEBUG
      printf("simplethreads: do_threaded_jobs joined thread %d.\n", i);
#endif
    }
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_threaded_jobs destroying mutex.\n");
#endif
  if ((err=pthread_mutex_destroy(&trace_admin_mutex))!=0)
    CRIT_ERR(printf("pthread_mutex_destroy error: %d\n", err));
#ifdef DEBUG_HEAVY
  printf("simplethreads: do_threaded_jobs destroyed mutex.\n");
#endif

  //
  // And that was it.
  //
  // puts("All threads now joined and terminated.");
}






//======================================================================
//
// 050620: Starting a new version which will use "condition variables".
//         Don't know how I could have missed this. Or was it added to 
//         pthreads after I made the first version many years ago?
//         Or maybe it was not in 'linuxthreads'???
//
//         I think my first implementation does exactly what it should,
//         and what usage of condition variables will do, but maybe
//         with a small (large?!) overhead due to the constant polling
//         in the "look for job" loop... Should be interesting to test
//         after the new version is working. (Maybe with the raytracing
//         code?)
//
//         Are the three classic examples, "producer/consumer",
//         "reader/writer" and "dining philosophers" really the same?
//         If not, which one (if any) do we have here?
//         "Reader/writer"?
//
//         Hmm... When thinking more about this, maybe there is no use
//         for condition variables in just this case... Think maybe not...
//
//======================================================================
