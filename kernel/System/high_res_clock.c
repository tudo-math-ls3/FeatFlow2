/*
 * High-Resolution Clock for FeatFlow 2
 * ------------------------------------
 * This file contains wrapper functions for the high-resolution timer
 * functions offered by different operating systems.
 *
 *
 * There are 2 functions implemented in this file:
 *
 * -> hrc_stamp
 * This function retrieves the current time stamp.
 *
 * -> hrc_diff
 * This function calculates the elapsed time between two time stamps,
 * measured in seconds.
 *
 *
 * Corresponding Fortran 90 interfaces:
 *
 *   subroutine hrc_stamp(_time)
 *     integer*8, intent(out) :: _time
 *   end subroutine
 *
 *   subroutine hrc_diff(_elapsed, _time1, _time2)
 *     real*8, intent(out) :: _elapsed
 *     integer*8, intent(in) :: _time1, _time2
 *   end subroutine
 */

#if defined(_WIN32)

/* Windows API */

/* prototype declarations */
int __stdcall QueryPerformanceCounter(long long int*);
int __stdcall QueryPerformanceFrequency(long long int*);

void hrc_stamp(long long int* _time)
{
  if(QueryPerformanceCounter(_time) == 0)
    *_time = 0; /* no high-resolution clock available */
}

void hrc_diff(double* _elapsed, long long int* _time1, long long int* _time2)
{
  long long int freq;

  /* query clock frequency and calculate elapsed time */
  if(QueryPerformanceFrequency(&freq) != 0)
    *_elapsed = (double)(*_time2 - *_time1) / (double)freq;
  else
    *_elapsed = 0.0; /* no high-resolution clock available */
}

#elif defined(unix) || defined (__unix) || defined(__unix__)

/* Unix API */

#include <sys/time.h>
#include <stdlib.h>

void hrc_stamp(long long int* _time)
{
  struct timeval stamp;
  gettimeofday(&stamp, NULL);
  *_time = 1000000ll*((long long int)stamp.tv_sec) + (long long int)stamp.tv_usec;
}

void hrc_diff(double* _elapsed, long long int* _time1, long long int* _time2)
{
  *_elapsed = 1E-6 * (double)(*_time2 - *_time1);
}

/* Wrapper functions with trailing underscores following */
void hrc_stamp_(long long int* _time)
{
  hrc_stamp(_time);
}

void hrc_stamp__(long long int* _time)
{
  hrc_stamp(_time);
}

void hrc_diff_(double* _elapsed, long long int* _time1, long long int* _time2)
{
  hrc_diff(_elapsed, _time1, _time2);
}

void hrc_diff__(double* _elapsed, long long int* _time1, long long int* _time2)
{
  hrc_diff(_elapsed, _time1, _time2);
}

#endif
