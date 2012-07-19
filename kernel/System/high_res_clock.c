/*
 * High-Resolution Clock for FeatFlow 2
 * -----------------------------------------------
 * This file contains wrapper functions for the high-resolution timer
 * functions offered by different operating systems.
 */

#if defined(_WIN32)

/* Windows API */

/* prototype declarations */
int __stdcall QueryPerformanceCounter(long long int*);
int __stdcall QueryPerformanceFrequency(long long int*);

/* Corresponding Fortran interface:
 *   subroutine hrc_stamp(_time)
 *     integer*8, intent(out) :: _time
 *   end subroutine
 */
void hrc_stamp(long long int* _time)
{
  QueryPerformanceCounter(_time);
}

/* Corresponding Fortran interface:
 *   subroutine hrc_diff(_elapsed, _time1, _time2)
 *     real*8, intent(out) :: _elapsed
 *     integer*8, intent(in) :: _time1, _time2
 *   end subroutine
 */
void hrc_diff(double* _elapsed, long long int* _time1, long long int* _time2)
{
  long long int freq;

  /* query clock frequency and calculate elapsed time */
  if(QueryPerformanceFrequency(&freq) != 0)
    *_elapsed = (double)(*_time2 - *_time1) / (double)freq;
  else
    *_elapsed = 0.0;
}

#elif defined(__unix__)

/* Unix API */

#include <sys/time.h>
#include <stdlib.h>

/* Corresponding Fortran interface:
 *   subroutine hrc_stamp(_time)
 *     integer*8, intent(out) :: _time
 *   end subroutine
 */
void hrc_stamp(long long int* _time)
{
  struct timeval stamp;
  gettimeofday(&stamp, NULL);
  *_time = 1000000*((long long int)stamp.tv_sec) + (long long int)stamp.tv_usec;
}

/* Corresponding Fortran interface:
 *   subroutine hrc_stamp(_time)
 *     integer*8, intent(out) :: _time
 *   end subroutine
 */
void hrc_stamp_(long long int* _time)
{
  struct timeval stamp;
  gettimeofday(&stamp, NULL);
  *_time = 1000000*((long long int)stamp.tv_sec) + (long long int)stamp.tv_usec;
}

/* Corresponding Fortran interface:
 *   subroutine hrc_stamp(_time)
 *     integer*8, intent(out) :: _time
 *   end subroutine
 */
void hrc_stamp__(long long int* _time)
{
  struct timeval stamp;
  gettimeofday(&stamp, NULL);
  *_time = 1000000*((long long int)stamp.tv_sec) + (long long int)stamp.tv_usec;
}

/* Corresponding Fortran interface:
 *   subroutine hrc_diff(_elapsed, _time1, _time2)
 *     real*8, intent(out) :: _elapsed
 *     integer*8, intent(in) :: _time1, _time2
 *   end subroutine
 */
void hrc_diff(double* _elapsed, long long int* _time1, long long int* _time2)
{
  *_elapsed = (double)(*_time2 - *_time1);
}

/* Corresponding Fortran interface:
 *   subroutine hrc_diff(_elapsed, _time1, _time2)
 *     real*8, intent(out) :: _elapsed
 *     integer*8, intent(in) :: _time1, _time2
 *   end subroutine
 */
void hrc_diff_(double* _elapsed, long long int* _time1, long long int* _time2)
{
  *_elapsed = (double)(*_time2 - *_time1);
}

/* Corresponding Fortran interface:
 *   subroutine hrc_diff(_elapsed, _time1, _time2)
 *     real*8, intent(out) :: _elapsed
 *     integer*8, intent(in) :: _time1, _time2
 *   end subroutine
 */
void hrc_diff__(double* _elapsed, long long int* _time1, long long int* _time2)
{
  *_elapsed = (double)(*_time2 - *_time1);
}

#endif 
