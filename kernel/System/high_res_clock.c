#ifdef _WIN32
/*
 * High-Resolution Clock for FeatFlow 2 on Windows
 * -----------------------------------------------
 * This file contains wrapper functions for the high-resolution timer
 * functions offered by the Windows API.
 */


/* prototype declarations */
int __stdcall QueryPerformanceCounter(long long*);
int __stdcall QueryPerformanceFrequency(long long*);

/* Corresponding Fortran interface:
 *   subroutine hrc_stamp(_time)
 *     integer*8, intent(out) :: _time
 *   end subroutine
 */
void hrc_stamp(long long* _time)
{
  QueryPerformanceCounter(_time);
}

/* Corresponding Fortran interface:
 *   subroutine hrc_stamp2(_low, _high)
 *     integer*4, intent(out) :: _low, _high
 *   end subroutine
 */
void hrc_stamp2(int* _low, int* _high)
{
  long long stamp = 0ll;
  if(QueryPerformanceCounter(&stamp) != 0)
  {
    *_low  = (int)( stamp        & 0xFFFFFFFFll);
    *_high = (int)((stamp >> 32) & 0xFFFFFFFFll);
  }
  else
  {
    /* query failed: no high-resolution clock available on platform (?) */
    *_low = *_high = 0;
  }
}

/* Corresponding Fortran interface:
 *   subroutine hrc_diff(_elapsed, _time1, _time2)
 *     real*8, intent(out) :: _elapsed
 *     integer*8, intent(in) :: _time1, _time2
 *   end subroutine
 */
void hrc_diff(double* _elapsed, long long* _time1, long long* _time2)
{
  long long freq;

  /* query clock frequency and calculate elapsed time */
  if(QueryPerformanceFrequency(&freq) != 0)
    *_elapsed = (double)(*_time2 - *_time1) / (double)freq;
  else
    *_elapsed = 0.0;
}


/* Corresponding Fortran interface:
 *   subroutine hrc_diff(_elapsed, _time1_l, _time1_h, _time2_l, _time2_h)
 *     real*8, intent(out) :: _elapsed
 *     integer*4, intent(in) :: _time1_l, _time1_h, _time2_l, _time2_h
 *   end subroutine
 */
void hrc_diff2(double* _elapsed, int* _time1_l, int* _time1_h, int* _time2_l, int* _time2_h)
{
  long long time1, time2, freq;

  /* recombine low and high order ints to long long */
  time1 = (((long long)*((unsigned int*)_time1_h)) << 32) | *((unsigned int*)_time1_l);
  time2 = (((long long)*((unsigned int*)_time2_h)) << 32) | *((unsigned int*)_time2_l);

  /* query clock frequency and calculate elapsed time */
  if(QueryPerformanceFrequency(&freq) != 0)
    *_elapsed = (double)(time2 - time1) / (double)freq;
  else
    *_elapsed = 0.0;
}

#endif /* _WIN32 */
