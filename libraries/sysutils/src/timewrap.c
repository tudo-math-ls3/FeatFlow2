#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

int getrusage (int, struct rusage *);             

long twrap_( long* c, long* r )
{
  struct rusage usage;
    
  if( getrusage(  RUSAGE_SELF, &usage ) == -1 )
    perror("Time");

  *r =  usage.ru_utime.tv_sec*1000 + usage.ru_utime.tv_usec/1000;
  *c =  1000;

  return 0;
}
