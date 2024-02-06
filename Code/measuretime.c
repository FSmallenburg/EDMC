#include <sys/time.h>


double getelapsedtime()
{
  static struct timeval t0;
  static int first = 1;
  
  struct timeval t1;
  if (first) 
  {
    first = 0;
    gettimeofday(&t0, NULL);
  }
  gettimeofday(&t1, NULL);
  double elapsedTime = (t1.tv_sec  - t0.tv_sec) ;      // sec to ms
  elapsedTime       += (t1.tv_usec - t0.tv_usec) / 1000000.0;   // us to ms 
  return elapsedTime;
}