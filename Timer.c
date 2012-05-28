#include <stdlib.h>
#include <sys/time.h>

struct itimerval now;

void resettimer()
{
  now.it_interval.tv_sec=100000000L;
  now.it_interval.tv_usec=0;
  now.it_value.   tv_sec=100000000L;
  now.it_value.tv_usec=0;
  setitimer(ITIMER_REAL,&now,NULL);
}

/*msec timer*/
int deltatime()
{
  int delta;
  struct itimerval now2;
  getitimer(ITIMER_REAL,&now2);
  delta = now.it_value.tv_usec - now2.it_value.tv_usec;
  delta = delta/1000 + (now.it_value.tv_sec - now2.it_value.tv_sec)*1000;
  /*itimer is a decremental timer*/
  /*reset interval*/
  now=now2;
  return delta;
}
