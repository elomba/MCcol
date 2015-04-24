#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <sys/times.h>
#include <signal.h>
void gethostname_(char* nom, int *ilen){
  int iret, lnt;
  size_t len;
  len = 20;
  iret = gethostname(nom,len) ;
  if ( iret != 0) {
    printf("Error return from hostname with code %i \n", iret);}
  for (lnt = 0 ; *nom != '\0'; nom++)
    lnt++;
  *ilen = lnt;
}

/*
    Use the following include in AIX 3.1
*/
/*#include <sys/m_param.h>*/
/*
   For Ultrix 4.2 uncomment the following definition instead
   checking the proper clock rate.
*/
#define HZ 100
float second_(void)
{
  struct tms *Buffer;
  struct tms fecha;
  float tiempo, tt1,tt2;
  Buffer = &fecha;
  tiempo = times(Buffer);
  tt1 = Buffer->tms_utime;
  tt2 = Buffer->tms_stime;
  tt1 = tt1/HZ;
  tt2 = tt2/HZ;
  return tt1;
}

/* This subroutine catches a SIGTERM kill signal. */
/* Taken from "Using C on the UNIX System", D.A. Curry, O'Reilly &Assoc */

void catch_ ()

{

extern int handler () ;
 
 signal(SIGTERM, handler) ;
 signal(SIGINT, handler) ;
}

handler ()

{ 
  extern void cierra_(int *clean) ;
  int clean ;
  clean = 0;
  printf (" *** kill signal intercepted \n");
  signal (SIGTERM, handler) ;
  signal (SIGINT, handler) ;
  cierra_(&clean) ;
}

