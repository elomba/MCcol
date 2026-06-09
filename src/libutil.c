#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <sys/times.h>
#include <signal.h>

void gethostname_(char *nom, int *ilen)
{
    int iret, lnt;
    size_t len = 20;

    iret = gethostname(nom, len);

    if (iret != 0) {
        printf("Error return from hostname with code %i\n", iret);
    }

    for (lnt = 0; *nom != '\0'; nom++)
        lnt++;

    *ilen = lnt;
}

#define HZ 100

float second_(void)
{
    struct tms fecha;
    float tt1;

    times(&fecha);

    tt1 = (float) fecha.tms_utime / HZ;

    return tt1;
}

/* Fortran routine */
extern void cierra_(int *clean);

/* Signal handler */
static void handler(int sig)
{
    int clean = 0;

    printf(" *** signal %d intercepted\n", sig);

    signal(SIGTERM, handler);
    signal(SIGINT, handler);

    cierra_(&clean);
}

void catch_(void)
{
    signal(SIGTERM, handler);
    signal(SIGINT, handler);
}
