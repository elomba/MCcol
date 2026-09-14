/*==============================================================================
 * File: libutil.c
 *
 * Purpose:
 *   Low-level C utility functions and POSIX system interfaces called by
 *   the Fortran Monte Carlo code gpMC.
 *
 * Capabilities:
 *   - Hostname inquiry: gethostname_()
 *   - Process user CPU time: second_()
 *   - Signal interception (SIGTERM, SIGINT): catch_() and handler()
 *     Interceps termination signals to trigger an orderly checkpoint dump
 *     via the Fortran subroutine cierra_().
 *
 * Notes on Fortran-C Calling Convention:
 *   Subroutine names include a trailing underscore (e.g. catch_(), gethostname_())
 *   which corresponds to the standard GNU/Intel Fortran name-mangling convention.
 *============================================================================*/

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <sys/times.h>
#include <signal.h>

/**
 * gethostname_ - Retrieves the system hostname.
 * @nom  : Character array pointer where hostname string will be written.
 * @ilen : Pointer to integer where hostname length will be stored.
 */
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

/**
 * second_ - Returns process user CPU time in seconds using times().
 * Returns: User CPU time as float.
 */
float second_(void)
{
    struct tms fecha;
    float tt1;

    times(&fecha);

    tt1 = (float) fecha.tms_utime / HZ;

    return tt1;
}

/* External Fortran subroutine: cierra(clean) defined in Dump.f90 */
extern void cierra_(int *clean);

/**
 * handler - Signal handler for SIGTERM and SIGINT.
 * Intercepts kill/interrupt signals and calls the Fortran checkpoint routine
 * cierra_(0) to dump the simulation state before exiting.
 */
static void handler(int sig)
{
    int clean = 0;

    printf(" *** signal %d intercepted: triggering emergency checkpoint...\n", sig);

    /* Re-establish signal handlers */
    signal(SIGTERM, handler);
    signal(SIGINT, handler);

    /* Call Fortran orderly shutdown routine */
    cierra_(&clean);
}

/**
 * catch_ - Registers POSIX signal handlers for SIGTERM and SIGINT.
 * Called from gpMC main initialization (Main.f90).
 */
void catch_(void)
{
    signal(SIGTERM, handler);
    signal(SIGINT, handler);
}
