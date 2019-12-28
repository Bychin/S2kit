/**
 * @file csecond.c
 * @brief Timing function.
 * @author Mark Taylor
 */

#include "csecond.h"

#include <limits.h>
#include <stdio.h>
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>

#ifdef CLK_TCK
#define DIVIDER CLK_TCK
#else
#define DIVIDER HZ // for old BSD systems
#endif

// define this to return wallclock instead of cpu time
// #define WALLCLOCK

#if (defined CRAY)
double CSECOND() {
#elif (defined T3D)
double CSECOND() {
#elif (defined IBM)
double csecond() {
#else // works for SUN, LINUX, DECALPHA, SGI
double csecond() {
#endif
    static int firstcall = 1;
    static struct tms buf0; // times structure
    static clock_t rv0;

    if (firstcall) {
        firstcall = 0;
        rv0 = times(&buf0);
    }

    struct tms buf;
    clock_t rv = times(&buf);
#ifdef WALLCLOCK
    return ((double)(rv - rv0) / (double)DIVIDER);
#else
    return (
        (double)(
            (buf.tms_utime + buf.tms_stime + buf.tms_cutime + buf.tms_cstime) -
            (buf0.tms_utime + buf0.tms_stime + buf0.tms_cutime + buf0.tms_cstime)
        ) / (double)DIVIDER
    );
#endif
}
