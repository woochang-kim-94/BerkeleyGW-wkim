/******************************************************************************
*
*  os_c                 Originally by FHJ     Last Modified 12/04/2012 (FHJ)
*
*    BerkeleyGW wrapper for some libc and operating-system routines.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

/* Enable POSIX extensions. Only required because we use gcc -std=c99 */
#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <sys/types.h>
#include <time.h>

int sleep_c(double s)
{
	struct timespec rqtp;
	rqtp.tv_sec = (time_t) s;
	rqtp.tv_nsec = (long)((double)(s - (time_t)s) * 1.0e9);
	return nanosleep(&rqtp, NULL);
}
