/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* BLAS header filename. */
#define ATLAS_BLAS_H "cblas.h"

/* CLAPACK header filename. */
#define ATLAS_LAPACK_H "clapack.h"

/* BLAS header filename. */
/* #undef BLAS_H */

/* Archive creation date */
#define DATE "2022-04-27"

/* FFTW header filename. */
#define FFTW_H "fftw3.h"

/* Define to 1 if you have the `atexit' function. */
#define HAVE_ATEXIT 1

/* Define if you have the ATLAS libraries and header files. */
#define HAVE_ATLAS 1

/* Define if you have the BLAS library and header files. */
/* #undef HAVE_BLAS */

/* Define to 1 if you have the <cblas.h> header file. */
#define HAVE_CBLAS_H 1

/* Define to 1 if you have the <clapack.h> header file. */
#define HAVE_CLAPACK_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define if you have the FFTW double precision libraries and header files. */
/* #undef HAVE_FFTW */

/* Define if you have the FFTW single precision libraries and header files. */
#define HAVE_FFTWF 1

/* Define if you have the FFTW single precision multithreaded libraries and
   header files. */
/* #undef HAVE_FFTWF_MP */

/* Define if you have the FFTW double precision multithreaded libraries and
   header files. */
/* #undef HAVE_FFTW_MP */

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Define to 1 if you have the `getenv' function. */
#define HAVE_GETENV 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `isinf' function. */
#define HAVE_ISINF 1

/* Define to 1 if you have the `isnan' function. */
#define HAVE_ISNAN 1

/* Define if you have the LAPACK libraries. */
/* #undef HAVE_LAPACK */

/* Define if you have the LAPACKe library and header files. */
/* #undef HAVE_LAPACKE */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `logf' function. */
#define HAVE_LOGF 1

/* Define to 1 if the system has the type `long long int'. */
#define HAVE_LONG_LONG_INT 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <mathimf.h> header file. */
/* #undef HAVE_MATHIMF_H */

/* Define to 1 if you have the `memcpy' function. */
#define HAVE_MEMCPY 1

/* Define to 1 if you have the `memmove' function. */
#define HAVE_MEMMOVE 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have the `mkdir' function. */
#define HAVE_MKDIR 1

/* Define if you have the MKL libraries. */
/* #undef HAVE_MKL */

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the `munmap' function. */
#define HAVE_MUNMAP 1

/* Define if you have the OpenBLAS library and header files. */
/* #undef HAVE_OPENBLAS */

/* Define if you have the OpenBLAS parallel library and header files. */
/* #undef HAVE_OPENBLASP */

/* Define to 1 if you have the `posix_memalign' function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define if you have POSIX threads libraries and header files. */
/* #undef HAVE_PTHREAD */

/* Define to 1 if you have the `setlinebuf' function. */
#define HAVE_SETLINEBUF 1

/* Define to 1 if you have the `sincosf' function. */
#define HAVE_SINCOSF 1

/* Define to 1 if `stat' has the bug that it succeeds when given the
   zero-length file name argument. */
/* #undef HAVE_STAT_EMPTY_STRING_BUG */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strftime' function. */
#define HAVE_STRFTIME 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strstr' function. */
#define HAVE_STRSTR 1

/* Define to 1 if you have the `sysconf' function. */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type `unsigned long long int'. */
#define HAVE_UNSIGNED_LONG_LONG_INT 1

/* LAPACKe header filename. */
/* #undef LAPACKE_H */

/* LAPACK header filename. */
/* #undef LAPACK_H */

/* Define to 1 if `lstat' dereferences a symlink specified with a trailing
   slash. */
#define LSTAT_FOLLOWS_SLASHED_SYMLINK 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* MKL header filename. */
/* #undef MKL_H */

/* Name of package */
#define PACKAGE "sextractor"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "astromatic@astromatic.net"

/* Define to the full name of this package. */
#define PACKAGE_NAME "SExtractor"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "SExtractor 2.24.4"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "sextractor"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://astromatic.net/software/sextractor"

/* Define to the version of this package. */
#define PACKAGE_VERSION "2.24.4"

/* Define to the necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define as the return type of signal handlers (`int' or `void'). */
#define RETSIGTYPE void

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Maximum number of POSIX threads */
#define THREADS_NMAX 1024

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* Triggers model-fitting */
#define USE_MODEL 1

/* Triggers multhreading */
/* #undef USE_THREADS */

/* Version number of package */
#define VERSION "2.24.4"

/* Default URL of the XSLT filter */
#define XSL_URL "file:///usr/local/share/sextractor/sextractor.xsl"

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
/* #undef _LARGEFILE_SOURCE */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `int' if <sys/types.h> doesn't define. */
/* #undef gid_t */

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> doesn't define. */
/* #undef uid_t */
