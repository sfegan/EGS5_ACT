// egs5_futils.h - Helper routines for interface with FORTRAN
// Stephen Fegan - sfegan@llr.in2p3.fr - July 2012
// $Id: egs5_futils.h 4530 2012-09-27 08:42:49Z sfegan $

#ifndef EGS5_FUTILS
#define EGS5_FUTILS

#ifdef  __cplusplus
extern "C" {
#endif

void fort_open_old_(int* theunit, const char* thename, int thenamelen);
void fort_open_new_(int* theunit, const char* thename, int thenamelen);
void fort_open_unknown_(int* theunit, const char* thename, int thenamelen);
void fort_close_(int* theunit);
void sjf_test_(void);

#ifdef  __cplusplus
} // extern "C"
#endif

#endif // !defined EGS5_FUTILS
