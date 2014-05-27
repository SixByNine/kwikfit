//  Copyright (C) 2014 Michael J. Keith, University of Manchester

#ifndef __kwikfit_h
#define __kwikfit_h
#include <string.h>
#ifndef STREQ
#define STREQ(a,b) (strcmp(a,b)==0)
#endif

#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

   extern int debugFlag; // tempo2 goodness
   // begin header

   typedef struct kwikfit_component {
	  double phase;
	  double height;
	  double concentration;
	  char scale_free;
   } kwikfit_component_t;

   typedef struct kwikfit_profile {
	  uint64_t ncomp;
	  kwikfit_component_t *comps;
	  char name[128];
   } kwikfit_profile_t;

   typedef struct kwikfit_template {
	  uint64_t nprof;
	  kwikfit_profile_t *profs;
	  char name[128];
   } kwikfit_template_t;


kwikfit_template_t *kwikfit_read_template(char* filename);
void kwikfit_write(kwikfit_template_t* template, FILE* out);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // __kwikfit_h
