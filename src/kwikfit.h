//  Copyright (C) 2014 Michael J. Keith, University of Manchester

#ifndef __kwikfit_h
#define __kwikfit_h
#include <string.h>
#include <stdio.h>
#ifndef STREQ
#define STREQ(a,b) (strcmp(a,b)==0)
#endif

#include <inttypes.h>

extern int debugFlag; // tempo2 goodness
#ifdef __cplusplus
extern "C" {
#endif

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

   typedef struct kwikfit_result {
	  double phase;
	  double error;
	  double chisq;
	  uint64_t nbins;
	  double *data;
	  double *fit;
	  double *residual;
	  double *data_cov;
	  kwikfit_template_t *tmpl;
	  double *amplitudes;
	  double **amp_cov;
   } kwikfit_result_t;

kwikfit_template_t *kwikfit_read_template(char* filename);
void kwikfit_write(kwikfit_template_t* tmpl, FILE* out);
double vonMises(double phase, kwikfit_component_t *component);
kwikfit_result_t *kwikfit_doFit(uint64_t nbins, double *profile, kwikfit_template_t *tmpl, uint64_t subbin_res);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif // __kwikfit_h
