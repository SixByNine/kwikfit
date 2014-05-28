#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester

#include <math.h>
#include <stdio.h>
#include <mjklog.h>
#include "kwikfit.h"

double vonMises(double phase, kwikfit_component_t *component){

   return component->height * exp((cos(2.0*M_PI*(phase - component->phase)) -1)*component->concentration);

   //logmsg("cos: %lf cos-1 %lf",cos(2.0*M_PI*(phase-component->phase)),(cos(2.0*M_PI*(phase - component->phase)) -1));
   //logmsg("conc: %lf",component->concentration);
   //logmsg("p %lf c %lf h %lf r %lf",phase,component->phase,component->height,ret);
}

