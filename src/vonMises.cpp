#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester

#include <math.h>
#include "kwikfit.h"

double vonMises(double phase, kwikfit_component_t *component){
   return component->height * exp((cos(phase - component->phase) -1)*component->concentration);
}

