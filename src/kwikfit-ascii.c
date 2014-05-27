#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester


#include <stdio.h>
#include <string.h>
#include "mjk_cmd.h"
#include "mjklog.h"
#include "kwikfit.h"

int main (int argc, char** argv){
   debugFlag=1;
   logmsg("TEST");
   logdbg("DEBUG");
   char *template_fname = getS("--template","-t",argc,argv,"example.kwtemplate");
   logmsg("reading '%s'",template_fname);
   kwikfit_template_t *template = kwikfit_read_template(template_fname);
   kwikfit_write(template,stdout);

return 0;
}
