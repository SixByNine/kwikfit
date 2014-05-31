#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2014 Michael J. Keith, University of Manchester

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <T2toolkit.h>
#include "mjklog.h"
#include "kwikfit.h"

/*
 *
   typedef struct kwikfit_result {
	  double phase;
	  double error;
	  double chisq;
	  double *phase_plot;
	  double *chisq_plot;
	  uint64_t nbins;
	  uint64_t nfree;
	  double *data;
	  double *fit;
	  double *residual;
	  double *data_cov;
	  kwikfit_template_t *tmpl;
	  double *amplitudes;
	  double *amp_err;
	  double **amp_cvm;
   } kwikfit_result_t;
*/


void kwikfit_draw_hist(const uint64_t nbins, double *x, double *data){
   uint64_t i;
   cpgmove(x[0],data[0]);
   for (i=1; i< nbins; i++){
	  cpgdraw(x[i-1],data[i]);
	  cpgdraw(x[i],data[i]);
   }
   cpgdraw(x[i-1],data[0]);
   cpgdraw(1.0,data[0]);

}

void kwikfit_draw_points(const uint64_t nbins, double *x, double *data,int ptype){
   uint64_t i;
   for (i=0; i< nbins; i++){
	  cpgpt1(x[i],data[i],ptype);
   }
}

void kwikfit_draw_line(const uint64_t nbins, double *x, double *data){
   uint64_t i;
   cpgmove(x[0],data[0]);
   for (i=1; i< nbins; i++){
	  cpgdraw(x[i],data[i]);
   }
}


void kwikfit_draw_hist_logY(const uint64_t nbins, double *x, double *data){
   uint64_t i;
   double y[nbins];
   for(i=0;i<nbins;i++){
	  y[i]=log10(data[i]);
   }
   kwikfit_draw_hist(nbins,x,y);
}



void kwikfit_plot_result(kwikfit_result_t* result, const char* device) {
   const uint64_t nbins = result->nbins;
   uint64_t i;
   double min=TKfindMin_d(result->residual,nbins);
   double max=TKfindMax_d(result->data,nbins);
   double y[nbins];
   cpgopen(device);
   cpgsvp(0.08,0.95,0.09,0.6);

   double phases[nbins];
   for (i=0; i< nbins; i++){
	  phases[i]=((double)i) /(double)(nbins)-0.5;
   }

   //   int64_t centre = ((int64_t)floor(result->phase * nbins + 0.5) + nbins/2)%nbins;
   int64_t centre = -result->phase*nbins;
   logmsg("phase %lf c %ld",result->phase,centre);

   centre=0;
   logmsg("Min: %lf, Max: %lf",min,max);
   logmsg("Centre: %ld",centre);
   cpgswin(-0.5,0.5,min,max);
   cpgbox("C",0.1,10,"",0,0);
   cpgbox("ABTSN",0.1,10,"BCTSN",0,0);
   cpglab("Phase","Amplitude","");
   cpgsci(4);
   kwikfit_rotate_array(result->data,y,nbins,centre);
   kwikfit_draw_hist(nbins,phases,y);
   cpgsci(3);
   kwikfit_rotate_array(result->fit,y,nbins,centre);
   kwikfit_draw_hist(nbins,phases,y);
   cpgsci(2);
   kwikfit_rotate_array(result->residual,y,nbins,centre);
   kwikfit_draw_hist(nbins,phases,y);

   cpgsci(1);
   double chisq_plot[result->nplot];
   double chifit[result->nplot];
   double chiphase[result->nplot];

   centre=0;
   for(i=1;i<result->nplot;i++){
	  if(result->phase_plot[i] < result->phase_plot[i-1]){
		 centre=i-1;
		 break;
	  }
   }
   kwikfit_rotate_array(result->chisq_plot,chisq_plot,result->nplot,-centre);
   kwikfit_rotate_array(result->chifit,chifit,result->nplot,-centre);
   kwikfit_rotate_array(result->phase_plot,chiphase,result->nplot,-centre);
   logmsg("ph %lf %lf ==> %lf %lf",result->phase_plot[0],result->phase_plot[result->nplot-1],chiphase[0],chiphase[result->nplot-1]);
   for(i=0;i<result->nplot;i++){
	  chisq_plot[i] = chisq_plot[i]/result->nfree;
	  chifit[i] = chifit[i]/result->nfree;
   }

   min=TKfindMin_d(chisq_plot,result->nplot);
   max = min*20;

   double cmin = TKfindMin_d(result->data_cov,nbins);
   double cmax = TKfindMax_d(result->data_cov,nbins);


   logmsg("Min: %lf, Max: %lf",min,max);
   cpgsvp(0.08,0.95,0.6,0.75);
   cpgswin(-0.5,0.5,0,max);
   cpgbox("C",0.1,10,"",0,0);
   cpgbox("BTS",0.1,10,"BCTSN",0,0);
   cpglab("","Chisq","");

   /*for(i=0; i < result->nplot; i++){
	 logmsg("%d %lg %lg",i,result->phase_plot[i],result->chisq_plot[i]);
	 }*/

   for(i=0;i<nbins;i++){
	  y[i]=max*(result->data_cov[i]-cmin)/(cmax-cmin);
   }
   cpgsci(4);
   kwikfit_draw_line(nbins,phases,y);

   double * cov = kwikfit_get_cov(result->fit,nbins);
   for(i=0;i<nbins;i++){
	  y[i]=max*(cov[i]-cmin)/(cmax-cmin);
   }
   cpgsci(5);
   kwikfit_draw_line(nbins,phases,y);



   cpgsci(2);
   kwikfit_draw_line(result->nplot,chiphase,chisq_plot);
   cpgsci(3);
   cpgsls(2);
   kwikfit_draw_line(result->nplot,chiphase,chifit);




   logmsg("top plot");
   max/=4;
   centre=result->nplot/2;
   centre=0;
   kwikfit_rotate_array(result->chisq_plot,chisq_plot,result->nplot,centre);
   kwikfit_rotate_array(result->chifit,chifit,result->nplot,centre);
   kwikfit_rotate_array(result->phase_plot,chiphase,result->nplot,centre);
   fflush(stdout);
   char txt[88];
   sprintf(txt,"%.4f",result->phase);
   double ref_phase = atof(txt);
   for(i=0;i<result->nplot;i++){
	  chisq_plot[i] = chisq_plot[i]/result->nfree;
	  chifit[i] = chifit[i]/result->nfree;
	  chiphase[i]= chiphase[i] - ref_phase;
	  logmsg("%d %lf %lf %lf",i,chiphase[i],chisq_plot[i],chifit[i]);
   }
   logmsg("find l/r");
   double l = -result->error*5;
   double r = result->error*5;
   /*for(i=0;i<result->nplot;i++){
	 if(chisq_plot[i] < max){
	 l=result->phase_plot[i];
	 break;
	 }
	 }
	 for(i=result->nplot-1; i>=0;i--){
	 if(chisq_plot[i] < max){
	 r=result->phase_plot[i];
	 break;
	 }
	 }

	 if (r < -l)r=-l;
	 else l=-r;
	 */
   if(r-l < 1e-3){
	  r+=1e-3;
	  l-=1e-3;
   }

   if (r>0.5)r=0.5;
   if (l<-0.5)l=0.5;

   logmsg("Min: %lf, Max: %lf",min,max);
   logmsg("L: %lf, R: %lf",l,r);

   cpgsci(1);
   cpgsvp(0.08,0.95,0.75,0.9);
   cpgswin(l,r,0,max);
   //cpgbox("",0,0,"A",0,0);
   cpgbox("BCTSM",0,0,"BCTSN",0,0);
   sprintf(txt,"Zoomed Phase (around %.4f)",ref_phase);
   cpglab("","Chisq",txt);
   cpgsls(2);
   cpgmove(result->phase-ref_phase,max);
   cpgdraw(result->phase-ref_phase,0);
   cpgsls(1);

   cpgsci(2);
   kwikfit_draw_line(result->nplot,chiphase,chisq_plot);

   cpgsls(2);
   cpgsci(3);
   kwikfit_draw_line(result->nplot,chiphase,chifit);

   cpgsls(1);
   cpgsci(2);
   kwikfit_draw_points(result->nplot,chiphase,chisq_plot,4);
   cpgclos();

}

