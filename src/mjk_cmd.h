#include <time.h>
#ifndef _mjk_cmd_h
#define _mjk_cmd_h
#ifndef STREQ
#define STREQ(a,b) (strcmp(a,b)==0)
#endif


#ifdef __MACH__
#define clockid_t int
#endif

char* getS(char *lo, char* so, int argc, char** argv, char* val);

char getB(char *lo, char* so, int argc, char** argv, char val);

double getF(char *lo, char* so, int argc, char** argv, double val);
int getI(char *lo, char* so, int argc, char** argv, int val);

void getArgs(int *argc, char** argv);

typedef struct mjk_clock{
   time_t s;
   long n;
   struct timespec last_time;
   char running;
   clockid_t spec;
} mjk_clock_t;

mjk_clock_t *init_clock();
void start_clock(mjk_clock_t* clock);
void stop_clock(mjk_clock_t* clock);
void reset_clock(mjk_clock_t* clock);
double read_clock(mjk_clock_t* clock);



#endif
