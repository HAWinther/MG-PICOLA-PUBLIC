#include <sys/time.h>
#include "msg.h"
#include "timer.h"

#define nCategory 4
static const char * CatName[]= {"Init        ", 
                                "GenerateIC  ", 
                                "TimeStepping",
                                "Output      "
                               };

#define nSubCategory 17
static const char * SubName[]= {"",                      
                                "Kick                 ",            
                                "Drift                ",         
                                "MoveParticles        ", 
                                "Forces               ",                
                                "ModifiedGravity      ", 
                                "PtoMesh              ",
                                "MtoParticles         ", 
                                "ReadParticlesFromFile",
                                "WriteOutput          ",
                                "FFTs                 ",
                                "NonGaussianIC        ", 
                                "DisplacementFields   ",
                                "OutputLightcone      ",
                                "DriftLightcone       ",
                                "PofkComputation      ",
                                "HaloFinding          "
                               };
static int initialized = 0;
static enum Category Cat;
static double Time[nCategory][nSubCategory], tBegin[nCategory][nSubCategory];

static double now(){
  struct timeval tp;
  gettimeofday(&tp, 0);
  return (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6;
}

void timer_init(){
  for(int i = 0; i < nCategory; i++) {
    for(int j = 0; j < nCategory; j++) {
      Time[i][j] = 0.0;
      tBegin[i][j] = 0.0;
    }
  }
  initialized = 1;
}

void timer_set_category(enum Category new_cat){
  double nw = now();

  if(initialized)
    Time[Cat][all] += nw - tBegin[Cat][all];
  else 
    timer_init();

  Cat= new_cat;
  tBegin[Cat][all] = nw;
}

void timer_start(enum SubCategory sub){
  tBegin[Cat][sub] = now();
}

void timer_stop(enum SubCategory sub){
  Time[Cat][sub] += now() - tBegin[Cat][sub];
}

void timer_print(){
  timer_set_category(0);
  double total = 0.0;

  for(int i = 0; i < nCategory; i++)
    total += Time[i][0];

  for(int icat = 0; icat < nCategory; icat++) {
    msg_printf(info, "%-16s %7.2f   %4.1f%%\n", 
        CatName[icat], Time[icat][0], 100.0*Time[icat][0]/total);

    for(int isub = 1; isub < nSubCategory; isub++) {
      if(Time[icat][isub] > 0.0)
        msg_printf(info, "  %-14s %7.2f   %4.1f%%\n", 
            SubName[isub], Time[icat][isub], 100*Time[icat][isub]/total);
    }
  }
  msg_printf(info, "----------------------------------\n");
  msg_printf(info, "%-16s %7.2f sec\n", "Total", total);
}

