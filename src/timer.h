#ifndef TIMER_H
#define TIMER_H 1

enum Category {
               _Init, 
               _GenerateIC, 
               _TimeStepping,
               _Output
              };

enum SubCategory {
                  all, 
                  _Kick, 
                  _Drift, 
                  _MoveParticles, 
                  _Forces, 
                  _ComputeFifthForce, 
                  _PtoMesh, 
                  _MtoParticles, 
                  _ReadParticlesFromFile, 
                  _WriteOutput, 
                  _FFT, 
                  _NonGaussianIC, 
                  _DisplacementFields, 
                  _OutputLightcone, 
                  _DriftLightcone,
                  _PofkComputation,
                  _HaloFinding
                 };

void timer_set_category(enum Category new_cat);
void timer_start(enum SubCategory sub);
void timer_stop(enum SubCategory sub);
void timer_print();
#endif
