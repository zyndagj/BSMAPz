#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include<time.h>
#include "param.h"

using namespace std;

time_t Initial_Time();
//time used during the past step
time_t Cal_StepTime();
//total time exhaust
time_t Cal_AllTime();
//current time on string format
char * Curr_Time();

bit32_t myrand(int i, bit32_t * rseed);
bool HitComp(gHit a, gHit b);
bool HitComp2(gHit a, gHit b);
bool HitCompChr(Hit a, Hit b);
void disp_bfa(bit32_t a, int len=16);
void disp_bfa64(bit64_t a, int len=32);

#endif //_UTILITIES_H_
