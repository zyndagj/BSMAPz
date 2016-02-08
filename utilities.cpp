#include "utilities.h"

using namespace std;

extern Param param;

static time_t time_begin;
static time_t time_last;

time_t Initial_Time() {
	time_begin = time(NULL);
	time_last = time_begin;
	return time_begin;
};

//time used during the past step
time_t Cal_StepTime() {
	time_t tused = time(NULL)-time_last;
	time_last = time(NULL);
	return tused;
};

//total time exhaust
time_t Cal_AllTime() {
	return time(NULL)-time_begin;
};

//current time on string format
char * Curr_Time() {
	time_t t=time(NULL);
	char *atime=ctime(&t);
	atime[strlen(atime)-1]=0;
	return atime;
}

//generate randomness

bit32_t myrand(int i, bit32_t* rseed) {
	if(param.randseed == 0) return rand_r(rseed);
	else{
		bit64_t v;
    	v = ((bit64_t) i+ param.randseed*1000000) * 3935559000370003845LL + 2691343689449507681LL;
        v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
        v *= 4768777513237032717LL;
        v ^= v << 20; v ^= v >> 41; v ^= v << 5;
        return (bit32_t) (v&0xffffffffUL);
	}
}


bool HitComp(gHit a, gHit b) {
    return ((a.chr<b.chr)||(a.chr==b.chr)&&(a.loc<b.loc));
}

bool HitComp2(gHit a, gHit b) {
    return ((a.chr/2<b.chr/2)||(a.chr/2==b.chr/2)&&(a.loc<b.loc));
}

bool HitCompChr(Hit a, Hit b) {return a.chr<b.chr;}

void disp_bfa(bit32_t a, int len) {
    for(int i=len-1; i>=0; i--) {
        cout<<param.useful_nt[(a>>(i*2))&0x3];
    }
    cout<<" ";
}

void disp_bfa64(bit64_t a, int len) {
    for(int i=len-1; i>=0; i--) {
        cout<<param.useful_nt[(a>>(i*2))&0x3];
    }
    cout<<" ";
}

