#ifndef _PAIRS_H_
#define _PAIRS_H_

#include "dbseq.h"
#include "reads.h"
#include "align.h"

using namespace std;


extern char chain_flag[];

struct PairHit
{
	bit16_t chain;
	bit8_t na, nb;   //# of snps
	bit32_t insert;
	gHit a;
	gHit b;
};

typedef PairHit PairArray[MAXHITS+1];

class PairAlign
{
public:
	PairAlign();
    ~PairAlign();
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a1, vector<ReadInf> &a2);
	int RunAlign(RefSeq &ref);
	void Do_Batch(RefSeq &ref);
	void StringAlign(RefSeq &ref, string &os);
	void StringAlign_ClosestUnpair(RefSeq &ref, string &os);
	
	//added by yxi
	int GetPairs(bit32_t na, bit32_t nb);
    int StringAlignPair(RefSeq &ref, string &os);
	void StringAlignUnpair(int fa, int fb, RefSeq &ref, string &os);
    void s_OutHitPair(PairHit pp, int n, RefSeq &ref, string &os);
    void s_OutHitUnpair(int readinpair, int chain_a, int chain_b, int ma, bit32_t na, gHit ha, int mb, gHit hb, RefSeq &ref, string &os);
    int TrimAdapter();
    void FixPairReadName();

public:	
	SingleAlign _sa;
	SingleAlign _sb;
	bit32_t num_reads;
	bit32_t n_aligned_pairs, n_aligned_a, n_aligned_b;	
	bit32_t n_unique_pairs, n_unique_a, n_unique_b;
	bit32_t n_multiple_pairs, n_multiple_a, n_multiple_b;
	string _str_align;
protected:
	bit32_t _cur_n_hits[2*MAXSNPS+1];
	//PairHit pairhits[2*MAXSNPS+1][MAXHITS+1];
    PairArray *pairhits;
    bit32_t rand_rSeed;	//thread safe RNG seed
	//by yxi
    char _mapseq[256];
    char _ch[1024];  
    char cigar[16];
   	SingleAlign * _stmp;  
    int checked_pair_mismatch[MAXSNPS+1][MAXSNPS+1];
};

#endif //_PAIR_ALIGH_H_
