#ifndef _DBSEQ_H_
#define _DBSEQ_H_

#include<vector>
#include<set>
#include<string>
#include<iostream>
#include<fstream>
#include "param.h"
#include "utilities.h"
#include "gzstream.h"

#define PREFETCH_CAL_UNIT 8
#define PREFETCH_CRT_UNIT 6
#define PREFETCH_LOOP 10U
#define REF_MARGIN 400
#define BINSEQPAD 2

using namespace std;

struct OneBfa
{
	bit32_t n;  //count
	bit64_t *s;
};
struct RefTitle
{
	string name;
	bit32_t size;
	//added by yxi
	bit32_t rc_offset;
};
struct Block
{
	bit32_t id;
	bit32_t begin;
	bit32_t end;
};

struct KmerLoc
{
	bit32_t n1; //ab, ac, ad seed
	Hit *loc1;
};

struct KmerLoc2
{
    bit32_t n[2]; //ab, ac, ad seed
    bit32_t *loc1; 
};

struct shortHit{
    ref_id_t chr;
    ref_loc_t loc;
};

struct shorthitcompclass{
    bool operator()(shortHit a, shortHit b) {
    	if(a.loc<b.loc) return 1;
    	else if(a.loc>b.loc) return 0;
    	else if(a.chr<b.chr) return 1;
    	else return 0;
    }
};



class RefSeq
{
public:
	RefSeq();
	ref_loc_t LoadNextSeq(ifstream &fin);
	ref_loc_t LoadNextSeq(igzstream &fin);
	void BinSeq(OneBfa &a);
	void cBinSeq(OneBfa &a);
	void UnmaskRegion();
	void Run_ConvertBinseq(ifstream &fin, igzstream &gzfin);
	inline bit32_t s_MakeSeed_1(bit64_t *_m, int _a);
	inline bit32_t s_MakeSeed_2(bit32_t *_m, int _a);
	//inline bit32_t s_MakeSeed_2(bit24_t *_m, int _a, int _b, int _c, int _d);
	
	void InitialIndex();
	void CalKmerFreq();
	void t_CalKmerFreq(bit32_t ref_chain);
	void AllocIndex();
	void FillIndex();
	void t_FillIndex(bit32_t ref_chain);
	void CreateIndex();
	void ReleaseIndex();
	void FinishIndex();
	void find_CCGG();
    pair<ref_loc_t,bit32_t> CCGG_seglen(ref_id_t chr, ref_loc_t pos, int readlen);
    ref_loc_t hit2int(Hit h);

public:
	bit32_t total_num;
	bit64_t sum_length;
	vector<OneBfa> bfa;
	bit32_t total_kmers;
	KmerLoc *index;
	vector<RefTitle> title;	
protected:
	ref_id_t _count;
	string _name;
	string _seq;
	ref_loc_t _length;
	shortHit tmploc;
public:	
	vector<Block> _blocks;  //unmasked ref region
	//map<shortHit,bit32_t,shorthitcompclass> ccgg_seglen;

    //by yxi
    //int max_seedseg_num;
	vector<vector<ref_loc_t> > CCGG_index[50];
	vector<vector<pair<ref_loc_t,bit32_t> > > CCGG_sites;
    vector<pair<ref_loc_t,bit32_t> > *CCGG_sites_chr;
    int n_CCGG;
    KmerLoc2 *index2;
    bit64_t *refcat, *crefcat, *xref[2];
    vector <bit32_t> ref_anchor, cref_anchor;
    vector <bit32_t*> mem_pool;
};

#endif //_DBSEQ_H_
