#ifndef _ALIGN_H_
#define _ALIGN_H_

#include<vector>
#include<string>
#include<cstdio>
#include<set>
#include<algorithm>
#include<cmath>

#include "param.h"
#include "reads.h"
#include "dbseq.h"

using namespace std;

const unsigned int FIXSIZE=SEGLEN*FIXELEMENT; 
typedef gHit HitArray[MAXHITS+1];
typedef HitArray HitMatrix[MAXSNPS+1];

typedef bit64_t SegArray[FIXELEMENT*2];
typedef bit32_t SeedArray[16];

extern Param param;
extern char rev_char[];
extern char chain_flag[];

class SingleAlign {
public:
	SingleAlign();
	~SingleAlign();
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a);
	bit32_t CountNs();
	void set_RRBS_start(); //by yxi
	int TrimLowQual();
	void ConvertBinaySeq();
    int CountSeeds(RefSeq &ref, int n, bit32_t start);
	
	inline void CountMismatch(bit64_t *q, bit32_t offset, bit64_t *s);
	void ClearHits();
	int RunAlign(RefSeq &ref);
	int FilterReads();
	void Do_Batch(RefSeq &ref);
	void StringAlign(RefSeq &ref, string &os);
	void Reverse_Seq();
	void s_OutHit(int chain, int n, bit8_t nspsn, gHit *hit, int insert_size, RefSeq &ref, string &os);
    
    //by yxi
    void SnpAlign(RefSeq &ref, bit32_t mode);
    int TrimAdapter();
    void SortHits4PE(int n);
    void ReorderSeed(RefSeq &ref);
    bit32_t GetTotalSeedLoc(RefSeq &ref, bit32_t start);
    void AdjustSeedStartArray(RefSeq &ref);
    inline bit32_t MismatchPattern0(bit64_t *q, bit64_t *s, bit32_t seglen_offset);
    inline bit32_t MismatchPattern1(bit64_t *q, bit64_t *s, bit32_t gap_index, bit32_t seglen_offset);
    bit32_t GapAlign(RefSeq &ref, bit32_t mode, bit32_t seed_pos);
    //int MatchGap(bit32_t mmi1[], bit32_t mmi2[], bit32_t shift);    
    inline bit32_t AddHit(RefSeq &ref, bit32_t w, bit32_t mode);
    gHit int2hit(RefSeq &ref, gHit gh, int gap_size, bit32_t gap_pos);

public:
    //Hot data section
	//cache line 1
	bit32_t tmp_snp __attribute__((aligned(64))); //keep in the same cache line, 4 byte
	bit32_t snp_thres, map_readlen, raw_readlen; // 16 byte
	bit32_t read_chain_index, ref_chain_index, _seed, cseed_offset; //32 byte
	bit32_t rand_rSeed, n_aligned, read_max_snp_num, seedseg_num; // 48 byte
	bit32_t xflag_chain[2], xseed_start_offset[2]; // 64 byte

	//cache line 2
	gHit _hit, _ghit; // 16 byte
	bit32_t *_pro, *_refloc2; // 32 byte
	Hit *_refloc; //40 byte
    HitMatrix *xhits; // 48 byte
    HitArray *hits, *chits; //64 byte

    //cache line 3
    bit32_t N_count, tmp_snp0, n_unique, n_multiple;

    SegArray xseq[2] __attribute__((aligned(64)));
    bit64_t tmp_seq, tmp_reg;
    bit32_t end_element, end_offset;
	SeedArray xseeds[2][MAXSNPS+1];
	bit32_t xseed_array[2][FIXSIZE-SEGLEN], xseedreg_array[2][FIXSIZE-SEGLEN];	
	bit32_t x_cur_n_hit[2][MAXSNPS+1], *_cur_n_hit, *_cur_n_chit;
	bit32_t xseed_start_array[2][MAXSNPS+1];

    set<ref_loc_t> *hitset, *ghitset; //, *chitset; 
    pair<int,int> xseedindex[2][MAXSNPS+1];	
	bit32_t mm_index[2*MAXGAPS+1][MAXSNPS+1];
	char cigar[16];
	bit64_t total_candidates, total_reads, total_seeds;

    //cold data section
	vector<ReadInf>::iterator _pread;
	string _outseq[2], _outqual[2];
	bit32_t num_reads;
	vector<ReadInf> mreads;
	string _str_align;   //align results, prepare for output
    vector<string>::iterator read_motif_iter;

	//local variables
	string::iterator _sp;
	string::reverse_iterator _sq;
	
	char _ch[1024];
    pair<ref_loc_t,bit32_t> seg_info;
    char _mapseq[256];
    string::iterator _readnt, _adapternt;
};

inline void SingleAlign::CountMismatch(bit64_t *q, bit32_t offset, bit64_t *s) {
    if(param.nt3) {
    	if((tmp_snp=N_count+param.XM64((((*q)>>offset)^param.XT64(*s))&((q[FIXELEMENT])>>offset)))>snp_thres) return;
        for(bit32_t i=1;i<FIXELEMENT;i++) {
            if((tmp_snp+=param.XM64(((((q[i-1]<<1)<<(63-offset))|q[i]>>offset)^param.XT64(s[i]))&(((q[FIXELEMENT-1+i]<<1)<<(63-offset))|q[i+FIXELEMENT]>>offset)))>snp_thres) return;
        }
    }
    else {
    	if((tmp_snp=N_count+param.XM64((((*q)>>offset)&param.XC64(*s)^(*s))&((q[FIXELEMENT])>>offset)))>snp_thres) return;
    	for(bit32_t i=1;i<FIXELEMENT;i++) {
    		if((tmp_snp+=param.XM64(((((q[i-1]<<1)<<(63-offset))|q[i]>>offset)&param.XC64(s[i])^s[i])&(((q[FIXELEMENT-1+i]<<1)<<(63-offset))|q[i+FIXELEMENT]>>offset)))>snp_thres) return; 
    	}
    	/*
        tmp_snp=N_count;
        for(bit32_t i=0;i<5;i++) { 
        	tmp=(s[i]<<offset)|((s[i+1]>>(63-offset))>>1);
            if((tmp_snp+=param.XM64((q[i]&param.XC64(tmp)^tmp)&q[i+5]))>snp_thres) return;
		}
        */
    }
}

inline bit32_t SingleAlign::MismatchPattern0(bit64_t *q, bit64_t *s, bit32_t offset) {
    register bit64_t tmp;
    bit32_t *mm_array; int i, j,jj, ss=0;
    //cout<<"in MMpattern()  gap_index:"<<gap_index<<endl;
    //cout<<" offset="<<seglen_offset<<" end_element="<<end_element<<" end_offset="<<end_offset<<endl;
    mm_array=mm_index[0]; 
    //cout<<"ss="<<ss<<" tmp_snp="<<tmp_snp<<" snp_thres="<<snp_thres<<endl;
    for(i=0;i<end_element;i++) {
		tmp=(s[i]<<offset)|((s[i+1]>>(63-offset))>>1);
        if(param.nt3) tmp=(q[i]^param.XT64(tmp));
        else tmp^=q[i]&param.XC64(tmp);
        j=(int)(i*SEGLEN)-1;
        while(tmp) {
            jj=__builtin_clzll(tmp)>>1; j+=jj+1;
            //param.disp_bin64(tmp);
            //cout<<"\ti="<<i<<" j="<<j<<" jj="<<jj<< "ss="<<ss+1<<endl;
            mm_array[ss++]=j; 
            if(ss>(int)snp_thres-2) return j;
            tmp<<=2; tmp<<=(jj<<1);
        }
        //cout<<"ss="<<ss<<" tmp_snp="<<tmp_snp<<" snp_thres="<<snp_thres<<endl;
    }
    tmp=(s[i]<<offset)|((s[i+1]>>(63-offset))>>1);
    if(param.nt3) tmp=(q[i]^param.XT64(tmp));
    else tmp^=q[i]&param.XC64(tmp);
	tmp>>=end_offset; tmp<<=end_offset;
    j=(int)(i*SEGLEN)-1;
    while(tmp) {
        jj=__builtin_clzll(tmp)>>1; j+=jj+1;
        mm_array[ss++]=j; 
        if(ss>(int)snp_thres-2) return j;
        tmp<<=2; tmp<<=(jj<<1);
    }
    for(;ss<=(int)snp_thres-2;ss++) mm_array[ss]=map_readlen;
    return map_readlen;
}

inline bit32_t SingleAlign::MismatchPattern1(bit64_t *q, bit64_t *s, bit32_t gap_index, bit32_t offset) {
    register bit64_t tmp; register bit32_t shift;
    bit32_t *mm_array; int i,ii, j,jj, ss=0;
    //cout<<"in MMpattern()  gap_index:"<<gap_index<<endl;             
    //cout<<"shift="<<(1-((int)gap_index%2)*2)*((int)gap_index+1)/2<<" offset="<<seglen_offset<<" end_element="<<end_element;
    //cout<<" end_offset="<<end_offset<<endl;  
    //for(i=0;i<5;i++) disp_bfa64(param.swap_endian64(q[i])); cout<<endl;
    //for(i=0;i<5;i++) disp_bfa64(param.swap_endian64(s[i])); cout<<endl;
    mm_array=mm_index[gap_index]; 
    for(i=end_element,ii=0;i>=0;i--,ii+=SEGLEN) { 
		tmp=(s[i]<<offset)|((s[i+1]>>(63-offset))>>1);
        if(param.nt3) tmp=q[i]^param.XT64(tmp);
        else tmp^=q[i]&param.XC64(tmp);
        shift=end_offset*(i==end_element); tmp=(tmp>>shift)<<shift;
        j=ii-(int) (end_offset>>1)-1; 
        while(tmp) {
            jj=__builtin_ctzll(tmp)>>1; j+=jj+1;
            //cout<<"jj="<<jj<<" j="<<j<<" "; param.disp_bin64(tmp); 
            //cout<<" mm_array["<<ss<<"]="<<j<<endl;
            mm_array[ss++]=j; 
            if(ss>(int)snp_thres-2) return 1;
            tmp>>=2; tmp>>=(jj<<1);
        }
    }
    for(;ss<=(int)snp_thres-2;ss++) mm_array[ss]=map_readlen;
    return map_readlen;
}   

/*
inline bit32_t SingleAlign::MismatchPattern1(bit64_t *q, bit64_t *s, bit32_t gap_index, bit32_t seglen_offset) {
    register bit64_t tmp; register bit32_t shift;
    bit32_t *mm_array; int i,ii, j,jj, ss=0, end_element, right, end_offset, max0, offset=seglen_offset<<1;
    end_element=(map_readlen+seglen_offset-1)/SEGLEN; 
    right=(map_readlen+seglen_offset-1)%SEGLEN+1;  end_offset=SEGLEN-right;
    //cout<<"in MMpattern()  gap_index:"<<gap_index<<endl;             
    //cout<<"shift="<<(1-((int)gap_index%2)*2)*((int)gap_index+1)/2<<" offset="<<seglen_offset<<" end_element="<<end_element;
    //cout<<" end_offset="<<end_offset<<endl;  
    //for(i=0;i<5;i++) disp_bfa64(param.swap_endian64(q[i])); cout<<endl;
    //for(i=0;i<5;i++) disp_bfa64(param.swap_endian64(s[i])); cout<<endl;
    mm_array=mm_index[gap_index]; 
	for(i=end_element,ii=0;i>=0;i--,ii+=SEGLEN) {
        tmp_seq=(i>0)*(q[i-(i>0)]<<63-offset)<<1; 
        if(param.nt3) tmp=(tmp_seq|q[i]>>offset)^param.XT64(s[i]);
		else tmp=(tmp_seq|q[i]>>offset)&param.XC64(s[i])^s[i];
		shift=seglen_offset*2*(i==0); tmp=(tmp<<shift)>>shift;
		shift=end_offset*2*(i==end_element); tmp=(tmp>>shift)<<shift;
		j=ii-(int) end_offset-1; 
		while(tmp) {
			jj=__builtin_ctzll(tmp)>>1; j+=jj+1;
			//cout<<"jj="<<jj<<" j="<<j<<" "; param.disp_bin64(tmp); 
			//cout<<" mm_array["<<ss<<"]="<<j<<endl;
			mm_array[ss++]=j; 
            if(ss>(int)snp_thres-2) return 1;
            tmp>>=2; tmp>>=(jj<<1);
        }
	}
    for(;ss<=(int)snp_thres-2;ss++) mm_array[ss]=map_readlen;
    return map_readlen;
}	
*/
inline bit32_t SingleAlign::AddHit(RefSeq &ref, bit32_t w, bit32_t mode) {
	if((int)_ghit.loc<0) return 0; //underflow the start of refseq
  	if(_ghit.loc+map_readlen>ref.title[_ghit.chr].size) return 0; //overflow the end of refseq
    if(_ghit.gap_size) {
    	if(!ghitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0;
    }
    else {
    	if(!hitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0; //hit already exist
	}
	//cout<<"###HIT### "<<ref.title[_ghit.chr].name<<":"<<_ghit.loc<<" mis:"<<tmp_snp<<" snp_thres:"<<snp_thres<<endl;
    xhits[read_chain_index][w][x_cur_n_hit[read_chain_index][w]++]=_ghit;
    //if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return 1;
    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) {
        if(w==0) return 1; 
        else snp_thres=w-1;
	}
    return 0;
}

inline void SingleAlign::Reverse_Seq() {
	_outseq[0]=_pread->seq; _outseq[1]=_pread->seq;
	reverse(_outseq[1].begin(), _outseq[1].end());
	for(string::iterator p=_outseq[1].begin(); p!=_outseq[1].end(); ++p) *p=rev_char[(unsigned char)*p];
	_outqual[0]=_pread->qual; _outqual[1]=_pread->qual;
	reverse(_outqual[1].begin(), _outqual[1].end());
}
#endif //_ALIGN_H_
