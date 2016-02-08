#ifndef _PARAM_H_
#define _PARAM_H_

#define SEGLEN 32

#include <unistd.h>
#include<string>
#include<cstring>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include<vector>
#include<algorithm>

using namespace std;

const unsigned int FIXELEMENT=6; // 160/32+1
const unsigned int MAXSNPS=15;
const unsigned int MAXGAPS=3;

typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;

typedef bit32_t ref_id_t;
typedef bit32_t ref_loc_t;
typedef ref_loc_t* NewIndex;

struct Hit {
	ref_id_t chr;       //index of chr
	ref_loc_t loc;     //location of first bp on reference seq, count from 0
};

struct gHit {
    ref_loc_t loc;
    bit32_t chr: 18;
    bit32_t strand: 2; // strand infomation 00: '++', 01: '+-', 10: '-+', 11: '--'
    bit16_t gap_size:4;  // positive for insert on read, negative for insert on ref, 0 for no-gap
    bit16_t gap_pos:8;
};

class Param {
public:
	Param();
	void SetSeedSize(int n);
    void InitMapping();
	void SetAdaptors(int n);
    void SetDigestionSite(const char *a);
    void SetAlign(char readnt, char refnt);

public:

    //hot data section
    
    //cache line 1
    bit32_t seed_size __attribute__((aligned(64))); // 4 byte
    bit32_t seed_bits, max_snp_num, max_num_hits; //16 byte
    bit32_t min_insert, max_insert, max_seedseg_num, chains; //32 byte
    bit32_t RRBS_flag, pairend, index_interval, randseed; //48 byte
    bit32_t gap, gap_edge, trim_lowQ, max_ns; //64 byte

    //cache line 2
    bit32_t min_read_size, n_adapter, read_start, read_end; // 16 byte
    bit32_t out_sam, out_ref, out_unmap, report_repeat_hits; //32 byte   
    bit32_t sam_header, max_readlen; bit16_t *_T; //48 byte
    int input_format, output_format, gz_input, gz_ref; // 64 byte    

    //cache line 3
    bit32_t max_dbseq_size, append_dbseq_size, total_ref_seq;
    vector <bit32_t> digest_pos; // 16 byte
    vector <string> digest_site; 
    string useful_nt; // 32 byte
    string nx_nt; int num_procs; //48 byte
    bit32_t stdout; bit8_t zero_qual, qual_threshold, default_qual; //51 byte
    bit8_t read_nt, ref_nt, N_mis;
    bit32_t max_kmer_num, verbose_level; 
    float max_kmer_ratio;
    bit32_t nt3, pipe_out, seed_bits_lz; 
    
    //cache line 4+    
    bit32_t profile[MAXSNPS+1][16] __attribute__((aligned(64)));
    
 //int chains;   //0: forward strands only ; 1: forward and reverse strands	
	string adapter[10];
    //int input_format; // 0: fasta, 1:fastq, 2: SAM, 3: BAM, -1: auto detect
    //int output_format; // 0: bsp, 1:SAM, 2: BAM, -1: auto detect
    //int gz_input, gz_ref; // 0: no, 1: yes, -1: auto detect

    void disp_bin64(bit64_t t, int len=32) {
        for(int i=len*2-1;i>=0;i--) {
            if(t&(1ULL<<i)) cout<<'1';
            else cout<<'0';
            if(i%4==0) cout<<' ';
        }
        cout<<endl;
    }

    inline bit64_t XT64(register bit64_t tt) {tt-=(tt<<1)&tt&0xAAAAAAAAAAAAAAAAULL; return tt;}
    inline bit32_t XT32(register bit32_t tt) {tt-=(tt<<1)&tt&0xAAAAAAAAUL; return tt;}

    inline bit32_t XT(register bit32_t tt) {
        register bit32_t ss;
        tt-=(tt<<1)&tt&0xAAAAAAAAUL;
        tt-=(tt>>2)&0x33333333UL; // 4 bit transform
        ss=(tt&0xF0F0F0F0UL)>>1;
        tt-=ss-(ss>>3); //8 bit transform, ss+(ss<<3) = s*9
        ss=(tt&0xFF00FF00UL)>>2; 
        tt=(tt&0x00FF00FFUL)+ss+(ss>>2)+(ss>>6);
        return (tt&0xFFFFUL)+(tt>>16)*6561;
        /*
        bit64_t t=(bit64_t) ((~(tt<<1&tt))|0x55555555UL)&tt;
        return (bit32_t) ((bit64_t) ((((t&0xFFCULL)*0x100401004010ULL)&0xC00C00C00C00C0ULL)*0x300901B0510F3ULL)>>54)+((bit32_t) t&3UL)  //low 6nt
              +(bit32_t) ((bit64_t) ((((t>>8&0x3FF0ULL)*0x10040100401ULL)&0x30030030030030ULL)*0x100300901B051ULL)>>52)*729UL // mid 5nt
              +(bit32_t) ((bit64_t) ((((t>>22)*0x10040100401ULL)&0x3003003003003ULL)*0x100300901B051ULL)>>48&0xFFFULL)*177147UL;// high 5nt
        */
    }                                                                                                       

    inline bit32_t XC(bit32_t tt) {return ((~tt)<<1)|tt|0x55555555UL;}  // generate T2C mask according to C locations
    inline bit64_t XC64(register bit64_t tt) {return ((~tt)<<1)|tt|0x5555555555555555ULL;}
    inline bit32_t swap_endian(bit32_t tt) {return ((tt>>16)|(tt<<16));}
    inline bit64_t swap_endian64(bit64_t tt) {return ((tt>>32)|(tt<<32));}

    inline bit32_t XM(bit32_t tt) {
        tt=(tt|(tt>>1))&0x55555555;
        tt=(tt+(tt>>2))&0x33333333;
        return (((bit32_t)(tt*0x1111111))>>28)+(tt&0x3);
    }

    inline bit32_t XM64(register bit64_t tt) {
        //tt=(tt|(tt>>1))&0x5555555555555555ULL;
        tt|=tt>>1;
        tt&=0x5555555555555555ULL;
        /*
        register bit64_t tt1;
        tt1=tt>>2;
        tt+=tt1;
        tt1=(tt+(tt1>>2))&0x0C30C30C30C30C30ULL;
        return ((tt1*0x0041041041041041ULL)>>58)+(tt&0x3ULL);    
        */
        tt+=tt>>2;
        tt&=0x3333333333333333ULL;
        tt+=tt>>4;
        tt&=0x0F0F0F0F0F0F0F0FULL;
        tt*=0x0101010101010101ULL;
        tt>>=56;
        return tt;
    }

    /*
    inline bit32_t XM64X2(bit64_t tt1, bit64_t tt2) {
        tt1=((tt1|(tt1>>1))&0x5555555555555555ULL)+((tt2|(tt2>>1))&0x5555555555555555ULL);
        tt1=(tt1&0x3333333333333333ULL)+((tt1>>2)&0x3333333333333333ULL);
        return (((tt1+(tt1>>4))&0x0F0F0F0F0F0F0F0FULL)*0x0101010101010101ULL)>>56;
    }
    */
    
    bit32_t map3to4(bit32_t tt){
      	int s=0, i; for(i=0;i<16;i++) {s|=(tt%3)<<i*2; tt/=3;};
	return s;
    }
};

#endif //_PARAM_H_
