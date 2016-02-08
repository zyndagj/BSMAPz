#include "param.h"
#include<iostream>

using namespace std;

Param::Param()
{
	num_procs=sysconf(_SC_NPROCESSORS_ONLN);
	if(num_procs>8) num_procs=8;

	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
	
	max_ns = 5;
	trim_lowQ=0;
	
	zero_qual= '!';
	qual_threshold= 0;
	default_qual=40;
	
	min_insert= 28;
	max_insert= 1000;

	SetSeedSize(16);	
	//seed_size= 16;
	//seed_bits=(1<<(seed_size*2))-1;
	
	max_snp_num = 108;
	max_num_hits = MAXHITS>100?100:MAXHITS;
	max_kmer_ratio = 5e-7;
	
	min_read_size=seed_size;
	input_format=gz_input=gz_ref=-1;	
    n_adapter=0;

	report_repeat_hits = 1;

	useful_nt="ACGTacgt";
	nx_nt="NXnx";
	
	for(bit32_t i=0; i<seed_size; i++) seed_bits|=0x3<<i*2;

	out_sam=0;
	read_start=1;
	read_end=~0;

    out_ref=0;
    out_unmap=0;
    RRBS_flag=0;
    index_interval=4;
    randseed=0;
    chains=0;
    pairend=0;
    gap=0;
    gap_edge=6;
    max_readlen=(FIXELEMENT-1)*SEGLEN;
    SetAlign('T','C');
    sam_header=1;
    stdout=1;
    N_mis=0;
    nt3=0;
    pipe_out=0;
    verbose_level=1;

};

void Param::InitMapping(){
	for(bit32_t i=0; i<index_interval; i++) for(bit32_t j=0; j<=MAXSNPS; j++) {  //for 4 binary seqs
	    profile[j][i]=((j*seed_size+i+index_interval-1)/index_interval)*index_interval;
	}
}

void Param::SetDigestionSite(const char *a) {
    string ds=a, ds1; bit32_t dp;
    int IUPAC_count[256], i, j, ii[256], l;
    char IUPAC[256][4];
    for(i=0;i<256;++i) {IUPAC_count[i]=0; for(j=0;j<4;++j) IUPAC[i][j]=(char) 0;}
    IUPAC['A'][0]='A'; IUPAC['C'][0]='C'; IUPAC['G'][0]='G'; IUPAC['T'][0]='T'; 
    IUPAC['N'][0]='A'; IUPAC['N'][1]='C'; IUPAC['N'][2]='G'; IUPAC['N'][3]='T';
    strcpy(IUPAC['R'], "AG"); strcpy(IUPAC['Y'], "CT"); strcpy(IUPAC['S'], "CG");
    strcpy(IUPAC['W'], "AT"); strcpy(IUPAC['K'], "GT"); strcpy(IUPAC['M'], "AC");
    strcpy(IUPAC['B'], "CGT"); strcpy(IUPAC['D'], "AGT"); strcpy(IUPAC['H'], "ACT"); strcpy(IUPAC['V'], "ACG");
    for(i=0;i<256;++i) for(j=0;j<4;++j) IUPAC_count[i]+=(IUPAC[i][j]>0);
    //cout<<digest_site<<endl;
    if((dp=ds.find('-'))<0) {
        cout<<"Digestion position not marked, use \'-\' to mark. example: \'C-CGG\'\n";  
        exit(1);
    }
    ds.erase(dp,1); l=ds.size(); ds1=ds;
    for(i=0;i<l;++i) ii[i]=0; j=0;
    while(j<l) {
        for(i=0;i<l;++i) ds1[i]=IUPAC[ds[i]][ii[i]];
        digest_site.push_back(ds1);
        digest_pos.push_back(dp);
        j=0; ii[j]++; 
        while(ii[j]>=IUPAC_count[ds[j]]&&j<l) {
            ii[j]=0; j++; ii[j]++;
        }
    }
    RRBS_flag=1;
    index_interval=1;
    //SetSeedSize(12);
}

void Param::SetSeedSize(int n) {
    //if(RRBS_flag&&n!=12) {      
    //    cerr<<"seed size fixed to 12 in RRBS mode, -s option ignored\n";
    //    return;
    //}
    if(n>16||n<10) {cerr<<"seed size must be between 10 and 16\n"; exit(1);}
	seed_size=n;
	seed_bits_lz=(SEGLEN-seed_size)*2;
	min_read_size=seed_size+index_interval-1;	
	seed_bits=0;
	for(bit32_t i=0; i<seed_size; i++) seed_bits|=0x3<<i*2;
}

bit8_t alphabet[256];
bit8_t rev_alphabet[256];
bit8_t alphabet0[256]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'A' */
0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'a' */ 
0,0,1,0,0,0,2,0,0,0,0,0,0,0,0, /* next is 'p' */
0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

bit8_t reg_alphabet[256]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'A' */ 
3,0,3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'a' */ 
3,0,3,0,0,0,3,0,0,0,0,0,0,0,0, /* next is 'p' */
0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

bit8_t nv3='N'; //set any unknown char as 'N'
char rev_char[256] = {
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'A' */ 
'T',nv3,'G',nv3,nv3,nv3,'C',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'P' */
nv3,nv3,nv3,nv3,'A',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'a' */ 
't',nv3,'g',nv3,nv3,nv3,'c',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'p' */
nv3,nv3,nv3,nv3,'a',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3
};

bit8_t bit_nt[4];

char nt_code[4] = {'A', 'C', 'G', 'T'};

void Param::SetAlign(char readnt, char refnt){
    int i, j; bit8_t tmp=0;
    read_nt=toupper(readnt); ref_nt=toupper(refnt);
    if((!reg_alphabet[read_nt])||(!reg_alphabet[ref_nt])) {
        cerr<<"Unknown nucleotide."<<endl;
        exit(1);
    }
    if(read_nt==ref_nt){
        cerr<<"Must specify different nucleotides for additional alignment."<<endl;
        exit(1);
    }

    for(i=0;i<4;i++) bit_nt[i]=100;    
    bit_nt[alphabet0[read_nt]]=3;
    bit_nt[alphabet0[ref_nt]]=1;
    for(i=0;i<4;i++)
        if(nt_code[i]!=ref_nt&&nt_code[i]!=read_nt) {
            for(j=0;j<4;j++) if(bit_nt[j]==100) break;
            bit_nt[i]=tmp; tmp=2;
        }

    //for(i=0;i<4;i++) cout<<" "<<nt_code[i]<<":"<<(int) bit_nt[i]; cout<<endl;

    for(i=0;i<256;i++) alphabet[i]=bit_nt[0];
    alphabet[(unsigned char)'c']=bit_nt[1]; alphabet[(unsigned char)'C']=bit_nt[1]; 
    alphabet[(unsigned char)'g']=bit_nt[2]; alphabet[(unsigned char)'G']=bit_nt[2];
    alphabet[(unsigned char)'t']=bit_nt[3]; alphabet[(unsigned char)'T']=bit_nt[3];

    for(i=0;i<256;i++) rev_alphabet[i]=bit_nt[3];
    rev_alphabet[(unsigned char)'c']=bit_nt[2]; rev_alphabet[(unsigned char)'C']=bit_nt[2]; 
    rev_alphabet[(unsigned char)'g']=bit_nt[1]; rev_alphabet[(unsigned char)'G']=bit_nt[1];
    rev_alphabet[(unsigned char)'t']=bit_nt[0]; rev_alphabet[(unsigned char)'T']=bit_nt[0];

    for(i=0;i<4;i++) useful_nt[bit_nt[i]]=nt_code[i];
    for(i=0;i<4;i++) useful_nt[bit_nt[i]+4]=tolower(nt_code[i]);

    //cout<<useful_nt<<endl;

    //for(i=0;i<4;i++) cout<<" "<<(char) nt_code[i]<<":"<<(int) alphabet[nt_code[i]]; cout<<endl;
    //for(i=0;i<4;i++) cout<<" "<<(char) tolower(nt_code[i])<<":"<<(int) tolower(alphabet[nt_code[i]]); cout<<endl;

    //for(i=0;i<4;i++) cout<<" "<<(char) nt_code[i]<<":"<<(int) rev_alphabet[nt_code[i]]; cout<<endl;
    //for(i=0;i<4;i++) cout<<" "<<(char) tolower(nt_code[i])<<":"<<(int) tolower(rev_alphabet[nt_code[i]]); cout<<endl;
}

char chain_flag[2] = {'+', '-'};

char revnt_code[4] = {'T', 'G', 'C', 'A'};



