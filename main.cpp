#include <iomanip>
#include <limits.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "reads.h"
#include "dbseq.h"
#include "align.h"
#include "param.h"
#include "pairs.h"
#include "utilities.h"

#ifdef THREAD
#include <pthread.h>
#endif

using namespace std;

//global variables
Param param;
string query_a_file;
string query_b_file;
string ref_file;
string out_align_file;
string command_line;

ifstream fin_db; igzstream gzfin_db;
ifstream fin_a; igzstream gzfin_a;
ifstream fin_b; igzstream gzfin_b;
ofstream fout;
FILE *pout;
ReadClass read_a, read_b;
RefSeq refseq_ref;

bit32_t n_aligned=0, n_unique=0, n_multiple=0;   //number of reads aligned
bit32_t n_aligned_pairs=0, n_unique_pairs=0, n_multiple_pairs=0;  //number of pairs aligned
bit32_t n_aligned_a=0, n_unique_a=0, n_multiple_a=0;  //number of a reads aligned
bit32_t n_aligned_b=0, n_unique_b=0, n_multiple_b=0;  //number of b reads aligned
bit32_t ref_time, read_time;
bit16_t tid[64];
char version[] = "1.1";
ostringstream message;

void info(int level) {
	if(level<=param.verbose_level) cerr<<message.str();
	message.str("");
} 

#ifdef THREAD
pthread_mutex_t mutex_fin=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_fout=PTHREAD_MUTEX_INITIALIZER;

void *t_SingleAlign(void *tid) {
	SingleAlign a;
	int n;
	bit32_t cur_at;
	while(1) {
		pthread_mutex_lock(&mutex_fin);
		n=read_a.LoadBatchReads(fin_a,gzfin_a,0);
		cur_at=read_a._index;
		a.ImportBatchReads(read_a.num, read_a.mreads);
		pthread_mutex_unlock(&mutex_fin);
		if(!n) break;
		a.Do_Batch(refseq_ref);
		pthread_mutex_lock(&mutex_fout);
		if(param.stdout) cout<<a._str_align; 
		else if(param.pipe_out) {fwrite(a._str_align.c_str(),1,a._str_align.size(),pout); fflush(pout);}
		else fout<<a._str_align;
		message<<"[bsmapz] @"<<Curr_Time()<<" \t"<<cur_at-param.read_start+1<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl; info(2);
		pthread_mutex_unlock(&mutex_fout);
	}
	pthread_mutex_lock(&mutex_fout);
	n_aligned+=a.n_aligned; n_unique+=a.n_unique; n_multiple+=a.n_multiple;
	pthread_mutex_unlock(&mutex_fout);
	read_time+=Cal_AllTime()-ref_time;
	return NULL;
};

void Do_SingleAlign() {
	vector<pthread_t> pthread_ids(param.num_procs);
	for(int i=0; i<param.num_procs; i++) {
		tid[i]=i; pthread_create(&pthread_ids[i], NULL, t_SingleAlign, (void*)&tid[i]);
	}
	for(int i=0; i<param.num_procs; i++) pthread_join(pthread_ids[i], NULL);
};


void *t_PairAlign(void *tid) {
	PairAlign a;
	int n1, n2;
	bit32_t cur_at;
	while(1) {
		pthread_mutex_lock(&mutex_fin);
		n1=read_a.LoadBatchReads(fin_a,gzfin_a,1);
		n2=read_b.LoadBatchReads(fin_b,gzfin_b,2);
		cur_at=read_a._index;
		a.ImportBatchReads(n1, read_a.mreads, read_b.mreads);
		pthread_mutex_unlock(&mutex_fin);
		if(!n1||(n1!=n2)) break;
		a.Do_Batch(refseq_ref);
		pthread_mutex_lock(&mutex_fout);
		if(param.stdout) cout<<a._str_align; 
		else if(param.pipe_out) {fwrite(a._str_align.c_str(),1,a._str_align.size(),pout); fflush(pout);}
		else fout<<a._str_align;
		message<<"[bsmapz] @"<<Curr_Time()<<" \t"<<cur_at-param.read_start+1<<" read pairs finished. "<<Cal_AllTime()<<" secs passed"<<endl;info(2);
		pthread_mutex_unlock(&mutex_fout);		
	}
	pthread_mutex_lock(&mutex_fout);
	n_aligned_pairs+=a.n_aligned_pairs; n_unique_pairs+=a.n_unique_pairs; n_multiple_pairs+=a.n_multiple_pairs;
	n_aligned_a+=a.n_aligned_a; n_unique_a+=a.n_unique_a; n_multiple_a+=a.n_multiple_a;
	n_aligned_b+=a.n_aligned_b; n_unique_b+=a.n_unique_b; n_multiple_b+=a.n_multiple_b;
	pthread_mutex_unlock(&mutex_fout);
	read_time+=Cal_AllTime()-ref_time;
	return NULL;		
};

void Do_PairAlign() {
	//if(param.max_snp_num>0) param.max_snp_num=2;
	vector<pthread_t> pthread_ids(param.num_procs);
	//create
	//cout <<param.num_procs<<"num_procs\n";
	for(int i=0; i<param.num_procs; i++) {
		tid[i]=i; pthread_create(&pthread_ids[i], NULL, t_PairAlign, (void*)&tid[i]);
	}
	for(int i=0; i<param.num_procs; i++) pthread_join(pthread_ids[i], NULL);
};

void* wrapper_CalKmerFreq0(void*) {refseq_ref.t_CalKmerFreq(0); return NULL;}
void* wrapper_CalKmerFreq1(void*) {refseq_ref.t_CalKmerFreq(1); return NULL;}
void* wrapper_FillIndex0(void*) {refseq_ref.t_FillIndex(0); return NULL;}
void* wrapper_FillIndex1(void*) {refseq_ref.t_FillIndex(1); return NULL;}

void Do_Formatdb() {
	if(param.RRBS_flag) refseq_ref.CreateIndex();
	else {
		pthread_t t0, t1;
		refseq_ref.InitialIndex();
		pthread_create(&t0, NULL, wrapper_CalKmerFreq0, NULL);
		pthread_create(&t1, NULL, wrapper_CalKmerFreq1, NULL);
		pthread_join(t0, NULL); pthread_join(t1, NULL);
		refseq_ref.AllocIndex();
		pthread_create(&t0, NULL, wrapper_FillIndex0, NULL);
		pthread_create(&t1, NULL, wrapper_FillIndex1, NULL);
		pthread_join(t0, NULL); pthread_join(t1, NULL);
		refseq_ref.FinishIndex();
	}
	message<<"[bsmapz] @"<<Curr_Time()<<" \tcreate seed table. "<<Cal_AllTime()<<" secs passed\n"; info(1);
};

#else
void Do_SingleAlign() {
	SingleAlign a;
	while(read_a.LoadBatchReads(fin_a,gzfin_a,0)) {
		a.ImportBatchReads(read_a.num, read_a.mreads);
		a.Do_Batch(refseq_ref);
		if(param.stdout) cout<<a._str_align; 
		else if(param.pipe_out) {fwrite(a._str_align.c_str(),1,a._str_align.size(),pout); fflush(pout);}
		else fout<<a._str_align;
		message<<"[bsmapz] @"<<Curr_Time()<<" \t"<<read_a._index-param.read_start+1<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl; info(2);
	}
	n_aligned=a.n_aligned; n_unique=a.n_unique; n_multiple=a.n_multiple;	
	read_time+=Cal_AllTime()-ref_time;
};

void Do_PairAlign() {
	PairAlign a;
	int n1, n2;
	while(1) {
		n1=read_a.LoadBatchReads(fin_a,gzfin_a,1);
		n2=read_b.LoadBatchReads(fin_b,gzfin_b,2);
		if(!n1||(n1!=n2))
			break;
		a.ImportBatchReads(n1, read_a.mreads, read_b.mreads);
		a.Do_Batch(refseq_ref);		
		if(param.stdout) cout<<a._str_align; 
		else if(param.pipe_out) {fwrite(a._str_align.c_str(),1,a._str_align.size(),pout); fflush(pout);}
		else fout<<a._str_align;
		message<<"[bsmapz] @"<<Curr_Time()<<" \t"<<read_a._index-param.read_start+1<<" read pairs finished. "<<Cal_AllTime()<<" secs passed"<<endl; info(2);
	}	

    n_aligned_pairs+=a.n_aligned_pairs; n_unique_pairs+=a.n_unique_pairs; n_multiple_pairs+=a.n_multiple_pairs;
    n_aligned_a+=a.n_aligned_a; n_unique_a+=a.n_unique_a; n_multiple_a+=a.n_multiple_a;
    n_aligned_b+=a.n_aligned_b; n_unique_b+=a.n_unique_b; n_multiple_b+=a.n_multiple_b;
	read_time+=Cal_AllTime()-ref_time;
};

void Do_Formatdb() {
	refseq_ref.CreateIndex();
	message<<"[bsmapz] @"<<Curr_Time()<<" \tcreate seed table. "<<Cal_AllTime()<<" secs passed\n"; info(1);
};

#endif

//usage
void usage(void)
{
cerr<<"Usage:	bsmapz [options]\n"
		<<"       -a  <str>   query a file, FASTA/FASTQ/BAM format\n"
		<<"       -d  <str>   reference sequences file, FASTA format\n"
		<<"       -o  <str>   output alignment file, BSP/SAM/BAM format, if omitted, the output will be written to STDOUT in SAM format.\n"
		<<"\n  Options for alignment:\n"
		<<"       -s  <int>   seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.\n"
		<<"       -v  <float> if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.\n"
		<<"                   otherwise it's interpreted as the maximum number of mismatches allowed on a read, <="<<MAXSNPS<<".\n"
		<<"                   example: -v 5 (max #mismatches = 5), -v 0.1 (max #mismatches = read_length * 10%)\n" 
		<<"                   default="<<(param.max_snp_num-100)/100.0<<".\n"
		<<"       -g  <int>   gap size, BSMAPz only allows 1 continuous gap (insert or deletion) with up to "<<MAXGAPS<<" nucleotides\n"
		<<"                   default="<<param.gap<<"\n"
		<<"       -w  <int>   maximum number of equal best hits to count, <="<<MAXHITS<<"\n"
		<<"       -3          using 3-nucleotide mapping approach\n" 
        <<"       -B  <int>   start from the Nth read or read pair, default: 1\n"
        <<"       -E  <int>   end at the Nth read or read pair, default: 4,294,967,295\n"
       	<<"       -I  <int>   index interval, default="<<param.index_interval<<"\n"        
       	<<"       -k  <float> set the cut-off ratio for over-represented kmers, default="<<param.max_kmer_ratio<<"\n"
       	<<"                   example: -k 1e-6 means the top 0.0001\% over-represented kmer will be skipped in alignment\n"
#ifdef THREAD
		<<"       -p  <int>   number of processors to use, default="<<param.num_procs<<"\n"
#endif	
        <<"       -D  <str>   activating RRBS mapping mode and set restriction enzyme digestion sites. \n"
        <<"                   digestion position marked by \'-\', example: -D C-CGG for MspI digestion.\n"
        <<"                   default: none (whole genome shotgun bisulfite mapping mode)\n"        
    	<<"       -S  <int>   seed for random number generation used in selecting multiple hits\n"
        <<"                   other seed values generate pseudo random number based on read index number, to allow reproducible mapping results. \n"
        <<"                   default="<<param.randseed<<". (get seed from system clock, mapping results not resproducible.)\n"
        <<"       -n  [0,1]   set mapping strand information. default: "<<param.chains<<"\n"
        <<"                   -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+), \n" 
        <<"                   for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.\n"
        <<"                   -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, -- \n"
        <<"       -M  <str>   set alignment information for the additional nucleotide transition. \n"
        <<"                   <str> is in the form of two different nucleotides N1N2, \n"
        <<"                   indicating N1 in the reads could be mapped to N2 in the reference sequences.\n"
        <<"                   default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion. \n"
        <<"                   example: -M GA could be used to detect A=>I(G) transition in RNA editing. \n"
		<<"\n  Options for trimming:\n"
		<<"       -q  <int>   quality threshold in trimming, 0-40, default=0 (no trim)\n"
		<<"       -z  <int>   base quality, default="<<(int) param.zero_qual<<" [Illumina is using 64, Sanger Institute is using 33]\n"
		<<"       -f  <int>   filter low-quality reads containing >n Ns, default="<<param.max_ns<<"\n"
        <<"       -A  <str>   3-end adapter sequence, default: none (no trim)\n"
        <<"       -L  <int>   map the first N nucleotides of the read, default:"<<param.max_readlen<<" (map the whole read).\n"
		<<"\n  Options for reporting:\n"
        <<"       -r  [0,1,2] how to report repeat hits, 0=none(unique hit/pair); 1=random one; 2=all(slow), default:"<<param.report_repeat_hits<<".\n"
        <<"       -R          print corresponding reference sequences in SAM output, default=off\n"
        <<"       -u          report unmapped reads, default=off\n"
		<<"       -H          do not print header information in SAM format output\n"
		<<"       -V  [0,1,2] verbose level: 0=no message displayed (quiet mode); \n"
		<<"                   1=major message (default); 2=detailed message.\n" 
		<<"\n  Options for pair-end alignment:\n"
		<<"       -b  <str>   query b file\n"
		<<"       -m  <int>   minimal insert size allowed, default="<<param.min_insert<<"\n"
		<<"       -x  <int>   maximal insert size allowed, default="<<param.max_insert<<"\n\n"

		<<"       -h          help\n\n";
	exit(1);
};

int mGetOptions(int rgc, char *rgv[])
{
	//[options]
	int i;
	for(i=1,command_line=rgv[0];i<rgc;i++) command_line=command_line+" "+rgv[i];
	for(i=1; i<rgc; i++) {
		if(rgv[i][0]!='-') return i;		
		switch(rgv[i][1]) {
			case 'a': if(rgv[i][2]==0) query_a_file = rgv[++i]; else if(rgv[i][2]=='=') query_a_file=rgv[i]+3; else return i; break;
			case 'b': if(rgv[i][2]==0) query_b_file = rgv[++i]; else if(rgv[i][2]=='=') query_b_file=rgv[i]+3; else return i; 
			    param.pairend=1; break;
			case 'd': if(rgv[i][2]==0) ref_file = rgv[++i]; else if(rgv[i][2]=='=') ref_file=rgv[i]+3; else return i; break;
			case 's': if(rgv[i][2]==0) 
				param.SetSeedSize(atoi(rgv[++i])); else if(rgv[i][2]=='=') param.SetSeedSize(atoi(rgv[i]+3)); else return i;
				//if(param.RRBS_flag) param.SetSeedSize(12); 
				break;
			case 'o': 
				if(rgv[i][2]==0) out_align_file = rgv[++i];
				else if(rgv[i][2]=='=') out_align_file=rgv[i]+3; 
				else return i;
				param.stdout=0; 
				break;
			case 'm': if(rgv[i][2]==0) param.min_insert = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.min_insert=atoi(rgv[i]+3); else return i; break;
			//case 'n': if(rgv[i][2]==0) param.chains = atoi(rgv[++i])%4; else if(rgv[i][2]=='=') param.chains=atoi(rgv[i]+3)%4; else return i; break;
			case 'n': if(rgv[i][2]==0) param.chains=(atoi(rgv[++i])!=0); else if(rgv[i][2]=='=') param.chains=(atoi(rgv[i]+3)!=0); else return i; break;
			case 'g': if(rgv[i][2]==0) param.gap=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.gap=atoi(rgv[i]+3); else return i; 
			    if(param.gap>MAXGAPS) {
			        cerr<<"warning: gap length exceeds max value:"<<MAXGAPS<<endl;
			        param.gap=MAXGAPS;
			    }
			    break;
			case 'x': if(rgv[i][2]==0) param.max_insert = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_insert=atoi(rgv[i]+3); else return i; break;
			case 'r': if(rgv[i][2]==0) param.report_repeat_hits = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.report_repeat_hits=atoi(rgv[i]+3); else return i;
				if(param.report_repeat_hits>2||param.report_repeat_hits<0) {
					cerr<<"invalid -r value: "<<param.report_repeat_hits<<", must be 0, 1, or 2.\n";
					exit(1);
				} 
				break;
			case 'V': if(rgv[i][2]==0) param.verbose_level = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.verbose_level=atoi(rgv[i]+3); else return i; 
				if(param.verbose_level>2||param.verbose_level<0) {
                    cerr<<"invalid -V value: "<<param.verbose_level<<", must be 0, 1, or 2.\n"; 
                    exit(1);
                }
				break;
			case 'I': if(rgv[i][2]==0) param.index_interval = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.index_interval=atoi(rgv[i]+3); else return i; 
			    if(param.RRBS_flag) param.index_interval=1;
			    if(param.index_interval>16) {cerr<<"index interval exceeds max value:16\n"; exit(1);}
			    break;			
			case 'k': if(rgv[i][2]==0) param.max_kmer_ratio=atof(rgv[++i]); else if(rgv[i][2]=='=') param.max_kmer_ratio=atof(rgv[i]+3); else return i; break;
			case 'v': 
                double tmp_max_snp; 
                if(rgv[i][2]==0) tmp_max_snp = atof(rgv[++i]); else if(rgv[i][2]=='=') tmp_max_snp=atof(rgv[i]+3); else return i;
                if(tmp_max_snp<1.0) {
                    param.max_snp_num=(int)(tmp_max_snp*100+0.5)+100;
                    if(param.max_snp_num==100) param.max_snp_num=0;
                }
                else{
                    param.max_snp_num=(int)(tmp_max_snp+0.5);
                    if(param.max_snp_num>MAXSNPS) {
                        cerr<<"warning: number of mismatches exceeds max value:"<<MAXSNPS<<endl;
                        param.max_snp_num=MAXSNPS;
                    }
                }
                break;			
			case 'w': if(rgv[i][2]==0) param.max_num_hits = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_num_hits=atoi(rgv[i]+3); else return i;
			    if(param.max_num_hits>MAXHITS) {cerr<<"number of multi-hits exceeds max value:"<<MAXHITS<<endl; exit(1);}
			    break;
			case 'q': if(rgv[i][2]==0) param.qual_threshold = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.qual_threshold = atoi(rgv[i]+3); else return i; break;
			case 'f': if(rgv[i][2]==0) param.max_ns=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_ns=atoi(rgv[i]+3); else return i; break;
			case 'z': if(rgv[i][2]==0) param.zero_qual=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.zero_qual=atoi(rgv[i]+3); else return i; break;
			case 'p': if(rgv[i][2]==0) param.num_procs=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.num_procs=atoi(rgv[i]+3); else return i; break;
			case 'A': if(rgv[i][2]==0) param.adapter[param.n_adapter++]=rgv[++i]; else if(rgv[i][2]=='=') param.adapter[param.n_adapter++]=rgv[i]+3; else return i;   
                break;
            case 'R': if(rgv[i][2]==0) param.out_ref=1; else return i; break;
            case '3': if(rgv[i][2]==0) param.nt3=1; else return i; break;
            case 'H': if(rgv[i][2]==0) param.sam_header=0; else return i; break;
            case 'u': if(rgv[i][2]==0) param.out_unmap=1; else return i; break;
		    case 'B': if(rgv[i][2]==0) param.read_start = max(atoi(rgv[++i]),1); else if(rgv[i][2]=='=') param.read_start=max(atoi(rgv[i]+3),1); else return i; break;
		    case 'E': if(rgv[i][2]==0) param.read_end = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.read_end=atoi(rgv[i]+3); else return i; break;
	        case 'D': if(rgv[i][2]==0) param.SetDigestionSite(rgv[++i]); else if(rgv[i][2]=='=') param.SetDigestionSite(rgv[i]+3); else return i; break;
	        case 'M': if(rgv[i][2]==0) {i++; param.SetAlign(rgv[i][0], rgv[i][1]);}
                      else if(rgv[i][2]=='=') param.SetAlign(rgv[i][3], rgv[i][4]); 
                      else return i; 
                      break;
            case 'L': if(rgv[i][2]==0) param.max_readlen = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_readlen = atoi(rgv[i]+3); else return i; break;
            case 'N': if(rgv[i][2]==0) param.N_mis=1; else return i; break;
            case 'S': if(rgv[i][2]==0) param.randseed = atoi(rgv[++i]);  else if(rgv[i][2]=='=') param.randseed = atoi(rgv[i]+3); else return i; break;	        
			case 'h':usage();   //usage information
            default: return i;
		}
	}
    param.InitMapping();
	return 0;
}

void check_ofile(string &filename, string &err_msg) {
	ofstream ff;                       
	ff.open(filename.c_str()); cerr.flush();
	if(!ff) {   
		cerr<<endl<<err_msg<<filename<<endl;                                              
		exit(1);
	}
}

int check_ifile(string &filename, string &err_msg) {
	int byte1, byte2; ifstream ff;
    ff.open(filename.c_str()); cerr.flush();
	if(!ff) {
		cerr<<endl<<err_msg<<filename<<endl;
		exit(1);
    }                	
	byte1=ff.get(); byte2=ff.get();
	return ((byte1==0x1f)&&(byte2==0x8b));  // if input is gzip file or not       
}

int check_ifile_format(string &filename, int gz_flag) {
	ifstream ff; igzstream gg;
	string s1,s2,s3,s4; char ch[1000];
	if(gz_flag) {
		gg.open(filename.c_str());
        gg>>s1; gg.getline(ch, 1000);
		if(s1[0]=='>') {message<<" \t(format: gzipped FASTA)\n"; info(1); return 0;} //fasta
		if(s1[0]=='@') {message<<" \t(format: gzipped FASTQ)\n"; info(1); return 1;} //fastq
		if(samopen(filename.c_str(), "rb", 0)!=0) {message<<" \t(format: BAM)\n"; info(1); return 3;} //BAM         
        if(samopen(filename.c_str(), "r", 0)!=0) {message<<" \t(format: SAM)\n"; info(1); return 2;} //SAM
	}	 
	else {
		ff.open(filename.c_str());
        ff>>s1; ff.getline(ch, 1000);
		if(s1[0]=='>') {message<<" \t(format: FASTA)\n"; info(1); return 0;} //fasta
		if(s1[0]=='@') {message<<" \t(format: FASTQ)\n"; info(1); return 1;} //fastq
        if(samopen(filename.c_str(), "rb", 0)!=0) {message<<" \t(format: BAM)\n"; info(1); return 3;} //BAM   		
        if(samopen(filename.c_str(), "r", 0)!=0) {message<<" \t(format: SAM)\n"; info(1); return 2;} //SAM
	}
	cerr<<"\t(format: unknown)\nUnknown input format.\n";
	exit(1);
}

void RunProcess(void) {
	char _ch[256]; string _str="@HD\tVN:1.0\n";
	if(out_align_file.size()>4){
	    if(out_align_file.compare(out_align_file.size()-4,4,".sam")==0) param.out_sam=1;
	    else if (out_align_file.compare(out_align_file.size()-4,4,".bam")==0) param.out_sam=2;
    }
    else param.out_sam=1;

    if(param.max_snp_num<100) message<<"\tmax number of mismatches: "<<param.max_snp_num;
    else message<<"\tmax number of mismatches: read_length * "<<param.max_snp_num-100<<"% ";
    message<<" \tmax gap size: "<<param.gap<<" \tkmer cut-off ratio: "<<param.max_kmer_ratio<<" \tmax multi-hits: "<<param.max_num_hits<<endl;
    message<<"\tquality cutoff: "<<(int)param.qual_threshold<<" \tbase quality char: '"<<param.zero_qual<<"'"<<" \tmax Ns: "<<param.max_ns;
    if(param.max_readlen<(FIXELEMENT-1)*SEGLEN) message<<" \thard clip length: "<<param.max_readlen; message<<endl;
    if(param.pairend||param.RRBS_flag) message<<"\tmin fragment size: "<<param.min_insert<<" \tmax fragemt size: "<<param.max_insert<<endl;
    if(param.nt3) message<<"\t3-nucleotide alignment mode: "<<param.read_nt<<" <=> "<<param.ref_nt;
    else message<<"\twildcard alignment mode: "<<param.read_nt<<" in reads => "<<param.ref_nt<<" in reference";
    message<<" \tseed size: "<<param.seed_size<<" \tindex interval: "<<param.index_interval<<endl;
    message<<"\talignment report: "; 
    if(param.report_repeat_hits==0) message<<"unique alignments only";
    if(param.report_repeat_hits==1) message<<"unique + one random multiple alignments";
    if(param.report_repeat_hits==2) message<<"unique + all multiple alignments"; message<<endl;
    if(param.pairend){
        message<<"\tmapping strand (read_1): ++,-+";
        if(param.chains) message<< ",+-,--";
        message<<"\tmapping strand (read_2): +-,--";
        if(param.chains) message<< ",++,-+";
    }
    else {
        message<<"\tmapping strand: ++,-+";
        if(param.chains) message<< ",+-,--";
    }
    message<<endl;
    for(bit32_t i=0; i<param.n_adapter;i++) message<<"\tadapter sequence #"<<i+1<<": "<<param.adapter[i]<<endl;
    if(param.RRBS_flag) {
        message<<"\tRRBS digestion site: ";
        for(bit32_t i=0; i<param.digest_site.size(); i++) 
        	message<<param.digest_site[i].substr(0,param.digest_pos[i])<<'-'<<param.digest_site[i].substr(param.digest_pos[i])<<" ";
		message<<endl;
	}
	info(2);
	string err_msg;
	//pair-end alignment    
	if(param.pairend) {
		message<<"[bsmapz] @"<<Curr_Time()<<" \tPair-end alignment("<<param.num_procs<<" threads),"; info(1);
		message<<" \tstart from read pair #"<<param.read_start;
		if(~param.read_end) message<<" \tend at read pair #"<<param.read_end; info(2); message<<endl; info(1);
		err_msg="failed to open read file #1 (check -a option): ";
        message<<"\tInput read file #1: "<<query_a_file;
		param.gz_input=check_ifile(query_a_file, err_msg);
		param.input_format=check_ifile_format(query_a_file, param.gz_input);
	    read_a.InitIndex(fin_a, gzfin_a, query_a_file);
        
		err_msg="failed to open read file #2 (check -b option): ";
        message<<"\tInput read file #2: "<<query_b_file;
		if(param.gz_input!=check_ifile(query_b_file, err_msg)) {
			cerr<<"Input read file #1 and #2 should be both regular files or gzip files."<<endl;
			exit(1);
		}
		if(param.input_format!=check_ifile_format(query_b_file, param.gz_input)) {
			cerr<<"Input read file #1 and #2 should be in same format."<<endl;
			exit(1);
		}
		read_b.InitIndex(fin_b, gzfin_b, query_b_file);

		if(param.stdout) message<<"\tOutput: STDOUT\t (format: SAM)\n";
		else {
			cerr<<"\tOutput file: "<<out_align_file;
			switch(param.out_sam) {
				case 0: message<<"\t (format: BSP)\n"; break;
				case 1: message<<"\t (format: SAM)\n"; break;
				case 2: message<<"\t (format: SAM, automatically convert to BAM)\n"; break;
			}
			err_msg="failed to open output file (check -o option): ";
			check_ofile(out_align_file, err_msg);
			if(param.out_sam==2) {
				string samtools_command="samtools view -bS - >"+out_align_file;
				if((pout=popen(samtools_command.c_str(),"w"))==NULL) {
					cerr.flush();
					message<<"unable to creat samtools pipe, using post alignment conversion instead.\n";
					fout.open(out_align_file.c_str());
				}
				else param.pipe_out=1;
			}
			else fout.open(out_align_file.c_str());
		}

		if(param.out_sam&&param.sam_header) {
			for(bit32_t i=0;i<refseq_ref.total_num;i++){
				sprintf(_ch,"@SQ\tSN:%s\tLN:%u\n",refseq_ref.title[i<<1].name.c_str(),refseq_ref.title[i<<1].size);
				_str.append(_ch);
			}
			sprintf(_ch,"@PG\tID:BSMAP\tVN:%s\tCL:\"%s\"\n",version,command_line.c_str()); _str.append(_ch);
			if(param.stdout) cout<<_str;
			else if(param.pipe_out) fwrite(_str.c_str(),1,_str.size(),pout);
			else fout<<_str;
		}
		
		n_aligned_pairs=n_aligned_a=n_aligned_b=0; info(1);
		Do_PairAlign();
		if(param.pipe_out) pclose(pout); else fout.close();
		if(param.input_format==3) {
		    bam_destroy1(read_a.SAM_b); samclose(read_a.SAM_fp);
		    bam_destroy1(read_b.SAM_b); samclose(read_b.SAM_fp);
        }
    
        message<<fixed<<setw(4)<<setprecision(1);
        message<<"[bsmapz] @"<<Curr_Time()<<" \ttotal read pairs: "<<read_a._index-param.read_start+1<<" \ttotal time consumed:  "<<Cal_AllTime()<<" secs\n";
		message<<"\taligned pairs: "<<n_aligned_pairs<<" ("<<100.0*n_aligned_pairs/(read_a._index-param.read_start+1)<<"%), ";
		message<<"unique pairs: "<<n_unique_pairs<<" ("<<100.0*n_unique_pairs/(read_a._index-param.read_start+1)<<"%), ";
		if(param.report_repeat_hits==0) message<<"suppressed ";
		message<<"non-unique pairs: "<<n_multiple_pairs<<" ("<<100.0*n_multiple_pairs/(read_a._index-param.read_start+1)<<"%)\n";

		message<<"\tunpaired read #1: "<<n_aligned_a<<" ("<<100.0*n_aligned_a/(read_a._index-param.read_start+1)<<"%), ";
        message<<"unique reads: "<<n_unique_a<<" ("<<100.0*n_unique_a/(read_a._index-param.read_start+1)<<"%), ";
        if(param.report_repeat_hits==0) message<<"suppressed ";
        message<<"non-unique reads: "<<n_multiple_a<<" ("<<100.0*n_multiple_a/(read_a._index-param.read_start+1)<<"%)\n";

		message<<"\tunpaired read #2: "<<n_aligned_b<<" ("<<100.0*n_aligned_b/(read_b._index-param.read_start+1)<<"%), ";
        message<<"unique reads: "<<n_unique_b<<" ("<<100.0*n_unique_b/(read_b._index-param.read_start+1)<<"%), ";
        if(param.report_repeat_hits==0) message<<"suppressed ";
        message<<"non-unique reads: "<<n_multiple_b<<" ("<<100.0*n_multiple_b/(read_b._index-param.read_start+1)<<"%)\n";
        info(1);
	}
	//single-read alignment
	else {
		message<<"[bsmapz] @"<<Curr_Time()<<" \tSingle-end alignment("<<param.num_procs<<" threads),"; info(1);
        message<<" \tstart from read #"<<param.read_start;
		if(~param.read_end) message<<" \tend at read #"<<param.read_end; info(2); cerr<<endl; info(1);
        err_msg="failed to open read file (check -a option): ";
		message<<"\tInput read file: "<<query_a_file;
		param.gz_input=check_ifile(query_a_file, err_msg);
		param.input_format=check_ifile_format(query_a_file, param.gz_input);
		read_a.InitIndex(fin_a, gzfin_a, query_a_file);

        if(param.stdout) message<<"\tOutput: STDOUT\t (format: SAM)\n";
        else {
        	message<<"\tOutput file: "<<out_align_file;
			switch(param.out_sam) {
				case 0: message<<"\t (format: BSP)\n"; break;
				case 1: message<<"\t (format: SAM)\n"; break;
				case 2: message<<"\t (format: SAM, automatically convert to BAM)\n"; break;
			}
	        err_msg="failed to open output file (check -o option): ";
			check_ofile(out_align_file, err_msg);
            if(param.out_sam==2) {
                string samtools_command="samtools view -bS - >"+out_align_file;
                if((pout=popen(samtools_command.c_str(),"w"))==NULL) {
                    message.flush();
                    message<<"unable to creat samtools pipe, using post alignment conversion instead.\n";
                    fout.open(out_align_file.c_str());
                }
                else param.pipe_out=1;
            }
			else fout.open(out_align_file.c_str());
        }
                                                		
		if(param.out_sam&&param.sam_header) {
    		char _ch[1000];
	    	for(bit32_t i=0;i<refseq_ref.total_num;i++) {
	    	    sprintf(_ch,"@SQ\tSN:%s\tLN:%u\n",refseq_ref.title[i<<1].name.c_str(),refseq_ref.title[i<<1].size);
	    	    _str.append(_ch);
	    	}
			sprintf(_ch,"@PG\tID:BSMAP\tVN:%s\tCL:\"%s\"\n",version,command_line.c_str()); _str.append(_ch);
            if(param.stdout) cout<< _str;
            else if(param.pipe_out) fwrite(_str.c_str(),1,_str.size(),pout);
			else fout<<_str;
	    }
		n_aligned=0;
		info(1);
		Do_SingleAlign();
		if(param.pipe_out) pclose(pout); else fout.close();
		if(param.input_format==3) {
		    bam_destroy1(read_a.SAM_b);
            samclose(read_a.SAM_fp);
        }
        message<<fixed<<setw(4)<<setprecision(1);
        message<<"[bsmapz] @"<<Curr_Time()<<" \ttotal reads: "<<read_a._index-param.read_start+1<<" \ttotal time:  "<<Cal_AllTime()<<" secs\n";
		message<<"\taligned reads: "<<n_aligned<<" ("<<100.0*n_aligned/(read_a._index-param.read_start+1)<<"%), ";
		message<<"unique reads: "<<n_unique<<" ("<<100.0*n_unique/(read_a._index-param.read_start+1)<<"%), ";  
		if(param.report_repeat_hits==0) message<<"suppressed ";
		message<<"non-unique reads: "<<n_multiple<<" ("<<100.0*n_multiple/(read_a._index-param.read_start+1)<<"%)\n";
		info(1);
	}
};

int main(int argc, char *argv[]) {
	//print usage
	if (argc == 1) usage();
	Initial_Time();
	int noptions=mGetOptions(argc, argv);
	if(noptions) {
        cerr<<"unknown option: "<<argv[noptions]<<endl;
        exit(noptions);
    }
    message<<"\nBSMAPz v"<<version<<endl; info(2); srand(time(NULL));
	string err_msg="failed to open reference file (check -d option): ";
	message<<"[bsmapz] @"<<Curr_Time()<<" \tloading reference file: "<<ref_file; info(1);
	param.gz_ref=check_ifile(ref_file, err_msg);
	if(check_ifile_format(ref_file, param.gz_ref)) {
		cerr<<"reference must be in FASTA format.\n";
		exit(1);
	}
	if(param.gz_ref) gzfin_db.open(ref_file.c_str());
	else fin_db.open(ref_file.c_str());

	refseq_ref.Run_ConvertBinseq(fin_db, gzfin_db);
	message<<"[bsmapz] @"<<Curr_Time()<<" \t"<<refseq_ref.total_num<<" reference seqs loaded, total size "<<refseq_ref.sum_length<<" bp. "<<Cal_AllTime()<<" secs passed"<<endl;
	info(1);			
	Do_Formatdb(); ref_time=Cal_AllTime(); read_time=0;
	RunProcess(); 
	/*
	cerr<<(read_a._index-param.read_start+1)<<" reads in "<<read_time<<" sec."<<endl;
	cerr.precision(1);
	cerr<<fixed<<1.0*(read_a._index-param.read_start+1)/read_time<<" read";
	if(param.pairend) cerr<<" pairs per sec."<<endl;
	else cerr<<"s per sec."<<endl;
	*/
    if(param.out_sam==2&&param.pipe_out==0){
		char sys_cmd[PATH_MAX+20], abs_bam_file[PATH_MAX];
		char *res=realpath(out_align_file.c_str(), abs_bam_file);
		if(res) {
		    sprintf(sys_cmd,"sam2bam.sh %s",abs_bam_file);
		    system(sys_cmd);
		} 
		else {
		    message<<"[bsmapz] warning: cannot get absolute path for sam2bam.sh input. alignment file remains in SAM format.\n"; 
		}
	}
	info(1);
    refseq_ref.ReleaseIndex();
	return 0;
}
