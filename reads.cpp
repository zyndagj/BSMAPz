#include "reads.h"

using namespace std;

extern Param param;

ReadClass::ReadClass()
{
	_index=0;
	mreads.resize(BatchNum);
}

void ReadClass::InitIndex(ifstream &fin, igzstream &gzfin, string filename) {
	string s1,s2,s3,s4;
	char ch[1000];
	bit32_t i, maxi;

    //cout<<"file format:"<<_file_format<<endl;
	switch(param.input_format) {
	case 0: //fasta
	case 1: //fastq
		maxi=(param.read_start-1)*(2+2*param.input_format);
		if(param.gz_input) {
			gzfin.open(filename.c_str());
			for(i=0;!gzfin.eof()&&i<maxi;i++) gzfin.getline(ch,1000);
		}
		else {
			fin.open(filename.c_str());
			for(i=0;!fin.eof()&&i<maxi;i++) fin.getline(ch,1000);
		}
		break;
    case 2: //sam
	case 3: //bam
		SAM_fp=samopen(filename.c_str(), "rb", 0);
		SAM_b=bam_init1();
		maxi=(param.read_start-1)*(1+param.pairend);
		for(i=0;i<maxi;i++) if(samread(SAM_fp,SAM_b)<0) break;
	}
	_index=param.read_start-1;
}

int ReadClass::LoadBatchReads(ifstream &fin, igzstream &gzfin, int readset)  // readset {0: single-end, 1:pair-end set1, 2:pair-end set2)  
{
	char ch[1000];
	char c;
	vector<ReadInf>::iterator p=mreads.begin();

	size_t i,l_seq;
	char *s, *t;
	
	if (param.input_format<2) //.fa and .fq format
		if(param.gz_input){
            for(num=0; num<BatchNum && _index<param.read_end; p++,num++,_index++){
				gzfin>>c; if(gzfin.eof()) break;
				p->index=_index;
				gzfin>>p->name; gzfin.getline(ch,1000);
				gzfin>>p->seq; p->readset=readset;
                if(param.input_format) {//*.fq
					gzfin>>ch; fin.getline(ch, 1000);
					gzfin>>p->qual;
				}
                else p->qual=string(p->seq.size(), param.zero_qual+param.default_qual);
				if(p->seq.size()>param.max_readlen) {
					p->seq.erase(param.max_readlen); p->qual.erase(param.max_readlen);
				}
            }                                                                                                                                                                                                            		
		}
		else{
			for(num=0; num<BatchNum && _index<param.read_end; p++,num++,_index++){
    			fin>>c; if(fin.eof()) break;
    			p->index=_index;
	    		fin>>p->name; fin.getline(ch,1000);
	    		fin>>p->seq; p->readset=readset;
	    		if(param.input_format) {//*.fq
	    			fin>>ch; fin.getline(ch, 1000);
	    			fin>>p->qual;
	    		}
	    		else p->qual=string(p->seq.size(), param.zero_qual+param.default_qual);
            	if(p->seq.size()>param.max_readlen) {
                	p->seq.erase(param.max_readlen); p->qual.erase(param.max_readlen);
				}
			}
	    }
	else //SAM/BAM format
    	for(num=0; num<BatchNum && _index<param.read_end; num++,p++,_index++){
		//cout<< "index:"<<_index<<endl;
            //cout<<"num:"<<num<<"  mode:"<<mode<<endl;
            if(readset==2) if (samread(SAM_fp,SAM_b)<0) break;
    		if (samread(SAM_fp,SAM_b)<0) break;
      		p->index=_index;
       		p->name=string((char*)bam1_qname(SAM_b));
           //  
       		l_seq=min((bit32_t)SAM_b->core.l_qseq,(bit32_t)param.max_readlen);
       		p->seq.assign(l_seq,0); p->qual.assign(l_seq,0);
       		s=(char*) bam1_seq(SAM_b);  t=(char*) bam1_qual(SAM_b);
            if(readset) {
                if(SAM_b->core.flag&0x40) p->readset=1;
                else if(SAM_b->core.flag&0x80) p->readset=2;
                else p->readset=readset;
            }
            else p->readset=readset;
       		for(i=0;i<l_seq;i++){
       		    p->seq[i]=bam_nt16_rev_table[bam1_seqi(s,i)];
       		    p->qual[i]=t[i]+33;
       		}
			//cout<<p->name<<" "<<p->seq<<endl;
            if(readset==1) if (samread(SAM_fp,SAM_b)<0) break;
		}

	return num;
}

