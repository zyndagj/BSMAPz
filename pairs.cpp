#include<cmath>
#include "pairs.h"

using namespace std;

extern Param param;

PairAlign::PairAlign()
{
	n_aligned_pairs=n_aligned_a=n_aligned_b=0;
    n_unique_pairs=n_unique_a=n_unique_b=0;
    n_multiple_pairs=n_multiple_a=n_multiple_b=0;
    pairhits = new PairArray[2*MAXSNPS+1];
    _str_align.reserve(BatchNum*1024);
    rand_rSeed=getpid()*time(NULL);
}

PairAlign::~PairAlign() {
    delete [] pairhits;
}

void PairAlign::ImportBatchReads(bit32_t n, vector<ReadInf> &a1, vector<ReadInf> &a2)
{
	_sa.ImportBatchReads(n, a1);
	_sb.ImportBatchReads(n, a2);
	num_reads=n;
}

int PairAlign::GetPairs(bit32_t na, bit32_t nb)
{
	bit32_t i,j,insert_size,bstart,bend,npair=0;
    ref_loc_t seg_start,seg_end;
    if(na>_sa.read_max_snp_num||nb>_sb.read_max_snp_num) return 0;
	PairHit pp;
	pp.na=na;pp.nb=nb;
	//_sa.SortHits4PE(na);
	//_sb.SortHits4PE(nb);
    ref_id_t chra;

    //cout<<"na:"<<na<<"    ";
    //for(i=0;i<=_sa.read_max_snp_num;++i) cout<<_sa._cur_n_hit[i]<<" "<<_sa._cur_n_chit[i]<<"    "; cout<<endl;
    //cout<<"f:"<<_sa._cur_n_hit[na]<<endl;
    //for(i=0;i<_sa._cur_n_hit[na]; i++) cout<<"chr"<<(int)_sa.hits[na][i].chr<<":"<<_sa.hits[na][i].loc<<endl;
    //cout<<"r:"<<_sa._cur_n_chit[na]<<endl;
    //for(i=0;i<_sa._cur_n_chit[na]; i++) cout<<"chr"<<(int)_sa.chits[na][i].chr<<":"<<_sa.chits[na][i].loc<<endl;
        
    //cout<<"nb:"<<nb<<"    ";
    //for(i=0;i<=_sb.read_max_snp_num;++i) cout<<_sb._cur_n_hit[i]<<" "<<_sb._cur_n_chit[i]<<"    "; cout<<endl;
    //cout<<"f:"<<_sb._cur_n_hit[nb]<<endl;
    //for(i=0;i<_sb._cur_n_hit[nb]; i++) cout<<"chr"<<(int)_sb.hits[nb][i].chr<<":"<<_sb.hits[nb][i].loc<<endl;
    //cout<<"r:"<<_sb._cur_n_chit[nb]<<endl;
    //for(i=0;i<_sb._cur_n_chit[nb]; i++) cout<<"chr"<<(int)_sb.chits[nb][i].chr<<":"<<_sb.chits[nb][i].loc<<endl;

	//a+ vs b-
	pp.chain=0; checked_pair_mismatch[na][nb]=1;
    chra=~0; bstart=0; bend=0;
	//cout<<"a0b1\n";
	for(i=0;i<_sa._cur_n_hit[na]; i++) {
        if(chra!=_sa.hits[na][i].chr){
            chra=_sa.hits[na][i].chr;
            for(bstart=bend;bstart<_sb._cur_n_chit[nb];bstart++) if(_sb.chits[nb][bstart].chr>=chra) break;
            for(bend=bstart;bend<_sb._cur_n_chit[nb];bend++) if(_sb.chits[nb][bend].chr>chra) break;
        }
        //cout<<"\tchr:"<<(int)chra<<"\tbstart:"<<bstart<<"\tbend:"<<bend<<endl;
        for(j=bstart;j<bend;j++) { 
		    //cout<<"a:"<<i<<" b:"<<j<<endl;
            if(chra&1) {seg_start=_sb.chits[nb][j].loc; seg_end=_sa.hits[na][i].loc+_sa._pread->seq.size();}
            else {seg_start=_sa.hits[na][i].loc; seg_end=_sb.chits[nb][j].loc+_sb._pread->seq.size();}
            insert_size=seg_end-seg_start;
            //cout<<"seg_start:"<<seg_start<<" seg_end:"<<seg_end<<endl;
		    //cout << "a: "<< _sa.hits[na][i].loc << "\tb: "<<_sb.chits[nb][j].loc<<"\tinsert: "<<insert_size<<endl;
	        if((insert_size>=param.min_insert)&&(insert_size<=param.max_insert)) {
	            //cout << "a: "<< (int)_sa.hits[na][i].chr << "\tb: "<<(int)_sb.chits[nb][j].chr<<"\tinsert: "<<insert_size<<endl;
	            pp.a=_sa.hits[na][i]; pp.b=_sb.chits[nb][j];
    			pp.insert=insert_size;
	    		pairhits[na+nb][_cur_n_hits[na+nb]++]=pp;
	    		npair++;
                if(_cur_n_hits[na+nb]>=param.max_num_hits) return npair;
            }
		}
	}

	//a- vs b+
	pp.chain=1; chra=~0; bstart=0; bend=0;
	//cout<<"a1b0\n";
	for(i=0;i<_sa._cur_n_chit[na]; i++) {
        if(chra!=_sa.chits[na][i].chr){
            chra=_sa.chits[na][i].chr;
            for(bstart=bend;bstart<_sb._cur_n_hit[nb];bstart++) if(_sb.hits[nb][bstart].chr>=chra) break;
            for(bend=bstart;bend<_sb._cur_n_hit[nb];bend++) if(_sb.hits[nb][bend].chr>chra) break;            
        }
        //cout<<"\tchr:"<<chra<<"\tbstart:"<<bstart<<"\tbend:"<<bend<<endl;
        for(j=bstart;j<bend;j++) { 
		    //cout<<"a:"<<i<<" b:"<<j<<endl;
            if((chra&1)==0) {seg_start=_sb.hits[nb][j].loc; seg_end=_sa.chits[na][i].loc+_sa._pread->seq.size();}
            else {seg_start=_sa.chits[na][i].loc; seg_end=_sb.hits[nb][j].loc+_sb._pread->seq.size();}
            insert_size=seg_end-seg_start;
            //cout<<"seg_start:"<<seg_start<<" seg_end:"<<seg_end<<endl;
		    //cout << "a: "<< _sa.chits[na][i].loc << "\tb: "<<_sb.hits[nb][j].loc<<"\tinsert: "<<insert_size<<endl;
	        if((insert_size>=param.min_insert)&&(insert_size<=param.max_insert)) {
	            //cout << "a: "<< (int)_sa.chits[na][i].chr << "\tb: "<<(int)_sb.hits[nb][j].chr<<"\tinsert: "<<insert_size<<endl;
	            pp.a=_sa.chits[na][i]; pp.b=_sb.hits[nb][j];
    			pp.insert=insert_size;
	    		pairhits[na+nb][_cur_n_hits[na+nb]++]=pp;
	    		npair++;
                if(_cur_n_hits[na+nb]>=param.max_num_hits) return npair;
            }
		}
	}

/*	
	cout<<"\tpair:"<<_cur_n_hits[na+nb]<<" na:"<<na<<" nb:"<<nb<<endl;

	cout<<"ma0:"<< _sa._cur_n_hit[na]<<endl;
	for(i=0;i<_sa._cur_n_hit[na]; i++) 
		cout << "\t"<< (int)_sa.hits[na][i].chr<< ' '<<_sa.hits[na][i].loc<<endl;
    cout<<"ma1:"<< _sa._cur_n_chit[na]<<endl;
    for(i=0;i<_sa._cur_n_chit[na]; i++) 
        cout << "\t"<< (int)_sa.chits[na][i].chr<< ' '<<_sa.chits[na][i].loc<<endl; 

    cout<<"mb0:"<< _sb._cur_n_hit[nb]<<endl;
    for(i=0;i<_sb._cur_n_hit[nb]; i++)
        cout << "\t"<< (int)_sb.hits[nb][i].chr<< ' '<<_sb.hits[nb][i].loc<<endl;
    cout<<"mb1:"<< _sb._cur_n_chit[nb]<<endl;
    for(i=0;i<_sb._cur_n_chit[nb]; i++)
        cout << "\t"<< (int)_sb.chits[nb][i].chr<< ' '<<_sb.chits[nb][i].loc<<endl;
	
*/
	return npair;	
}    

int PairAlign::RunAlign(RefSeq &ref)
{
	bit32_t n,i,j,seedseg_num;

    for(i=0; i<=MAXSNPS*2; i++) _cur_n_hits[i]=0;
    for(i=0; i<=MAXSNPS; i++) for(j=0; j<=MAXSNPS; j++) checked_pair_mismatch[i][j]=0;

	_sa.ClearHits();
	_sb.ClearHits();

	_sa.seedseg_num=min((bit32_t)((_sa._pread->seq.size()-param.index_interval+1)/param.seed_size),_sa.read_max_snp_num+1);
	_sb.seedseg_num=min((bit32_t)((_sb._pread->seq.size()-param.index_interval+1)/param.seed_size),_sb.read_max_snp_num+1);

    seedseg_num=max(_sa.seedseg_num,_sb.seedseg_num);
	_sa.ConvertBinaySeq();
	_sb.ConvertBinaySeq();
	_sa.snp_thres=_sa.read_max_snp_num;
    _sb.snp_thres=_sb.read_max_snp_num;
    _sa.cseed_offset=_sa._pread->seq.size()%param.seed_size;
    _sb.cseed_offset=_sb._pread->seq.size()%param.seed_size;
        //cout<<"_sa "; for (int ii=0; ii<= param.max_snp_num;ii++) cout<<' '<<_sa._cur_n_hit[ii]<<' '<<_sa._cur_n_chit[ii]<<' '; cout<<endl;
        //cout<<"_sb "; for (int ii=0; ii<= param.max_snp_num;ii++) cout<<' '<<_sb._cur_n_hit[ii]<<' '<<_sb._cur_n_chit[ii]<<' '; cout<<endl;
        //cout<<endl;
    _sa.ReorderSeed(ref); _sb.ReorderSeed(ref);
    bit32_t maxi=max(_sa.read_max_snp_num,_sb.read_max_snp_num);                
	//for(i=0; i<=_sa.read_max_snp_num+_sb.read_max_snp_num; i++){
	n=0;
	for(i=0; i<=maxi; i++){
	    //cout << "i="<<i<<endl;
        if(i<_sa.seedseg_num) _sa.SnpAlign(ref,i); 
        if(i<_sb.seedseg_num) _sb.SnpAlign(ref,i);
        if(i<=_sa.read_max_snp_num) _sa.SortHits4PE(i);
        if(i<=_sb.read_max_snp_num) _sb.SortHits4PE(i);       
        n+=GetPairs(i,i); 
        for(j=0; j<i; j++) n+=GetPairs(i,j)+GetPairs(j,i); 
        if(param.nt3) continue;
        if(n>0) return 1;
	}

/*
    PairHit pp;
    cout<<_sa._pread->name.c_str()<<endl;
    for(i=0; i<=param.max_snp_num*2; i++) if(_cur_n_hits[i]) {
        cout<<"mis:"<<i<<"  pairs:"<<_cur_n_hits[i]<<endl;
        for(j=0;j<_cur_n_hits[i];j++) {
            pp=pairhits[i][j];
            cout<<ref.title[pp.a.chr].name<<":"<<pp.a.loc+1<<":"<<chain_flag[pp.a.chr%2]<<chain_flag[pp.chain]<<endl;
            cout<<ref.title[pp.b.chr].name<<":"<<pp.b.loc+1<<":"<<chain_flag[pp.b.chr%2]<<chain_flag[!pp.chain]<<endl;
        }
    }
    return 1;
*/

    //cout <<"nothing\n";
	return n;		
}

void PairAlign::Do_Batch(RefSeq &ref)
{
	_str_align.clear();
	bit32_t tt=0,pair_reported=0; int filter1, filter2, paired;
    for(_sa._pread=_sa.mreads.begin(), _sb._pread=_sb.mreads.begin(); tt<num_reads; _sa._pread++, _sb._pread++, tt++) {
        //if(_sa._pread->seq.size()==_sb._pread->seq.size()) TrimAdapter();
        //else 
        //    {_sa.TrimAdapter(); _sb.TrimAdapter();}
        //cout << _sa._pread->name.c_str()<<"\t"<<_sb._pread->name.c_str()<<endl;
        filter1=_sa.FilterReads(); filter2=_sb.FilterReads();
        FixPairReadName();
        //cout << _sa._pread->name.c_str()<<" "<<filter1<<"\t"<<_sb._pread->name.c_str()<<" "<<filter2<<endl;
        if(filter1==0 && filter2==0) paired=RunAlign(ref);
        else {
            paired=0; 
            if(filter1==0) _sa.RunAlign(ref);
            if(filter2==0) _sb.RunAlign(ref);
        }
        //cout<<"pair:"<<paired<<endl;
        if(paired) pair_reported=StringAlignPair(ref, _str_align);
        if(pair_reported==0||paired==0) StringAlignUnpair(filter1, filter2, ref, _str_align);
	}
//	cout<<_str_align<<endl;
}

int PairAlign::StringAlignPair(RefSeq &ref, string &os) {//return 1: pair reported, 0: pair not reported
	bit32_t i, j, sum;
	for(i=0;i<=MAXSNPS*2; i++) {
	    if((sum=_cur_n_hits[i])>0) break;
    }
    if(sum==0) return 0;
    if(sum==1) { //unique pair
        ++n_unique_pairs; ++n_aligned_pairs;
        s_OutHitPair(pairhits[i][0], 1, ref, os);\
        return 1;
    }
    else {
        ++n_multiple_pairs;
		if(1==param.report_repeat_hits) {   //randomly pick up one
		    ++n_aligned_pairs;
		    j=myrand(_sa._pread->index, &rand_rSeed)%_cur_n_hits[i];
		    s_OutHitPair(pairhits[i][j], _cur_n_hits[i], ref, os);
		    return 1;
        }
        else if(2==param.report_repeat_hits) {
            ++n_aligned_pairs;
            for(j=0;j<_cur_n_hits[i];++j) s_OutHitPair(pairhits[i][j], _cur_n_hits[i], ref, os); 
            return 1;
		}
		else return 0; //multiple pair supressed by -r 0 option
	}
}

void PairAlign::StringAlignUnpair(int fa, int fb, RefSeq &ref, string &os) {
	int ma, mb, ma1, mb1; bit32_t na, nb, ra, rb, ca, cb, i;
	gHit ha, hb;

    na=nb=ra=rb=ma=mb=ca=cb=0; ha=_sa.hits[0][0]; hb=_sb.hits[0][0];
    if(fa) ma=-1;
    else {
        for(na=0; na<=_sa.read_max_snp_num; na++) if((ma=_sa._cur_n_hit[na]+_sa._cur_n_chit[na])>0) break;
        if(ma>0) {
            ra=myrand(_sa._pread->index, &rand_rSeed)%ma;
            if((ca=(ra>=_sa._cur_n_hit[na]))==0) ha=_sa.hits[na][ra];
            else ha=_sa.chits[na][ra-_sa._cur_n_hit[na]];
        }
        na%=(_sa.read_max_snp_num+1);     
    }

    if(fb) mb=-1;
    else {
        for(nb=0; nb<=_sb.read_max_snp_num; nb++) if((mb=_sb._cur_n_hit[nb]+_sb._cur_n_chit[nb])>0) break;
        if(mb>0) {
            rb=myrand(_sb._pread->index, &rand_rSeed)%mb;
            if((cb=(rb>=_sb._cur_n_hit[nb]))==0) hb=_sb.hits[nb][rb];
            else hb=_sb.chits[nb][rb-_sb._cur_n_hit[nb]];
        }
        nb%=(_sb.read_max_snp_num+1);
    }

    if(ma>1&&param.report_repeat_hits==0) ma1=0; else ma1=ma;
    if(mb>1&&param.report_repeat_hits==0) mb1=0; else mb1=mb;
    if(ma<=0) { // no hit
        if(param.out_unmap) s_OutHitUnpair(0, 0, cb, ma, 0, ha, mb1, hb, ref, os);
    }
    else if(ma==1) { // unique hit
        n_aligned_a++; n_unique_a++;
        s_OutHitUnpair(0, ca, cb, 1, na, ha, mb1, hb, ref, os);
    }
    else { // multiple hit
        n_multiple_a++;
        if(param.report_repeat_hits==1) {
            n_aligned_a++;
            s_OutHitUnpair(0, ca, cb, ma, na, ha, mb1, hb, ref, os);
        }
        else if(param.report_repeat_hits==2) {
            n_aligned_a++;
            for(i=0;i<_sa._cur_n_hit[na];++i) s_OutHitUnpair(0, 0, cb, ma, na, _sa.hits[na][i], mb1, hb, ref, os);
            for(i=0;i<_sa._cur_n_chit[na];++i) s_OutHitUnpair(0, 1, cb, ma, na, _sa.chits[na][i], mb1, hb, ref, os);
        }
        else if(param.out_unmap) s_OutHitUnpair(0, 0, cb, 0, 0, ha, mb1, hb, ref, os);                         
    }

    if(mb<=0) {
        if(param.out_unmap) s_OutHitUnpair(1, 0, ca, mb, 0, hb, ma1, ha, ref, os);
    }
    else if(mb==1) {
        n_aligned_b++; n_unique_b++;
        s_OutHitUnpair(1, cb, ca, 1, nb, hb, ma1, ha, ref, os); 
    }
    else {
        n_multiple_b++;
        if(param.report_repeat_hits==1) {
            n_aligned_b++;
            s_OutHitUnpair(1, cb, ca, mb, nb, hb, ma1, ha, ref, os);
        }
        else if(param.report_repeat_hits==2) {
            n_aligned_b++;
            for(i=0;i<_sb._cur_n_hit[nb];++i) s_OutHitUnpair(1, 0, cb, mb, nb, _sb.hits[nb][i], ma1, ha, ref, os);
            for(i=0;i<_sb._cur_n_chit[nb];++i) s_OutHitUnpair(1, 1, cb, mb, nb, _sb.chits[nb][i], ma1, ha, ref, os);
        }
        else if (param.out_unmap) s_OutHitUnpair(1, 0, ca, 0, 0, hb, ma1, ha, ref, os);
    }

    //s_OutHitUnpair(0, (ra>=_sa._cur_n_hit[na]), (rb>=_sb._cur_n_hit[nb]), ma, na, ha, mb, hb, ref, os);
    //s_OutHitUnpair(1, (rb>=_sb._cur_n_hit[nb]), (ra>=_sa._cur_n_hit[na]), mb, nb, hb, ma, ha, ref, os);
}

void PairAlign::s_OutHitPair(PairHit pp, int n, RefSeq &ref, string &os) {
    bit32_t ii, flag, rev_seq_a, rev_seq_b; int pp_insert;
    pair<ref_loc_t,int> seg_info;
    ref_loc_t seg_start;

    //if fraglen < readlen then we need to remove adapter sequence
    rev_seq_a=pp.chain^(pp.a.chr%2); rev_seq_b=(!pp.chain)^(pp.b.chr%2);
    /*
    if(pp.insert<_sa.map_readlen){ 
        if(pp.chain^(pp.a.chr%2)) {
            pp.a.loc+=_sa.map_readlen-pp.insert;
            if(pp.a.gap_pos>_sa.map_readlen-pp.insert) pp.a.gap_pos-=_sa.map_readlen-pp.insert;
            else pp.a.gap_size=0;
        }
        _sa._pread->seq.erase(pp.insert);
        if(_sa._pread->qual.size()>pp.insert) _sa._pread->qual.erase(pp.insert);
        _sa.map_readlen=_sa._pread->seq.size();
    }

    if(pp.insert<_sb.map_readlen){
        if((!pp.chain)^(pp.b.chr%2)) {
            pp.b.loc+=_sb._pread->seq.size()-pp.insert;
            if(pp.b.gap_pos>_sb.map_readlen-pp.insert) pp.b.gap_pos-=_sb.map_readlen-pp.insert;
            else pp.b.gap_size=0;
        }
        _sb._pread->seq.erase(pp.insert);
        if(_sb._pread->qual.size()>pp.insert) _sb._pread->qual.erase(pp.insert);
        _sb.map_readlen=_sb._pread->seq.size();
    }
    */

    if(param.out_sam) {
        flag=0x3;  //read_a
        if(n>1) flag|=0x100; 

        if(rev_seq_a){
            flag|=0x10; //reverse read sequence
            seg_start=pp.b.loc+1; pp_insert=-pp.insert;
        }
        else {
            flag|=0x20; //reverse mate read sequence
            seg_start=pp.a.loc+1; pp_insert=pp.insert;
        }
        flag|=0x40*_sa._pread->readset;

        if(pp.a.gap_size==0) sprintf(cigar, "%uM", _sa.map_readlen);
        else if(pp.a.gap_size>0) sprintf(cigar, "%dM%dD%dM", (int)pp.a.gap_pos, (int)pp.a.gap_size, (int)_sa.map_readlen-(int)pp.a.gap_pos);
        else sprintf(cigar, "%dM%dI%dM", (int)pp.a.gap_pos, -(int)pp.a.gap_size, (int)_sa.map_readlen-(int)pp.a.gap_pos+(int)pp.a.gap_size);


        sprintf(_ch,"%s\t%d\t%s\t%u\t255\t%s\t=\t%u\t%d\t%s\t%s\tNM:i:%d",
            _sa._pread->name.c_str(),flag,ref.title[pp.a.chr].name.c_str(),pp.a.loc+1,cigar,pp.b.loc+1,pp_insert,_sa._outseq[rev_seq_a].c_str(),_sa._outqual[rev_seq_a].c_str(),pp.na);
            
        os.append(_ch);

        if(param.out_ref) {
            int ptr=0;
            for(ii=2;ii>0;ii--) {
                if(pp.a.loc<ii) continue;
                _mapseq[ptr++]=param.useful_nt[*(ref.bfa[(pp.a.chr>>1)<<1].s+(pp.a.loc-ii)/SEGLEN)>>(SEGLEN*2-2-((pp.a.loc-ii)%SEGLEN)*2)&0x3]+32;
            }
            for(ii=0;ii<_sa._pread->seq.size()+2;ii++) {
                _mapseq[ptr++]=param.useful_nt[*(ref.bfa[(pp.a.chr>>1)<<1].s+(pp.a.loc+ii)/SEGLEN)>>(SEGLEN*2-2-((pp.a.loc+ii)%SEGLEN)*2)&0x3];
            }
            _mapseq[ptr]=0;_mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;

            sprintf(_ch, "\tXR:Z:%s",_mapseq);
            os.append(_ch);
        }

        //seg_info=ref.CCGG_seglen(pp.a.chr, pp.a.loc); 
        if(param.RRBS_flag){       
            sprintf(_ch,"\tZP:i:%d\tZL:i:%d",seg_start, pp.insert); 
            os.append(_ch);
        }

        sprintf(_ch,"\tZS:Z:%c%c\n", chain_flag[pp.a.chr%2], chain_flag[pp.chain]); 
        os.append(_ch);        

        //sprintf(_ch,"\tZP:i:%d\tZL:i:%d\n",seg_start,pp.insert); os.append(_ch);


        flag=0x3;  //read_b
        if(n>1) flag|=0x100;
        if(rev_seq_b) {
            flag|=0x10; //reverse read sequence
            seg_start=pp.a.loc+1; pp_insert=-pp.insert;
        }
        else {
            flag|=0x20; //reverse mate read sequence
            seg_start=pp.b.loc+1; pp_insert=pp.insert;
        }

        flag|=0x40*_sb._pread->readset;

        if(pp.b.gap_size==0) sprintf(cigar, "%uM", _sb.map_readlen);
        else if(pp.b.gap_size>0) sprintf(cigar, "%dM%dD%dM", (int)pp.b.gap_pos, (int)pp.b.gap_size, (int)_sb.map_readlen-(int)pp.b.gap_pos);
        else sprintf(cigar, "%dM%dI%dM", (int)pp.b.gap_pos, -(int)pp.b.gap_size, (int)_sb.map_readlen-(int)pp.b.gap_pos+(int)pp.b.gap_size);

        sprintf(_ch,"%s\t%d\t%s\t%u\t255\t%s\t=\t%u\t%d\t%s\t%s\tNM:i:%d",
            _sb._pread->name.c_str(),flag,ref.title[pp.b.chr].name.c_str(),pp.b.loc+1,cigar,pp.a.loc+1,pp_insert,_sb._outseq[rev_seq_b].c_str(),_sb._outqual[rev_seq_b].c_str(),pp.nb);

        os.append(_ch);        

        if(param.out_ref) {
            int ptr=0;
            for(ii=2;ii>0;ii--) {
                if(pp.b.loc<ii) continue;
                _mapseq[ptr++]=param.useful_nt[*(ref.bfa[(pp.b.chr>>1)<<1].s+(pp.b.loc-ii)/SEGLEN)>>(SEGLEN*2-2-((pp.b.loc-ii)%SEGLEN)*2)&0x3]+32;
            }
            for(ii=0;ii<_sb._pread->seq.size()+2;ii++) {
                _mapseq[ptr++]=param.useful_nt[*(ref.bfa[(pp.b.chr>>1)<<1].s+(pp.b.loc+ii)/SEGLEN)>>(SEGLEN*2-2-((pp.b.loc+ii)%SEGLEN)*2)&0x3];
            }
            _mapseq[ptr]=0;_mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;
            sprintf(_ch, "\tXR:Z:%s",_mapseq);
            os.append(_ch);
        }

        //seg_info=ref.CCGG_seglen(pp.a.chr, pp.a.loc); 
        if(param.RRBS_flag){       
            sprintf(_ch,"\tZP:i:%d\tZL:i:%d",seg_start, pp.insert); 
            os.append(_ch);
        }

        sprintf(_ch,"\tZS:Z:%c%c\n", chain_flag[pp.b.chr%2], chain_flag[!pp.chain]); 
        os.append(_ch);        

    }
    else {
    	_sa.s_OutHit(pp.chain, n, pp.na, &pp.a, pp.insert, ref, os);
    	_sb.s_OutHit(!pp.chain, n, pp.nb, &pp.b, pp.insert, ref, os);			
	}
}

void PairAlign::s_OutHitUnpair(int readinpair, int chain_a, int chain_b, int ma, bit32_t na, gHit ha, int mb, gHit hb, RefSeq &ref, string &os) {
    bit32_t ii, flag=1, rev_seq=chain_a^(ha.chr%2);
    pair<ref_loc_t,int> seg_info;
    if(readinpair==0) _stmp=&_sa; else _stmp=&_sb;
    if(param.out_sam){
        flag|=0x40*_stmp->_pread->readset;
        if(ma<=0){ //QC, NM
            if(ma<0) flag|=0x204; 
            if(ma==0) flag|=0x004; 
            if(ma>1) flag|=0x104;
            if(mb<=0) {
                flag|=0x008;
           	    sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_stmp->_pread->name.c_str(),flag,_stmp->_pread->seq.c_str(),_stmp->_pread->qual.c_str()); 
           	}
           	else {
           	    if(chain_b^(hb.chr%2)) flag|=0x020;
           	    sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t%s\t%u\t0\t%s\t%s\n",_stmp->_pread->name.c_str(),flag,ref.title[hb.chr].name.c_str(),hb.loc+1,_stmp->_pread->seq.c_str(),_stmp->_pread->qual.c_str()); 
           	}
            os.append(_ch);                
        }           	    
        else{
            if(ma>1) flag|=0x100; //MA, OF
            if(rev_seq) flag|=0x010; //reverse read sequence

            if(ha.gap_size==0) sprintf(cigar, "%uM", _stmp->map_readlen);
            else if(ha.gap_size>0) sprintf(cigar, "%dM%dD%dM", (int)ha.gap_pos, (int)ha.gap_size, (int)_stmp->map_readlen-(int)ha.gap_pos);
            else sprintf(cigar, "%dM%dI%dM", (int)ha.gap_pos, -(int)ha.gap_size, (int)_stmp->map_readlen-(int)ha.gap_pos+(int)ha.gap_size);
            if(mb<=0) {
                flag|=0x008;
                sprintf(_ch,"%s\t%d\t%s\t%u\t255\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d",
                    _stmp->_pread->name.c_str(),flag,ref.title[ha.chr].name.c_str(),ha.loc+1,cigar,_stmp->_outseq[rev_seq].c_str(),_stmp->_outqual[rev_seq].c_str(),na);
            }
            else {
           	    if(chain_b^(hb.chr%2)) flag|=0x020;
                sprintf(_ch,"%s\t%d\t%s\t%u\t255\t%s\t%s\t%u\t0\t%s\t%s\tNM:i:%d",
                    _stmp->_pread->name.c_str(),flag,ref.title[ha.chr].name.c_str(),ha.loc+1,cigar,ref.title[hb.chr].name.c_str(),hb.loc+1,_stmp->_outseq[rev_seq].c_str(),_stmp->_outqual[rev_seq].c_str(),na);
           	}

            os.append(_ch);                
           	
            if(param.out_ref) {
                int ptr=0;
                for(ii=2;ii>0;ii--) {
                    if(ha.loc<ii) continue;
                    _mapseq[ptr++]=param.useful_nt[*(ref.bfa[(ha.chr>>1)<<1].s+(ha.loc-ii)/SEGLEN)>>(SEGLEN*2-2-((ha.loc-ii)%SEGLEN)*2)&0x3]+32;
                }
                for(ii=0;ii<_stmp->_pread->seq.size()+2;ii++) {
                    _mapseq[ptr++]=param.useful_nt[*(ref.bfa[(ha.chr>>1)<<1].s+(ha.loc+ii)/SEGLEN)>>(SEGLEN*2-2-((ha.loc+ii)%SEGLEN)*2)&0x3];
                }
                _mapseq[ptr]=0;_mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;

                sprintf(_ch, "\tXR:Z:%s",_mapseq);
                os.append(_ch);
            }

            if(param.RRBS_flag) {
                seg_info=ref.CCGG_seglen(ha.chr, ha.loc, _stmp->_pread->seq.size());
                sprintf(_ch,"\tZP:i:%d\tZL:i:%d",seg_info.first,seg_info.second); 
                os.append(_ch);
            }

            sprintf(_ch,"\tZS:Z:%c%c\n", chain_flag[ha.chr%2], chain_flag[chain_a]); 
            os.append(_ch);

        }
    }
    else _stmp->s_OutHit(chain_a, ma, na, &ha, 0, ref, os);
}


/*    
int PairAlign::TrimAdapter(){
    int i,j,k,l,h, ma,mb,ma0, mb0,pos;
    char tmpnt;
    for(pos=_sa._pread->seq.size()-4;pos>=param.seed_size;pos--){
        for(i=0; i<param.n_adapter; i++){
            ma0=mb0=0;
            for(k=0; k<param.adapter[i].size(); k++) {
                if((l=pos+k)>=_sa._pread->seq.size()) break;
                if((tmpnt=param.adapter[i][k])=='N') continue;
                if(_sa._pread->seq[l]!=tmpnt&&_sa._pread->seq[l]!='N') ma0++;
                if(_sb._pread->seq[l]!=tmpnt&&_sb._pread->seq[l]!='N') mb0++;
            }
            if((2*k+6)<(ma0+mb0)*5) continue;
            for(j=0;j<4;j++) {
                ma=ma0; mb=mb0;
                for(h=0,l=pos-3; h<3; h++,l++) {
                    if(_sa._pread->seq[l]!=param.adapter_anchor[j][h]&&_sa._pread->seq[l]!='N') ma++;
                    if(_sb._pread->seq[l]!=param.adapter_anchor[j][h]&&_sb._pread->seq[l]!='N') mb++;
                }
                if((2*k+6)>=(ma+mb)*5) {
                    _sa._pread->seq.erase(pos); 
                    _sa._pread->qual.erase(pos);
                    _sb._pread->seq.erase(pos); 
                    _sb._pread->qual.erase(pos);
                    return 1;
                }
            }
        }
    }
    return 0;
}
*/

void PairAlign::FixPairReadName() { //make sure the name of the two reads are same
    //cout<<_sa._pread->name<<endl;
    //cout<<_sb._pread->name<<endl;
    if(!param.out_sam) return;
    if(_sa._pread->name==_sb._pread->name) return;
    int i,d=-1, i0=min(_sa._pread->name.size(),_sb._pread->name.size());
    for(i=0;i<i0;i++) {
        if(_sa._pread->name[i]!=_sb._pread->name[i]) break;
        else if(isdigit(_sa._pread->name[i])) d=i;
    }
    if(i>0) {
        if(d<0) d=i-1;
        _sa._pread->name.erase(d+1);
        _sb._pread->name.erase(d+1);
    }
    else{
        cerr<<"Error: Paired reads name not match:\n" << _sa._pread->name.c_str()<<endl<< _sb._pread->name.c_str()<<endl;
        exit(1);
    }
    //cout<<_sa._pread->name<<endl;
}
