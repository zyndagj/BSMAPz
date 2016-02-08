#ifndef _READS_H_
#define _READS_H_

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include "param.h"
#include "sam.h"
#include "gzstream.h"

using namespace std;

const bit32_t BatchNum=50000;

struct ReadInf
{
	bit32_t index;
	bit32_t readset; //added by yxi
	string name;
	string seq;
	string qual;
};

class ReadClass
{
public:
	ReadClass();
	void InitIndex(ifstream &fin, igzstream &gzfin, string filename);
	int LoadBatchReads(ifstream &fin, igzstream &gzfin, int mode);
public:
	vector<ReadInf> mreads;
	bit32_t num;
	bit32_t _index;
    
    //added by yxi, for BAM input support
    samfile_t *SAM_fp;
    bam1_t *SAM_b;
};

#endif //_READS_H_
