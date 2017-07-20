#ifndef STATISTICSLIBSVM_H
#define STATISTICSLIBSVM_H

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include "svm.h"

using namespace seqan;

struct Info_sep
{
// example: SC=-1.000000;VT=1;SE=0;PE=1;CE=0;RE=0.019990;RD=53.071429;GC=0.392500;CP=1.966275;SVLEN=-548;SVTYPE=DEL
    //SC
    std::string sc;
    //VT
    std::string vt;
    //SE
    std::string se;
    //PE
    std::string pe;
    //CE
    std::string ce;
    //RE
    std::string re;
    //RD
    std::string rd;
    //GC
    std::string gc;
    //CP
    std::string cp;
    //SVLEN
    std::string svlen;
    //SVTYPE
    std::string svtype;
    //ID of record(not in info)
    std::string recordID;
    //position in record file
    int entryNumber;
};

void clearInfo(Info_sep & infoInStruct);
void feedInfo(int entryNumber, seqan::VcfRecord record, Info_sep & infoInStruct);

#endif