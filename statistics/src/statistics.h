#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

using namespace std;

struct Info_sep
{
// example: SC=-1.000000;VT=1;SE=0;PE=1;CE=0;RE=0.019990;RD=53.071429;GC=0.392500;CP=1.966275;SVLEN=-548;SVTYPE=DEL
    //SC
    float sc;
    //VT
    unsigned vt;
    //SE
    unsigned se;
    //PE
    unsigned pe;
    //CE
    unsigned ce;
    //RE
    float re;
    //RD
    float rd;
    //GC
    float gc;
    //CP
    float cp;
    //SVLEN
    int svlen;
    //SVTYPE
    string svtype;
};

void feedInfo(seqan::VcfRecord record, Info_sep & infoInStruct);

#endif