#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

using namespace seqan;

struct Info_sep
{
// example: SC=-1.000000;VT=1;SE=0;PE=1;CE=0;RE=0.019990;RD=53.071429;GC=0.392500;CP=1.966275;SVLEN=-548;SVTYPE=DEL
    //SC
    CharString sc;
    //VT
    CharString vt;
    //SE
    CharString se;
    //PE
    CharString pe;
    //CE
    CharString ce;
    //RE
    CharString re;
    //RD
    CharString rd;
    //GC
    CharString gc;
    //CP
    CharString cp;
    //SVLEN
    CharString svlen;
    //SVTYPE
    CharString svtype;
};

void feedInfo(seqan::VcfRecord record, Info_sep & infoInStruct);
void clearInfo(Info_sep & infoInStruct);

#endif