#ifndef RESULTCOMPARE_H
#define RESULTCOMPARE_H

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace std;

struct hit
{
// chromosome and position of true positiv structural variation
    int chr;
    int start;
    int end;
    double score;
};

void clearHit(hit & thisHit);
// 
int feedFromVCF(seqan::VcfRecord record, hit & thisHit);

void quickSort(std::vector<hit> & unknown, unsigned left, unsigned right);
#endif