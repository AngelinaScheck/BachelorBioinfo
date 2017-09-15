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
    //std::string svtype;                       // type
};

void clearHit(hit & thisHit);
// 
void feedFromVCFrank(seqan::VcfRecord record, struct hit & thisHit);

void insertSort(std::vector<hit> & presorted, struct hit newHit);

// void quickSort(std::vector<hit> & unknown, unsigned left, unsigned right);

bool recepMatch(struct hit trueChecked, struct hit unknown);
#endif