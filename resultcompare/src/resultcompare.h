#ifndef RESULTCOMPARE_H
#define RESULTCOMPARE_H

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/bed_io.h>

using namespace std;

struct hit
{
// chromosome and position of true positiv structural variation
    int chr;                //x-chromosome is 23, y chromosome is 24
    long start;             //int was too short sometimes during testing
    long end;
    double score;
    std::string svtype;     // type
    long size;
    int vt;
    double se;
    double pe;
    double ce;
    double re;
};

struct info{
    double sc;
    int vt;
    double se;
    double pe;
    double ce;
    double re;
    double rd;
    double gc;
    double cp;
    long svlen;
    std::string svtype;
};


std::vector<struct hit> feedRank(std::string filename);

struct info infoExtract(std::string originalInfo);

std::vector<struct hit> readReference (std::string filename);

void clearHit(hit & thisHit);
// 
void feedFromVCFrank(seqan::VcfRecord record, struct hit & thisHit);

void insertSort(std::vector<hit> & presorted, struct hit newHit);

void quickSort(std::vector<hit> & unknown, unsigned left, unsigned right);

bool recepMatch(struct hit trueChecked, struct hit unknown);
#endif