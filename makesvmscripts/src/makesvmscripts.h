#ifndef MAKESVMSCRIPTS_H
#define MAKESVMSCRIPTS_H

struct infoStore{
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

struct infoStore infoExtract(std::string originalInfo);
seqan::CharString infoToCString(struct infoStore infoFromStruct);
#endif