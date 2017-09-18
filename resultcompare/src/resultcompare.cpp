#include <seqan/bed_io.h>
#include <seqan/vcf_io.h>
#include <seqan/find.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <algorithm> 

#include "resultcompare.h"

#define CharStringToStdString(x) std::string(toCString(x))

using namespace seqan;

int main()
{
    //__________________Read_and_process_list_of_true_positives_in_bed_file_format__________________________________
    std::cout<< "processing of list of true positives" <<std::endl;
    
    std::vector<struct hit> listoftrues=readReference ("GIAB_NA12878.bed");
//     for(unsigned i=0; i<listoftrues.size();i++){
//         std::cout<<"hit "<<i<<": chromosome "<<listoftrues[i].chr<<", start "<<listoftrues[i].start<<", end "<<listoftrues[i].end<<std::endl;
//     }
    
    
    
//     //__________________Read_and_process_positives_from_rankaggregation_________________________________________
//     //--prepare structures----
//     
     std::cout<< "processing and file reading for rank aggregation" <<std::endl;
     std::vector<struct hit> listfromrank=feedRank("47_out.vcf");
//     
//     
//     //----loop-over-records-in-original-vcf-file,-fill-structures-------------------
//     std::vector<struct hit> listfromrank;
//     listoftrues.reserve(1000000);
//     
//     struct hit forRank;
//     
//     VcfFileIn vcfInRank(toCString(getAbsolutePath("../../resultcompare-build/47_out.vcf")));
//     VcfRecord recordRank;
//     // read header.
//     VcfHeader headerRank;
//     readHeader(headerRank, vcfInRank);
//    // int recordNumber=0;
//     double sc;
//     int vt;
//     double se;
//     double pe;
//     double ce;
//     double re;
//     double rd;
//     double gc;
//     double cp;
//     long svlen;
//     std::string svtype;
// 
//     while (!atEnd(vcfInRank)){
//         clearHit(forRank);
//         readRecord(recordRank, vcfInRank);
//         //for testing only deletions
//         if(recordRank.alt=="<DEL>"){
//             infoExtract(CharStringToStdString(recordRank.info), sc, vt, se, pe, ce, re, rd, gc, cp, svlen, svtype);
//     
//             forRank.chr=to_string(recordRank.rID);
//             forRank.start= long(recordRank.beginPos);
//             forRank.end= forRank.start + abs(svlen);
//             forRank.score=sc;
//             
//             listfromrank.push_back(forRank);
//             
//             
//             
//             //sort while reading to save memory
// //             if(recordNumber>0){
// //                 //insertSort(listfromrank, forRank);
// //             }
// //             else{
// //                 listfromrank.push_back(forRank);
// //             }
// //             recordNumber++;
//         }
//      }
//      
//      listfromrank.shrink_to_fit();
    
     
//     for(unsigned i=0; i<listfromrank.size();i++){
//         std::cout<<"hit "<<i<<": chromosome "<<listfromrank[i].chr<<", start "<<listfromrank[i].start<<", end "<<listfromrank[i].end <<", score "<<listfromrank[i].score<<std::endl;
//     }
    
    //quickSort(listfromrank, 0, listfromrank.size()-1);
     
     //testing purpose only! is correctly sorted?
     for (unsigned i=0; i<listfromrank.size()-1;i++){
        if(listfromrank[i].score>listfromrank[i+1].score){
            std::cerr<<"sorting did not work. "<<listfromrank[i].score<<" not smaller "<<listfromrank[i+1].score<< std::endl;
        }
    }
    
    std::cout<<"out of function size: "<<listfromrank.size()<<std::endl;
//      
//      
//     
    //------sort into quantiles-------------------------------------------------------------------
    int numberRank=listfromrank.size();
    int marginQ1= int (0.25*numberRank);
    int marginQ2= int (0.5*numberRank);
    int marginQ3= int (0.75*numberRank);
    
    std::vector<struct hit> rankQ1;
    rankQ1.reserve(marginQ1 + 1);
    std::vector<struct hit> rankQ2;
    rankQ2.reserve(marginQ1 + 1);
    std::vector<struct hit> rankQ3;
    rankQ3.reserve(marginQ1 + 1);
    std::vector<struct hit> rankQ4;
    rankQ4.reserve(marginQ1 + 1);
//     
    //std::cout<< "quick sort records" <<std::endl;
    //struct hit temp;
     //quickSort(listfromrank, 0, numberRank-1);
    
    
    std::cout<< "make quantiles" <<std::endl;
    for (unsigned i=0; i<numberRank; i++){
        if(i<marginQ1){rankQ1.push_back(listfromrank[i]);}
        else if(i<marginQ2){rankQ2.push_back(listfromrank[i]);}
        else if(i<marginQ3){rankQ3.push_back(listfromrank[i]);}
        else {rankQ4.push_back(listfromrank[i]);}
    }
    
    rankQ1.shrink_to_fit();
    rankQ2.shrink_to_fit();
    rankQ3.shrink_to_fit();
    rankQ4.shrink_to_fit();
    
    std::cout<<rankQ1.size()<<" "<<rankQ2.size()<<" "<<rankQ3.size()<<" "<<rankQ4.size()<<std::endl;
    //_______________________Compare__________________________________________________________________
    std::cout<< "comparison" <<std::endl;
    
    int numberTrues=listoftrues.size();

    
    //-----comparison-for-rank-aggregation------------------------------------------------------
    int totalRank=listfromrank.size();
    int truePosRankQ1 = 0;
    int truePosRankQ2 = 0;
    int truePosRankQ3 = 0;
    int truePosRankQ4 = 0;
    int truePosRankTotal = 0;
    
    int falseNegRank = 0;
    int falsePosRank= 0;
    
    for(unsigned i=0; i<rankQ1.size();i++){
        for(unsigned j=0;j<numberTrues;j++){
            if(recepMatch(rankQ1[i], listoftrues[j])){
                truePosRankQ1++;
                j=numberTrues;
            }
        }
    }
    
    for(unsigned i=0; i<rankQ2.size();i++){
        for(unsigned j=0;j<numberTrues;j++){
            if(recepMatch(rankQ2[i], listoftrues[j])){
                truePosRankQ2++;
                j=numberTrues;
            }
        }
    }
    
    for(unsigned i=0; i<rankQ3.size();i++){
        for(unsigned j=0;j<numberTrues;j++){
            if(recepMatch(rankQ3[i], listoftrues[j])){
                truePosRankQ3++;
                j=numberTrues;
            }
        }
    }
    
    for(unsigned i=0; i<rankQ4.size();i++){
        for(unsigned j=0;j<numberTrues;j++){
            if(recepMatch(rankQ4[i], listoftrues[j])){
                truePosRankQ4++;
                j=numberTrues;
            }
        }
    }
    
//     for(int i=0;i<numberTrues;i++){
//         for(int j=0; j<rankQ1.size();j++){
//             if(recepMatch(listoftrues[i], rankQ1[j])){
//                 truePosRankQ1++;
//                 break;
//             }
//         }
//         
//         for(int k=0; k<rankQ2.size();k++){ //!!!improve here for speed
//             if(recepMatch(listoftrues[i], rankQ2[k])){
//                 truePosRankQ2++;
//                 break;
//             }
//         }
//         
//         for(int l=0; l<rankQ3.size();l++){ //!!!improve here for speed
//             if(recepMatch(listoftrues[i], rankQ3[l])){
//                 truePosRankQ3++;
//                 break;
//             }
//         }
//         
//         
//         for(int m=0; m<rankQ4.size();m++){ //!!!improve here for speed
//             if(recepMatch(listoftrues[i], rankQ4[m])){
//                 truePosRankQ4++;
//                 break;
//             }
//         }
//         
//     }
    
    truePosRankTotal=truePosRankQ1+truePosRankQ2+truePosRankQ3+truePosRankQ4;
    
    falseNegRank= numberTrues-truePosRankTotal;
    falsePosRank=totalRank-truePosRankTotal;
    
    //-------------------------print-results-------------------------------------------------
    std::cout<<numberTrues<<" true positives"<<std::endl;
    
    std::cout<<"rank aggregation identified "<<truePosRankTotal<<" true positives, " <<std::endl;
    std::cout<<truePosRankQ1<<" in Q1 "<<truePosRankQ2<<" in Q2 "<<truePosRankQ3<<" in Q3 "<<truePosRankQ4<<" in Q4 "<<std::endl;
    
    std::cout<<"in total "<<falseNegRank<<" false negatives and "<<falsePosRank<<" false positives "<<std::endl;

    return 0;
}


//=======================================================================================================================

std::vector<struct hit> feedRank(std::string filename){
    std::vector<struct hit> listfromrank;
    struct hit forRank;
    
    listfromrank.reserve(1000000);
    
    std::string chrUnprocessed;
    int start;
    std::string kauderwelsch1;
    std::string kauderwelsch2;
    std::string svType;
    std::string kauderwelsch3;
    std::string kauderwelsch4;
    std::string infoFromStream;
    
    std::ifstream rankin(filename);
    std::string line;
    
    while (std::getline(rankin, line)){
        
        std::stringstream rankstream(line);
        if (rankstream >> chrUnprocessed >> start >> kauderwelsch1 >> kauderwelsch2 >> svType >> kauderwelsch3 >> kauderwelsch4 >> infoFromStream){
            
            clearHit(forRank);
            struct info infoForRank;
            infoForRank= infoExtract(infoFromStream);
            
            chrUnprocessed=chrUnprocessed.erase(0, 3);
            if(chrUnprocessed.length()>1 && isdigit(chrUnprocessed[1]) ){
                chrUnprocessed=chrUnprocessed.substr(0,2);
                forRank.chr=std::stoi(chrUnprocessed);
            }
            else if (chrUnprocessed=="X" || chrUnprocessed=="Y"){
                chrUnprocessed=chrUnprocessed.substr(0,1);
                 if (chrUnprocessed=="X"){
                    forRank.chr=23;
                }
                else {
                    forRank.chr=24;
                }
            }
            
            forRank.start= start;
            forRank.size=abs(infoForRank.svlen);
            forRank.end= start + forRank.size;
            forRank.score=infoForRank.sc;
            forRank.svtype=infoForRank.svtype;
            forRank.vt=infoForRank.vt;
            forRank.se=infoForRank.se;
            forRank.pe=infoForRank.pe;
            
            if(listfromrank.size()>=1 && (forRank.chr>=10 && forRank.chr<=24) && forRank.size>50 && ((forRank.vt==3)) && ((forRank.se + forRank.pe)>=4) && forRank.svtype=="DEL"){
                insertSort(listfromrank, forRank);
            }
            else if (forRank.chr>=10 && forRank.chr<=24 && forRank.size>50 && ((forRank.vt==3)) && ((forRank.se + forRank.pe)>=4) && forRank.svtype=="DEL" ) {
                listfromrank.push_back(forRank);
            }
        }
    }
    listfromrank.shrink_to_fit();
//     for(unsigned i=0; i<3;i++){
//         std::cout<<"hit "<<i<<": chromosome "<<listfromrank[i].chr<<", start "<<listfromrank[i].start<<", end "<<listfromrank[i].end<<std::endl;
//     }
    return listfromrank;
}

struct info infoExtract(std::string originalInfo){
    //Example:    SC=-1.000000;VT=2;SE=0;PE=1;CE=0;RE=39.848095;RD=110.195238;GC=0.414500;CP=1.790252;SVLEN=-576;SVTYPE=DEL
//     sc=0;
//     vt=0;
//     se=0;
//     pe=0;
//     ce=0;
//     re=0;
//     rd=0;
//     gc=0;
//     cp=0;
//     svlen=0;
//     svtype="";
    
    struct info fillInHere;
    
    int endOfInfo=originalInfo.size();
    for (unsigned i=0; i<endOfInfo; i++){
        if(originalInfo[i]=='S'){
            if(originalInfo[i+1]=='C'){
                int k=i+3;
                std::string subSC;
                while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subSC += originalInfo[k];
                    k++;
                }
                fillInHere.sc=std::stod(subSC);
            }
            else if(originalInfo[i+1]=='E'){
                int k=i+3;
                std::string subSE;
                while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subSE += originalInfo[k];
                    k++;
                }
                fillInHere.se=std::stod(subSE);
            }
            else if(originalInfo[i+1]=='V'){
                if(originalInfo.substr(i+2, 4)=="LEN="){
                    int k=i+6;
                    std::string subSVLEN;
                    while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subSVLEN += originalInfo[k];
                    k++;
                    }
                    fillInHere.svlen=std::stol(subSVLEN);
                }
                else if(originalInfo.substr(i+2, 4)=="TYPE"){
                    int k=i+7;
                    std::string subSVTYPE;
                    while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subSVTYPE += originalInfo[k];
                    k++;
                    }
                    fillInHere.svtype=subSVTYPE;
                }
            }
        }
        else if(originalInfo.substr(i-1, 3)==";VT"){
            int k=i+3;
            std::string subVT;
            while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subVT += originalInfo[k];
                    k++;
            }
            fillInHere.vt=std::stoi(subVT);
        }
        else if(originalInfo.substr(i-1, 3)==";PE"){
            int k=i+3;
            std::string subPE;
            while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subPE += originalInfo[k];
                    k++;
            }
            fillInHere.pe=std::stod(subPE);
        }
        else if(originalInfo[i]=='R'){
            if(originalInfo[i+1]=='E'){
                int k=i+3;
                std::string subRE;
                while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subRE += originalInfo[k];
                    k++;
                }
                fillInHere.re=std::stod(subRE);
            }
            else if(originalInfo[i+1]=='D'){
                int k=i+3;
                std::string subRD;
                while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subRD += originalInfo[k];
                    k++;
                }
                fillInHere.rd=std::stod(subRD);
            }
        }
        else if(originalInfo.substr(i, 2)=="GC"){
            int k=i+3;
            std::string subGC;
            while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subGC += originalInfo[k];
                    k++;
            }
            fillInHere.gc=std::stod(subGC);
        }
        else if(originalInfo[i]=='C'){
            if(originalInfo[i+1]=='E'){
                int k=i+3;
                std::string subCE;
                while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subCE += originalInfo[k];
                    k++;
                }
                fillInHere.ce=std::stod(subCE);
            }
            else if(originalInfo[i+1]=='P'){
                int k=i+3;
                std::string subCP;
                while((originalInfo[k]!=';') && (k<endOfInfo)){
                    subCP += originalInfo[k];
                    k++;
                }
                fillInHere.cp=std::stod(subCP);
            }
        }
    }
    return fillInHere;
}



std::vector<struct hit> readReference (std::string filename){
    std::vector<struct hit> listoftrues;
    listoftrues.reserve(1000000);
    struct hit truepositive;
    std::ifstream bedin(filename);
    std::string bedchromosome;
    long bedstart, bedend;
    std::string line;
    
    while (std::getline(bedin, line)){
        std::stringstream bedstream(line);
        if (bedstream >> bedchromosome >> bedstart >> bedend){
            clearHit(truepositive);
            if(bedchromosome=="X"){truepositive.chr=23;
            }
            else if(bedchromosome=="Y"){truepositive.chr=24;
            }
            else{truepositive.chr=std::stoi(bedchromosome);
            }
            truepositive.start= bedstart;
            truepositive.end=bedend;
            if(truepositive.end-truepositive.start>=50){
                listoftrues.push_back(truepositive);
            }
        }
    }
    listoftrues.shrink_to_fit();
    return listoftrues;
    
//     for(unsigned i=0; i<listoftrues.size();i++){
//         std::cout<<"hit "<<i<<": chromosome "<<listoftrues[i].chr<<", start "<<listoftrues[i].start<<", end "<<listoftrues[i].end<<std::endl;
//     }
}


void clearHit(hit & thisHit){
    thisHit.chr=0;
    thisHit.start=0;
    thisHit.end=0;
    thisHit.score=0;
    thisHit.svtype="";
    thisHit.size=0;
}


void feedFromVCFrank(seqan::VcfRecord record, struct hit & thisHit){
    //loop over Info CharString, pattern matching for keywords like SE, copy relevant info into Hit structure
    double sc;
    double svlen;
    
   // ______Processing of the Recordinfo________-
    Finder<CharString> finder(record.info);
   // ____Prepare Patterns(keywords)____-
    String<CharString> needles;
    appendValue(needles, "SC=");
    appendValue(needles, ";SVLEN=");
    //Find keywords, value starts at position Matchposition+Patternlength+2 and ends before semicolon
    Pattern<String<CharString>, WuManber> pattern(needles);
    unsigned endOfValue;
    while (find(finder, pattern)){
        for(unsigned i=endPosition(finder); i<=length(record.info); i++){
            if(record.info[i]==';' || i==length(record.info)) {
                //std::cout<< "end of hit at position " << i<< std::endl;
                endOfValue=i;
                break;
            }
        }
        //save info 
         CharString value = infix(record.info, endPosition(finder), endOfValue);
         
        if(infix(finder)=="SC="){
            if(empty(value)){
                 return;
             }
             else {
                sc=std::stod(CharStringToStdString(value));
             }
        }
        if(infix(finder)==";SVLEN="){
            if(empty(value)){
                 return;
             }
             else {
            svlen=std::stod(CharStringToStdString(value));
             }
        }

        else {
            thisHit.chr=record.rID;
            thisHit.start=record.beginPos;
            thisHit.end=record.beginPos + abs(svlen);
            thisHit.score=sc;
        }
    }
}

//--------------quickSort--------------------------------------------------------------------
void quickSort(std::vector<hit> & unknown, unsigned left, unsigned right){
    //intialize counter for left and right sublists
    unsigned i = left; 
    unsigned j = right;
    //initialize pivot element (middle)
    unsigned mid= left + (right - left) / 2;
    double pivot = unknown[mid].score;
    struct hit temp;
    
    //divide and conq.
    while (i<=j) {
        while (unknown[i].score <= pivot){
            //all elements in left sublist smaller than the pivot element remain untouched
            if (i>=j){
                break;
            }
            i++;
        }
        while (unknown[j].score > pivot){
            //all elements in right sublist bigger/equal than the pivot remain untouched  
            if (j<=i){
                break;
            }
            j--;
        }
        //swap the entries if dif. to pivot
        if (i <= j) {
            //assign the temporary variables
            temp=unknown[i];
            
            //swap
            unknown[i]=unknown[j];
            unknown[j]=temp;
            
            //continue
            i++;
            j--;
            }
    
        //if sublist already sorted go one level deeper    
        if(left<j){
            quickSort(unknown, left, j);
        }
        if(i<right){
            quickSort(unknown,i,right);
        }
    }
}

void insertSort(std::vector<hit> & presorted, struct hit newHit){
    //if value smaller than current biggest value, move left, if equal or smaller insert and move all bigger elements to the right

    if(newHit.score>presorted[presorted.size()-1].score){
        presorted.push_back(newHit);
        return;
    }
    else if(newHit.score<=presorted[0].score){
        presorted.insert(presorted.begin(), newHit);
        return;
    }

    else{
        int i=presorted.size()-1;
        while(i>0 && newHit.score<presorted[i].score){
            i--;
        }
        //std::cout<<newHit.score<<" i:"<<presorted[i+1].score<<std::endl;
        presorted.insert(presorted.begin()+i+1, newHit);
    }
}
// 
bool recepMatch(struct hit trueChecked, struct hit unknown){
    int sizeTrueChecked = abs(trueChecked.end-trueChecked.start);
    int sizeUnknown = abs(unknown.end-unknown.start);
    int sizeOfMatch = std::min(unknown.end, trueChecked.end) - std::max(unknown.start, trueChecked.start);
    
    if (((0.8* sizeTrueChecked) <= sizeOfMatch)  && ((0.8* sizeUnknown) <= sizeOfMatch) && (unknown.chr==trueChecked.chr) ){
        return true;
    }
    else {
        return false;
    }
}
