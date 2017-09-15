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
    
    std::vector<struct hit> listoftrues;
    listoftrues.reserve(1000000);
    struct hit truepositive;
    
    std::ifstream bedin("../GIAB_NA12878.bed");
    int bedchromosome;
    int bedstart, bedend;
    std::string line;
    
    while (std::getline(bedin, line)){
        std::stringstream bedstream(line);
        if (bedstream >> bedchromosome >> bedstart >> bedend){
            clearHit(truepositive);
            truepositive.chr=bedchromosome;
            truepositive.start= bedstart;
            truepositive.end=bedend;
            listoftrues.push_back(truepositive);
        }
    }
    listoftrues.shrink_to_fit();
    
    //__________________Read_and_process_positives_from_rankaggregation_________________________________________
    //--prepare structures----
    
    std::cout<< "processing and file reading for rank aggregation (deletions only)" <<std::endl;
    std::vector<struct hit> listfromrank;
    struct hit forRank;
    
    listfromrank.reserve(1000000);
    
    //----loop-over-records-in-original-vcf-file,-fill-structures-------------------
    VcfFileIn vcfInRank(toCString(getAbsolutePath("../../resultcompare-build/47_out.vcf")));
    VcfRecord recordRank;
    // read header.
    VcfHeader headerRank;
    readHeader(headerRank, vcfInRank);
    int vtRank;

    while (!atEnd(vcfInRank)){
        clearHit(forRank);
        readRecord(recordRank, vcfInRank);
        //for testing only deletions
        if(recordRank.alt=="<DEL>"){
            feedFromVCFrank(recordRank, forRank);
            listfromrank.push_back(forRank);
        }
     }
     listfromrank.shrink_to_fit();
    
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
    
    std::cout<< "quick sort records into quantiles" <<std::endl;
    quickSort(listfromrank, 0, numberRank-1);
    
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
    
    for(int i=0;i<numberTrues;i++){
        for(int j=0; j<rankQ1.size();j++){
            if(recepMatch(listoftrues[i], rankQ1[j])){
                truePosRankQ1++;
                break;
            }
        }
        
        for(int k=0; k<rankQ2.size();k++){ //!!!improve here for speed
            if(recepMatch(listoftrues[i], rankQ2[k])){
                truePosRankQ2++;
                break;
            }
        }
        
        for(int l=0; l<rankQ3.size();l++){ //!!!improve here for speed
            if(recepMatch(listoftrues[i], rankQ3[l])){
                truePosRankQ3++;
                break;
            }
        }
        
        
        for(int m=0; m<rankQ4.size();m++){ //!!!improve here for speed
            if(recepMatch(listoftrues[i], rankQ4[m])){
                truePosRankQ4++;
                break;
            }
        }
        
    }
    
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
void clearHit(hit & thisHit){
    thisHit.chr=0;
    thisHit.start=0;
    thisHit.end=0;
    thisHit.score=0;
}


void feedFromVCFrank(seqan::VcfRecord record, struct hit & thisHit){
    //loop over Info CharString, pattern matching for keywords like SE, copy relevant info into Hit structure
    double sc;
    double svlen;
    
    //______Processing of the Recordinfo________-
    Finder<CharString> finder(record.info);
    //____Prepare Patterns(keywords)____-
    String<CharString> needles;
    appendValue(needles, ";SC=");
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
         
        if(infix(finder)==";SC="){
            sc=std::stod(CharStringToStdString(value));
        }
        if(infix(finder)==";SVLEN="){
            svlen=std::stod(CharStringToStdString(value));
        }
        
        if(empty(sc)){
                 return;
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
    //temporary values
    struct hit temp;
    
    //divide and conq.
    while (i<=j) {
        while (unknown[i].score <= pivot){
            //all elements in left sublist smaller/equal than the pivot element remain untouched
            if (i>=j){
                break;
            }
            i++;
        }
        while (unknown[j].score > pivot){
            //all elements in right sublist bigger than the pivot remain untouched  
            if (j<=i){
                break;
            }
            j--;
        }
        //if an element in the left sublist bigger pivot or an element in the right sublist smaller/equal the pivot element was find, swap the entries
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
    
        //if sublist already sorted got one level deeper    
        if(left<j){
            quickSort(unknown, left, j);
        }
        if(i<right){
            quickSort(unknown,i,right);
        }
    }
}

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