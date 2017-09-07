#include <seqan/bed_io.h>
#include <seqan/vcf_io.h>
#include <seqan/find.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

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
    
    //seqan read record for bed files does not work on the given file
    // Open input bed file, store in struct
//     BedFileIn bedIn(toCString(getAbsolutePath("../../GIAB_NA12878.bed")));
//     BedRecord<Bed3> record;
// 
//     while (!atEnd(bedIn))
//     {
//          readRecord(record, bedIn);
//         truepositive.chr=record.ref;
//         truepositive.start=record.beginPos;
//         truepositive.end=record.endPos;
//         truepositive.score=1;
//         listoftrues.push_back(truepositive);
//         clearHit(truepositive);
//     }
    
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
    
    //__________________Read_and_process_positives_from_libsvm-prediction__________________________________
    //positives with vt=3 (training dataset) will be read later in parallel to rank aggregation
    //47_out_ListOfStrucVar.vcf is list of vt=2 that were identified as true by libsvm
    std::cout<< "processing and file reading for svm results" <<std::endl;
    int vt;
    std::vector<struct hit> listfromsvm;
    listfromsvm.reserve(1000000);
    struct hit positivesvm;
    VcfFileIn vcfIn(toCString(getAbsolutePath("../../resultcompare-build/47_out_ListOfStrucVar.vcf")));
    VcfRecord record;
    
    // Copy over header.
    VcfHeader header;
    readHeader(header, vcfIn);
    //copy over record
    while (!atEnd(vcfIn)){
        clearHit(positivesvm);
        readRecord(record, vcfIn);
        //for testing only deletions
        if(record.alt=="<DEL>"){
            vt=feedFromVCF(record, positivesvm);
            if(vt==2 && positivesvm.chr!=0 && positivesvm.start!=0){
                listfromsvm.push_back(positivesvm);
            }
        }
    }

    
    
    //__________________Read_and_process_positives_from_rankaggregation_________________________________________
    //--prepare structures----
    
    std::cout<< "processing and file reading for rank aggregation" <<std::endl;
    std::vector<struct hit> listfromrank;
    std::vector<struct hit> unknownForRank;
    std::vector<struct hit> negativeForRank;
    struct hit forRank;
    
    listfromrank.reserve(1000000);
    negativeForRank.reserve(1000000);
    unknownForRank.reserve(1000000);
    
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
            vtRank=feedFromVCF(recordRank, forRank);
            if(vtRank==3 && forRank.chr!=0 && forRank.start!=0){
                listfromsvm.push_back(forRank); //training set is not included in resultlist that was read before
                listfromrank.push_back(forRank);
            }
            if(vtRank==2 && forRank.chr!=0 && forRank.start!=0){
                unknownForRank.push_back(forRank);
            }
            if(vtRank==1 && forRank.chr!=0 && forRank.start!=0){
                negativeForRank.push_back(forRank);
            }
        }
     }
     listfromsvm.shrink_to_fit();
    
    //--------------rank aggregation: top 25% of unknown are considered to be positives----------------------    
    unknownForRank.shrink_to_fit();
    negativeForRank.shrink_to_fit();
    
    int numberOfUnknown=unknownForRank.size();

    std::cout<< "quick sort records with vt=2" <<std::endl;
    //sort unknown according to score
    quickSort(unknownForRank, 0, numberOfUnknown-1);
    std::cout<< "rank aggregation" <<std::endl;
    //top 25 % get written into positive list (round by type cast)
    int margin= int (0.25*numberOfUnknown);
    for (unsigned i=0; i<margin; i++){
        listfromrank.push_back(unknownForRank[i]); 
    }
    listfromrank.shrink_to_fit();
    
    //_______________________Compare__________________________________________________________________
    std::cout<< "comparison" <<std::endl;
    
    int numberTrues=listoftrues.size();
    
    int totalSVM=listfromsvm.size();
    int truePosSVM = 0;
    int falseNegSVM = 0;
    int falsePosSVM= 0;
    
    
    for(unsigned i=0;i<numberTrues;i++){
        for(unsigned j=0; j<totalSVM;j++){
            if(listoftrues[i].chr==listfromsvm[j].chr && listoftrues[i].start<=listfromsvm[j].start && listoftrues[i].end>listfromsvm[j].start){
                truePosSVM++;
                break;
            }
        }
    }
    
    //false negatives are the ones contained in listoftrues but not in listfromsvm
    falseNegSVM= numberTrues-truePosSVM;
    //false positives are contained in listfromsvm but not in listoftrues
    falsePosSVM=totalSVM-truePosSVM;
    
    //-----same-for-rank-aggregation------------------------------------------------------
    int totalRank=listfromrank.size();
    int truePosRank = 0;
    int falseNegRank = 0;
    int falsePosRank= 0;
    
    for(int i=0;i<numberTrues;i++){
        for(int j=0; j<totalRank;j++){
            if(listoftrues[i].chr==listfromrank[j].chr && listoftrues[i].start<=listfromrank[j].start && listoftrues[i].end>listfromrank[j].start){
                truePosRank++;
                break;
            }
        }
    }
    
    falseNegRank= numberTrues-truePosRank;
    falsePosRank=totalRank-truePosRank;
    
    //-------------------------print-results-------------------------------------------------
    std::cout<<numberTrues<<" true positives"<<std::endl;
    std::cout<<"svm identified "<<totalSVM<<" hits"<<std::endl;
    std::cout<<"of which "<<truePosSVM<<" true positives, "<<falseNegSVM<<" false negatives and "<<falsePosSVM<<" false positives "<<std::endl;
    std::cout<<"rank aggregation identified "<<totalRank<<" hits"<<std::endl;
    std::cout<<"of which "<<truePosRank<<" true positives, " <<falseNegRank<<" false negatives and "<<falsePosRank<<" false positives "<<std::endl;
   

//     std::ofstream result("47_checkedIfTrue.txt");
//     result<<numberTrues<<"TotalTrues "<<"totalHits "<<"truePositive "<<"falseNegative "<<"falsePositive"<<"\n";
//     result<<"svm: "<<totalSVM<<" "<<truePosSVM<<" "<<falseNegSVM<<" "<<falsePosSVM<<"\n";
//     result<<"rankAgg: "<<totalRank<<" "<<truePosRank<<" "<<falseNegRank<<" "<<falsePosRank<<"\n";
//     result.close();

    return 0;
}


//=======================================================================================================================
void clearHit(hit & thisHit){
    thisHit.chr=0;
    thisHit.start=0;
    thisHit.end=0;
    thisHit.score=0;
}


int feedFromVCF(seqan::VcfRecord record, struct hit & thisHit){
    //loop over Info CharString, pattern matching for keywords like SE, copy relevant info into Hit structure
    int vt=999;
    double se;
    double pe;
    double re;
    
    //______Processing of the Recordinfo________-
    Finder<CharString> finder(record.info);
    //____Prepare Patterns(keywords)____-
    String<CharString> needles;
    appendValue(needles, ";VT=");
    appendValue(needles, ";SE=");
    appendValue(needles, ";PE=");
    appendValue(needles, ";RE=");
    
    //Find keywords, value starts at position Matchposition+Patternlength+2 and ends before semicolon
    Pattern<String<CharString>, WuManber> pattern(needles);
    unsigned endOfValue;
    while (find(finder, pattern))
    {
        for(unsigned i=endPosition(finder); i<=length(record.info); i++){
            if(record.info[i]==';' || i==length(record.info)) {
                //std::cout<< "end of hit at position " << i<< std::endl;
                endOfValue=i;
                break;
            }
        }
        //save info 
         CharString value = infix(record.info, endPosition(finder), endOfValue);
         
        if(infix(finder)==";VT="){
            vt=std::stoi(CharStringToStdString(value));
            //std::cout<<vt<<std::endl;
        }
        if(infix(finder)==";SE="){
                se=std::stod(CharStringToStdString(value));
        }
        if(infix(finder)==";PE="){
            pe=std::stod(CharStringToStdString(value));
        }
        if(infix(finder)==";RE="){
            re=std::stod(CharStringToStdString(value));
        }
            
        if(empty(se)&& empty(pe) && empty(re)){
                 return 0;
             }
        else{
            //thisHit.svtype=CharStringToStdString(record.alt);
            thisHit.chr=record.rID;
            thisHit.start=record.beginPos;
            thisHit.score=se + pe + re;
        }
    }
    return vt;
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

// void sortByType(){
//     
// }