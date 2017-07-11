#include <seqan/vcf_io.h>
#include <seqan/find.h>

#include <iostream>
#include <fstream>

#include "statistics.h"

using namespace seqan;

int main()
{
//     Cheat Sheet for VCfFile record names in seqan class
//     class VcfRecord
//     {
//     public:
//     int32_t rID;                          // CHROM
//     int32_t beginPos;                     // POS
//     CharString id;                        // ID
//     CharString ref;                       // REF
//     CharString alt;                       // ALT
//     float qual;                           // QUAL
//     CharString filter;                    // FILTER
//     CharString info;                      // INFO
//     CharString format;                    // FORMAT
//     StringSet<CharString> genotypeInfos;  // <individual1> <individual2> ..
// 
//     // Constants for marking reference id and position as invalid.
//     static const int32_t INVALID_REFID = -1;
//     static const int32_t INVALID_POS = -1;
//     // This function returns the float value for "invalid quality".
//     static float MISSING_QUAL();
//     };
    
    //Format for svm:
    //
    
    //Current strategy: 
    //Make my own class derived from seqan vcf class but info string additionally seperated
    
    // Open input file.
    VcfFileIn vcfIn(toCString(getAbsolutePath("../../47_out.vcf")));

    // Attach to standard output.
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, Vcf());

    // Copy over header.
    VcfHeader header;
    readHeader(header, vcfIn);

    // Copy the file record by record into positive, negative and test-StringSet
    //Positive: VT=3, negative: VT<2, Test: VT=2
    VcfRecord record;
    Info_sep infoInStruct;
    int entryNumber=0;
    
    //Output File
    //Libsvm format <label> <index1>:<value1> <index2>:<value2> ... '\n'
    std::ofstream libsvmformat;
    libsvmformat.open ("statistics_scale");
    
    while (!atEnd(vcfIn))
    {
        entryNumber++;
        readRecord(record, vcfIn);
        feedInfo(record, infoInStruct);
        std::cout<< "we are in record number " << entryNumber << std::endl;
        //write output
        if(empty(infoInStruct.se)&& empty(infoInStruct.pe) && empty(infoInStruct.re)){
                continue;
            }
        else if(infoInStruct.vt=="3"){
            libsvmformat<<"1";
            
            if(!empty(infoInStruct.se)){
                libsvmformat << " 1:" << infoInStruct.se;
            }
            if(!empty(infoInStruct.pe)){
                libsvmformat << " 2:" <<infoInStruct.pe;
            }
            if(!empty(infoInStruct.re)){
                libsvmformat << " 3:" <<infoInStruct.re;
            }
            libsvmformat<<"\n";
        }
        else if(infoInStruct.vt=="1"){
            libsvmformat<<"-1";
            
            if(!empty(infoInStruct.se)){
                libsvmformat << " 1:" << infoInStruct.se;
            }
            if(!empty(infoInStruct.pe)){
                libsvmformat << " 2:" <<infoInStruct.pe;
            }
            if(!empty(infoInStruct.re)){
                libsvmformat << " 3:" <<infoInStruct.re;
            }
            libsvmformat<<"\n";
        }
        
        clearInfo(infoInStruct);
    }
    return 0;
}

//==========================================================================================================================
//save seperated record info
void feedInfo(seqan::VcfRecord record, Info_sep & infoInStruct){
    //loop over Info CharString, pattern matching for keywords like SC, copy everything between keyword= and ; intro entry in Info_sep structure and write to outputfile for libsvm
    
    //Processing of the Recordinfo
    Finder<CharString> finder(record.info);
    //Prepare Patterns(keywords)
    String<CharString> needles;
    appendValue(needles, "SC=");
    appendValue(needles, ";VT=");
    appendValue(needles, ";SE=");
    appendValue(needles, ";PE=");
    appendValue(needles, ";CE=");
    appendValue(needles, ";RE=");
    appendValue(needles, ";RD=");
    appendValue(needles, ";GC=");
    appendValue(needles, ";CP=");
    appendValue(needles, ";SVLEN=");
    appendValue(needles, ";SVTYPE=");
    
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
        if(infix(finder)=="SC"){
        infoInStruct.sc=value;
        }
        if(infix(finder)==";VT="){
        infoInStruct.vt=value;
        }
        if(infix(finder)==";SE="){
        infoInStruct.se=value;
        }
        if(infix(finder)==";PE="){
        infoInStruct.pe=value;
        }
        if(infix(finder)==";CE="){
        infoInStruct.ce=value;
        }
        if(infix(finder)==";RE="){
        infoInStruct.re=value;
        }
        if(infix(finder)==";RD="){
        infoInStruct.rd=value;
        }
        if(infix(finder)==";GC="){
        infoInStruct.gc=value;
        }
        if(infix(finder)==";CP="){
        infoInStruct.cp=value;
        }
        if(infix(finder)==";SVLEN="){
        infoInStruct.svlen=value;
        }
        if(infix(finder)==";SVTYPE="){
        infoInStruct.svtype=value;
        }
        
        
        //std::cout << position(pattern) << '\t' << infix(finder) << '\t' << infix(record.info, endPosition(finder), endOfValue) << std::endl;
    }
//     libsvmformat.close();
}

//---------------------------------------------------------------
//clear up saved separated record info
void clearInfo(Info_sep & infoInStruct){
    clear(infoInStruct.sc);
    clear(infoInStruct.vt);
    clear(infoInStruct.se);
    clear(infoInStruct.pe);
    clear(infoInStruct.ce);
    clear(infoInStruct.re);
    clear(infoInStruct.rd);
    clear(infoInStruct.gc);
    clear(infoInStruct.cp);
    clear(infoInStruct.svlen);
    clear(infoInStruct.svtype);
}

