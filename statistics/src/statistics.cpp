#include <seqan/vcf_io.h>
#include <seqan/find.h>
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
//     writeHeader(vcfOut, header);

    // Copy the file record by record into positive, negative and test-StringSet
        //Positive: VT=3, negative: VT<2, Test: VT=2
    VcfRecord record;
    Info_sep infoInStruct;
    int entryNumber=0;
    while (!atEnd(vcfIn))
    {
        //entryNumber++;
        //std::cout<< "we are in record number " << entryNumber << std::endl;
        readRecord(record, vcfIn);
        feedInfo(record, infoInStruct);
    }

    return 0;
}


void feedInfo(seqan::VcfRecord record, Info_sep & infoInStruct){
    //loop over Info CharString, pattern matching for keywords like SC, copy everything between keyword= and ; intro entry in Info_sep structure
    
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
        //get value
        CharString value = infix(record.info, endPosition(finder), endOfValue);
        
        //std::cout << position(pattern) << '\t' << infix(finder) << '\t' << infix(record.info, endPosition(finder), endOfValue) << std::endl;
    }

    
}