#include <seqan/vcf_io.h>
#include <seqan/find.h>

#include <iostream>
#include <fstream>
#include <limits>

#include "statistics.h"

using namespace seqan;

int main()
{

   
/*    
    std::fstream originalData("47_out.vcf");
    std::ifstream results("47TestResults");
    std::ofstream combined("47_out_ListOfStrucVar.txt");
    int isVar, lineNumber;
    
std::string line;
std::string line2;
while (std::getline(results, line))
{
    std::stringstream resultstream(line);
    
    if (resultstream >> isVar >> lineNumber )
    {
        if(isVar==1){
            GotoLine(originalData, (lineNumber+18));
            std::string lineFromOriginal;
            getline(originalData, lineFromOriginal);
            combined << lineFromOriginal <<"\n";
            
        }
    }
}*/
    
    
    //GotoLine(file, 10);

//     std::string line8;
//     file >> line8;
// 
//     std::cout << line8;
//     std::cin.get();
    
    
    
    
//     // Open input file.
//     VcfFileIn vcfIn(toCString(getAbsolutePath("../../47_out.vcf")));
// 
//     // Attach to standard output.
//     VcfFileOut vcfOut(vcfIn);
//     open(vcfOut, std::cout, Vcf());
// 
//     // Copy over header.
//     VcfHeader header;
//     readHeader(header, vcfIn);
// 
//     // Copy the file record by record into positive, negative and test-StringSet
//     //Positive: VT=3, negative: VT<2, Test: VT=2
//     VcfRecord record;
//     Info_sep infoInStruct;
//     int entryNumber=0;
//     
     //Output File

    
    
//     //Libsvm format <label> <index1>:<value1> <index2>:<value2> ... '\n'
//     std::ofstream libsvmformat;
//     libsvmformat.open ("47_out_realtrain");
//     std::ofstream svmtest;
//     svmtest.open("47_out_realtest");
//     
//     while (!atEnd(vcfIn))
//     {
//         entryNumber++;
//         std::cout<< "we are in record number " << entryNumber << std::endl;
//         //std::cout<< infix(record.id, 0, 3);
//         
//         readRecord(record, vcfIn);
//         feedInfo(record, infoInStruct);
//         if(empty(infoInStruct.se)&& empty(infoInStruct.pe) && empty(infoInStruct.re)){
//                 continue;
//             }
//         writeLibSvm(svmtest, libsvmformat, infoInStruct, entryNumber);
//         clearInfo(infoInStruct);
//     }
//     libsvmformat.close();
//     svmtest.close();
//     //split training file into two: every third for testing of algorithm, rest remains in training
//     std::ifstream originaltraining;
//     originaltraining.open("47_out_realtrain");
//     std::ofstream testing;
//     testing.open("47testSplit3");
//     std::ofstream training;
//     training.open("47trainSplit3");
//     int linecount=0;
//     
//     std::string line;
//     while (std::getline(originaltraining, line)){
//         linecount++;
//         if((linecount % 3)-2 ==0){
//             testing<<line <<"\n";
//         }
//         else{
//             training<<line<<"\n";
//         }
//     }
//     originaltraining.close();
//     testing.close();
//     training.close();
    
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


//----------------------Write outputs into extra files in libsvm format -----------------------------------------------
void writeLibSvm(std::ofstream & svmtest, std::ofstream & libsvmformat, Info_sep & infoInStruct, int entryNumber){
            if(infoInStruct.vt=="3"){
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
        //write testing file, label in first column is record ID
        else if(infoInStruct.vt=="2"){
             svmtest<<entryNumber;
            
            if(!empty(infoInStruct.se)){
                svmtest << " 1:" << infoInStruct.se;
            }
            if(!empty(infoInStruct.pe)){
                svmtest << " 2:" <<infoInStruct.pe;
            }
            if(!empty(infoInStruct.re)){
                svmtest << " 3:" <<infoInStruct.re;
            }
            svmtest<<"\n";
        }
        
}
//---------------------combine libsvm output with record info and ID --------------------------------------------------
void writeOutput(std::ifstream & testfile, std::ifstream & libsvmresult, std::ofstream & combindedResult){
    int linecount=0;
    linecount++;
}

//--------------go to line_--------
std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}