#include <seqan/bed_io.h>
#include <seqan/vcf_io.h>
#include <seqan/find.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <algorithm> 

#include "makesvmscripts.h"

#define CharStringToStdString(x) std::string(toCString(x))

using namespace seqan;

int main()
{
    //====read original vcf file, copy header into train and test vcf file, split into train and test set by vt, write into seperate files
//     std::vector<VcfRecord> trainRecordList;
//     std::vector<VcfRecord> testRecordList;
//     VcfRecord record;
//     struct infoStore thisInfo;
//     
//     int entryNumber=0;
//     
//     // Open input file.
//      VcfFileIn vcfIn(toCString(getAbsolutePath("../47_out.vcf")));
// //     // Attach outputs
// //     //VcfFileOut vcfOut(toCString(getAbsolutePath("../47_Train.vcf")));
// //     VcfFileOut vcfOut(toCString(getAbsolutePath("../47_Test.vcf")));
//     //VcfFileOut vcfOut(vcfIn);
//     //open(vcfOut, std::cout, Vcf());
//      
//     // Copy over header.
//     VcfHeader header;
//     readHeader(header, vcfIn);
//     //writeHeader(vcfOut, header);
//     
//     //sort and copy
//     while (!atEnd(vcfIn)){
//         entryNumber++;
//         readRecord(record, vcfIn);
//         thisInfo= infoExtract(CharStringToStdString(record.info));
//         
//         if(thisInfo.vt==2){
//             testRecordList.push_back(record);
//             //writeRecord(vcfOutTest, record);
//         }
//         else if(thisInfo.vt==1 || thisInfo.vt==3){
//             trainRecordList.push_back(record);
//             //writeRecord(vcfOut, record);
//         }
//         
//     }
    
//     ======================Read Training/Testting Data Translate to libsvm format and write==================
//     VcfRecord record;
//     struct infoStore thisInfo;
//     std::ofstream libsvmformat;
//     libsvmformat.open ("47_Test_SE_CE_PE_RE_CP");
//     
//     
//     int entryNumber=0;
//     
//     // Open input file.
//      VcfFileIn vcfIn(toCString(getAbsolutePath("../../makesvmscripts-build/47_Test.vcf")));
//      VcfHeader header;
//      readHeader(header, vcfIn);
//      while (!atEnd(vcfIn)){
//         entryNumber++;
//         readRecord(record, vcfIn);
//         thisInfo= infoExtract(CharStringToStdString(record.info));
//         int sece = thisInfo.se + thisInfo.ce;
//         if(thisInfo.vt==2){
//             libsvmformat << "1" << " 1:" <<sece<< " 2:" << thisInfo.pe << " 3:" << thisInfo.re << " 4:" << thisInfo.cp;
//             libsvmformat<<"\n";
//             
//         }
// //         else if(thisInfo.vt==1){
// //             libsvmformat << "-1" << " 1:" <<sece<< " 2:" << thisInfo.pe << " 3:" << thisInfo.re << " 4:" << thisInfo.cp;
// //             libsvmformat<<"\n";
// //         }
//      }
    
    //===========================Merge Back======================================================================
    
    std::vector<VcfRecord> trainRecordList;
    std::vector<VcfRecord> testRecordList;
    VcfRecord record;
    struct infoStore thisInfo;
    
    int entryNumber=0;
    
    // Open input file.
    VcfFileIn vcfIn(toCString(getAbsolutePath("../47_out.vcf")));
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, Vcf());
     
    // Copy over header.
    VcfHeader header;
    readHeader(header, vcfIn);
    writeHeader(vcfOut, header);
    
    //sort and copy
    while (!atEnd(vcfIn)){
        entryNumber++;
        readRecord(record, vcfIn);
        thisInfo= infoExtract(CharStringToStdString(record.info));
        
        if(thisInfo.vt==2){
            testRecordList.push_back(record);
        }
        else if(thisInfo.vt==1 || thisInfo.vt==3){
            trainRecordList.push_back(record);
        }
    }
        
    //---------------process probability scores vt=2------------------------
    std::vector<double> testScores;
    std::ifstream libsvmin("47_SE_CE_PE_RE.out");
    std::string line;
    int label;
    double positiveScore;
    double negativeScore;
    
    while (std::getline(libsvmin, line)){
        std::stringstream libsvminstream(line);
        if (libsvminstream >> label >> positiveScore >> negativeScore){
            testScores.push_back(positiveScore);
        }
    }
    
    //------combine and write-------------------
    VcfRecord outrecord;
    CharString writeInfoString;
    struct infoStore outinfo;
    
    for(unsigned i=0; i<testRecordList.size();i++){
        outinfo=infoExtract(CharStringToStdString(testRecordList[i].info));
        outinfo.sc=testScores[i];
        writeInfoString=infoToCString(outinfo);
        
        outrecord.info=writeInfoString;
        outrecord.rID=testRecordList[i].rID;
        outrecord.beginPos=testRecordList[i].beginPos;
        outrecord.id=testRecordList[i].id;
        outrecord.ref=testRecordList[i].ref;
        outrecord.alt=testRecordList[i].alt;
        outrecord.qual=testRecordList[i].qual;
        outrecord.filter=testRecordList[i].filter;
        outrecord.format=testRecordList[i].format;
        outrecord.genotypeInfos=testRecordList[i].genotypeInfos;
        
        writeRecord(vcfOut, outrecord);
    }
    
    for(unsigned j=0; j<trainRecordList.size(); j++){
        outinfo=infoExtract(CharStringToStdString(trainRecordList[j].info));
        if(outinfo.vt==3){
           outinfo.sc=1; 
        }
        else if(outinfo.vt==1){
           outinfo.sc=0; 
        }
        writeInfoString=infoToCString(outinfo);
        
        outrecord.info=writeInfoString;
        outrecord.rID=trainRecordList[j].rID;
        outrecord.beginPos=trainRecordList[j].beginPos;
        outrecord.id=trainRecordList[j].id;
        outrecord.ref=trainRecordList[j].ref;
        outrecord.alt=trainRecordList[j].alt;
        outrecord.qual=trainRecordList[j].qual;
        outrecord.filter=trainRecordList[j].filter;
        outrecord.format=trainRecordList[j].format;
        outrecord.genotypeInfos=trainRecordList[j].genotypeInfos;
        
        writeRecord(vcfOut, outrecord);
    }
    
    return 0;
}

struct infoStore infoExtract(std::string originalInfo){
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
    
    struct infoStore fillInHere;
    
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

CharString infoToCString(struct infoStore infoFromStruct){
    //SC=-1.000000;VT=1;SE=0;PE=1;CE=0;RE=0.134058;RD=64.228571;GC=0.375500;CP=1.952462;SVLEN=-671;SVTYPE=DEL
    CharString outputString;
    outputString += "SC=";
    outputString += toCString(std::to_string(infoFromStruct.sc));
    outputString += ";VT=";
    outputString += toCString(std::to_string(infoFromStruct.vt));
    outputString += ";SE=";
    outputString += toCString(std::to_string(infoFromStruct.se));
    outputString += ";PE=";
    outputString += toCString(std::to_string(infoFromStruct.pe));
    outputString += ";CE=";
    outputString += toCString(std::to_string(infoFromStruct.ce));
    outputString += ";RE=";
    outputString += toCString(std::to_string(infoFromStruct.re));
    outputString += ";RD=";
    outputString += toCString(std::to_string(infoFromStruct.rd));
    outputString += ";GC=";
    outputString += toCString(std::to_string(infoFromStruct.gc));
    outputString += ";CP=";
    outputString += toCString(std::to_string(infoFromStruct.cp));
    outputString += ";SVLEN=";
    outputString += toCString(std::to_string(infoFromStruct.svlen));
    outputString += ";SVTYPE=";
    outputString += toCString(infoFromStruct.svtype);
    
    return outputString;
}