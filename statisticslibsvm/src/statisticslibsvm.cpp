#include <seqan/vcf_io.h>
#include <seqan/find.h>

#include <iostream>
#include <fstream>
#include <limits>

#include <ctype.h>
#include <stdlib.h>
#include <vector>

#include "statisticslibsvm.h"
#include "svm.h"



#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define CharStringToStdString(x) std::string(toCString(x))


//initialize lib svm structures
struct svm_parameter param;		
struct svm_problem prob;		
struct svm_node *node; //(index,value) pair pointers
struct svm_model *model;

//=============================================================================================
int main(){
    
//------define parameters for libsvm (now default)------------------------------------------------
    struct svm_parameter param;
    param.svm_type = C_SVC;
    param.kernel_type = RBF;
    param.degree = 3;
    param.gamma = 0;	// 1/num_features
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 100;
    param.C = 1;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    
//------------------storage and variables for training and testing data-------------------------------------
    int attributeNumber =3; //now SE, PE, RE, evtl. change or parse
    std::vector<Info_sep> testInfo;
    std::vector<Info_sep> trainInfo;
    Info_sep infoInStruct;
    VcfRecord record;
    int entryNumber=0;
    std::vector<int> labelsTrain;
    //labelsTrain.reserve(1000000);
    std::vector<std::vector<double> > attributesTrain;
    //attributesTrain.reserve(1000000);
    std::vector<int> labelsTest;
    //labelsTest.reserve(1000000);
    std::vector<std::vector<double> > attributesTest;
    //attributesTest.reserve(1000000);
    std::vector<double> attributes(attributeNumber);
    for(unsigned i=0; i < attributeNumber; ++i){
        attributes.push_back(0);
    }
//-------read vcf file, while reading need to get: number of positives+negatives, number of unknown, generate int vector for labels, generate matrix for values-----
    // Open input file.
    VcfFileIn vcfIn(toCString(getAbsolutePath("../../47_out.vcf")));
    // Attach to standard output.
    VcfFileOut vcfOut(vcfIn);
    open(vcfOut, std::cout, Vcf());
    // Copy over header.
    VcfHeader header;
    readHeader(header, vcfIn);
    while (!atEnd(vcfIn)){
        entryNumber++;
        readRecord(record, vcfIn);
        feedInfo(entryNumber, record, infoInStruct);
        if(empty(infoInStruct.se)&& empty(infoInStruct.pe) && empty(infoInStruct.re)){
                 continue;
             }
        //negative case
        if(infoInStruct.vt=="1"){
            labelsTrain.push_back(-1);
            attributes[0]=(std::stod(infoInStruct.se));
            attributes[1]=(std::stod(infoInStruct.pe));
            attributes[2]=(std::stod(infoInStruct.re));
            
            attributesTrain.push_back(attributes);
        }
         //positive case
        else if(infoInStruct.vt=="3"){
            labelsTrain.push_back(1);
            attributes[0]=(std::stod(infoInStruct.se));
            attributes[1]=(std::stod(infoInStruct.pe));
            attributes[2]=(std::stod(infoInStruct.re));
            
            attributesTrain.push_back(attributes);
        }
         //unknown case
        else if(infoInStruct.vt=="2"){
            labelsTest.push_back(0);
            attributes[0]=(std::stod(infoInStruct.se));
            attributes[1]=(std::stod(infoInStruct.pe));
            attributes[2]=(std::stod(infoInStruct.re));
            
            attributesTest.push_back(attributes);
        }
        
        clearInfo(infoInStruct);
    }
    
//---------initialization for libsvm-------------------------------------------------------------------------
    struct svm_problem prob;
    //struct svm_node *nodePointer; //(index,value) pair pointers
    //struct svm_model *model;
    //get size
    attributesTrain.shrink_to_fit();
    attributesTest.shrink_to_fit();
    labelsTrain.shrink_to_fit();
    labelsTest.shrink_to_fit();
    
//     std::vector<double> forprint;
//         for(unsigned i=0; i<attributesTrain.size(); i++){
//             forprint=attributesTrain[i];
//             std::cout<<labelsTrain[i]<<" 1:"<<forprint[0]<<" 2:"<<forprint[1]<<" 3:"<<forprint[2]<<std::endl;
//         }
        
    int numberOfLinesTest=attributesTest.size();
    int numberOfLinesTrain=attributesTrain.size();
    std::cout<<"vcf file read, "<< numberOfLinesTrain<<" records for training"<<std::endl;
// //------fill in the libsvm problem structure and the libsvm node structure, synchronise the addresses-----------
    prob.l=numberOfLinesTrain;
    
    prob.y = Malloc(double,prob.l); //target values
    for (unsigned i=0; i < prob.l; i++){
        prob.y[i] = labelsTrain[i];
        }
        
    prob.x = Malloc(struct svm_node *, prob.l);
    //nodePointer= Malloc(struct svm_node, (attributeNumber+1));    
    
    for (unsigned entryRow=0; entryRow< prob.l; entryRow++){
        svm_node* nodePointer = Malloc(svm_node,(attributeNumber+1));
        for(unsigned attributeColumn=0; attributeColumn<attributeNumber; attributeColumn++){
            nodePointer[attributeColumn].index=attributeColumn+1; //columns are calculated from 0 but indices from 1 except for precomputed kernel
            nodePointer[attributeColumn].value=attributesTrain[entryRow][attributeColumn];
        }
        //in libsvm index = -1 indicates the end of one vector (readme libsvm)
        nodePointer[attributeNumber].index=-1;
        nodePointer[attributeNumber].value=0;
        
        //hand over address of the nodes-vector to prob.x
        prob.x[entryRow] = nodePointer;
        free(nodePointer);
    }
    
    //print
//     	for (int i = 0; i < prob.l; i++)
// 	{
// 		std::cout<<prob.y[i]<<" ";
// 		for (int k = 0; k < attributeNumber; k++)
// 		{
// 			int index = (prob.x[i][k].index);
// 			double value = (prob.x[i][k].value); 
// 			std::cout<<index<<":"<<value<<" ";
// 		}
// 		std::cout<<std::endl;
// 	}
//         std::cout<<"all ok"<<std::endl;
    
//-------make model-------------------------------------------------------------------------------
    std::cout<<"making model"<<std::endl;
    svm_model *model = svm_train(&prob,&param); //train
    svm_save_model("modelFromImplem.model", model); //save
    int check= svm_check_probability_model(model);
    std::cout<<"model made, beginning prediction"<<std::endl;
    
//-------predict---------------------------------------------------------------------------------
    //make node and predict for every entry in testing data (attributesTest), write result to outputfile
    std::ofstream svmout;
    svmout.open("predictFromImplem.out");
    
    for(unsigned testInstance=0; testInstance<numberOfLinesTest; testInstance++){
        svm_node* testNode = Malloc(svm_node,(attributeNumber+1));
        for(unsigned attributeColumn=0; attributeColumn<attributeNumber; attributeColumn++){
            testNode[attributeColumn].index=attributeColumn+1; //columns are calculated from 0 but indices from 1 except for precomputed kernel
            testNode[attributeColumn].value=attributesTest[testInstance][attributeColumn];
        }
        //in libsvm index = -1 indicates the end of one vector (readme libsvm)
        testNode[attributeNumber].index=-1;
        testNode[attributeNumber].value=0;
        double result=svm_predict(model, testNode);
        svmout<<result<<std::endl;
        free(testNode);
    }
    svmout.close();
    
    std::cout<<"predictions written to predictFromImplem.out" <<std::endl;
    
//---------clean up------------------------------------------------------------------------------------
    svm_destroy_param(&param);
    free(prob.y);
    free(prob.x);
    svm_free_model_content(model);
    
    return 0;
}


//==========================================================================================
//clear up saved separated record info
 void clearInfo(Info_sep & infoInStruct){
     infoInStruct.sc.clear();
    infoInStruct.vt.clear();
    infoInStruct.se.clear();
    infoInStruct.pe.clear();
    infoInStruct.ce.clear();
    infoInStruct.re.clear();
    infoInStruct.rd.clear();
    infoInStruct.gc.clear();
    infoInStruct.cp.clear();
    infoInStruct.svlen.clear();
    infoInStruct.svtype.clear();
    infoInStruct.recordID.clear();
    infoInStruct.entryNumber=0;
 }

//--------------save record info from vcf as struct------------------------------------------
void feedInfo(int entryNumber, seqan::VcfRecord record, Info_sep & infoInStruct){
    //save recordID and entryNumber
    infoInStruct.entryNumber= entryNumber;
    infoInStruct.recordID=CharStringToStdString(record.id);
    
    
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
        infoInStruct.sc=CharStringToStdString(value);
        }
        if(infix(finder)==";VT="){
        infoInStruct.vt=CharStringToStdString(value);
        }
        if(infix(finder)==";SE="){
        infoInStruct.se=CharStringToStdString(value);
        }
        if(infix(finder)==";PE="){
        infoInStruct.pe=CharStringToStdString(value);
        }
        if(infix(finder)==";CE="){
        infoInStruct.ce=CharStringToStdString(value);
        }
        if(infix(finder)==";RE="){
        infoInStruct.re=CharStringToStdString(value);
        }
        if(infix(finder)==";RD="){
        infoInStruct.rd=CharStringToStdString(value);
        }
        if(infix(finder)==";GC="){
        infoInStruct.gc=CharStringToStdString(value);
        }
        if(infix(finder)==";CP="){
        infoInStruct.cp=CharStringToStdString(value);
        }
        if(infix(finder)==";SVLEN="){
        infoInStruct.svlen=CharStringToStdString(value);
        }
        if(infix(finder)==";SVTYPE="){
        infoInStruct.svtype=CharStringToStdString(value);
        }
    }
    
}