//
//  individual.c
//  HCV Model
//
//  Created by Neil Bretana on 1/10/2014.
//  Copyright (c) 2014 Neil Bretana. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include "individual.h"

#define COLCTR 12  //Array of locations containing:
//[row][0] IDU+; HCV+; ATSI
//[row][1] IDU+; HCV+; NON-ATSI
//[row][2] IDU+; HCV-; ATSI; PREV. EXPOSED
//[row][3] IDU+; HCV-; ATSI; SUSCEPTIBLE
//[row][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED
//[row][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE
//[row][6] IDU-; HCV+; ATSI
//[row][7] IDU-; HCV+; NON-ATSI
//[row][8] IDU-; HCV-; ATSI; PREV. EXPOSED
//[row][9] IDU-; HCV-; ATSI; SUSCEPTIBLE
//[row][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED
//[row][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE
#define ROWPRIS 4 //# prisons +1
#define RGROUPS 7 //Risk groups: 0: non-injecting; 1: injecting less than daily, not sharing; 2: injecting daily or more, not sharing; 3: injecting less than daily, sharing less than daily; 4: injecting less than daily, sharing daily or more; 5: injecting daily or more, sharing less than daily; 6: injecting daily or more, sharing daily or more

int generateRand(){
    int rNum;
    
    rNum= rand();
    //rNum = rand()%10 +1;
    return rNum;
}

int generateAge(int minimum_number, int max_number){
    int rNum;
    rNum = rand() % (max_number + 1 - minimum_number) + minimum_number;
    return rNum;
}

int draw_multinom(gsl_rng **r, int nEvents, double probsInput[]){
    int nIndex=999, iFlag=0, iCrawl=0;
    unsigned outArray[nEvents];

    gsl_ran_multinomial(*r, nEvents, 1, probsInput, outArray);

    while(iFlag==0){
        if(outArray[iCrawl]==1){
            nIndex=iCrawl;
            iFlag=1;
        }
        iCrawl++;
    }
    //printf("draw: %d\n", nIndex);
    return nIndex;
}

int distributePop(int input){
    int n;
    do
        n=generateRand();
    while(n<=0||n>input);
    return n;
}

void newIndiv (int *idGlobal, int *COVIDentry, int *COVIDentryAb, int currDay, sIndiv **pHeadCopy, sIndiv **pTailCopy, gsl_rng **r, int nEvents, int (*pLocArray2)[ROWPRIS][COLCTR], int prison, int category, int seedPop, int seedZone, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *povercap){
    sIndiv *obj;
    int normI;
    int objAge, objLoc, objMeta, binAge, binCM, binATSI, binCOVID, binPREV, binInjFrq, binOpd, binSha, binShaFrq, binEverIDU, binZone, categoryFlag=0;
    int cellFlag=0;
    int cellCap=0;
    int mCtr1, mCtr2, mCtr3, mCtr4;
    float rMin, rMed, rMax, lb, ub;
    double normSum;
    double historyIDU;
    double probsMove[ROWPRIS];//array to move people from community to prison
    double probsMeta[6];
    double probsRisk[2]; //for proportion of IDUs
    double probsR[2];
    double probsAge[7];
    double probsMaxZone[18], probsMedZone[11], probsMinZone[27], probsMinArea[2], probsMedArea[2], probsMaxArea[4], probsMinUnit[6], probsMedUnit[4], probsMaxUnit[5], probsMinCell[13], probsMedCell[19], probsMaxCell[20];
    //Age 0-19
    //Age 20-44
    //Age 45-54
    //Age 55-64
    //Age 65-74
    //Age 75-84
    //Age >= 85
    gsl_rng *rCopy;
    rCopy=*r;
    
    obj=(sIndiv *)malloc(sizeof(sIndiv));
    
    //Fill in characteristics
    obj->ID=*idGlobal;
    (*idGlobal)++;
    
    //Age using Jan 2020 distribution from MRRC
    probsAge[0]=0.02;//Age 0-19
    probsAge[1]=0.74;//Age 20-44
    probsAge[2]=0.15;//Age 45-54
    probsAge[3]=0.05;//Age 55-64
    probsAge[4]=0.01;//Age 65-74
    probsAge[5]=0.003;//Age 75-84
    probsAge[6]=0.0009;//Age >= 85
    binAge=draw_multinom(&rCopy, 7, probsAge);
    printf("Age group selected: %d\n", binAge);
    switch(binAge){
        case 0:
            printf("Selecting between 0-19\n");
            objAge=generateAge(1, 19);
            printf("Age selected: %d\n", objAge);
            break;
        case 1:
            printf("Selecting between 20-44\n");
            objAge=generateAge(20, 44);
            printf("Age selected: %d\n", objAge);
            break;
        case 2:
            printf("Selecting between 45-54\n");
            objAge=generateAge(45, 54);
            printf("Age selected: %d\n", objAge);
            break;
        case 3:
            printf("Selecting between 55-64\n");
            objAge=generateAge(55, 64);
            printf("Age selected: %d\n", objAge);
            break;
        case 4:
            printf("Selecting between 65-74\n");
            objAge=generateAge(65, 74);
            printf("Age selected: %d\n", objAge);
            break;
        case 5:
            printf("Selecting between 75-84\n");
            objAge=generateAge(75, 84);
            printf("Age selected: %d\n", objAge);
            break;
        case 6:
            printf("Selecting between 85-99\n");
            objAge=generateAge(85, 99);
            printf("Age selected: %d\n", objAge);
            break;
    }
    
    //do
    //    objAge=generateRand(); //REPLACE WITH AGE DISTRIBUTION
    
    //while(objAge<18||objAge>100);
    obj->age=objAge;
    
    //set DAA OST NSP
    //obj->DAA=0;
    //obj->OST=0;
    
    //SET GENDER
    
    //BASIC ATTRBUTES
    if(category==0||category==2||category==3||category==6||category==8||category==9){
        obj->atsi=1;
        obj->timeOfImprisonment=0;
    }else if(category==1||category==4||category==5||category==7||category==10||category==11){
        obj->atsi=0;
        obj->timeOfImprisonment=0;
    }else if(category==99){//new individuals from the community
        //assign injecting category
        //check first if enough people in this category
        //choose if injector or not p=0.16 (dummy: NSW Inmate Health Survey 2009)
        //loop for if category count is 0
        if(prison!=99&&seedZone==99){ //seed pop
            obj->timeOfImprisonment=0;
        }else{ //not seed pop
            obj->timeOfImprisonment=currDay;
        }
        
        printf("category99\n");
        
        while(categoryFlag==0){
            //binIDU=gsl_ran_binomial(*r, 0.16, 1);
            //binIDU=gsl_ran_binomial(*r, 0, 0.16);
            probsRisk[0]=1-0.30;//.20//1-0.164; //does not have co-morbidity
            probsRisk[1]=0.30;//0.20//0.164; //probability of having co-morbidity
            binCM=draw_multinom(&rCopy, 2, probsRisk);
            printf("draw attributes\n");
            if(binCM==1){ //has co-morbidity
                obj->comorbidity=1; //comorb
                //probability of being ATSI; p=.204 (NSW Inmate Health Survey 2009)
                //binATSI=gsl_ran_binomial(*r, 0.204, 1);
                probsRisk[0]=1-0.20; //not ATSI
                probsRisk[1]=0.20; //probability of being ATSI
                binATSI=draw_multinom(&rCopy, 2, probsRisk);
                printf("co-morbidity Y\n");
                if(binATSI==1){
                    obj->atsi=1; //IDU ATSI
                    //probability of being infected with HCV upon entry IDU ATSI p= 36% antibody prev
                    //binHCV=gsl_ran_binomial(*r, 0.36, 1);
                    probsRisk[0]=1-0.68;//0.70//0.60//1-0.36; //no COVID
                    probsRisk[1]=0.68;//0.70//0.60//0.36; //probability of being infected with COVID
                    printf("ATSI Y\n");
                    if(seedPop==1){
                        binCOVID=0;
                    }else{
                        binCOVID=1;
                        //binCOVID=draw_multinom(&rCopy, 2, probsRisk);
                    }
                    
                    if(binCOVID==1){ //IDU ATSI HCV
                        //(*COVIDentry)++;
                        //probability of being viremic (as opposed to prev. exposed) IDU ATSI p= 80% antibody prev
                        //binPREV=gsl_ran_binomial(*r, 0.80, 1);
                        probsRisk[0]=1-1.0;//1-0.75;//1-0.80; //probability of not being viraemic (Has antibodies)
                        probsRisk[1]=1.0;//0.75;//0.80; //probability of being viraemic (no antibodies)
                        printf("COVID Y\n");
                        binPREV=draw_multinom(&rCopy, 2, probsRisk);
                        if(binPREV==1){
                            //Co-Morbidity+; COVID+; ATSI
                            obj->COVID=1;
                            obj->COVIDAb=0;
                            category=0;
                            (*COVIDentry)++;
                        }else{
                            //Co-Morbidity+; COVID-; ATSI; Antibody+
                            obj->COVID=0;
                            obj->COVIDAb=1;
                            category=2;
                            (*COVIDentryAb)++;
                        }
                    }else{ //Co-Morbidity+; COVID-; ATSI; Antibody-
                        obj->COVID=0;
                        obj->COVIDAb=0;
                        category=3;
                    }
                }else{
                    obj->atsi=0; //IDU non-ATSI
                    //probability of being infected with HCV upon entry IDU non-ATSI p= 23.9% antibody prev
                    //binHCV=gsl_ran_binomial(*r, 0.239, 1);
                    probsRisk[0]=1-0.50;//0.524;//0.424//1-0.239; //not IDU
                    probsRisk[1]=0.50;//0.524;//0.474//0.239; //probability of being an IDU
                    if(seedPop==1){
                        binCOVID=0;
                    }else{
                        binCOVID=1;
                        //binCOVID=draw_multinom(&rCopy, 2, probsRisk);
                    }
                    printf("ATSI N\n");
                    if(binCOVID==1){ //IDU non-ATSI HCV
                        //(*HCVentryAb)++;
                        //probability of being viremic (as opposed to prev. exposed) IDU ATSI p= 80% antibody prev
                        //binPREV=gsl_ran_binomial(*r, 0.80, 1);
                        probsRisk[0]=1-1.0;//1-0.80; //probability of not being viraemic (Has antibodies)
                        probsRisk[1]=1.0;//0.75;//0.80; //probability of being viraemic (Has NO antibodies)
                        binPREV=draw_multinom(&rCopy, 2, probsRisk);
                        printf("COVID N\n");
                        if(binPREV==1){
                            //Co-Morbidity+; COVID+; Non-ATSI
                            obj->COVID=1;
                            obj->COVIDAb=0;
                            category=1;
                            (*COVIDentry)++;
                        }else{
                            //Co-Morbidity+; COVID-; Non-ATSI; Antibody+
                            obj->COVID=0;
                            obj->COVIDAb=1;
                            category=4;
                            (*COVIDentryAb)++;
                        }
                    }else{ //IDU ATSI HCV- NEVER BEEN EXPOSED
                        obj->COVID=0;
                        obj->COVIDAb=0;
                        category=5;
                    }
                }
            }else{ //has NO co-morbidity
                obj->comorbidity=0; //comorb
                //probability of being ATSI for non-injectors; p=.204 (NSW Inmate Health Survey 2009)
                //binATSI=gsl_ran_binomial(*r, 0.16, 1);
                probsRisk[0]=1-0.20; //not IDU
                probsRisk[1]=0.20; //probability of being an IDU
                binATSI=draw_multinom(&rCopy, 2, probsRisk);
                printf("Co-morbidity N\n");
                if(binATSI==1){
                    obj->atsi=1; //non-IDU ATSI
                    //probability of being infected with HCV upon entry IDU ATSI p= dummy
                    //binHCV=gsl_ran_binomial(*r, 0.05, 1);
                    probsRisk[0]=1-0.37;//0.39;//0.29//1-0.05; //not IDU
                    probsRisk[1]=0.37;//0.39;//0.29//0.05; //probability of being an IDU
                    printf("ATSI Y\n");
                    if(seedPop==1){
                        binCOVID=0;
                    }else{
                        binCOVID=1;
                        //binCOVID=draw_multinom(&rCopy, 2, probsRisk);
                    }
                    if(binCOVID==1){ //non-IDU ATSI HCV
                        //(*HCVentryAb)++;
                        //probability of being viremic (as opposed to prev. exposed) IDU ATSI p= 80% antibody prev
                        //binPREV=gsl_ran_binomial(*r, 0.80, 1);
                        probsRisk[0]=1-1.0;//1-0.75;//1-0.80; //probability of not being viraemic (Has antibodies)
                        probsRisk[1]=1.0;//0.75;//0.80;  //probability of being viraemic (Has NO antibodies)
                        binPREV=draw_multinom(&rCopy, 2, probsRisk);
                        printf("COVID Y\n");
                        if(binPREV==1){
                            //Co-Morbidity-; COVID+; ATSI;
                            obj->COVID=1;
                            obj->COVIDAb=0;
                            category=6;
                            (*COVIDentry)++;
                        }else{
                            //Co-Morbidity-; COVID-; ATSI; Antibody+
                            obj->COVID=0;
                            obj->COVIDAb=1;
                            category=8;
                            (*COVIDentryAb)++;
                        }
                    }else{ //non-IDU ATSI HCV- NEVER BEEN EXPOSED
                        obj->COVID=0;
                        obj->COVIDAb=0;
                        category=9;
                    }
                }else{
                    obj->atsi=0; //non-IDU non-ATSI
                    //probability of being infected with HCV upon entry IDU ATSI p= dummy
                    //binHCV=gsl_ran_binomial(*r, 0.01, 1);
                    probsRisk[0]=1-0.33;//0.35;//0.30//1-0.01; //not IDU
                    probsRisk[1]=0.33;//0.35;//0.30//0.01; //probability of being an IDU
                    printf("ATSI N\n");
                    if(seedPop==1){
                        binCOVID=0;
                    }else{
                        binCOVID=1;
                        //binCOVID=draw_multinom(&rCopy, 2, probsRisk);
                    }
                    if(binCOVID==1){ //non-IDU ATSI HCV
                        //(*HCVentryAb)++;
                        //probability of being viremic (as opposed to prev. exposed) IDU ATSI p= 80% antibody prev
                        //binPREV=gsl_ran_binomial(*r, 0.80, 1);
                        probsRisk[0]=1-1.0;//1-0.75;//1-0.80;  //probability of not being viraemic (Has antibodies)
                        probsRisk[1]=1.0;//0.75;//0.80;  //probability of being viraemic (Has NO antibodies)
                        binPREV=draw_multinom(&rCopy, 2, probsRisk);
                        printf("COVID Y\n");
                        if(binPREV==1){
                            //Co-Morbidity-; COVID+; Non-ATSI;
                            obj->COVID=1;
                            obj->COVIDAb=0;
                            category=7;
                            (*COVIDentry)++;
                        }else{
                            //Co-Morbidity-; COVID-; Non-ATSI; Antibody+
                            obj->COVID=0;
                            obj->COVIDAb=1;
                            category=10;
                            (*COVIDentryAb)++;
                        }
                    }else{ //non-IDU ATSI HCV- NEVER BEEN EXPOSED
                        obj->COVID=0;
                        obj->COVIDAb=0;
                        category=11;
                    }
                }
                //probability of being infected with HCV upon entry
            }
            
            //if((*pLocArray2)[0][category]>0){
                categoryFlag=1;
            //}
        }//while
    }
    
    obj->tested=0;
    obj->testDate=0;
    obj->testDaysBeforeResult=0;
    obj->testResultTrue=0;
    obj->testAccuracy=0;
    obj->TP=0;
    obj->FP=0;
    obj->TN=0;
    obj->FN=0;
    
    printf("category selected\n");
    //Infected?
    if(category==0||category==1||category==6||category==7){
        probsMeta[0]=1.0; //exposed
        probsMeta[1]=0.0; //pre-clinical
        probsMeta[2]=0.0; //asymptomatic
        probsMeta[3]=0.0; //mild
        probsMeta[4]=0.0; //moderate
        probsMeta[5]=0.0; //severe
        objMeta=draw_multinom(&rCopy, 6, probsMeta);
        obj->severity=objMeta;
        
        obj->COVID=1;
        obj->COVIDAb=0;
        obj->timeOfInfection=0;//set to 0 everytime cleared.
        obj->timeOfProgression=0;
        obj->placeInfected=0; //NA for now
        
        //obj->infectionNumber=0; //counting only reinfected in prison
    }else{
        if(category==2||category==4||category==8||category==10){
            obj->severity=6; //previously infected or clearer
            //obj->infectionNumber=0; //counting only reinfections in prison
            obj->COVID=0;
            obj->COVIDAb=1;
            obj->timeOfInfection=0;//set to 0 everytime cleared.
            obj->timeOfProgression=0;
            obj->placeInfected=0; //NA for now
        }
        else if(category==3||category==5||category==9||category==11){
            obj->severity=7; //susceptible
            //obj->infectionNumber=0; //counting only reinfections in prison
            obj->COVID=0;
            obj->COVIDAb=0;
            obj->timeOfInfection=0;//set to 0 everytime cleared.
            obj->timeOfProgression=0;
            obj->placeInfected=0; //NA for now
        }
    }
    
    if(prison==99){ //random injector from community to be converted into individual
        printf("prison99\n");

        //ASSUMING B option selected
        //(*pLocArray2)[0][category]--; //subtract this individual from community population count
        //End assumption
        
        objLoc=11; //Maximum prison reception first, inmates
        //if(objLoc==2){
            //objLoc=5;
        //}else if(objLoc==3){
            //objLoc=9;
        //}
        printf("prison set\n");
        //set where individual got infected
        if(category==0||category==1||category==6||category==7){
            obj->placeInfected=4; //infected in community
        }else{
            obj->placeInfected=0; //not infected
        }
        
    }else{ //input from user
        objLoc=prison;
        
        if(currDay==0){
            obj->placeInfected=0;
        }else{
            //set where individual got infected
            if(category==0||category==1||category==6||category==7){
                obj->placeInfected=0; //infected in community
            }else{
                obj->placeInfected=0; //not infected
            }
        }
    }
    obj->group=category; //printf("new indiv group: %d . new indiv loc: %d\n", category, objLoc);
    
    if(objLoc==1||objLoc==2||objLoc==3||objLoc==4||objLoc==5){
        obj->location=1;
    }else if(objLoc==6||objLoc==7||objLoc==8||objLoc==9||objLoc==10){
        obj->location=2;
    }else if(objLoc==11||objLoc==12||objLoc==13||objLoc==14||objLoc==15){
        obj->location=3;
    }
    
    obj->locType=0;
    obj->locArea=0;
    obj->locPrison=0;
    obj->locUnit=0;
    obj->locPod=0;
    obj->locCell=0;
    
    obj->cohorting=0;
    obj->cohortNumber=0;
    obj->cohortMaxDays=0;
    obj->timeStartCohort=0;
    
    //printf("ID: %d, ATSI: %d, Location: %d, metavir: %d, injFreq: %d, sharing: %d, shaFreq: %d, injOpd: %d\n", obj->ID, obj->atsi, obj->location, obj->metavir, obj->injFreq, obj->sharing, obj->shaFreq, obj->injOpd);
    (*pLocArray2)[objLoc][category]++;//add count of at risk or not-at-risk individuals
    
    //SET TYPE according to "prison"
    if(objLoc==1||objLoc==6||objLoc==11){
        obj->indivType=0;//inmate
    }else if(objLoc==2||objLoc==7||objLoc==12){
        obj->indivType=1;//prison staff
    }else if(objLoc==3||objLoc==8||objLoc==13){
        obj->indivType=2;//healthcare staff
    }else if(objLoc==4||objLoc==9||objLoc==14){
        obj->indivType=3;//essential visitor
    }else if(objLoc==5||objLoc==10||objLoc==15){
        obj->indivType=4;//family visitor
    }

    if(obj->location==1){
        if(obj->indivType==0){
            probsMinZone[0]=0.04;
            probsMinZone[1]=0.004;
            probsMinZone[2]=0.004;
            probsMinZone[3]=0.12;
            probsMinZone[4]=0.014;
            probsMinZone[5]=0.010;
            probsMinZone[6]=0.11;
            probsMinZone[7]=0.003;
            probsMinZone[8]=0.03;
            probsMinZone[9]=0.04;
            probsMinZone[10]=0.03;
            probsMinZone[11]=0.02;
            probsMinZone[12]=0.01;
            probsMinZone[13]=0.006;
            probsMinZone[14]=0.03;
            probsMinZone[15]=0.06;
            probsMinZone[16]=0.04;
            probsMinZone[17]=0.06;
            probsMinZone[18]=0.08;
            probsMinZone[19]=0.02;
            probsMinZone[20]=0.03;
            probsMinZone[21]=0.06;
            probsMinZone[22]=0.03;
            probsMinZone[23]=0.04;
            probsMinZone[24]=0.058;
            probsMinZone[25]=0.005;
            probsMinZone[26]=0.05;
            obj->zone=draw_multinom(&rCopy, 27, probsMinZone);
        
            /*
            //set units
            //obj->locUnit=generateAge(1, 10);
            probsMinArea[0]=0.5;
            probsMinArea[1]=0.5;
            obj->locArea=draw_multinom(&rCopy, 2, probsMinArea);
            //obj->locArea=generateAge(1, 2);
            
            probsMinUnit[0]=0.16666666666;
            probsMinUnit[1]=0.16666666666;
            probsMinUnit[2]=0.16666666666;
            probsMinUnit[3]=0.16666666666;
            probsMinUnit[4]=0.16666666666;
            probsMinUnit[5]=0.16666666666;
            obj->locUnit=draw_multinom(&rCopy, 6, probsMinUnit);
            //obj->locUnit=generateAge(1, 6);
            
            probsMinCell[0]=0.07692307692;
            probsMinCell[1]=0.07692307692;
            probsMinCell[2]=0.07692307692;
            probsMinCell[3]=0.07692307692;
            probsMinCell[4]=0.07692307692;
            probsMinCell[5]=0.07692307692;
            probsMinCell[6]=0.07692307692;
            probsMinCell[7]=0.07692307692;
            probsMinCell[8]=0.07692307692;
            probsMinCell[9]=0.07692307692;
            probsMinCell[10]=0.07692307692;
            probsMinCell[11]=0.07692307692;
            probsMinCell[12]=0.07692307692;
             
             cellFlag=0;
             cellCap=0;
             while(cellFlag==0){
                 obj->locCell=draw_multinom(&rCopy, 13, probsMinCell);
                 
                 cellCap=(*pMinCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell];
                 printf("cellCap: %d\n", cellCap);
                 
                 if(cellCap<2){
                     cellFlag=1;
                 }
                 printf("cellFlag: %d\n", cellFlag);
             }
            */
            
            cellFlag=0;
            cellCap=0;
            mCtr2=0;
            //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                    mCtr3=0;
                    while(cellFlag==0||mCtr3<6){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                        mCtr4=0;
                        while(cellFlag==0||mCtr4<13){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                            cellCap=(*pMinCellArray)[obj->zone][mCtr2][mCtr3][mCtr4];
                            printf("cellCap: %d\n", cellCap);
                            
                            if(cellCap<2){
                                cellFlag=1;
                                obj->locArea=mCtr2;
                                obj->locUnit=mCtr3;
                                obj->locCell=mCtr4;
                            }
                            printf("cellFlag: %d\n", cellFlag);
                            mCtr4++;
                        }
                        mCtr3++;
                    }
                    mCtr2++;
                }
            //}
            if(cellFlag==0){
                obj->locArea=mCtr2-1;
                obj->locUnit=mCtr3-1;
                obj->locCell=mCtr4-1;
            }
            
            
            //obj->locCell=generateAge(1, 13);
            
            //(*pMinCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
        
        }else if(obj->indivType==1){
            obj->zone=seedZone;
        }else if(obj->indivType==2){
            obj->zone=seedZone;
        }else if(obj->indivType==3){
            obj->zone=seedZone;
        }else if(obj->indivType==4){
            obj->zone=seedZone;
        }
        
    }else if(obj->location==2){
        if(obj->indivType==0){
            probsMedZone[0]=0.19;
            probsMedZone[1]=0.02;
            probsMedZone[2]=0.02;
            probsMedZone[3]=0.07;
            probsMedZone[4]=0.09;
            probsMedZone[5]=0.06;
            probsMedZone[6]=0.17;
            probsMedZone[7]=0.30;
            probsMedZone[8]=0.04;
            probsMedZone[9]=0.02;
            probsMedZone[10]=0.03;
            obj->zone=draw_multinom(&rCopy, 11, probsMedZone);
            
            //set units
            //obj->locUnit=generateAge(1, 8);
            
            /*
            probsMedArea[0]=0.5;
            probsMedArea[1]=0.5;
            obj->locArea=draw_multinom(&rCopy, 2, probsMedArea);
            //obj->locArea=generateAge(1, 2);
            
            probsMedUnit[0]=0.16666666666;
            probsMedUnit[1]=0.16666666666;
            probsMedUnit[2]=0.16666666666;
            probsMedUnit[3]=0.16666666666;
            obj->locUnit=draw_multinom(&rCopy, 4, probsMedUnit);
            //obj->locUnit=generateAge(1, 6);
            
            probsMedCell[0]=0.05263157894;
            probsMedCell[1]=0.05263157894;
            probsMedCell[2]=0.05263157894;
            probsMedCell[3]=0.05263157894;
            probsMedCell[4]=0.05263157894;
            probsMedCell[5]=0.05263157894;
            probsMedCell[6]=0.05263157894;
            probsMedCell[7]=0.05263157894;
            probsMedCell[8]=0.05263157894;
            probsMedCell[9]=0.05263157894;
            probsMedCell[10]=0.05263157894;
            probsMedCell[11]=0.05263157894;
            probsMedCell[12]=0.05263157894;
            probsMedCell[13]=0.05263157894;
            probsMedCell[14]=0.05263157894;
            probsMedCell[15]=0.05263157894;
            probsMedCell[16]=0.05263157894;
            probsMedCell[17]=0.05263157894;
            probsMedCell[18]=0.05263157894;
            
            //obj->locArea=generateAge(1, 2);
            //obj->locUnit=generateAge(1, 4);
            //obj->locCell=generateAge(1, 19);
            
            cellFlag=0;
            cellCap=0;
            while(cellFlag==0){
                obj->locCell=draw_multinom(&rCopy, 19, probsMedCell);
                
                cellCap=(*pMedCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell];
                printf("cellCap: %d\n", cellCap);
                
                if(cellCap<2){
                    cellFlag=1;
                }
                printf("cellFlag: %d\n", cellFlag);
            }
            //(*pMedCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
             */
            
            cellFlag=0;
            cellCap=0;
            mCtr2=0;
            //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                    mCtr3=0;
                    while(cellFlag==0||mCtr3<4){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                        mCtr4=0;
                        while(cellFlag==0||mCtr4<19){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                            cellCap=(*pMedCellArray)[obj->zone][mCtr2][mCtr3][mCtr4];
                            printf("cellCap: %d\n", cellCap);
                            
                            if(cellCap<2){
                                cellFlag=1;
                                obj->locArea=mCtr2;
                                obj->locUnit=mCtr3;
                                obj->locCell=mCtr4;
                            }
                            printf("cellFlag: %d\n", cellFlag);
                            mCtr4++;
                        }
                        mCtr3++;
                    }
                    mCtr2++;
                }
            
            if(cellFlag==0){
                obj->locArea=mCtr2-1;
                obj->locUnit=mCtr3-1;
                obj->locCell=mCtr4-1;
            }
            
        }else if(obj->indivType==1){
            obj->zone=seedZone;
        }else if(obj->indivType==2){
            obj->zone=seedZone;
        }else if(obj->indivType==3){
            obj->zone=seedZone;
        }else if(obj->indivType==4){
            obj->zone=seedZone;
        }
        
        //obj->zone=binZone; //put individual in an area
    }else if(obj->location==3){
        if(obj->indivType==0){
            probsMaxZone[0]=0.06;
            probsMaxZone[1]=0.01;
            probsMaxZone[2]=0.06;
            probsMaxZone[3]=0.06;
            probsMaxZone[4]=0.01;
            probsMaxZone[5]=0.05;
            probsMaxZone[6]=0.06;
            probsMaxZone[7]=0.01;
            probsMaxZone[8]=0.01;
            probsMaxZone[9]=0.05;
            probsMaxZone[10]=0.01;
            probsMaxZone[11]=0.15;
            probsMaxZone[12]=0.07;
            probsMaxZone[13]=0.13;
            probsMaxZone[14]=0.09;
            probsMaxZone[15]=0.05;
            probsMaxZone[16]=0.07;
            probsMaxZone[17]=0.07;
            obj->zone=draw_multinom(&rCopy, 18, probsMaxZone);
            
            /*
            //set units
            //obj->locUnit=generateAge(1, 14);
            probsMaxArea[0]=0.25;
            probsMaxArea[1]=0.25;
            probsMaxArea[2]=0.25;
            probsMaxArea[3]=0.25;
            obj->locArea=draw_multinom(&rCopy, 4, probsMaxArea);
            //obj->locArea=generateAge(1, 2);
            
            probsMaxUnit[0]=0.2;
            probsMaxUnit[1]=0.2;
            probsMaxUnit[2]=0.2;
            probsMaxUnit[3]=0.2;
            probsMaxUnit[4]=0.2;
            obj->locUnit=draw_multinom(&rCopy, 5, probsMaxUnit);
            //obj->locUnit=generateAge(1, 6);
            
            probsMaxCell[0]=0.05;
            probsMaxCell[1]=0.05;
            probsMaxCell[2]=0.05;
            probsMaxCell[3]=0.05;
            probsMaxCell[4]=0.05;
            probsMaxCell[5]=0.05;
            probsMaxCell[6]=0.05;
            probsMaxCell[7]=0.05;
            probsMaxCell[8]=0.05;
            probsMaxCell[9]=0.05;
            probsMaxCell[10]=0.05;
            probsMaxCell[11]=0.05;
            probsMaxCell[12]=0.05;
            probsMaxCell[13]=0.05;
            probsMaxCell[14]=0.05;
            probsMaxCell[15]=0.05;
            probsMaxCell[16]=0.05;
            probsMaxCell[17]=0.05;
            probsMaxCell[18]=0.05;
            probsMaxCell[19]=0.05;
            
            cellFlag=0;
            cellCap=0;
            while(cellFlag==0){
                obj->locCell=draw_multinom(&rCopy, 20, probsMaxCell);
                
                cellCap=(*pMaxCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell];
                printf("cellCap: %d\n", cellCap);
                
                if(cellCap<2){
                    cellFlag=1;
                }
                printf("cellFlag: %d\n", cellFlag);
            }
            //(*pMaxCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
             */
            
            cellFlag=0;
            cellCap=0;
            mCtr2=0;
            //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                while(cellFlag==0||mCtr2<4){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                    mCtr3=0;
                    while(cellFlag==0||mCtr3<5){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                        mCtr4=0;
                        while(cellFlag==0||mCtr4<20){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                            cellCap=(*pMaxCellArray)[obj->zone][mCtr2][mCtr3][mCtr4];
                            printf("cellCap: %d\n", cellCap);
                            
                            if(cellCap<2){
                                cellFlag=1;
                                obj->locArea=mCtr2;
                                obj->locUnit=mCtr3;
                                obj->locCell=mCtr4;
                            }
                            printf("cellFlag: %d\n", cellFlag);
                            mCtr4++;
                        }
                        mCtr3++;
                    }
                    mCtr2++;
                }
            
            if(cellFlag==0){
                obj->locArea=mCtr2-1;
                obj->locUnit=mCtr3-1;
                obj->locCell=mCtr4-1;
            }
            
        }else if(obj->indivType==1){
            obj->zone=seedZone;
        }else if(obj->indivType==2){
            obj->zone=seedZone;
        }else if(obj->indivType==3){
            obj->zone=seedZone;
        }else if(obj->indivType==4){
            obj->zone=seedZone;
        }
        
    }
    
    if(prison==99){ // new entry, set zone to max reception
        probsMaxZone[0]=0.011;
        probsMaxZone[1]=0.0;
        probsMaxZone[2]=0.0;
        probsMaxZone[3]=0.0;
        probsMaxZone[4]=0.001;
        probsMaxZone[5]=0.002;
        probsMaxZone[6]=0.0;
        probsMaxZone[7]=0.0;
        probsMaxZone[8]=0.0;
        probsMaxZone[9]=0.0;
        probsMaxZone[10]=0.0;
        probsMaxZone[11]=0.286;
        probsMaxZone[12]=0.029;
        probsMaxZone[13]=0.211;
        probsMaxZone[14]=0.076;
        probsMaxZone[15]=0.103;
        probsMaxZone[16]=0.009;
        probsMaxZone[17]=0.028;
        
        /*
        //set units
        //obj->locUnit=generateAge(1, 13);
        probsMaxArea[0]=0.25;
        probsMaxArea[1]=0.25;
        probsMaxArea[2]=0.25;
        probsMaxArea[3]=0.25;
        obj->locArea=draw_multinom(&rCopy, 4, probsMaxArea);
        //obj->locArea=generateAge(1, 2);
        
        probsMaxUnit[0]=0.2;
        probsMaxUnit[1]=0.2;
        probsMaxUnit[2]=0.2;
        probsMaxUnit[3]=0.2;
        probsMaxUnit[4]=0.2;
        obj->locUnit=draw_multinom(&rCopy, 5, probsMaxUnit);
        //obj->locUnit=generateAge(1, 6);
        
        probsMaxCell[0]=0.05;
        probsMaxCell[1]=0.05;
        probsMaxCell[2]=0.05;
        probsMaxCell[3]=0.05;
        probsMaxCell[4]=0.05;
        probsMaxCell[5]=0.05;
        probsMaxCell[6]=0.05;
        probsMaxCell[7]=0.05;
        probsMaxCell[8]=0.05;
        probsMaxCell[9]=0.05;
        probsMaxCell[10]=0.05;
        probsMaxCell[11]=0.05;
        probsMaxCell[12]=0.05;
        probsMaxCell[13]=0.05;
        probsMaxCell[14]=0.05;
        probsMaxCell[15]=0.05;
        probsMaxCell[16]=0.05;
        probsMaxCell[17]=0.05;
        probsMaxCell[18]=0.05;
        probsMaxCell[19]=0.05;
        //obj->locCell=draw_multinom(&rCopy, 20, probsMaxCell);
        
        cellFlag=0;
        cellCap=0;
        while(cellFlag==0){
            obj->locCell=draw_multinom(&rCopy, 20, probsMaxCell);
            
            cellCap=(*pMaxCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell];
            printf("cellCap: %d\n", cellCap);
            
            if(cellCap<2){
                cellFlag=1;
            }
            printf("cellFlag: %d\n", cellFlag);
        }
         */
        
        printf("set reception prison\n");
        
        normSum=0.0;
        for(normI=0;normI<18;normI++){
            normSum+=probsMaxZone[normI];
            printf("normSum %d!!!\n", normSum);
        }
        
        for(normI=0;normI<18;normI++){
            probsMaxZone[normI]=probsMaxZone[normI]/normSum;
            printf("max zone %d prob: %f!!!\n", normI, probsMaxZone[normI]);
        }
        
        obj->zone=draw_multinom(&rCopy, 18, probsMaxZone);
        
        if(seedPop==0){
            obj->zone=11;
        }
        
        cellFlag=0;
        cellCap=0;
        mCtr2=0;
        //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
            while(cellFlag==0||mCtr2<4){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                mCtr3=0;
                while(cellFlag==0||mCtr3<5){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                    mCtr4=0;
                    while(cellFlag==0||mCtr4<20){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                        cellCap=(*pMaxCellArray)[obj->zone][mCtr2][mCtr3][mCtr4];
                        printf("cellCap: %d\n", cellCap);
                        
                        //if((*pIsolateMaxCell)[obj->zone][mCtr2][mCtr3][mCtr4]==2){
                            if(cellCap<2){ //new cells only
                                cellFlag=1;
                                obj->locArea=mCtr2;
                                obj->locUnit=mCtr3;
                                obj->locCell=mCtr4;
                            }
                            printf("cellFlag: %d\n", cellFlag);
                        //}
                        mCtr4++;
                    }
                    mCtr3++;
                }
                mCtr2++;
            }
        
        if(cellFlag==0){
            obj->locArea=mCtr2-1;
            obj->locUnit=mCtr3-1;
            obj->locCell=mCtr4-1;
            *povercap++;
        }
        
        //(*pMaxCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
        
        printf("reception prison drawn!!!\n");
    }
    
    obj->timeOfInfection=0;//set to 0 everytime cleared.
    obj->timeOfProgression=0;
    obj->placeInfected=0; //NA for now
    obj->isolated=0;
    obj->timeIsolated=0;
    
    /*if(obj->indivType==0){
        if(prison==99){
            obj->isolated=2; // 1 for 1 out isolation, 2 for 2 out
            obj->timeIsolated=currDay;
            (*pIsolateMaxCell)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]=2;
        }else{
            obj->isolated=0;
            obj->timeIsolated=0;
        }
    }*/
    
    obj->hospitalised=0;
    obj->court=0;
    obj->courtNumber=99;
    obj->moving=0;
    obj->truckNumber=99;
    
    obj->ppe=1; //apply face mask for all
    obj->ppeprotection=10.70; //reduction by 10.70%
    
    obj->PCRtested=0;
    obj->PCRdate=0;
    obj->PCRresult=0; //reduction by 10.70%
    
    obj->lockdown=0;
    obj->timeLockdown=0;
    
    printf("individual type set\n");
    
    //Set Head and Tails
    if(*pHeadCopy==NULL){
        obj->prevIndiv=NULL;
        obj->nextIndiv=NULL;
        *pHeadCopy=obj;
    }else{
        (*pTailCopy)->nextIndiv=obj;
        obj->prevIndiv=*pTailCopy;
        obj->nextIndiv=NULL;
    }
    *pTailCopy=obj;
    
    if(obj->indivType==0){
        if(obj->location==1){
            (*pMinCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
        }else if(obj->location==2){
            (*pMedCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
        }else if(obj->location==3){
            (*pMaxCellArray)[obj->zone][obj->locArea][obj->locUnit][obj->locCell]++;
        }
    }
    
    printf("individual initialised\n");
}

void traverse (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    while(current!=NULL){
        //printf("traversing: ID: %d, loc: %d, group: %d\n",current->ID, current->location, current->group);
        //printf("traversing: ID: %d, Age: %d, ATSI: %d, loc: %d, metavir: %d, injFreq: %d, sharing: %d, shaFreq: %d, injOpd: %d, group: %d\n",current->ID, current->age, current->atsi, current->location, current->metavir, current->injFreq, current->sharing, current->shaFreq, current->injOpd, current->group);
        current=current->nextIndiv;
    }
}

void staffUnisolate (sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay){
    sIndiv *current;
    current=*pHeadCopy;
    int testDate=current->testDate;
    
    printf("unisolating staff\n");
    while(current!=NULL){
        //printf("traversing: ID: %d, loc: %d, group: %d\n",current->ID, current->location, current->group);
        //printf("traversing: ID: %d, Age: %d, ATSI: %d, loc: %d, metavir: %d, injFreq: %d, sharing: %d, shaFreq: %d, injOpd: %d, group: %d\n",current->ID, current->age, current->atsi, current->location, current->metavir, current->injFreq, current->sharing, current->shaFreq, current->injOpd, current->group);
        if(current->indivType==1||current->indivType==2){ //pstaff
            if(current->isolated==1){
                if(currDay>(testDate)){
                    if(current->COVID==0){
                        current->isolated=0;
                    }
                }
            }
        }
        
        current=current->nextIndiv;
    }
}

int infectRedZone (gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*zoneMinCIArray)[30][7], int (*zoneMedCIArray)[30][7], int (*zoneMaxCIArray)[30][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int currDay){
    
    printf("infecting in red zone\n");
    sIndiv *target, *tempCurrent;
    target=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rStatus, binInf;
    //int zone=target->zone;
    int targetLocation=target->location, aCtr, iTarget, currLoc, currGroup, newGroup=99, eFlag=0;
    //int totalAtRisk=(*pLocArray2)[target->location][3]+(*pLocArray2)[target->location][5]+(*pLocArray2)[target->location][9]+(*pLocArray2)[target->location][11]; //add all at risk except those with immunity
    //sIndiv *refInfectList[totalAtRisk];
    //double probInfectList[totalAtRisk];
    double eventDraw[2]; //proability of recipient to get infected
    double probInf[2]; //probability of source infecting recipient using per event probability
    int success, eventDecision, susCtr;
    int targetArrayLoc;
    int contactIMax=0, contactIMed=0, contactIMin=0, contactPSMax=0, contactPSMed=0, contactPSMin=0, contactHSMax=0, contactHSMed=0, contactHSMin=0, contactEVMax=0, contactEVMed=0, contactEVMin=0, contactFVMax=0, contactFVMed=0, contactFVMin=0;
    int contactCapacityI=0, contactCapacityPS=0, contactCapacityHS=0, contactCapacityEV=0, contactCapacityFV=0;
    int targetAge, inmateFlag=0;
    int infectedCountI=0, infectedCountPS=0, infectedCountHS=0, infectedCountEV=0, infectedCountFV=0;
    double probsRisk[2];
    int isolateStatus;
    double ppeEffect;
    
    printf("going through list of individuals\n");
    success=0;
    //probInf[0]=0.0; probInf[1]=0.0;
    eventDraw[0]=0; eventDraw[1]=0;
    //aCtr=0;
    //Go through list of individuals
    
    if(target->indivType==0){ //Inmate
        contactIMax=generateAge(175, 388)+generateAge(30, 76);//gsl_ran_flat(*r, 36, 81);//Interaction within Units/Pods only
        contactIMed=generateAge(107, 238)+generateAge(14, 35);//gsl_ran_flat(*r, 107, 130);//Interaction within Units/Pods
        contactIMin=generateAge(88, 195)+generateAge(15, 38);//gsl_ran_flat(*r, 87, 107);//Interaction within Units/Pods
        contactPSMax=2;
        contactPSMed=2;
        contactPSMin=2;
        probsRisk[0]=1-0.17; //not a patient
        probsRisk[1]=0.17;
        contactHSMax=draw_multinom(&rCopy, 2, probsRisk);
        probsRisk[0]=1-0.07; //not a patient
        probsRisk[1]=0.07;
        contactHSMed=draw_multinom(&rCopy, 2, probsRisk);
        probsRisk[0]=1-0.09; //not a patient
        probsRisk[1]=0.09;
        contactHSMin=draw_multinom(&rCopy, 2, probsRisk);
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=1;
        }else if(targetLocation==2){
            targetArrayLoc=6;
        }else if(targetLocation==3){
            targetArrayLoc=11;
        }
        inmateFlag=1;
    }else if(target->indivType==1){//ps
        contactIMax=5;
        contactIMed=5;
        contactIMin=5;
        contactPSMax=39;
        contactPSMed=14;
        contactPSMin=12;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=2;
        }else if(targetLocation==2){
            targetArrayLoc=7;
        }else if(targetLocation==3){
            targetArrayLoc=12;
        }
    }else if(target->indivType==2){ //hs
        contactIMax=8; //8
        contactIMed=6; //6
        contactIMin=7; //7
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=9; //9
        contactHSMed=7; //7
        contactHSMin=5; //5
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=3;
        }else if(targetLocation==2){
            targetArrayLoc=8;
        }else if(targetLocation==3){
            targetArrayLoc=13;
        }
    }else if(target->indivType==3){ //ev
        contactIMax=1;
        contactIMed=1;
        contactIMin=1;
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=4;
        }else if(targetLocation==2){
            targetArrayLoc=9;
        }else if(targetLocation==3){
            targetArrayLoc=14;
        }
    }else if(target->indivType==4){ //fv
        contactIMax=1;
        contactIMed=1;
        contactIMin=1;
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=5;
        }else if(targetLocation==2){
            targetArrayLoc=10;
        }else if(targetLocation==3){
            targetArrayLoc=15;
        }
    }
    
    if(targetLocation==1){
        contactCapacityI=contactIMin;
        contactCapacityPS=contactPSMin;
        contactCapacityHS=contactHSMin;
        contactCapacityEV=contactEVMin;
        contactCapacityFV=contactFVMin;
    }else if(targetLocation==2){
        contactCapacityI=contactIMed;
        contactCapacityPS=contactPSMed;
        contactCapacityHS=contactHSMed;
        contactCapacityEV=contactEVMed;
        contactCapacityFV=contactFVMed;
    }else if(targetLocation==3){
        contactCapacityI=contactIMax;
        contactCapacityPS=contactPSMax;
        contactCapacityHS=contactHSMax;
        contactCapacityEV=contactEVMax;
        contactCapacityFV=contactFVMax;
    }
    
    tempCurrent=*pHeadCopy; //Target points to the individual at the beginning of the list
    //printf("total at-risk agents in prison %d: %d\n", targetLocation, totalAtRisk);
    
    while(tempCurrent!=NULL){
        //add if injecting
        //FIX location matching to infect other populations
        //
        isolateStatus=tempCurrent->isolated;
        
        if((targetLocation==tempCurrent->location)&&(tempCurrent!=target)){//susceptible if same location and not same indiv



            if(tempCurrent->zone==target->zone){  //tempCurrent->zone==target->zone
                if(tempCurrent->hospitalised==0&&tempCurrent->isolated==0&&tempCurrent->court==0){ //not in any hospital
                    //if(tempCurrent->COVID==0){ //not yet been infected and not currently infected
                    if(tempCurrent->indivType==0&&isolateStatus!=0&&tempCurrent->isolated==2){ //if recepient is inmate
                        if(infectedCountI<contactCapacityI){
                            if(target->indivType==0){ //if source is inmate
                                if(tempCurrent->locArea==target->locArea&&tempCurrent->locUnit==target->locUnit&&tempCurrent->locCell==target->locCell){
                                    //herd immunity
                                    if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                        probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[0]=1-probInf[1];
                                    }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                        ppeEffect=0.67; //N95
                                        if(target->severity==4||target->severity==5){ //mod or severe
                                            probInf[1]=0.01*2;//0.001;//0.05;//target = probability of getting infected with COVID-19
                                            probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                        }else{ //mild and asymptomatic
                                            probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//0.05;//probability of getting infected with COVID-19
                                            probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                        }
                                        probInf[0]=1-probInf[1];
                                    }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                        probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[0]=1-probInf[1];
                                    }
                                    
                                    printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                                    
                                    binInf=draw_multinom(&rCopy, 2, probInf);
                                    //printf("Infection commencing\n");
                                    printf("Infect Y or N: %d\n", binInf);
                                    
                                    if(binInf==1){ //Target is to be infected
                                        printf("Infection commencing\n");
                                        
                                        tempCurrent->COVID=1; // infected!!
                                        tempCurrent->severity=0;
                                        currLoc=tempCurrent->location;
                                        tempCurrent->placeInfected=currLoc;
                                        targetAge=tempCurrent->age;
                                        
                                        if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                            (*newInfectedInmatesMin)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=1;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMinCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMinCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMinCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMinCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMinCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMinCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMinCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=2;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=3;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=4;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=5;
                                            }
                                        }else if(currLoc==2){
                                            (*newInfectedInmatesMed)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=6;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMedCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMedCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMedCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMedCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMedCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMedCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMedCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=7;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=8;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=9;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=10;
                                            }
                                        }else if(currLoc==3){
                                            (*newInfectedInmatesMax)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=11;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=12;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=13;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=14;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=15;
                                            }
                                        }
                                        
                                        if(tempCurrent->age<=19){ //means infected in a prison
                                            (*exp19)++;
                                        }else if(tempCurrent->age>=20&&tempCurrent->age<=44){
                                            (*exp20to44)++;
                                        }else if(tempCurrent->age>=45&&tempCurrent->age<=54){
                                            (*exp45to54)++;
                                        }else if(tempCurrent->age>=55&&tempCurrent->age<=64){
                                            (*exp55to64)++;
                                        }else if(tempCurrent->age>=65&&tempCurrent->age<=74){
                                            (*exp65to74)++;
                                        }else if(tempCurrent->age>=75&&tempCurrent->age<=84){
                                            (*exp75to84)++;
                                        }else if(tempCurrent->age>=85){
                                            (*exp85)++;
                                        }
                                        
                                        tempCurrent->timeOfInfection=currDay; //Record date of infection
                                        ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                        currGroup=tempCurrent->group;
                                        
                                        (*pLocArray2)[currLoc][currGroup]--;
                                        
                                        //system("pause");
                                        switch(currGroup){
                                            case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                                break;
                                            case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                                break;
                                            case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                                break;
                                            case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=0;
                                                break;
                                            case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                                break;
                                            case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=1;
                                                break;
                                            case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                                break;
                                            case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                                break;
                                            case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                                break;
                                            case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=6;
                                                break;
                                            case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                                break;
                                            case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=7;
                                                break;
                                                
                                        }
                                        
                                        //if(currGroup!=newGroup){
                                        tempCurrent->group=newGroup;
                                        //}
                                        
                                        printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                        success++;
                                        
                                        (*pLocArray2)[currLoc][newGroup]++;
                                        printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                        
                                        (*newInfectedInmates)++;
                                    }
                                    
                                    infectedCountI++;
                                }
                            }else{ //if source not inmate
                                if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                    probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[0]=1-probInf[1];
                                }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                    ppeEffect=0.67;
                                    if(target->severity==4||target->severity==5){
                                        probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                    }else{
                                        probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                    }
                                    probInf[0]=1-probInf[1];
                                }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                    probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[0]=1-probInf[1];
                                }
                                
                                printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                                
                                binInf=draw_multinom(&rCopy, 2, probInf);
                                //printf("Infection commencing\n");
                                printf("Infect Y or N: %d\n", binInf);
                                
                                if(binInf==1){ //Target is to be infected
                                    printf("Infection commencing\n");
                                    
                                    tempCurrent->COVID=1; // infected!!
                                    tempCurrent->severity=0;
                                    currLoc=tempCurrent->location;
                                    tempCurrent->placeInfected=currLoc;
                                    targetAge=tempCurrent->age;
                                    
                                    if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                        (*newInfectedInmatesMin)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=1;
                                            
                                            if(targetAge<=19){
                                                (*zoneMinCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMinCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMinCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMinCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMinCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMinCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMinCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=2;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=3;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=4;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=5;
                                        }
                                    }else if(currLoc==2){
                                        (*newInfectedInmatesMed)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=6;
                                            
                                            if(targetAge<=19){
                                                (*zoneMedCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMedCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMedCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMedCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMedCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMedCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMedCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=7;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=8;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=9;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=10;
                                        }
                                    }else if(currLoc==3){
                                        (*newInfectedInmatesMax)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=11;
                                            
                                            if(targetAge<=19){
                                                (*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=12;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=13;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=14;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=15;
                                        }
                                    }
                                    
                                    tempCurrent->timeOfInfection=currDay; //Record date of infection
                                    ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                    currGroup=tempCurrent->group;
                                    
                                    (*pLocArray2)[currLoc][currGroup]--;
                                    
                                    //system("pause");
                                    switch(currGroup){
                                        case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                            break;
                                        case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                            break;
                                        case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                            break;
                                        case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=0;
                                            break;
                                        case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                            break;
                                        case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=1;
                                            break;
                                        case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                            break;
                                        case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                            break;
                                        case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                            break;
                                        case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=6;
                                            break;
                                        case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                            break;
                                        case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=7;
                                            break;
                                            
                                    }
                                    
                                    //if(currGroup!=newGroup){
                                    tempCurrent->group=newGroup;
                                    //}
                                    
                                    printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                    success++;
                                    
                                    (*pLocArray2)[currLoc][newGroup]++;
                                    printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                    
                                    (*newInfectedInmates)++;
                                }
                                
                                infectedCountI++;
                            }
                                

                        }
                    }else if(tempCurrent->indivType==1){
                        if(infectedCountPS<contactCapacityPS){
                            
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedPSMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedPSMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedPSMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                
                                (*newInfectedPS)++;
                            }
                            
                            infectedCountPS++;
                        }
                    }else if(tempCurrent->indivType==2){
                        if(infectedCountHS<contactCapacityHS){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedHSMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedHSMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedHSMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedHS)++;
                            }
                            
                            infectedCountHS++;
                        }
                    }else if(tempCurrent->indivType==3){
                        if(infectedCountEV<contactCapacityEV){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedEVMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedEVMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedEVMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedEV)++;
                            }
                            
                            infectedCountEV++;
                        }
                    }else if(tempCurrent->indivType==4){
                        if(infectedCountFV<contactCapacityFV){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedFVMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedFVMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedFVMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedFV)++;
                            }
                            
                            infectedCountFV++;
                        }
                    }
                    
                    //}//covid
                }//hospitalised
            }//same zone
        }//same location
        tempCurrent=tempCurrent->nextIndiv;
    }
    printf("Number of new infections for this instance: %d", success);
    return success;
}

int StaffInfectRedZone (gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*zoneMinCIArray)[30][7], int (*zoneMedCIArray)[30][7], int (*zoneMaxCIArray)[30][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int currDay){
    
    printf("infecting in red zone\n");
    sIndiv *target, *tempCurrent;
    target=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rStatus, binInf;
    //int zone=target->zone;
    int targetLocation=target->location, aCtr, iTarget, currLoc, currGroup, newGroup=99, eFlag=0;
    //int totalAtRisk=(*pLocArray2)[target->location][3]+(*pLocArray2)[target->location][5]+(*pLocArray2)[target->location][9]+(*pLocArray2)[target->location][11]; //add all at risk except those with immunity
    //sIndiv *refInfectList[totalAtRisk];
    //double probInfectList[totalAtRisk];
    double eventDraw[2]; //proability of recipient to get infected
    double probInf[2]; //probability of source infecting recipient using per event probability
    int success, eventDecision, susCtr;
    int targetArrayLoc;
    int contactIMax=0, contactIMed=0, contactIMin=0, contactPSMax=0, contactPSMed=0, contactPSMin=0, contactHSMax=0, contactHSMed=0, contactHSMin=0, contactEVMax=0, contactEVMed=0, contactEVMin=0, contactFVMax=0, contactFVMed=0, contactFVMin=0;
    int contactCapacityI=0, contactCapacityPS=0, contactCapacityHS=0, contactCapacityEV=0, contactCapacityFV=0;
    int targetAge, inmateFlag=0;
    int infectedCountI=0, infectedCountPS=0, infectedCountHS=0, infectedCountEV=0, infectedCountFV=0;
    double probsRisk[2];
    int isolateStatus;
    double ppeEffect;
    
    printf("going through list of individuals\n");
    success=0;
    //probInf[0]=0.0; probInf[1]=0.0;
    eventDraw[0]=0; eventDraw[1]=0;
    //aCtr=0;
    //Go through list of individuals
    
    if(target->indivType==0){ //Inmate
        contactIMax=generateAge(175, 388)+generateAge(30, 76);//gsl_ran_flat(*r, 36, 81);//Interaction within Units/Pods only
        contactIMed=generateAge(107, 238)+generateAge(14, 35);//gsl_ran_flat(*r, 107, 130);//Interaction within Units/Pods
        contactIMin=generateAge(88, 195)+generateAge(15, 38);//gsl_ran_flat(*r, 87, 107);//Interaction within Units/Pods
        contactPSMax=2;
        contactPSMed=2;
        contactPSMin=2;
        probsRisk[0]=1-0.17; //not a patient
        probsRisk[1]=0.17;
        contactHSMax=draw_multinom(&rCopy, 2, probsRisk);
        probsRisk[0]=1-0.07; //not a patient
        probsRisk[1]=0.07;
        contactHSMed=draw_multinom(&rCopy, 2, probsRisk);
        probsRisk[0]=1-0.09; //not a patient
        probsRisk[1]=0.09;
        contactHSMin=draw_multinom(&rCopy, 2, probsRisk);
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=1;
        }else if(targetLocation==2){
            targetArrayLoc=6;
        }else if(targetLocation==3){
            targetArrayLoc=11;
        }
        inmateFlag=1;
    }else if(target->indivType==1){//ps
        contactIMax=5;
        contactIMed=5;
        contactIMin=5;
        contactPSMax=1;//39;
        contactPSMed=1;//14;
        contactPSMin=1;//12;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=2;
        }else if(targetLocation==2){
            targetArrayLoc=7;
        }else if(targetLocation==3){
            targetArrayLoc=12;
        }
    }else if(target->indivType==2){ //hs
        contactIMax=8; //8
        contactIMed=6; //6
        contactIMin=7; //7
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=1;//9; //9
        contactHSMed=1;//7; //7
        contactHSMin=1;//5; //5
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=3;
        }else if(targetLocation==2){
            targetArrayLoc=8;
        }else if(targetLocation==3){
            targetArrayLoc=13;
        }
    }else if(target->indivType==3){ //ev
        contactIMax=1;
        contactIMed=1;
        contactIMin=1;
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=4;
        }else if(targetLocation==2){
            targetArrayLoc=9;
        }else if(targetLocation==3){
            targetArrayLoc=14;
        }
    }else if(target->indivType==4){ //fv
        contactIMax=1;
        contactIMed=1;
        contactIMin=1;
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=5;
        }else if(targetLocation==2){
            targetArrayLoc=10;
        }else if(targetLocation==3){
            targetArrayLoc=15;
        }
    }
    
    if(targetLocation==1){
        contactCapacityI=contactIMin;
        contactCapacityPS=contactPSMin;
        contactCapacityHS=contactHSMin;
        contactCapacityEV=contactEVMin;
        contactCapacityFV=contactFVMin;
    }else if(targetLocation==2){
        contactCapacityI=contactIMed;
        contactCapacityPS=contactPSMed;
        contactCapacityHS=contactHSMed;
        contactCapacityEV=contactEVMed;
        contactCapacityFV=contactFVMed;
    }else if(targetLocation==3){
        contactCapacityI=contactIMax;
        contactCapacityPS=contactPSMax;
        contactCapacityHS=contactHSMax;
        contactCapacityEV=contactEVMax;
        contactCapacityFV=contactFVMax;
    }
    
    tempCurrent=*pHeadCopy; //Target points to the individual at the beginning of the list
    //printf("total at-risk agents in prison %d: %d\n", targetLocation, totalAtRisk);
    
    while(tempCurrent!=NULL){
        //add if injecting
        //FIX location matching to infect other populations
        //
        isolateStatus=tempCurrent->isolated;
        
        if((targetLocation==tempCurrent->location)&&(tempCurrent!=target)){//susceptible if same location and not same indiv



            if(tempCurrent->zone==target->zone){  //tempCurrent->zone==target->zone
                if(tempCurrent->hospitalised==0&&tempCurrent->isolated==0&&tempCurrent->court==0){ //not in any hospital
                    //if(tempCurrent->COVID==0){ //not yet been infected and not currently infected
                    if(tempCurrent->indivType==0&&isolateStatus!=0&&tempCurrent->isolated==2){ //if recepient is inmate
                        if(infectedCountI<contactCapacityI){
                            if(target->indivType==0){ //if source is inmate
                                if(tempCurrent->locArea==target->locArea&&tempCurrent->locUnit==target->locUnit&&tempCurrent->locCell==target->locCell){
                                    //herd immunity
                                    if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                        probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[0]=1-probInf[1];
                                    }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                        ppeEffect=0.85;
                                        if(target->severity==4||target->severity==5){ //mod or severe
                                            probInf[1]=0.01*2;//0.001;//0.05;//target = probability of getting infected with COVID-19
                                            probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                        }else{ //mild and asymptomatic
                                            probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//0.05;//probability of getting infected with COVID-19
                                            probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                        }
                                        probInf[0]=1-probInf[1];
                                    }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                        probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[0]=1-probInf[1];
                                    }
                                    
                                    printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                                    
                                    binInf=draw_multinom(&rCopy, 2, probInf);
                                    //printf("Infection commencing\n");
                                    printf("Infect Y or N: %d\n", binInf);
                                    
                                    if(binInf==1){ //Target is to be infected
                                        printf("Infection commencing\n");
                                        
                                        tempCurrent->COVID=1; // infected!!
                                        tempCurrent->severity=0;
                                        currLoc=tempCurrent->location;
                                        tempCurrent->placeInfected=currLoc;
                                        targetAge=tempCurrent->age;
                                        
                                        if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                            (*newInfectedInmatesMin)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=1;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMinCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMinCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMinCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMinCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMinCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMinCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMinCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=2;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=3;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=4;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=5;
                                            }
                                        }else if(currLoc==2){
                                            (*newInfectedInmatesMed)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=6;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMedCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMedCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMedCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMedCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMedCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMedCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMedCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=7;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=8;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=9;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=10;
                                            }
                                        }else if(currLoc==3){
                                            (*newInfectedInmatesMax)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=11;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=12;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=13;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=14;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=15;
                                            }
                                        }
                                        
                                        if(tempCurrent->age<=19){ //means infected in a prison
                                            (*exp19)++;
                                        }else if(tempCurrent->age>=20&&tempCurrent->age<=44){
                                            (*exp20to44)++;
                                        }else if(tempCurrent->age>=45&&tempCurrent->age<=54){
                                            (*exp45to54)++;
                                        }else if(tempCurrent->age>=55&&tempCurrent->age<=64){
                                            (*exp55to64)++;
                                        }else if(tempCurrent->age>=65&&tempCurrent->age<=74){
                                            (*exp65to74)++;
                                        }else if(tempCurrent->age>=75&&tempCurrent->age<=84){
                                            (*exp75to84)++;
                                        }else if(tempCurrent->age>=85){
                                            (*exp85)++;
                                        }
                                        
                                        tempCurrent->timeOfInfection=currDay; //Record date of infection
                                        ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                        currGroup=tempCurrent->group;
                                        
                                        (*pLocArray2)[currLoc][currGroup]--;
                                        
                                        //system("pause");
                                        switch(currGroup){
                                            case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                                break;
                                            case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                                break;
                                            case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                                break;
                                            case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=0;
                                                break;
                                            case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                                break;
                                            case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=1;
                                                break;
                                            case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                                break;
                                            case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                                break;
                                            case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                                break;
                                            case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=6;
                                                break;
                                            case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                                break;
                                            case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=7;
                                                break;
                                                
                                        }
                                        
                                        //if(currGroup!=newGroup){
                                        tempCurrent->group=newGroup;
                                        //}
                                        
                                        printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                        success++;
                                        
                                        (*pLocArray2)[currLoc][newGroup]++;
                                        printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                        
                                        (*newInfectedInmates)++;
                                    }
                                    
                                    infectedCountI++;
                                }
                            }else{ //if source not inmate
                                if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                    probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[0]=1-probInf[1];
                                }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                    ppeEffect=0.85;
                                    if(target->severity==4||target->severity==5){
                                        probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                    }else{
                                        probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                    }
                                    probInf[0]=1-probInf[1];
                                }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                    probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[0]=1-probInf[1];
                                }
                                
                                printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                                
                                binInf=draw_multinom(&rCopy, 2, probInf);
                                //printf("Infection commencing\n");
                                printf("Infect Y or N: %d\n", binInf);
                                
                                if(binInf==1){ //Target is to be infected
                                    printf("Infection commencing\n");
                                    
                                    tempCurrent->COVID=1; // infected!!
                                    tempCurrent->severity=0;
                                    currLoc=tempCurrent->location;
                                    tempCurrent->placeInfected=currLoc;
                                    targetAge=tempCurrent->age;
                                    
                                    if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                        (*newInfectedInmatesMin)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=1;
                                            
                                            if(targetAge<=19){
                                                (*zoneMinCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMinCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMinCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMinCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMinCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMinCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMinCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=2;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=3;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=4;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=5;
                                        }
                                    }else if(currLoc==2){
                                        (*newInfectedInmatesMed)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=6;
                                            
                                            if(targetAge<=19){
                                                (*zoneMedCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMedCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMedCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMedCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMedCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMedCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMedCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=7;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=8;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=9;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=10;
                                        }
                                    }else if(currLoc==3){
                                        (*newInfectedInmatesMax)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=11;
                                            
                                            if(targetAge<=19){
                                                (*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=12;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=13;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=14;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=15;
                                        }
                                    }
                                    
                                    tempCurrent->timeOfInfection=currDay; //Record date of infection
                                    ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                    currGroup=tempCurrent->group;
                                    
                                    (*pLocArray2)[currLoc][currGroup]--;
                                    
                                    //system("pause");
                                    switch(currGroup){
                                        case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                            break;
                                        case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                            break;
                                        case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                            break;
                                        case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=0;
                                            break;
                                        case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                            break;
                                        case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=1;
                                            break;
                                        case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                            break;
                                        case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                            break;
                                        case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                            break;
                                        case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=6;
                                            break;
                                        case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                            break;
                                        case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=7;
                                            break;
                                            
                                    }
                                    
                                    //if(currGroup!=newGroup){
                                    tempCurrent->group=newGroup;
                                    //}
                                    
                                    printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                    success++;
                                    
                                    (*pLocArray2)[currLoc][newGroup]++;
                                    printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                    
                                    (*newInfectedInmates)++;
                                }
                                
                                infectedCountI++;
                            }
                                

                        }
                    }else if(tempCurrent->indivType==1){
                        if(infectedCountPS<contactCapacityPS){
                            
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedPSMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedPSMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedPSMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                
                                (*newInfectedPS)++;
                            }
                            
                            infectedCountPS++;
                        }
                    }else if(tempCurrent->indivType==2){
                        if(infectedCountHS<contactCapacityHS){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedHSMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedHSMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedHSMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedHS)++;
                            }
                            
                            infectedCountHS++;
                        }
                    }else if(tempCurrent->indivType==3){
                        if(infectedCountEV<contactCapacityEV){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedEVMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedEVMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedEVMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedEV)++;
                            }
                            
                            infectedCountEV++;
                        }
                    }else if(tempCurrent->indivType==4){
                        if(infectedCountFV<contactCapacityFV){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.85;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedFVMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedFVMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedFVMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedFV)++;
                            }
                            
                            infectedCountFV++;
                        }
                    }
                    
                    //}//covid
                }//hospitalised
            }//same zone
        }//same location
        tempCurrent=tempCurrent->nextIndiv;
    }
    printf("Number of new infections for this instance: %d", success);
    return success;
}

void checkIsolate (sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *povercap){

    printf("updating inmates in isolation\n");
    //if COVID++ keep in isolation
    sIndiv *current;
    current=*pHeadCopy;
    int cellFlag=0;
    int cellCap=0;
    int deltaTime=0;
    int mCtr2, mCtr3, mCtr4;
    int prevArea, prevUnit, prevCell;
    int mUnits, mAreas, mCells; //int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20]);
    
    while(current!=NULL){
        printf("going through inmates\n");
        if(current->indivType==0){ //inmate
            if(current->isolated==1){ //isolated 1out
                if(current->COVID==0){
                    //if no more covid, release
                    
                    if(current->location==1){
                        (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]--;
                        
                        if((*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]==0){
                            (*pIsolateMinCell)[current->zone][current->locArea][current->locUnit][current->locCell]=0;
                            //if room is empty, set isolate array == 0
                        }
                        
                        ////move to non-isolated cell
                        cellFlag=0;
                        cellCap=0;
                        mCtr2=0;
                        //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                            while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                                mCtr3=0;
                                while(cellFlag==0||mCtr3<6){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                    mCtr4=0;
                                    while(cellFlag==0||mCtr4<13){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                        cellCap=(*pMinCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                        printf("cellCap: %d\n", cellCap);
                                        
                                        //update count
                                        //if((*pIsolateMaxCell)[obj->zone][mCtr2][mCtr3][mCtr4]==2){
                                            if(cellCap<2 && (*pIsolateMinCell)[current->zone][mCtr2][mCtr3][mCtr4]==0){ //cells that are not full and not a red zone
                                                
                                                cellFlag=1;
                                                current->locArea=mCtr2;
                                                current->locUnit=mCtr3;
                                                current->locCell=mCtr4;
                                            }
                                            printf("cellFlag: %d\n", cellFlag);
                                        //}
                                        mCtr4++;
                                    }
                                    mCtr3++;
                                }
                                mCtr2++;
                            }
                        
                        if(cellFlag==0){
                            current->locArea=mCtr2-1;
                            current->locUnit=mCtr3-1;
                            current->locCell=mCtr4-1;
                            *povercap++;
                        }
                        
                        current->isolated==0;
                        current->timeIsolated==0;
                        
                        (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        
                    }else if(current->location==2){
                        (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]--;
                        
                        if((*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]==0){
                            (*pIsolateMedCell)[current->zone][current->locArea][current->locUnit][current->locCell]=0;
                            //if room is empty, set isolate array == 0
                        }
                        
                        ////move to non-isolated cell
                        cellFlag=0;
                        cellCap=0;
                        mCtr2=0;
                        //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                            while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                                mCtr3=0;
                                while(cellFlag==0||mCtr3<4){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                    mCtr4=0;
                                    while(cellFlag==0||mCtr4<19){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                        cellCap=(*pMedCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                        printf("cellCap: %d\n", cellCap);
                                        
                                        //update count
                                        //if((*pIsolateMaxCell)[obj->zone][mCtr2][mCtr3][mCtr4]==2){
                                            if(cellCap<2 && (*pIsolateMedCell)[current->zone][mCtr2][mCtr3][mCtr4]==0){ //cells that are not full and not a red zone
                                                
                                                cellFlag=1;
                                                current->locArea=mCtr2;
                                                current->locUnit=mCtr3;
                                                current->locCell=mCtr4;
                                            }
                                            printf("cellFlag: %d\n", cellFlag);
                                        //}
                                        mCtr4++;
                                    }
                                    mCtr3++;
                                }
                                mCtr2++;
                            }
                        
                        if(cellFlag==0){
                            current->locArea=mCtr2-1;
                            current->locUnit=mCtr3-1;
                            current->locCell=mCtr4-1;
                            *povercap++;
                        }
                        
                        current->isolated==0;
                        current->timeIsolated==0;
                        
                        (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        
                    }else if(current->location==3){
                        (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]--;
                        
                        if((*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]==0){
                            (*pIsolateMaxCell)[current->zone][current->locArea][current->locUnit][current->locCell]=0;
                            //if room is empty, set isolate array == 0
                        }
                        
                        ////move to non-isolated cell
                        cellFlag=0;
                        cellCap=0;
                        mCtr2=0;
                        //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                            while(cellFlag==0||mCtr2<4){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                                mCtr3=0;
                                while(cellFlag==0||mCtr3<5){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                    mCtr4=0;
                                    while(cellFlag==0||mCtr4<20){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                        cellCap=(*pMaxCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                        printf("cellCap: %d\n", cellCap);
                                        
                                        //update count
                                        //if((*pIsolateMaxCell)[obj->zone][mCtr2][mCtr3][mCtr4]==2){
                                            if(cellCap<2 && (*pIsolateMaxCell)[current->zone][mCtr2][mCtr3][mCtr4]==0){ //cells that are not full and not a red zone
                                                
                                                cellFlag=1;
                                                current->locArea=mCtr2;
                                                current->locUnit=mCtr3;
                                                current->locCell=mCtr4;
                                            }
                                            printf("cellFlag: %d\n", cellFlag);
                                        //}
                                        mCtr4++;
                                    }
                                    mCtr3++;
                                }
                                mCtr2++;
                            }
                        
                        if(cellFlag==0){
                            current->locArea=mCtr2-1;
                            current->locUnit=mCtr3-1;
                            current->locCell=mCtr4-1;
                            *povercap++;
                        }
                        
                        current->isolated==0;
                        current->timeIsolated==0;
                        
                        (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                    }
                    
                }
            }
        }
        current=current->nextIndiv;
    }
    
    printf("done updating inmates in isolation\n");
}

void checkQuarantine(sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *povercap){
    
    printf("updating inmates in quarantine\n");
    sIndiv *current;
    current=*pHeadCopy;
    int cellFlag=0;
    int cellCap=0;
    int deltaTime=0;
    int mCtr2, mCtr3, mCtr4;
    int location, prevZone, prevArea, prevUnit, prevCell;
    int mUnits, mAreas, mCells; //int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20]);
    int check=0;
    
    while(current!=NULL){
        //printf("traversing: ID: %d, loc: %d, group: %d\n",current->ID, current->location, current->group);
        //printf("traversing: ID: %d, Age: %d, ATSI: %d, loc: %d, metavir: %d, injFreq: %d, sharing: %d, shaFreq: %d, injOpd: %d, group: %d\n",current->ID, current->age, current->atsi, current->location, current->metavir, current->injFreq, current->sharing, current->shaFreq, current->injOpd, current->group);
        if(current->indivType==0){ //inmate
            if(current->isolated==2){ //only for quarantined staff
                deltaTime=currDay-(current->timeIsolated);
                prevArea=current->locArea;
                prevUnit=current->locUnit;
                prevCell=current->locCell;
                location=current->location;
                
                if(location==1){
                    mAreas=2;
                    mUnits=6;
                    mCells=13;
                }else if(location==2){
                    mAreas=2;
                    mUnits=4;
                    mCells=19;
                }else if(location==3){
                    mAreas=4;
                    mUnits=5;
                    mCells=20;
                }
                
                if(deltaTime==2 || deltaTime==13){
                    ////if COVID, isolate
                    if(current->COVID==1){
                        //isolate
                        printf("set isolation prison\n");
                        
                        cellFlag=0;
                        cellCap=0;
                        mCtr2=0;
                        //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                            while(cellFlag==0||mCtr2<mAreas){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                                mCtr3=0;
                                while(cellFlag==0||mCtr3<mUnits){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                    mCtr4=0;
                                    while(cellFlag==0||mCtr4<mCells){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                        cellCap=(*pMaxCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                        printf("cellCap: %d\n", cellCap);
                                        
                                        //update count
                                        //if((*pIsolateMaxCell)[obj->zone][mCtr2][mCtr3][mCtr4]==2){
                                            if(cellCap==0){ //new cells only
                                                
                                                cellFlag=1;
                                                current->locArea=mCtr2;
                                                current->locUnit=mCtr3;
                                                current->locCell=mCtr4;
                                                (*pIsolateMaxCell)[current->zone][mCtr2][mCtr3][mCtr4]==1; //set as 1out cell
                                            }
                                            printf("cellFlag: %d\n", cellFlag);
                                        //}
                                        mCtr4++;
                                    }
                                    mCtr3++;
                                }
                                mCtr2++;
                            }
                        
                        if(cellFlag==0){
                            current->locArea=mCtr2-1;
                            current->locUnit=mCtr3-1;
                            current->locCell=mCtr4-1;
                            (*povercap)++;
                        }
                        
                        current->isolated=1;
                        current->timeIsolated=currDay;
                        
                        if(location==1){
                            (*pMinCellArray)[current->zone][prevArea][prevUnit][prevCell]--;
                            (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        }else if(location==2){
                            (*pMedCellArray)[current->zone][prevArea][prevUnit][prevCell]--;
                            (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        }else if(location==3){
                            (*pMaxCellArray)[current->zone][prevArea][prevUnit][prevCell]--;
                            (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        }
                    }
                }
                
                if(deltaTime>14){
                    //if delta time is == 15 and covid ==0 unQuarantine
                    if(current->COVID==0){
                        //if delta time is == 15 and covid ==0 unQuarantine
                        
                        //remove from current cell
                        if(location==1){
                            (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]--;
                            if((*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]==0){
                                (*pIsolateMinCell)[current->zone][current->locArea][current->locUnit][current->locCell]=0;
                                //if room is empty, set isolate array == 0
                            }
                        }else if(location==2){
                            (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]--;
                            if((*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]==0){
                                (*pIsolateMedCell)[current->zone][current->locArea][current->locUnit][current->locCell]=0;
                                //if room is empty, set isolate array == 0
                            }
                        }else if(location==3){
                            (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]--;
                            if((*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]==0){
                                (*pIsolateMaxCell)[current->zone][current->locArea][current->locUnit][current->locCell]=0;
                                //if room is empty, set isolate array == 0
                            }
                        }
                        
                        ////move to non-isolated cell
                        cellFlag=0;
                        cellCap=0;
                        mCtr2=0;
                        location=current->location;
                        
                        if(location==1){
                            mAreas=2;
                            mUnits=6;
                            mCells=13;
                        }else if(location==2){
                            mAreas=2;
                            mUnits=4;
                            mCells=19;
                        }else if(location==3){
                            mAreas=4;
                            mUnits=5;
                            mCells=20;
                        }
                        //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                            while(cellFlag==0||mCtr2<mAreas){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                                mCtr3=0;
                                while(cellFlag==0||mCtr3<mUnits){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                    mCtr4=0;
                                    while(cellFlag==0||mCtr4<mCells){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                        cellCap=(*pMaxCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                        printf("cellCap: %d\n", cellCap);
                                        
                                        //update count
                                        //if((*pIsolateMaxCell)[obj->zone][mCtr2][mCtr3][mCtr4]==2){
                                        if(location==1){
                                            check=(*pIsolateMinCell)[current->zone][mCtr2][mCtr3][mCtr4];
                                        }else if(location==2){
                                            check=(*pIsolateMedCell)[current->zone][mCtr2][mCtr3][mCtr4];
                                        }else if(location==3){
                                            check=(*pIsolateMaxCell)[current->zone][mCtr2][mCtr3][mCtr4];
                                        }
                                        
                                        if(cellCap<2 && check==0){ //cells that are not full and not a red zone
                                                
                                                cellFlag=1;
                                                current->locArea=mCtr2;
                                                current->locUnit=mCtr3;
                                                current->locCell=mCtr4;
                                            }
                                            printf("cellFlag: %d\n", cellFlag);
                                        //}
                                        mCtr4++;
                                    }
                                    mCtr3++;
                                }
                                mCtr2++;
                            }
                        
                        if(cellFlag==0){
                            current->locArea=mCtr2-1;
                            current->locUnit=mCtr3-1;
                            current->locCell=mCtr4-1;
                            *povercap++;
                        }
                        
                        current->isolated==0;
                        current->timeIsolated==0;
                        
                        if(location==1){
                            (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        }else if(location==2){
                            (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        }else if(location==3){
                            (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                        }
                        
                    }
                    
                    //if room is empty, set isolate array == 0
                }
            }
        }
        
        current=current->nextIndiv;
    }
}

//this function performs PCR tests and implements lockdown if infection is confirmed
void PCRinmate(gsl_rng **r, sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18], int (*pLockdownMinArea)[27][2], int (*pLockdownMedArea)[11][2], int (*pLockdownMaxArea)[18][4], int (*pLockdownMinUnit)[27][2][6], int (*pLockdownMedUnit)[11][2][4], int (*pLockdownMaxUnit)[18][4][5], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *pMoveFlag){
    sIndiv *current;
    sIndiv *lockCurrent;
    current=*pHeadCopy;
    lockCurrent=*pHeadCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rN=0;
    double probTest[2];
    int binTest;
    
    int lockLocation;
    int lockPrison;
    int lockArea;
    int lockUnit;
    int lockCell;
    
    printf("PCR test inmate\n");
    //probTest[0]=0.712; //detected
    //probTest[1]=1-0.712; //not detected
    
    while(current!=NULL){
        if(current->indivType==0&&current->isolated==0){ //inmate not isolated
            if(current->severity>=3){ //have symptoms
                if(current->PCRtested==1){
                    if((current->PCRdate)<currDay){//if PCR tested since yesterday
                        if(current->COVID==1){
                            *pMoveFlag=1; //set movement restrictions
                        }
                        current->PCRtested=0;
                        current->PCRdate=0;
                    }
                }else if(current->PCRtested==0){//if not PCR tested
                    current->PCRtested=1;
                    current->PCRdate=currDay;
                }
            }
        }
        current=current->nextIndiv;
    }
}

//unlockdown
void unlock(gsl_rng **r, sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18], int (*pLockdownMinArea)[27][2], int (*pLockdownMedArea)[11][2], int (*pLockdownMaxArea)[18][4], int (*pLockdownMinUnit)[27][2][6], int (*pLockdownMedUnit)[11][2][4], int (*pLockdownMaxUnit)[18][4][5], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *pMoveFlag){
    
    sIndiv *current;
    sIndiv *lockCurrent;
    current=*pHeadCopy;
    lockCurrent=*pHeadCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rN=0;
    double probTest[2];
    int binTest;
    
    int lockLocation;
    int lockPrison;
    int lockArea;
    int lockUnit;
    int lockCell;
    
    int posCtr=0;
    
    printf("check to stop move restrictions\n");
    //probTest[0]=0.712; //detected
    //probTest[1]=1-0.712; //not detected
    
    while(current!=NULL){
        if(current->indivType==0){ //inmate not isolated
            if(current->PCRtested==1){
                if((current->PCRdate)<currDay){//if PCR tested since yesterday
                    if(current->COVID==1){
                        posCtr++;
                    }
                }
            }
        }
        current=current->nextIndiv;
    }
    
    if(posCtr==0){
        *pMoveFlag=0;
    }
    
}

int rapidTestPStaff (gsl_rng **r, sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int *PSTP, int *PSTN, int *PSFP, int *PSFN){
    sIndiv *current;
    current=*pHeadCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rN=0;
    double probTest[2];
    int binTest;
    
    printf("rapid test pstaff\n");
    //probTest[0]=0.712; //detected
    //probTest[1]=1-0.712; //not detected
    
    //if indiv is positive, probtest0 TP or probTest1 FN
    //if indiv is negative, probtest0 TN or probTest1 FP
    
    while(current!=NULL){
        if(current->indivType==1&&current->isolated==0){ //pstaff
            rN++;
            current->testDate=currDay;
            if(current->COVID==1){//if indiv is positive, probtest0 TP or probTest1 FN
                probTest[0]=0.71; //TP
                probTest[1]=1-0.71; //FN
                
                binTest=draw_multinom(&rCopy, 2, probTest);
                if(binTest==0){ //test success
                    current->detected=1;
                    current->testResultTrue=1;
                    current->isolated=1;
                    current->timeIsolated=currDay;
                    (*PSTP)++;
                }else if(binTest==1){ //test fail
                    current->detected=0;
                    current->testResultTrue=0;
                    (*PSFN)++;
                }
            }else if(current->COVID==0){//indiv is negative
                probTest[0]=0.98; //TN
                probTest[1]=1-0.98; //FP
                
                binTest=draw_multinom(&rCopy, 2, probTest);
                if(binTest==0){ //test success
                    current->detected=0;
                    current->testResultTrue=1;
                    (*PSTN)++;
                }else if(binTest==1){ //test fail
                    current->detected=1;
                    current->testResultTrue=0;
                    current->isolated=1;
                    current->timeIsolated=currDay;
                    (*PSFP)++;
                }
                //FP
            }
            current->tested=1;
            current->testDate=currDay;
        }
        current=current->nextIndiv;
    }
    return rN;
}

int rapidTestHStaff (gsl_rng **r, sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int *HSTP, int *HSTN, int *HSFP, int *HSFN){
    sIndiv *current;
    current=*pHeadCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rN=0;
    double probTest[2];
    int binTest;
    
    //probTest[0]=1-0.712; //not detected
    //probTest[1]=0.712; //detected
    printf("rapid test hstaff\n");
    
    while(current!=NULL){
        if(current->indivType==2&&current->isolated==0){ //pstaff
            rN++;
            current->testDate=currDay;
            if(current->COVID==1){//if indiv is positive, probtest0 TP or probTest1 FN
                probTest[0]=0.71; //TP
                probTest[1]=1-0.71; //FN
                
                binTest=draw_multinom(&rCopy, 2, probTest);
                if(binTest==0){ //test success
                    current->detected=1;
                    current->testResultTrue=1;
                    current->isolated=1;
                    current->timeIsolated=currDay;
                    (*HSTP)++;
                }else if(binTest==1){ //test fail
                    current->detected=0;
                    current->testResultTrue=0;
                    (*HSFN)++;
                }
            }else if(current->COVID==0){//indiv is negative
                probTest[0]=0.98; //TN
                probTest[1]=1-0.98; //FP
                
                binTest=draw_multinom(&rCopy, 2, probTest);
                if(binTest==0){ //test success
                    current->detected=0;
                    current->testResultTrue=1;
                    (*HSTN)++;
                }else if(binTest==1){ //test fail
                    current->detected=1;
                    current->testResultTrue=0;
                    current->isolated=1;
                    current->timeIsolated=currDay;
                    (*HSFP)++;
                }
                //FP
            }
            current->tested=1;
            current->testDate=currDay;
        }
        current=current->nextIndiv;
    }
    return rN;
}

void countZonePop (sIndiv **pHeadCopy, int (*pZoneArray)[30][3], int (*pZoneCovidArray)[30][3]){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0, tZone;
    
    if(current->indivType==0){
        while(current!=NULL){
            tZone=current->zone;
            
            if(current->location==1){ // Min
                (*pZoneArray)[tZone][0]++;
                if(current->COVID==1){
                    (*pZoneCovidArray)[tZone][0]++;
                }
            }else if(current->location==2){// Med
                (*pZoneArray)[tZone][1]++;
                if(current->COVID==1){
                    (*pZoneCovidArray)[tZone][1]++;
                }
            }else if(current->location==3){// Max
                (*pZoneArray)[tZone][2]++;
                if(current->COVID==1){
                    (*pZoneCovidArray)[tZone][2]++;
                }
            }
            
            current=current->nextIndiv;
        }
    }
}

int countHCVpris (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){

        if(current->placeInfected==1||current->placeInfected==2||current->placeInfected==3){ //means infected in a prison
            num++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

void countMovingPop (sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*nMovingPop)[20]){
    sIndiv *current;
    current=*pHeadCopy;
    //int num=0;
    
    while(current!=NULL){
        
        if(current->moving==1){ //means infected in a prison
            if(current->truckNumber!=99){
                (*nMovingPop)[current->truckNumber]++;
            }
        }
        
        current=current->nextIndiv;
    }
}

void countSevereInmates (sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pAgeSeverityArray)[7][8]){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    int tAge, aGroup, sGroup;
    int tSeverity;
    
    while(current!=NULL){
        
        if(current->indivType==0){ //inmates only
            if(current->COVID==1){
                tSeverity=current->severity;
                tAge=current->age;
                
                if(tAge<=19){ //means infected in a prison
                    (*pAgeSeverityArray)[0][tSeverity]++;
                }else if(tAge>=20&&tAge<=44){
                    (*pAgeSeverityArray)[1][tSeverity]++;
                }else if(tAge>=45&&tAge<=54){
                    (*pAgeSeverityArray)[2][tSeverity]++;
                }else if(tAge>=55&&tAge<=64){
                    (*pAgeSeverityArray)[3][tSeverity]++;
                }else if(tAge>=65&&tAge<=74){
                    (*pAgeSeverityArray)[4][tSeverity]++;
                }else if(tAge>=75&&tAge<=84){
                    (*pAgeSeverityArray)[5][tSeverity]++;
                }else if(tAge>=85){
                    (*pAgeSeverityArray)[6][tSeverity]++;
                }
            }
            
        }
        
        current=current->nextIndiv;
    }

}

int countOpd (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){

        if(current->injOpd==1){ //means infected in a prison
            num++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countOpdNotOST (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==1||current->risk==2||current->risk==3||current->risk==7||current->risk==8||current->risk==9){ //injecting opioid
            if(current->OST==0){ //if injecting opioid and not on OST
                num++;
            }
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countEverIDU (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){

        if(current->everIDU==1){ //means infected in a prison
            num++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countHCVCom (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->placeInfected==4){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countHCVantibody (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->metavir<=5){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countHCVRNA (sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->metavir<5){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countAveStay (sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay){
    sIndiv *current;
    current=*pHeadCopy;
    int totalDays=0;
    int totalSubjects=0;
    int average=0;
    int deltaDays=0;
    
    while(current!=NULL){
        deltaDays=currDay-(current->timeOfImprisonment);
        totalDays=totalDays+deltaDays;
        totalSubjects++;
        
        current=current->nextIndiv;
    }
    average=totalDays/totalSubjects;
    return average;
}

//0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 This is for the injecting risk group
int count0(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==0){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

//0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 This is for the injecting risk group
int count1(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==1){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

//0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 This is for the injecting risk group
int count2(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==2){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count3(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==3){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count4(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==4){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count5(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==5){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count6(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==6){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count7(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==7){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count8(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==8){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count9(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==9){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count10(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==10){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count11(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==11){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int count12(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->risk==12){ //means infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

int countR(sIndiv **pHeadCopy, sIndiv **pTailCopy){
    sIndiv *current;
    current=*pHeadCopy;
    int num=0;
    
    while(current!=NULL){
        
        if(current->infectionNumber>0){ //means re-infected in a prison
            num ++;
        }
        
        current=current->nextIndiv;
    }
    return num;
}

//float probAge(sIndiv **pTarget){ //Probability to age //MOVE TO OUTSIDE
  //  return 0.1;
//}

/*
 //normal
 double normMean;
 double normSigma=1.581139; //standard deviation of normalized vector (0-1) and minus the mean to adjust to mean 0.
 xI=gsl_ran_gaussian(r, normSigma);
 xI=xI+normMean; //to adjust back to a number between 0 and 1, add the normalised mean to the random number
 //printf("%f\n", xI);

 */

//float probMoveCom(gsl_rng **r, sIndiv **pTarget, double probsInput[ROWPRIS][ROWPRIS]){ //Probability to move locations
double probMoveCom(gsl_rng **r, sIndiv **pTarget, int currDay, int *pMoveFlag){ //Probability to get released
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    rN=0.0;
    if(current->indivType==0){
        if(current->hospitalised==0&&current->court==0&&current->isolated==0&&current->moving==0&&current->lockdown==0){
            if(current->location==1){//prison 1 MIN to Community *from 9010parameters.xlsx
                //rN=gsl_ran_exponential(*r, 0.006178163);
                if(current->zone==0){
                    rN=0.0136;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.003;
                }else if(current->zone==3){
                    rN=0.048;
                }else if(current->zone==4){
                    rN=0.001;
                }else if(current->zone==5){
                    rN=0.00238;
                }else if(current->zone==6){
                    rN=0.072;
                }else if(current->zone==7){
                    rN=0.00102;
                }else if(current->zone==8){
                    rN=0.01827;
                }else if(current->zone==9){
                    rN=0.010;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.013;
                }else if(current->zone==12){
                    rN=0.004;
                }else if(current->zone==13){
                    rN=0.0004;
                }else if(current->zone==14){
                    rN=0.00904;
                }else if(current->zone==15){
                    rN=0.00831;
                }else if(current->zone==16){
                    rN=0.00614;
                }else if(current->zone==17){
                    rN=0.011;
                }else if(current->zone==18){
                    rN=0.015;
                }else if(current->zone==19){
                    rN=0.010;
                }else if(current->zone==20){
                    rN=0.0034;
                }else if(current->zone==21){
                    rN=0.01027;
                }else if(current->zone==22){
                    rN=0.015;
                }else if(current->zone==23){
                    rN=0.009;
                }else if(current->zone==24){
                    rN=0.014;
                }else if(current->zone==25){
                    rN=0.004;
                }else if(current->zone==26){
                    rN=0.01525;
                }
            }else if(current->location==2){//prison 2 MED to Community
                if(current->zone==0){
                    rN=0.035;
                }else if(current->zone==1){
                    rN=0.007;
                }else if(current->zone==2){
                    rN=0.007;
                }else if(current->zone==3){
                    rN=0.00904;
                }else if(current->zone==4){
                    rN=0.019;
                }else if(current->zone==5){
                    rN=0.0215;
                }else if(current->zone==6){
                    rN=0.0465;
                }else if(current->zone==7){
                    rN=0.046;
                }else if(current->zone==8){
                    rN=0.003;
                }else if(current->zone==9){
                    rN=0.035;
                }else if(current->zone==10){
                    rN=0.012;
                }
            }else if(current->location==3){//prison 3 MAX to Community
                //rN=gsl_ran_exponential(*r, 0.001700223);//0.001800229);
                //}else if(currDay>4380){
                //    rN=gsl_ran_exponential(*r, 0.002700223);//0.001800229);
                //}
                if(current->zone==0){
                    rN=0.019;
                }else if(current->zone==1){
                    rN=0.00005;
                }else if(current->zone==2){
                    rN=0.005;
                }else if(current->zone==3){
                    rN=0.011;
                }else if(current->zone==4){
                    rN=0.006;
                }else if(current->zone==5){
                    rN=0.01343;
                }else if(current->zone==6){
                    rN=0.018;
                }else if(current->zone==7){
                    rN=0.0024;
                }else if(current->zone==8){
                    rN=0.002;
                }else if(current->zone==9){
                    rN=0.0066;
                }else if(current->zone==10){
                    rN=0.007;
                }else if(current->zone==11){
                    rN=0.1041;
                }else if(current->zone==12){
                    rN=0.04134;
                }else if(current->zone==13){
                    rN=0.09;
                }else if(current->zone==14){
                    rN=0.03343;
                }else if(current->zone==15){
                    rN=0.056;
                }else if(current->zone==16){
                    rN=0.0244;
                }else if(current->zone==17){
                    rN=0.037;
                }
            }
        }
        
        rN=rN/5;
        
        if(current->timeOfImprisonment==currDay){ //if inmate just got it, don't move straight away
            rN=0.0;
        }
        
    }else if(current->indivType==1){
        
    }else if(current->indivType==2){
        
    }else if(current->indivType==3){

    }else if(current->indivType==4){

    }
    
    if((*pMoveFlag)==1){
        rN=0.0;
    }
    
    return rN;
    
}

double probMoveP1(gsl_rng **r, sIndiv **pTarget, int currDay, int *pMoveFlag){ //Probability to move to MIN
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    rN=0.0;
    if(current->indivType==0){
        if(current->hospitalised==0&&current->court==0&&current->isolated==0&&current->moving==0&&current->lockdown==0){
            /*if(current->location==2){//prison 2 MED to MIN
                rN=gsl_ran_exponential(*r, 0.003022);
            }else if(current->location==3){//prison 3 MAX to MIN
                rN=gsl_ran_exponential(*r, 0.0005001);
            }else{//ALREADY in MIN, STAY IN MIN
                rN=gsl_ran_exponential(*r, 7.94E-04);
            }*/
            if(current->location==1){//prison 1 MIN to Community *from 9010parameters.xlsx
                //rN=gsl_ran_exponential(*r, 0.006178163);
                if(current->zone==0){
                    rN=0.022;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.001;
                }else if(current->zone==3){
                    rN=0.063;
                }else if(current->zone==4){
                    rN=0.0;
                }else if(current->zone==5){
                    rN=0.001;
                }else if(current->zone==6){
                    rN=0.059;
                }else if(current->zone==7){
                    rN=0.0;
                }else if(current->zone==8){
                    rN=0.005;
                }else if(current->zone==9){
                    rN=0.006;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.015;
                }else if(current->zone==12){
                    rN=0.002;
                }else if(current->zone==13){
                    rN=0.001;
                }else if(current->zone==14){
                    rN=0.005;
                }else if(current->zone==15){
                    rN=0.005;
                }else if(current->zone==16){
                    rN=0.005;
                }else if(current->zone==17){
                    rN=0.006;
                }else if(current->zone==18){
                    rN=0.009;
                }else if(current->zone==19){
                    rN=0.005;
                }else if(current->zone==20){
                    rN=0.008;
                }else if(current->zone==21){
                    rN=0.004;
                }else if(current->zone==22){
                    rN=0.019;
                }else if(current->zone==23){
                    rN=0.007;
                }else if(current->zone==24){
                    rN=0.006;
                }else if(current->zone==25){
                    rN=0.004;
                }else if(current->zone==26){
                    rN=0.008;
                }
            }else if(current->location==2){//prison 2 MED to Community
                if(current->zone==0){
                    rN=0.060;
                }else if(current->zone==1){
                    rN=0.004;
                }else if(current->zone==2){
                    rN=0.002;
                }else if(current->zone==3){
                    rN=0.004;
                }else if(current->zone==4){
                    rN=0.005;
                }else if(current->zone==5){
                    rN=0.025;
                }else if(current->zone==6){
                    rN=0.025;
                }else if(current->zone==7){
                    rN=0.024;
                }else if(current->zone==8){
                    rN=0.005;
                }else if(current->zone==9){
                    rN=0.019;
                }else if(current->zone==10){
                    rN=0.013;
                }
            }else if(current->location==3){//prison 3 MAX to Community
                //rN=gsl_ran_exponential(*r, 0.001700223);//0.001800229);
                //}else if(currDay>4380){
                //    rN=gsl_ran_exponential(*r, 0.002700223);//0.001800229);
                //}
                if(current->zone==0){
                    rN=0.017;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.007;
                }else if(current->zone==3){
                    rN=0.014;
                }else if(current->zone==4){
                    rN=0.011;
                }else if(current->zone==5){
                    rN=0.027;
                }else if(current->zone==6){
                    rN=0.010;
                }else if(current->zone==7){
                    rN=0.001;
                }else if(current->zone==8){
                    rN=0.001;
                }else if(current->zone==9){
                    rN=0.011;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.174;
                }else if(current->zone==12){
                    rN=0.023;
                }else if(current->zone==13){
                    rN=0.112;
                }else if(current->zone==14){
                    rN=0.077;
                }else if(current->zone==15){
                    rN=0.032;
                }else if(current->zone==16){
                    rN=0.019;
                }else if(current->zone==17){
                    rN=0.020;
                }
            }
        }
        
        if(current->timeOfImprisonment==currDay){ //if inmate just got it, don't move straight away
            rN=0.0;
        }
        
    }else if(current->indivType==1){
        
    }else if(current->indivType==2){
        
    }else if(current->indivType==3){
        rN=0.0;
    }else if(current->indivType==4){
        rN=0.0;
    }
    
    if((*pMoveFlag)==1){
        rN=0.0;
    }
    
    return rN;
}

double probMoveP2(gsl_rng **r, sIndiv **pTarget){ //Probability to move to MED
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){
        if(current->hospitalised==0&&current->court==0&&current->isolated==0&&current->moving==0){
            if(current->location==1){//prison 1 MIN to MED
                rN=gsl_ran_exponential(*r, 3.82E-04);
            }else if(current->location==3){//prison 3 MAX to MED
                rN=gsl_ran_exponential(*r, 0.0005038);
            }else{//ALREADY IN MED, STAY IN MED
                rN=gsl_ran_exponential(*r, 0.0004631);
            }
            /*if(current->location==1){//prison 1 MIN to Community *from 9010parameters.xlsx
                //rN=gsl_ran_exponential(*r, 0.006178163);
                if(current->zone==0){
                    rN=0.022;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.001;
                }else if(current->zone==3){
                    rN=0.063;
                }else if(current->zone==4){
                    rN=0.0;
                }else if(current->zone==5){
                    rN=0.001;
                }else if(current->zone==6){
                    rN=0.059;
                }else if(current->zone==7){
                    rN=0.0;
                }else if(current->zone==8){
                    rN=0.005;
                }else if(current->zone==9){
                    rN=0.006;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.015;
                }else if(current->zone==12){
                    rN=0.002;
                }else if(current->zone==13){
                    rN=0.001;
                }else if(current->zone==14){
                    rN=0.005;
                }else if(current->zone==15){
                    rN=0.005;
                }else if(current->zone==16){
                    rN=0.005;
                }else if(current->zone==17){
                    rN=0.006;
                }else if(current->zone==18){
                    rN=0.009;
                }else if(current->zone==19){
                    rN=0.005;
                }else if(current->zone==20){
                    rN=0.008;
                }else if(current->zone==21){
                    rN=0.004;
                }else if(current->zone==22){
                    rN=0.019;
                }else if(current->zone==23){
                    rN=0.007;
                }else if(current->zone==24){
                    rN=0.006;
                }else if(current->zone==25){
                    rN=0.004;
                }else if(current->zone==26){
                    rN=0.008;
                }
            }else if(current->location==2){//prison 2 MED to Community
                if(current->zone==0){
                    rN=0.060;
                }else if(current->zone==1){
                    rN=0.004;
                }else if(current->zone==2){
                    rN=0.002;
                }else if(current->zone==3){
                    rN=0.004;
                }else if(current->zone==4){
                    rN=0.005;
                }else if(current->zone==5){
                    rN=0.025;
                }else if(current->zone==6){
                    rN=0.025;
                }else if(current->zone==7){
                    rN=0.024;
                }else if(current->zone==8){
                    rN=0.005;
                }else if(current->zone==9){
                    rN=0.019;
                }else if(current->zone==10){
                    rN=0.013;
                }
            }else if(current->location==3){//prison 3 MAX to Community
                //rN=gsl_ran_exponential(*r, 0.001700223);//0.001800229);
                //}else if(currDay>4380){
                //    rN=gsl_ran_exponential(*r, 0.002700223);//0.001800229);
                //}
                if(current->zone==0){
                    rN=0.017;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.007;
                }else if(current->zone==3){
                    rN=0.014;
                }else if(current->zone==4){
                    rN=0.011;
                }else if(current->zone==5){
                    rN=0.027;
                }else if(current->zone==6){
                    rN=0.010;
                }else if(current->zone==7){
                    rN=0.001;
                }else if(current->zone==8){
                    rN=0.001;
                }else if(current->zone==9){
                    rN=0.011;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.174;
                }else if(current->zone==12){
                    rN=0.023;
                }else if(current->zone==13){
                    rN=0.112;
                }else if(current->zone==14){
                    rN=0.077;
                }else if(current->zone==15){
                    rN=0.032;
                }else if(current->zone==16){
                    rN=0.019;
                }else if(current->zone==17){
                    rN=0.020;
                }
            }*/

        }
    }else if(current->indivType==1){
        
    }else if(current->indivType==2){
        
    }else if(current->indivType==3){
        rN=0.0;
    }else if(current->indivType==4){
        rN=0.0;
    }
    
    return rN;
}

double probMoveP3(gsl_rng **r, sIndiv **pTarget){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){
        if(current->hospitalised==0&&current->court==0&&current->isolated==0&&current->moving==0){
            if(current->location==1){//prison 1 MIN to MAX
                rN=gsl_ran_exponential(*r, 0.000225);
            }else if(current->location==2){//prison 2 MED to MAX
                rN=gsl_ran_exponential(*r, 2.88E-04);
            }else{//ALREADY IN MAX, STAY IN MAX
                rN=gsl_ran_exponential(*r, 0.0008053);
            }
            /*if(current->location==1){//prison 1 MIN to Community *from 9010parameters.xlsx
                //rN=gsl_ran_exponential(*r, 0.006178163);
                if(current->zone==0){
                    rN=0.022;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.001;
                }else if(current->zone==3){
                    rN=0.063;
                }else if(current->zone==4){
                    rN=0.0;
                }else if(current->zone==5){
                    rN=0.001;
                }else if(current->zone==6){
                    rN=0.059;
                }else if(current->zone==7){
                    rN=0.0;
                }else if(current->zone==8){
                    rN=0.005;
                }else if(current->zone==9){
                    rN=0.006;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.015;
                }else if(current->zone==12){
                    rN=0.002;
                }else if(current->zone==13){
                    rN=0.001;
                }else if(current->zone==14){
                    rN=0.005;
                }else if(current->zone==15){
                    rN=0.005;
                }else if(current->zone==16){
                    rN=0.005;
                }else if(current->zone==17){
                    rN=0.006;
                }else if(current->zone==18){
                    rN=0.009;
                }else if(current->zone==19){
                    rN=0.005;
                }else if(current->zone==20){
                    rN=0.008;
                }else if(current->zone==21){
                    rN=0.004;
                }else if(current->zone==22){
                    rN=0.019;
                }else if(current->zone==23){
                    rN=0.007;
                }else if(current->zone==24){
                    rN=0.006;
                }else if(current->zone==25){
                    rN=0.004;
                }else if(current->zone==26){
                    rN=0.008;
                }
            }else if(current->location==2){//prison 2 MED to Community
                if(current->zone==0){
                    rN=0.060;
                }else if(current->zone==1){
                    rN=0.004;
                }else if(current->zone==2){
                    rN=0.002;
                }else if(current->zone==3){
                    rN=0.004;
                }else if(current->zone==4){
                    rN=0.005;
                }else if(current->zone==5){
                    rN=0.025;
                }else if(current->zone==6){
                    rN=0.025;
                }else if(current->zone==7){
                    rN=0.024;
                }else if(current->zone==8){
                    rN=0.005;
                }else if(current->zone==9){
                    rN=0.019;
                }else if(current->zone==10){
                    rN=0.013;
                }
            }else if(current->location==3){//prison 3 MAX to Community
                //rN=gsl_ran_exponential(*r, 0.001700223);//0.001800229);
                //}else if(currDay>4380){
                //    rN=gsl_ran_exponential(*r, 0.002700223);//0.001800229);
                //}
                if(current->zone==0){
                    rN=0.017;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.007;
                }else if(current->zone==3){
                    rN=0.014;
                }else if(current->zone==4){
                    rN=0.011;
                }else if(current->zone==5){
                    rN=0.027;
                }else if(current->zone==6){
                    rN=0.010;
                }else if(current->zone==7){
                    rN=0.001;
                }else if(current->zone==8){
                    rN=0.001;
                }else if(current->zone==9){
                    rN=0.011;
                }else if(current->zone==10){
                    rN=0.005;
                }else if(current->zone==11){
                    rN=0.174;
                }else if(current->zone==12){
                    rN=0.023;
                }else if(current->zone==13){
                    rN=0.112;
                }else if(current->zone==14){
                    rN=0.077;
                }else if(current->zone==15){
                    rN=0.032;
                }else if(current->zone==16){
                    rN=0.019;
                }else if(current->zone==17){
                    rN=0.020;
                }
            }*/

        }
    }else if(current->indivType==1){
        
    }else if(current->indivType==2){
        
    }else if(current->indivType==3){
        rN=0.0;
    }else if(current->indivType==4){
        rN=0.0;
    }
    
    return rN;
}

double probOutsideVisit(gsl_rng **r, sIndiv **pTarget, int currDay, int *pMoveFlag){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub;
    
    if(current->indivType==0){ //individual
        if(current->hospitalised==0&&current->court==0&&current->isolated==0&&current->moving==0&&current->lockdown==0){
            /*if(current->location==1){//prison 1 MIN to MAX
                rN=0.008;
                //rN=gsl_ran_exponential(*r, 0.000225);
            }else if(current->location==2){//prison 2 MED to MAX
                //rN=gsl_ran_exponential(*r, 2.88E-04);
                rN=0.008;
            }else{//ALREADY IN MAX, STAY IN MAX
                //rN=gsl_ran_exponential(*r, 0.0008053);
                rN=0.008;
            }*/
            if(current->location==1){//prison 1 MIN to Community *from 9010parameters.xlsx
                //rN=gsl_ran_exponential(*r, 0.006178163);
                if(current->zone==0){
                    rN=0.007;
                }else if(current->zone==1){
                    rN=0.0001;
                }else if(current->zone==2){
                    rN=0.002;
                }else if(current->zone==3){
                    rN=0.02;
                }else if(current->zone==4){
                    rN=0.0;
                }else if(current->zone==5){
                    rN=0.0;
                }else if(current->zone==6){
                    rN=0.048;
                }else if(current->zone==7){
                    rN=0.0006;
                }else if(current->zone==8){
                    rN=0.009;
                }else if(current->zone==9){
                    rN=0.0001;
                }else if(current->zone==10){
                    rN=0.001;
                }else if(current->zone==11){
                    rN=0.014;
                }else if(current->zone==12){
                    rN=0.0;
                }else if(current->zone==13){
                    rN=0.0;
                }else if(current->zone==14){
                    rN=0.006;
                }else if(current->zone==15){
                    rN=0.001;
                }else if(current->zone==16){
                    rN=0.0003;
                }else if(current->zone==17){
                    rN=0.011;
                }else if(current->zone==18){
                    rN=0.015;
                }else if(current->zone==19){
                    rN=0.008;
                }else if(current->zone==20){
                    rN=0.0001;
                }else if(current->zone==21){
                    rN=0.002;
                }else if(current->zone==22){
                    rN=0.021;
                }else if(current->zone==23){
                    rN=0.003;
                }else if(current->zone==24){
                    rN=0.0001;
                }else if(current->zone==25){
                    rN=0.008;
                }else if(current->zone==26){
                    rN=0.0172;
                }
            }else if(current->location==2){//prison 2 MED to Community
                if(current->zone==0){
                    rN=0.019;
                }else if(current->zone==1){
                    rN=0.0002;
                }else if(current->zone==2){
                    rN=0.006;
                }else if(current->zone==3){
                    rN=0.0003;
                }else if(current->zone==4){
                    rN=0.012;
                }else if(current->zone==5){
                    rN=0.024;
                }else if(current->zone==6){
                    rN=0.06;
                }else if(current->zone==7){
                    rN=0.035;
                }else if(current->zone==8){
                    rN=0.001;
                }else if(current->zone==9){
                    rN=0.031;
                }else if(current->zone==10){
                    rN=0.027;
                }
            }else if(current->location==3){//prison 3 MAX to Community
                //rN=gsl_ran_exponential(*r, 0.001700223);//0.001800229);
                //}else if(currDay>4380){
                //    rN=gsl_ran_exponential(*r, 0.002700223);//0.001800229);
                //}
                if(current->zone==0){
                    rN=0.006;
                }else if(current->zone==1){
                    rN=0.001;
                }else if(current->zone==2){
                    rN=0.002;
                }else if(current->zone==3){
                    rN=0.002;
                }else if(current->zone==4){
                    rN=0.018;
                }else if(current->zone==5){
                    rN=0.044;
                }else if(current->zone==6){
                    rN=0.019;
                }else if(current->zone==7){
                    rN=0.002;
                }else if(current->zone==8){
                    rN=0.002;
                }else if(current->zone==9){
                    rN=0.002;
                }else if(current->zone==10){
                    rN=0.008;
                }else if(current->zone==11){
                    rN=0.19;
                }else if(current->zone==12){
                    rN=0.03;
                }else if(current->zone==13){
                    rN=0.12;
                }else if(current->zone==14){
                    rN=0.04;
                }else if(current->zone==15){
                    rN=0.038;
                }else if(current->zone==16){
                    rN=0.009;
                }else if(current->zone==17){
                    rN=0.041;
                }
            }
        }
        
        rN=rN-(rN*.2);
        
        if(current->timeOfImprisonment==currDay){ //if inmate just got it, don't move straight away
            rN=0.0;
        }
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){
        rN=0.0;
    }else if(current->indivType==4){
        rN=0.0;
    }
    
    if((*pMoveFlag)==1){
        rN=0.0;
    }
    
    return rN;
}

double probMoveHospitalF(gsl_rng **r, sIndiv **pTarget, int nMinHospitalF, int nMedHospitalF, int nMaxHospitalF){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    int capMinHospitalF=100, capMedHospitalF=100, capMaxHospitalF=100;
    double rN=0.0, lb, ub;
    
    if(current->indivType==0){
        if(current->hospitalised==0&&current->court==0){
            if(current->severity==4){//if moderate
                //if(current->location==1){//prison 1 MIN to MAX
                //  if(nMinHospitalF<capMinHospitalF){
                rN=1.0;
                //}
                //}else if(current->location==2){//prison 2 MED to MAX
                //  if(nMinHospitalF<capMinHospitalF){
                //    rN=1.0;
                //}
                //}else{//ALREADY IN MAX, STAY IN MAX
                //  if(nMinHospitalF<capMinHospitalF){
                //    rN=1.0;
                    //}
                //}
            }
        }
    }

    return rN;
}

double probMoveHospitalC(gsl_rng **r, sIndiv **pTarget, int nMinHospitalF, int nMedHospitalF, int nMaxHospitalF){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    int capMinHospitalF=100, capMedHospitalF=100, capMaxHospitalF=100;
    double rN=0.0, lb, ub;
    
    if(current->indivType==0){
        if(current->severity==5){//if severe
            if(current->hospitalised==0){
                if(current->location==1){//prison 1 MIN to MAX
                    if(nMinHospitalF>=capMinHospitalF){
                        rN=1.0;
                    }
                }else if(current->location==2){//prison 2 MED to MAX
                    if(nMinHospitalF>=capMinHospitalF){
                        rN=1.0;
                    }
                }else{//ALREADY IN MAX, STAY IN MAX
                    if(nMinHospitalF>=capMinHospitalF){
                        rN=1.0;
                    }
                }
            }
        }
    }
    
    return rN;
}

double probReturnPrison(gsl_rng **r, sIndiv **pTarget, int nMinHospitalF, int nMedHospitalF, int nMaxHospitalF){
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub;
    
    if(current->indivType==0){
        if(current->hospitalised!=0){
            if(current->COVID==0){//if moderate
                rN=1.0;
            }
        }else if(current->court!=0){
            rN=1.0;
        }else if(current->moving!=0){
            rN=1.0;
        }
    }
    
    return rN;
}

double probTest(gsl_rng **r, sIndiv **pTarget, int currDay, int nInmateTests, int nPStaffTests, int nHStaffTests, int nInmatesTested, int nPStaffTested, int nHStaffTested, int nEVisitorsTested, int nFVisitorsTested){ //Detect
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub, testDate, testDaysNeeded, testResultDay;
    
    //case 1: test if not tested yet,
    //case 2: if awaiting result
    if(current->indivType==0){ //individual
        if(current->tested==0){ //not yet tested
            if(current->detected==0){ //not yet detected
                //if(nInmatesTested<nInmateTests){ //enough tests available
                    rN=1.0;
                //}
            }
        }
    }else if(current->indivType==1){ //prison staff
        if(current->tested==0){ //not yet tested
            if(current->detected==0){ //not yet detected
                //if(nPStaffTested<nPStaffTests){ //enough tests available
                    rN=1.0;
                //}
            }
        }
    }else if(current->indivType==2){ //hospital staff
        if(current->tested==0){ //not yet tested
            if(current->detected==0){ //not yet detected
                //if(nHStaffTested<nHStaffTests){ //enough tests available
                    rN=1.0;
                //}
            }
        }
    }
    
    return rN;
}

double probTestResult(gsl_rng **r, sIndiv **pTarget, int currDay){ //Detect
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub, testDate, testDaysNeeded, testResultDay;
    
    if(current->tested==1){ //tested
        testDate=current->testDate;
        testDaysNeeded=current->testDaysBeforeResult;
        testResultDay=testDate+testDaysNeeded;
        //don't test those who are already detected
        if(testResultDay==currDay){
            rN=1.0;
        }//////what about for those who were tested but undetected?
    }
    
    return rN;
}

double probRapidTest(gsl_rng **r, sIndiv **pTarget, int currDay, int nTests){ //Detect
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub;
    
    if(current->tested==0){ //if already tested today or not in isolation
        if(current->indivType==0){ //individual
            rN=0.0;
        }else if(current->indivType==1){ //prison staff
            rN=1.0;
        }else if(current->indivType==2){ //hospital staff
            rN=1.0;
        }else if(current->indivType==3){ //essential visitor rapid teast
            rN=0.0;
        }else if(current->indivType==4){ //family visitor rapid test
            rN=0.0;
        }
    }
    
    return rN;
}

double probThermalTest(gsl_rng **r, sIndiv **pTarget, int currDay, int nTests){ //Detect
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->tested==0){
        if(current->indivType==0){ //individual
            rN=0.0;
        }else if(current->indivType==1){ //prison staff
            rN=0.0;
        }else if(current->indivType==2){ //hospital staff
            rN=0.0;
        }else if(current->indivType==3){ //essential visitor
            rN=0.0;
        }else if(current->indivType==4){ //family visitor
            rN=0.0;
        }
    }
    
    return rN;
}

double probInfectIndiv(sIndiv **pTarget, int currDay, int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18], int (*pLockdownMinArea)[27][2], int (*pLockdownMedArea)[11][2], int (*pLockdownMaxArea)[18][4], int (*pLockdownMinUnit)[27][2][6], int (*pLockdownMedUnit)[11][2][4], int (*pLockdownMaxUnit)[18][4][5], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20]){ //Infected can infect everyone so 1 if symptomatic
    sIndiv *current;
    current=*pTarget;
    double rN=0.0;
    int infTime;
    infTime=current->timeOfInfection;
    
    int lockLocation;
    int lockPrison;
 //   int lockArea;
 //   int lockUnit;

    //reduce depending on PPE
    if(current->COVID==1){
        if(current->hospitalised==0){ //&&current->isolated==0
            if(current->severity>=2&&current->severity<=5){ //asymptomatic to severe
                if(current->indivType==0){ //individual
                    //infTime=infTime+3;
                    //if(infTime>=currDay){
                    rN=1.0-0.0504; //assuming standard facemask worn at all times
                    //}else{
                    //    rN=0.0;
                    //}
                }else if(current->indivType==1){ //prison staff
                    lockLocation=current->location;
                    lockPrison=current->zone;
                    //lockArea=current->locArea; staff doesn't have assigned area and unit
                    //lockUnit=current->locUnit;
                    
                    rN=1.0-0.0504;//assuming standard facemask worn at all times
                    if(lockLocation==1){ //min
                        if((*pLockdownMinZone)[lockPrison]==1){
                            rN=1.0-0.18;//0.0504; //assuming standard facemask worn at all times
                        }
                    }else if(lockLocation==2){ //med
                        if((*pLockdownMedZone)[lockPrison]==1){
                            rN=1.0-0.18;//0.0504; //assuming standard facemask worn at all times
                        }
                    }else if(lockLocation==3){ //max
                        if((*pLockdownMaxZone)[lockPrison]==1){
                            rN=1.0-0.18;//0.0504; //assuming standard facemask worn at all times
                        }
                    }
                }else if(current->indivType==2){ //hospital staff
                    lockLocation=current->location;
                    lockPrison=current->zone;
                    //lockArea=current->locArea;
                    //lockUnit=current->locUnit;
                    
                    rN=1.0-0.0504;//assuming standard facemask worn at all times
                    if(lockLocation==1){ //min
                        if((*pLockdownMinZone)[lockPrison]==1){
                            rN=1.0-0.18;//0.0504; //assuming standard facemask worn at all times
                        }
                    }else if(lockLocation==2){ //med
                        if((*pLockdownMedZone)[lockPrison]==1){
                            rN=1.0-0.18;//0.0504; //assuming standard facemask worn at all times
                        }
                    }else if(lockLocation==3){ //max
                        if((*pLockdownMaxZone)[lockPrison]==1){
                            rN=1.0-0.18;//0.0504; //assuming standard facemask worn at all times
                        }
                    }
                }else if(current->indivType==3){ //essential visitor
                    rN=1.0-0.0504; //assuming standard facemask
                }else if(current->indivType==4){ //family visitor
                    rN=1.0-0.0504; //assuming standard facemask
                }
            }
        }
    }
    
    return rN;
}

double probNtrClear(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to clear (only if infected and depending on severity and delta time? HCV: metavir scale low && infection < 6 mths)
    sIndiv *current;
    current=*pTarget;
    int deltaTime;
    double rN=0.0;
    double normMean;
    int peakTime, startTime;
    startTime=current->timeOfInfection;
    deltaTime=currDay-startTime;
    double factor;
    
    if(current->COVID==1){
        if(current->severity==0){ //exposed
            //peakTime=startTime+7;
            //if(deltaTime>=3){ //can clear only after 5 days
            //factor=deltaTime-1;
            //rN=0.10*(factor*0.5);
            //if(rN>=1.0){
            //    rN=1.0;
            //}
            rN=0.0; //no clearance during exposed
            //}
        }else if(current->severity==1){ //pre-clinical
            //if(deltaTime>=5&&deltaTime<=16){ //can clear only after 16 days
            //factor=deltaTime-7;
            //rN=0.12*(factor*0.5);
            //if(rN>=1.0){
            //    rN=1.0;
            //}
            rN=0.0;
            //}
        }else if(current->severity==2){ //asympotamtic
            //if(deltaTime>=5&&deltaTime<=16){ //can clear only after 5 days
            //factor=deltaTime-10;
            //rN=0.05*(factor*0.5);
            //if(rN>=1.0){
            //    rN=1.0;
            //}
            //rN=0.125;
            rN=1-exp(-0.125);//gsl_ran_lognormal(*r, 2.0794, 0.2394);
            //}
        }else if(current->severity==3){ //mild
            //if(deltaTime>=5&&deltaTime<=16){ //can clear only after 5 days
            //factor=deltaTime-14;
            //rN=0.05*(factor*0.5);
            //if(rN>=1.0){
            //    rN=1.0;
            //}
            //rN=0.125;
            rN=1-exp(-0.125);//gsl_ran_lognormal(*r, 2.0794, 0.2394);
            //}
        }else if(current->severity==4){ //moderate
            //if(deltaTime>=17&&deltaTime<=23){ //can clear only after 5 days
            //factor=deltaTime-14;
            //rN=0.05*(factor*0.5);
            //if(rN>=1.0){
            //    rN=1.0;
            //}
            //rN=0.071428571;
            rN=1-exp(-0.071);//gsl_ran_lognormal(*r, 2.6391, 0.1678);
            //}
        }else if(current->severity==5){ //severe
            //if(deltaTime>=24){ //can clear only after 5 days
            //factor=deltaTime-14;
            //rN=0.05*(factor*0.5);
            //if(rN>=1.0){
            //    rN=1.0;
            //}
            //rN=0.071428571;
            rN=1-exp(-0.026);//gsl_ran_lognormal(*r, 3.627, 0.0636);
            //}
        }
    }
    
    /*if(current->COVID==1){
     if(current->hospitalised==0){
     infTime=infTime+7; //within the first 7 days, there is chance to recover 1/7 * no chance of recovering after
     if(currDay>=infTime){ //delta days and severity condition?
     rN=0.5;
     }else{
     rN=0.1;
     }
     }else{
     rN=0.5;
     }
     }
     */
    
    return rN;
    
    /*deltaTime=currDay-(current->timeOfInfection);
     if(current->metavir<3&&deltaTime<=180){ //infected and less than 180 days
     //exponential
     normMean=0.02/365;//0.05/365;
     rN=gsl_ran_exponential(*r, normMean);
     //rN=0.001150685;
     ////rN=0.115068493;
     //rN=0.139726;
     return rN; //In function for clearing add condition, if lower than threshold (mean) then clear.
     }else if(current->metavir<3&&(deltaTime>180&&deltaTime<360)){ //infected and over 180 days, different rate
     normMean=0.009/365; //0.01
     rN=gsl_ran_exponential(*r, normMean);
     //rN=0.000273973;
     ////rN=0.02739726;
     //rN=0.03561643836;
     return rN;
     }else //not infected
     return 0.0;*/
    
}


double probProgress(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to progress metavir scale
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub;
    int infTime, age, deltaTime;
    
    infTime=current->timeOfInfection;
    age=current->age;
    
    deltaTime=currDay-infTime;
    
    if(current->COVID==1){
        rN=1.0;
        if(current->severity==0){ //exposed -> pre-clinical
            //infTime=infTime+2;
            /*if(deltaTime>=4&&deltaTime<=7){ //delta days and severity condition?
             if(age>=0&&age<=19){
             rN=gsl_ran_flat(*r, 0.625, 0.85);
             }else if(age>=20&&age<=44){
             rN=gsl_ran_flat(*r, 0.56, 0.77);
             }else if(age>=45&&age<=54){
             rN=gsl_ran_flat(*r, 0.44, 0.675);
             }else if(age>=55&&age<=64){
             rN=gsl_ran_flat(*r, 0.32, 0.57);
             }else if(age>=65&&age<=74){
             rN=gsl_ran_flat(*r, 0.21, 0.47);
             }else if(age>=75&&age<=84){
             rN=gsl_ran_flat(*r, 0.18, 0.43);
             }else{
             rN=gsl_ran_flat(*r, 0.18, 0.43);
             }
             }*/
            //rN=0.30;
            rN=1-exp(-0.30);//gsl_ran_gamma(*r, 0.0481, 6.2322);
        }else if(current->severity==1){ //pre-clinical -> asymptomatic or mild
            //infTime=infTime+2;
            //currDay>=infTime
            /*if(deltaTime>=4&&deltaTime<=7){ //delta days and severity condition?
             if(age>=0&&age<=19){
             rN=gsl_ran_flat(*r, 0.15, 0.375);
             }else if(age>=20&&age<=44){
             rN=gsl_ran_flat(*r, 0.23, 0.44);
             }else if(age>=45&&age<=54){
             rN=gsl_ran_flat(*r, 0.325, 0.56);
             }else if(age>=55&&age<=64){
             rN=gsl_ran_flat(*r, 0.43, 0.68);
             }else if(age>=65&&age<=74){
             rN=gsl_ran_flat(*r, 0.53, 0.79);
             }else if(age>=75&&age<=84){
             rN=gsl_ran_flat(*r, 0.57, 0.82);
             }else{
             rN=gsl_ran_flat(*r, 0.57, 0.82);
             }
             }*/
            rN=(1-exp(-0.40))*0.33;//0.40; //shorten time to infectiousness from 3 to 2 days
        }else if(current->severity==3){ //mild to moderate
            //infTime=infTime+5;
            /*if(deltaTime>=23){ //delta days and severity condition?
             if(age>=0&&age<=23){
             rN=gsl_ran_flat(*r, 0.0243, 0.0832);
             }else if(age>=20&&age<=44){
             rN=gsl_ran_flat(*r, 0.03927, 0.1347);
             }else if(age>=45&&age<=54){
             rN=gsl_ran_flat(*r, 0.03695, 0.1269);
             }else if(age>=55&&age<=64){
             rN=gsl_ran_flat(*r, 0.0998, 0.2035);
             }else if(age>=65&&age<=74){
             rN=gsl_ran_flat(*r, 0.0844, 0.289);
             }else if(age>=75&&age<=84){
             rN=gsl_ran_flat(*r, 0.10435, 0.357);
             }else{
             rN=gsl_ran_flat(*r, 0.11, 0.376);
             }
             }*/
            //rN=0.19;
            rN=1-exp(-0.19);//gsl_ran_gamma(*r, 0.1633, 1.0412);
        }else if(current->severity==4){ //moderate to severe
            //infTime=infTime+5;
            //if(deltaTime>=23){ //delta days and severity condition?
            if(age>=0&&age<=19){
                //rN=0.0001;//gsl_ran_flat(*r, 0.0243, 0.0832);
                rN=gsl_ran_gamma(*r, 0.0075, 0.0145);
            }else if(age>=20&&age<=44){
                //rN=0.0157;//gsl_ran_flat(*r, 0.03927, 0.1347);
                rN=gsl_ran_gamma(*r, 0.0082, 1.913);
            }else if(age>=45&&age<=54){
                //rN=0.0336;//gsl_ran_flat(*r, 0.03695, 0.1269);
                rN=gsl_ran_gamma(*r, 0.0035, 9.5657);
            }else if(age>=55&&age<=64){
                //rN=0.0541;//gsl_ran_flat(*r, 0.0998, 0.2035);
                rN=gsl_ran_gamma(*r, 0.0084, 6.4473);
            }else if(age>=65&&age<=74){
                //rN=0.0769;//gsl_ran_flat(*r, 0.0844, 0.289);
                rN=gsl_ran_gamma(*r, 0.0084, 9.1154);
            }else if(age>=75&&age<=84){
                //rN=0.0948;//gsl_ran_flat(*r, 0.10435, 0.357);
                rN=gsl_ran_gamma(*r, 0.0084, 11.3455);
            }else{
                //rN=0.0997;//gsl_ran_flat(*r, 0.11, 0.376);
                rN=gsl_ran_gamma(*r, 0.0083, 12.0019);
            }
            //}
            //rN=1.0;
        }
    }
    
    return rN;
}


double probCOVIDdeath(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to die via HCV (only if infected and metavir scale up)
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub;
    int infTime;
    
    infTime=current->timeOfInfection;
    int age;
    age=current->age;
    
    //
    if(current->COVID==1){
        if(current->severity==5){ //delta days and severity condition?
            if(age>=0&&age<=19){
                rN=0.0000092;//gsl_ran_flat(*r, 0.00001, 0.00038); //0.00001-0.00038
                //rN=gsl_ran_gamma(*r, 0.000125108, 0.073536);
            }else if(age>=20&&age<=44){
                rN=0.002;//0.0455;//0.01490//gsl_ran_flat(*r, 0.00044, 0.002); //0.00044-0.002
                //rN=gsl_ran_gamma(*r, 0.01475821, 3.0830);
            }else if(age>=45&&age<=54){
                rN=0.00398;//0.00611//gsl_ran_flat(*r, 0.002102, 0.008015); //0.002102-0.008015
                //rN=gsl_ran_gamma(*r, 0.007346, 0.541736);
            }else if(age>=55&&age<=64){
                rN=0.00495;//0.02027;//gsl_ran_flat(*r, 0.00727, 0.02585);
                //rN=gsl_ran_gamma(*r, 0.0008173, 0.60558);
            }else if(age>=65&&age<=74){
                rN=0.0275;//0.05003;//gsl_ran_flat(*r, 0.0128, 0.06165);
                //rN=gsl_ran_gamma(*r, 0.0093674, 2.93569);
            }else if(age>=75&&age<=84){
                rN=0.05;//0.09716;//gsl_ran_flat(*r, 0.02625, 0.1087);
                //rN=gsl_ran_gamma(*r, 0.014135, 3.537261);
            }else{
                rN=0.10807;//0.12547;//gsl_ran_flat(*r, 0.038, 0.133);
                //rN=gsl_ran_gamma(*r, 0.0181797, 5.944538);
            }
            
            if(current->atsi==1){
                rN=rN*1.6;
            }
        }else{
            rN=0.0;
        }
    }else{
        rN=0.0;
    }
    
    return rN;
    
    /*if(current->metavir==4){ //infected
     lb=4.23625E-06/365;
     ub=5.35469E-06/365;
     rN=gsl_ran_flat(*r, lb, ub);
     return rN;
     }else //not infected
     return 0.0;*/
}

double probNtrDeath(gsl_rng **r, sIndiv **pTarget){ //Probability to die of natural causes (age)
    //uniform distribution
    double rN; double ub=5.35469E-06; double lb=4.23625E-06;
    rN=gsl_ran_flat(*r, lb, ub);
    return rN;
}

double probIsolate(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    int tSeverity=current->severity;
    double rN=0.0, lb, ub;
    
    if(current->isolated==0){
        if(current->indivType==0){ //individual
            if(current->detected==1){
                if(tSeverity==3){ //isolate only mild cases
                    rN=0.0;
                }
            }
        }else if(current->indivType==1){ //prison staff
            if(tSeverity==3){
                //isolate only mild cases
                rN=0.0;
            }
        }else if(current->indivType==2){ //hospital staff
            if(tSeverity==3){
                //isolate only mild cases
                rN=0.0;
            }
        }else if(current->indivType==3){ //essential visitor
            rN=0.0;
        }else if(current->indivType==4){ //family visitor
            rN=0.0;
        }
    }
    
    return rN;
}

double probEndIsolate(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN=0.0, lb, ub;

    if(current->isolated==1){
        if(current->indivType==0){ //individual
            if(current->COVID==0){
                rN=1.0;
            }
        }else if(current->indivType==1){ //prison staff
            if(current->COVID==0){
                rN=1.0;
            }
        }else if(current->indivType==2){ //hospital staff
            if(current->COVID==0){
                rN=1.0;
            }
        }else if(current->indivType==3){ //essential visitor
            rN=0.0;
        }else if(current->indivType==4){ //family visitor
            rN=0.0;
        }
    }
    
    return rN;
}

double probQuarantine(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probEndQuarantine(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probCohorting(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probEndCohorting(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probPPE(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probEndPPE(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probGroupVulnerable(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probEndGroupVulnerable(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to move to MAX
    sIndiv *current;
    current=*pTarget;
    double rN, lb, ub;
    
    if(current->indivType==0){ //individual
        
    }else if(current->indivType==1){ //prison staff
        
    }else if(current->indivType==2){ //hospital staff
        
    }else if(current->indivType==3){ //essential visitor
        
    }else if(current->indivType==4){ //family visitor
        
    }
    
    return rN;
}

double probOSTstart(gsl_rng **r, sIndiv **pTarget, int nOST, int nDailyOST, int nOSTDailyCap, int currDay){ //Probability to start DAA therapy (only if infected)
    sIndiv *current;
    current=*pTarget;
    
    if(currDay<=4380){
        if(current->risk==1||current->risk==2||current->risk==3||current->risk==7||current->risk==8||current->risk==9){ // only to opioid users
            if(nOST>0){
                if(nDailyOST<nOSTDailyCap){
                    if(current->OST==0){ //not yet on OST
                        return 0.9;
                    }else //not infected
                        return 0.0;
                }
            }
        }else
            return 0.0;
    }else if(currDay>4380){
        if(current->risk==1||current->risk==2||current->risk==3||current->risk==7||current->risk==8||current->risk==9){ // only to opioid users
            if(nOST>0){
                if(nDailyOST<nOSTDailyCap){
                    if(current->OST==0){ //not yet on OST
                        return 0.9;
                    }else //not infected
                        return 0.0;
                }
            }
        }else
            return 0.0;
    }
}

double probOSTstop(gsl_rng **r, sIndiv **pTarget){ //Probability to start DAA therapy (only if infected)
    sIndiv *current;
    current=*pTarget;
    
    if(current->OST==1){ //on OST
        return 0.049;
    }else //not infected
        return 0.0;
}

double probDAAstart(gsl_rng **r, sIndiv **pTarget, int nDAA, int currDay, int nDaily, int nDailyCap){ //Probability to start DAA therapy (only if infected)
    sIndiv *current;
    current=*pTarget;
    
    if(nDAA>0){
        if(nDaily<nDailyCap){
            if(current->metavir<5&&current->DAA==0){ //infected and is not yet on DAA
                return 0.9; //200/730
            }else //not infected
                return 0.0;
        }else if(nDaily>=nDailyCap)
            return 0.0;
    }else //no more DAA
        return 0.0;
}

double probDAAstop(gsl_rng **r, sIndiv **pTarget){ //Probability to start DAA therapy (only if infected)
    sIndiv *current;
    current=*pTarget;
    
    if(current->metavir<5&&current->DAA==1){ //infected and is not yet on DAA
        return 0.0;
    }else //not infected
        return 0.0;
}

double probDAAclear(gsl_rng **r, sIndiv **pTarget, int currDay){ //Probability to clear HCV via DAA (only if infected)
    sIndiv *current;
    int deltaDays, deltaWeeks;
    current=*pTarget;
    
    if(current->DAA==1){ //infected and is on DAA
        if(current->metavir<5){
        //compute for time undergoing DAA
        // 12 weeks=84 days
            deltaDays=currDay-(current->timeStartDAA);
            deltaWeeks=deltaDays/7;
            if(deltaWeeks<12){
                return 0.0792*deltaWeeks; //7.92 from 95% clearance in 12-week regimen (Martin, NK, et al.)
            }else
                return 0.0;
        }else //not infected
            return 0.0;
    }else
        return 0.0;
}

void age(sIndiv **pTargetCopy){
    sIndiv *current;
    current=*pTargetCopy;
    current->age++;
}

void removeIndiv(sIndiv **pTarget, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20]){
    //If target is the head
    sIndiv *temp, *current;
    int currLoc, currGroup, tempLoc;
    current=*pTarget;
    tempLoc=current->location;
    currGroup=current->group;
    
    int zone=current->zone;
    int area=current->locArea;
    int unit=current->locUnit;
    int cell=current->locCell;
    
    printf("removing individual.\n");
    if(tempLoc==1){ //objLoc==1||objLoc==6||objLoc==11
        if(current->indivType==0){
            currLoc=1;
            (*pMinCellArray)[zone][area][unit][cell]--;
        }else if(current->indivType==1){
            currLoc=2;
        }else if(current->indivType==2){
            currLoc=3;
        }else if(current->indivType==3){
            currLoc=4;
        }else if(current->indivType==4){
            currLoc=5;
        }
        
    }else if(tempLoc==2){
        if(current->indivType==0){
            currLoc=6;
            (*pMedCellArray)[zone][area][unit][cell]--;
        }else if(current->indivType==1){
            currLoc=7;
        }else if(current->indivType==2){
            currLoc=8;
        }else if(current->indivType==3){
            currLoc=9;
        }else if(current->indivType==4){
            currLoc=10;
        }
        
    }else if(tempLoc==3){
        if(current->indivType==0){
            currLoc=11;
            (*pMaxCellArray)[zone][area][unit][cell]--;
        }else if(current->indivType==1){
            currLoc=12;
        }else if(current->indivType==2){
            currLoc=13;
        }else if(current->indivType==3){
            currLoc=14;
        }else if(current->indivType==4){
            currLoc=15;
        }
        
    }
    //currLoc=current->location;
    printf("removing individual..\n");
    //subtract from current population
    //printf("Remove individual from loction %d group %d\n", currLoc, currGroup);
    (*pLocArray2)[currLoc][currGroup]--;
    printf("removing individual...\n");
    
    if(*pTarget==*pHeadCopy&&*pTarget!=*pTailCopy){
        *pHeadCopy=(*pTarget)->nextIndiv;
        temp=(*pTarget)->nextIndiv;
        temp->prevIndiv=NULL;
    }
    //If target is the tail
    else if(*pTarget==*pTailCopy&&*pTarget!=*pHeadCopy){
        temp=(*pTarget)->prevIndiv;
        temp->nextIndiv=NULL;
        *pTailCopy=(*pTarget)->prevIndiv;
    }
    //If target is both head and tails
    else if(*pTarget==*pTailCopy&&*pTarget==*pHeadCopy){
        *pHeadCopy=NULL;
        *pTailCopy=NULL;
    }
    //Otherwise, if target is just in the middle (neither head nor tail)
    else{
        temp=(*pTarget)->prevIndiv;
        temp->nextIndiv=(*pTarget)->nextIndiv;
        temp=(*pTarget)->nextIndiv;
        temp->prevIndiv=(*pTarget)->prevIndiv;
    }
    temp=*pTarget;
    *pTarget=(*pTarget)->nextIndiv;
    free(temp);
    
    printf("removed\n");
}

void changeRisk(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int newRisk){
    sIndiv *current;
    current=*pTargetCopy;
    int currGroup, currLoc, newGroup=99;
    
    currGroup=current->group; //Get current population group
    currLoc=current->location; //Current Location
    
        switch(newRisk){
            case 0: //stopInj
                current->injecting=0; //stop sharing and injecting opioid
                current->injFreq=0;
                current->sharing=0;
                current->shaFreq=0;
                current->injOpd=0;
                current->risk=0;
                //turn into non injecting
                if(currGroup<6){//if previously injecting, turn into non-injecting group
                    newGroup=currGroup+6;
                }else
                    newGroup=currGroup;
                break;
            case 1: //injLdOpNoSh
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=1;
                current->sharing=0;
                current->shaFreq=0;
                current->injOpd=1;
                current->risk=1;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 2: //injLdOpShLd
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=1;
                current->sharing=1;
                current->shaFreq=1;
                current->injOpd=1;
                current->risk=2;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                current->group=newGroup;
                break;
            case 3: //injLdOpShDm
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=1;
                current->sharing=1;
                current->shaFreq=2;
                current->injOpd=1;
                current->risk=3;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                
                break;
            case 4: //injLdNoNoSh
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=1;
                current->sharing=0;
                current->shaFreq=0;
                current->injOpd=0;
                current->risk=4;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 5: //injLdNoShLd
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=1;
                current->sharing=1;
                current->shaFreq=1;
                current->injOpd=0;
                current->risk=5;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 6: //injLdNoShDm
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=1;
                current->sharing=1;
                current->shaFreq=2;
                current->injOpd=0;
                current->risk=6;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 7: //injDmOpNoSh
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=2;
                current->sharing=0;
                current->shaFreq=0;
                current->injOpd=1;
                current->risk=7;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 8: //injDmOpShLd
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=2;
                current->sharing=1;
                current->shaFreq=1;
                current->injOpd=1;
                current->risk=8;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 9: //injDmOpShDm
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=2;
                current->sharing=1;
                current->shaFreq=2;
                current->injOpd=1;
                current->risk=9;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 10: //injDmNoNoSh
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=2;
                current->sharing=0;
                current->shaFreq=0;
                current->injOpd=0;
                current->risk=10;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 11: //injDmNoShLd
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=2;
                current->sharing=1;
                current->shaFreq=1;
                current->injOpd=0;
                current->risk=11;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
            case 12: //injDmNoShDm
                current->injecting=1; //stop sharing and injecting opioid
                current->injFreq=2;
                current->sharing=1;
                current->shaFreq=2;
                current->injOpd=0;
                current->risk=12;
                current->everIDU=1;
                //turn into injecting group if originally not injecting
                if(currGroup>5){//if previously not-injecting, turn into injecting group
                    newGroup=currGroup-6;
                }else //if already injecting, retain in current group
                    newGroup=currGroup;
                break;
        }

    
    current->group=newGroup;
    current->risk=newRisk;
    
    //update pLocArray
    (*pLocArray2)[currLoc][currGroup]--;
    (*pLocArray2)[currLoc][newGroup]++;
    printf("ex-group: %d . new group: %d\n",currGroup,newGroup);
    printf("Pop ex-group: %d . Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
}

void moveLocation(gsl_rng **r, sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int newLoc, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18]){ //loc=0,1,2
    sIndiv *current;
    current=*pTargetCopy;
    int currLoc, currGroup, tempLoc, targetZone, nZCtr, oldZone, oldZoneFlag=0;
    double normTotal, targetProb;
    gsl_rng *rCopy;
    rCopy=*r;
    
    int cellFlag=0;
    int cellCap=0;
    int zoneFlag=0;
    int nZoneLock=0;
    
    int mCtr1, mCtr2, mCtr3, mCtr4, mZoneCtr;
    
    double probsMaxZone[18], probsMedZone[11], probsMinZone[27], probsTruckNumber[20], probsSecZone[3];
    double probsMinArea[2], probsMedArea[2], probsMaxArea[4], probsMinUnit[6], probsMedUnit[4], probsMaxUnit[5], probsMinCell[13], probsMedCell[19], probsMaxCell[20];
    
    int pZone=current->zone;
    int pArea=current->locArea;
    int pUnit=current->locUnit;
    int pCell=current->locCell;
    
    //set new newLoc
    if(current->location==1){//prison 1 MIN *from 9010parameters.xlsx
        probsSecZone[0]=7.94E-04;
        probsSecZone[1]=3.82E-04;
        probsSecZone[2]=0.000225;
        
        normTotal=0.0;
        
        for(nZCtr=0; nZCtr<3; nZCtr++){
            normTotal=normTotal+probsSecZone[nZCtr];
        }
        
        for(nZCtr=0; nZCtr<3; nZCtr++){
            probsSecZone[nZCtr]=probsSecZone[nZCtr]/normTotal;
        }
        
        if(current->indivType==0){
            (*pMinCellArray)[pZone][pArea][pUnit][pCell]--;
        }
        
        newLoc=draw_multinom(&rCopy, 3, probsSecZone);
        newLoc++;
    }else if(current->location==2){//prison 2 MED
        probsSecZone[0]=0.003022;
        probsSecZone[1]=0.0004631;
        probsSecZone[2]=2.88E-04;
        
        normTotal=0.0;
        
        for(nZCtr=0; nZCtr<3; nZCtr++){
            normTotal=normTotal+probsSecZone[nZCtr];
        }
        
        for(nZCtr=0; nZCtr<3; nZCtr++){
            probsSecZone[nZCtr]=probsSecZone[nZCtr]/normTotal;
        }
        
        if(current->indivType==0){
            (*pMedCellArray)[pZone][pArea][pUnit][pCell]--;
        }
        
        newLoc=draw_multinom(&rCopy, 3, probsSecZone);
        newLoc++;
    }else if(current->location==3){//prison 3 MAX
        probsSecZone[0]=0.0005001;
        probsSecZone[1]=0.0005038;
        probsSecZone[2]=0.0008053;
        
        normTotal=0.0;
        
        for(nZCtr=0; nZCtr<3; nZCtr++){
            normTotal=normTotal+probsSecZone[nZCtr];
        }
        
        for(nZCtr=0; nZCtr<3; nZCtr++){
            probsSecZone[nZCtr]=probsSecZone[nZCtr]/normTotal;
        }
        
        if(current->indivType==0){
            (*pMaxCellArray)[pZone][pArea][pUnit][pCell]--;
        }
        
        newLoc=draw_multinom(&rCopy, 3, probsSecZone);
        newLoc++;
    }
    
    
    currGroup=current->group; //Current population group
    //currLoc=current->location; //Current Location
    tempLoc=current->location;
    
    if(tempLoc==1){ //objLoc==1||objLoc==6||objLoc==11
        if(current->indivType==0){
            currLoc=1;
        }else if(current->indivType==1){
            currLoc=2;
        }else if(current->indivType==2){
            currLoc=3;
        }else if(current->indivType==3){
            currLoc=4;
        }else if(current->indivType==4){
            currLoc=5;
        }
    }else if(tempLoc==2){
        if(current->indivType==0){
            currLoc=6;
        }else if(current->indivType==1){
            currLoc=7;
        }else if(current->indivType==2){
            currLoc=8;
        }else if(current->indivType==3){
            currLoc=9;
        }else if(current->indivType==4){
            currLoc=10;
        }
    }else if(tempLoc==3){
        if(current->indivType==0){
            currLoc=11;
        }else if(current->indivType==1){
            currLoc=12;
        }else if(current->indivType==2){
            currLoc=13;
        }else if(current->indivType==3){
            currLoc=14;
        }else if(current->indivType==4){
            currLoc=15;
        }
    }
    
    //printf("current location: %d -> new location: %d\n", current->location, newLoc);
    if(newLoc!=current->location){
        if(newLoc==0){
            if(current->indivType==0){ //don't count prison staff, health care staff, visitors, and family visitors coming out.
                printf("MOVE TO COMMUNITY!!!\n");
                (*pLocArray2)[0][currGroup]++; //move to at-risk Community, no need to subtract from current location as that is done in the function removeIndiv
            }//return 2;//individual went to community; raise dead flag
        }else{
            //printf("MOVE TO PRISON %d!!!\n", newLoc);
            (*pLocArray2)[currLoc][currGroup]--; //remove infected (at-risk) from current loc/group
            //printf("delete from %d\n", current->location);
            current->location=newLoc;
            if(newLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                if(current->indivType==0){
                    newLoc=1;
                }else if(current->indivType==1){
                    newLoc=2;
                }else if(current->indivType==2){
                    newLoc=3;
                }else if(current->indivType==3){
                    newLoc=4;
                }else if(current->indivType==4){
                    newLoc=5;
                }
                
                probsMinZone[0]=0.025;
                probsMinZone[1]=0.001;
                probsMinZone[2]=0.001;
                probsMinZone[3]=0.083;
                probsMinZone[4]=0.001;
                probsMinZone[5]=0.002;
                probsMinZone[6]=0.088;
                probsMinZone[7]=0.001;
                probsMinZone[8]=0.012;
                probsMinZone[9]=0.010;
                probsMinZone[10]=0.005;
                probsMinZone[11]=0.010;
                probsMinZone[12]=0.003;
                probsMinZone[13]=0.001;
                probsMinZone[14]=0.005;
                probsMinZone[15]=0.008;
                probsMinZone[16]=0.008;
                probsMinZone[17]=0.011;
                probsMinZone[18]=0.014;
                probsMinZone[19]=0.006;
                probsMinZone[20]=0.009;
                probsMinZone[21]=0.008;
                probsMinZone[22]=0.010;
                probsMinZone[23]=0.011;
                probsMinZone[24]=0.011;
                probsMinZone[25]=0.002;
                probsMinZone[26]=0.010;
                
                nZoneLock=0;
                for(mZoneCtr=0; mZoneCtr<27; mZoneCtr++){
                    if((*pLockdownMinZone)[mZoneCtr]==1){
                        probsMinZone[mZoneCtr]=0.0;
                        nZoneLock++;
                    }
                }
                
                if(nZoneLock!=27){
                    normTotal=0.0;
                    
                    for(nZCtr=0; nZCtr<27; nZCtr++){
                        normTotal=normTotal+probsMinZone[nZCtr];
                    }
                    
                    for(nZCtr=0; nZCtr<27; nZCtr++){
                        probsMinZone[nZCtr]=probsMinZone[nZCtr]/normTotal;
                    }
                    
                    current->zone=draw_multinom(&rCopy, 27, probsMinZone);
                    
                    cellFlag=0;
                    cellCap=0;
                    mCtr2=0;
                    //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                        while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                            mCtr3=0;
                            while(cellFlag==0||mCtr3<6){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                mCtr4=0;
                                while(cellFlag==0||mCtr4<13){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                    cellCap=(*pMinCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                    printf("cellCap: %d\n", cellCap);
                                    
                                    if(cellCap<2 && ((*pIsolateMinCell)[current->zone][mCtr2][mCtr3][mCtr4]==0) ){
                                        cellFlag=1;
                                        current->locArea=mCtr2;
                                        current->locUnit=mCtr3;
                                        current->locCell=mCtr4;
                                    }
                                    printf("cellFlag: %d\n", cellFlag);
                                    mCtr4++;
                                }
                                mCtr3++;
                            }
                            mCtr2++;
                        }
                    
                    if(cellFlag==0){
                        current->locArea=mCtr2-1;
                        current->locUnit=mCtr3-1;
                        current->locCell=mCtr4-1;
                    }
                    
                }
                
                if(current->indivType==0){
                    (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                }
                
                //current->zone=binZone; //put individual in an area
            }else if(newLoc==2){
                if(current->indivType==0){
                    newLoc=6;
                }else if(current->indivType==1){
                    newLoc=7;
                }else if(current->indivType==2){
                    newLoc=8;
                }else if(current->indivType==3){
                    newLoc=9;
                }else if(current->indivType==4){
                    newLoc=10;
                }
                
                probsMedZone[0]=0.067;
                probsMedZone[1]=0.006;
                probsMedZone[2]=0.002;
                probsMedZone[3]=0.007;
                probsMedZone[4]=0.013;
                probsMedZone[5]=0.018;
                probsMedZone[6]=0.043;
                probsMedZone[7]=0.025;
                probsMedZone[8]=0.006;
                probsMedZone[9]=0.023;
                probsMedZone[10]=0.007;
                
                nZoneLock=0;
                for(mZoneCtr=0; mZoneCtr<11; mZoneCtr++){
                    if((*pLockdownMedZone)[mZoneCtr]==1){
                        probsMedZone[mZoneCtr]=0.0;
                        nZoneLock++;
                    }
                }
                
                if(nZoneLock!=11){
                    normTotal=0.0;
                    
                    for(nZCtr=0; nZCtr<11; nZCtr++){
                        normTotal=normTotal+probsMedZone[nZCtr];
                    }
                    
                    for(nZCtr=0; nZCtr<11; nZCtr++){
                        probsMedZone[nZCtr]=probsMedZone[nZCtr]/normTotal;
                    }
                    
                    current->zone=draw_multinom(&rCopy, 11, probsMedZone);
                    
                    cellFlag=0;
                    cellCap=0;
                    mCtr2=0;
                    //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                        while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                            mCtr3=0;
                            while(cellFlag==0||mCtr3<4){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                mCtr4=0;
                                while(cellFlag==0||mCtr4<19){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                    cellCap=(*pMedCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                    printf("cellCap: %d\n", cellCap);
                                    
                                    if(cellCap<2 && ((*pIsolateMedCell)[current->zone][mCtr2][mCtr3][mCtr4]==0)){
                                        cellFlag=1;
                                        current->locArea=mCtr2;
                                        current->locUnit=mCtr3;
                                        current->locCell=mCtr4;
                                    }
                                    printf("cellFlag: %d\n", cellFlag);
                                    mCtr4++;
                                }
                                mCtr3++;
                            }
                            mCtr2++;
                        }
                    
                    if(cellFlag==0){
                        current->locArea=mCtr2-1;
                        current->locUnit=mCtr3-1;
                        current->locCell=mCtr4-1;
                    }
                }
                
                if(current->indivType==0){
                    (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                }
                
            }else if(newLoc==3){
                if(current->indivType==0){
                    newLoc=11;
                }else if(current->indivType==1){
                    newLoc=12;
                }else if(current->indivType==2){
                    newLoc=13;
                }else if(current->indivType==3){
                    newLoc=14;
                }else if(current->indivType==4){
                    newLoc=15;
                }
                
                probsMaxZone[0]=0.019;
                probsMaxZone[1]=0.001;
                probsMaxZone[2]=0.009;
                probsMaxZone[3]=0.018;
                probsMaxZone[4]=0.013;
                probsMaxZone[5]=0.032;
                probsMaxZone[6]=0.017;
                probsMaxZone[7]=0.002;
                probsMaxZone[8]=0.002;
                probsMaxZone[9]=0.013;
                probsMaxZone[10]=0.008;
                probsMaxZone[11]=0.098;
                probsMaxZone[12]=0.028;
                probsMaxZone[13]=0.063;
                probsMaxZone[14]=0.066;
                probsMaxZone[15]=0.012;
                probsMaxZone[16]=0.030;
                probsMaxZone[17]=0.023;
                
                nZoneLock=0;
                for(mZoneCtr=0; mZoneCtr<18; mZoneCtr++){
                    if((*pLockdownMaxZone)[mZoneCtr]==1){
                        probsMaxZone[mZoneCtr]=0.0;
                        nZoneLock++;
                    }
                }
                
                if(nZoneLock!=18){
                    normTotal=0.0;
                    
                    for(nZCtr=0; nZCtr<18; nZCtr++){
                        normTotal=normTotal+probsMaxZone[nZCtr];
                    }
                    
                    for(nZCtr=0; nZCtr<18; nZCtr++){
                        probsMedZone[nZCtr]=probsMaxZone[nZCtr]/normTotal;
                    }
                    
                    current->zone=draw_multinom(&rCopy, 18, probsMaxZone);
                    
                    cellFlag=0;
                    cellCap=0;
                    mCtr2=0;
                    //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                        while(cellFlag==0||mCtr2<4){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                            mCtr3=0;
                            while(cellFlag==0||mCtr3<5){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                                mCtr4=0;
                                while(cellFlag==0||mCtr4<20){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                    cellCap=(*pMaxCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                    printf("cellCap: %d\n", cellCap);
                                    
                                    if(cellCap<2 && ((*pIsolateMaxCell)[current->zone][mCtr2][mCtr3][mCtr4]==0) ){
                                        cellFlag=1;
                                        current->locArea=mCtr2;
                                        current->locUnit=mCtr3;
                                        current->locCell=mCtr4;
                                    }
                                    printf("cellFlag: %d\n", cellFlag);
                                    mCtr4++;
                                }
                                mCtr3++;
                            }
                            mCtr2++;
                        }
                    
                    if(cellFlag==0){
                        current->locArea=mCtr2-1;
                        current->locUnit=mCtr3-1;
                        current->locCell=mCtr4-1;
                    }
                }
                
                if(current->indivType==0){
                    (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
                }
    
            }
            (*pLocArray2)[newLoc][currGroup]++; //add infected (at-risk) from list
            //return 1; //changed location
        }
    }else{
        if(newLoc==1){
            
            oldZone=current->zone;
            
            probsMinZone[0]=0.025;
            probsMinZone[1]=0.001;
            probsMinZone[2]=0.001;
            probsMinZone[3]=0.083;
            probsMinZone[4]=0.001;
            probsMinZone[5]=0.002;
            probsMinZone[6]=0.088;
            probsMinZone[7]=0.001;
            probsMinZone[8]=0.012;
            probsMinZone[9]=0.010;
            probsMinZone[10]=0.005;
            probsMinZone[11]=0.010;
            probsMinZone[12]=0.003;
            probsMinZone[13]=0.001;
            probsMinZone[14]=0.005;
            probsMinZone[15]=0.008;
            probsMinZone[16]=0.008;
            probsMinZone[17]=0.011;
            probsMinZone[18]=0.014;
            probsMinZone[19]=0.006;
            probsMinZone[20]=0.009;
            probsMinZone[21]=0.008;
            probsMinZone[22]=0.010;
            probsMinZone[23]=0.011;
            probsMinZone[24]=0.011;
            probsMinZone[25]=0.002;
            probsMinZone[26]=0.010;
            
            nZoneLock=0;
            for(mZoneCtr=0; mZoneCtr<27; mZoneCtr++){
                if((*pLockdownMinZone)[mZoneCtr]==1){
                    probsMinZone[mZoneCtr]=0.0;
                    nZoneLock++;
                }
            }
            
            if(nZoneLock!=27){
                //targetZone=current->zone;
                normTotal=0.0;
                //probsMinZone[targetZone]=0.0;
                
                for(nZCtr=0; nZCtr<27; nZCtr++){
                    normTotal=normTotal+probsMinZone[nZCtr];
                }
                
                for(nZCtr=0; nZCtr<27; nZCtr++){
                    probsMinZone[nZCtr]=probsMinZone[nZCtr]/normTotal;
                }
                
                current->zone=draw_multinom(&rCopy, 27, probsMinZone);
                
                if(current->zone==oldZone){
                    oldZoneFlag=1;
                }
                
                cellFlag=0;
                cellCap=0;
                mCtr2=0;
                //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                    while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                        mCtr3=0;
                        while(cellFlag==0||mCtr3<6){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                            mCtr4=0;
                            while(cellFlag==0||mCtr4<13){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                cellCap=(*pMinCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                printf("cellCap: %d\n", cellCap);
                                
                                if(cellCap<2 && ((*pIsolateMinCell)[current->zone][mCtr2][mCtr3][mCtr4]==0) ){
                                    cellFlag=1;
                                    current->locArea=mCtr2;
                                    current->locUnit=mCtr3;
                                    current->locCell=mCtr4;
                                }
                                printf("cellFlag: %d\n", cellFlag);
                                mCtr4++;
                            }
                            mCtr3++;
                        }
                        mCtr2++;
                    }
                 
                if(cellFlag==0){
                    current->locArea=mCtr2-1;
                    current->locUnit=mCtr3-1;
                    current->locCell=mCtr4-1;
                }
            }
            
            if(current->indivType==0){
                (*pMinCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
            }
                        
        }else if(newLoc==2){
            
            oldZone=current->zone;
            
            probsMedZone[0]=0.067;
            probsMedZone[1]=0.006;
            probsMedZone[2]=0.002;
            probsMedZone[3]=0.007;
            probsMedZone[4]=0.013;
            probsMedZone[5]=0.018;
            probsMedZone[6]=0.043;
            probsMedZone[7]=0.025;
            probsMedZone[8]=0.006;
            probsMedZone[9]=0.023;
            probsMedZone[10]=0.007;
            
            nZoneLock=0;
            for(mZoneCtr=0; mZoneCtr<11; mZoneCtr++){
                if((*pLockdownMedZone)[mZoneCtr]==1){
                    probsMedZone[mZoneCtr]=0.0;
                    nZoneLock++;
                }
            }
            
            if(nZoneLock!=11){
                //targetZone=current->zone;
                normTotal=0.0;
                //probsMedZone[targetZone]=0.0;
                
                for(nZCtr=0; nZCtr<11; nZCtr++){
                    normTotal=normTotal+probsMedZone[nZCtr];
                }
                
                for(nZCtr=0; nZCtr<11; nZCtr++){
                    probsMedZone[nZCtr]=probsMedZone[nZCtr]/normTotal;
                }
                
                current->zone=draw_multinom(&rCopy, 11, probsMedZone);
                
                if(current->zone==oldZone){
                    oldZoneFlag=1;
                }
                
                cellFlag=0;
                cellCap=0;
                mCtr2=0;
                //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                    while(cellFlag==0||mCtr2<2){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                        mCtr3=0;
                        while(cellFlag==0||mCtr3<4){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                            mCtr4=0;
                            while(cellFlag==0||mCtr4<19){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                cellCap=(*pMedCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                printf("cellCap: %d\n", cellCap);
                                
                                if(cellCap<2 && ((*pIsolateMedCell)[current->zone][mCtr2][mCtr3][mCtr4]==0) ){
                                    cellFlag=1;
                                    current->locArea=mCtr2;
                                    current->locUnit=mCtr3;
                                    current->locCell=mCtr4;
                                }
                                printf("cellFlag: %d\n", cellFlag);
                                mCtr4++;
                            }
                            mCtr3++;
                        }
                        mCtr2++;
                    }
                
                //current->locCell=draw_multinom(&rCopy, 19, probsMedCell);
                
                if(cellFlag==0){
                    current->locArea=mCtr2-1;
                    current->locUnit=mCtr3-1;
                    current->locCell=mCtr4-1;
                }
            }
            
            if(current->indivType==0){
                (*pMedCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
            }
            
        }else if(newLoc==3){
            
            oldZone=current->zone;
            
            probsMaxZone[0]=0.019;
            probsMaxZone[1]=0.001;
            probsMaxZone[2]=0.009;
            probsMaxZone[3]=0.018;
            probsMaxZone[4]=0.013;
            probsMaxZone[5]=0.032;
            probsMaxZone[6]=0.017;
            probsMaxZone[7]=0.002;
            probsMaxZone[8]=0.002;
            probsMaxZone[9]=0.013;
            probsMaxZone[10]=0.008;
            probsMaxZone[11]=0.098;
            probsMaxZone[12]=0.028;
            probsMaxZone[13]=0.063;
            probsMaxZone[14]=0.066;
            probsMaxZone[15]=0.012;
            probsMaxZone[16]=0.030;
            probsMaxZone[17]=0.023;
            
            nZoneLock=0;
            for(mZoneCtr=0; mZoneCtr<18; mZoneCtr++){
                if((*pLockdownMinZone)[mZoneCtr]==1){
                    probsMinZone[mZoneCtr]=0.0;
                    nZoneLock++;
                }
            }
            
            if(nZoneLock!=18){
                //targetZone=current->zone;
                normTotal=0.0;
                //probsMaxZone[targetZone]=0.0;
                
                for(nZCtr=0; nZCtr<18; nZCtr++){
                    normTotal=normTotal+probsMaxZone[nZCtr];
                }
                
                for(nZCtr=0; nZCtr<18; nZCtr++){
                    probsMaxZone[nZCtr]=probsMaxZone[nZCtr]/normTotal;
                }
                
                current->zone=draw_multinom(&rCopy, 18, probsMaxZone);
                
                if(current->zone==oldZone){
                    oldZoneFlag=1;
                }
                
                cellFlag=0;
                cellCap=0;
                mCtr2=0;
                //for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
                    while(cellFlag==0||mCtr2<4){//(mcCtr2=0; mcCtr2<2; mcCtr2++)//min med max
                        mCtr3=0;
                        while(cellFlag==0||mCtr3<5){//}for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                            mCtr4=0;
                            while(cellFlag==0||mCtr4<20){ //for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                                cellCap=(*pMaxCellArray)[current->zone][mCtr2][mCtr3][mCtr4];
                                printf("cellCap: %d\n", cellCap);
                                
                                if(cellCap<2 && ((*pIsolateMaxCell)[current->zone][mCtr2][mCtr3][mCtr4]==0) ){
                                    cellFlag=1;
                                    current->locArea=mCtr2;
                                    current->locUnit=mCtr3;
                                    current->locCell=mCtr4;
                                }
                                printf("cellFlag: %d\n", cellFlag);
                                mCtr4++;
                            }
                            mCtr3++;
                        }
                        mCtr2++;
                    }
                
                if(cellFlag==0){
                    current->locArea=mCtr2-1;
                    current->locUnit=mCtr3-1;
                    current->locCell=mCtr4-1;
                }
            }
            
            if(current->indivType==0){
                (*pMaxCellArray)[current->zone][current->locArea][current->locUnit][current->locCell]++;
            }
        }
    }
    
    if(oldZoneFlag!=1){
    //set moving state
    current->moving=1;
    //set moving truck number
    probsTruckNumber[0]=0.05;
    probsTruckNumber[1]=0.05;
    probsTruckNumber[2]=0.05;
    probsTruckNumber[3]=0.05;
    probsTruckNumber[4]=0.05;
    probsTruckNumber[5]=0.05;
    probsTruckNumber[6]=0.05;
    probsTruckNumber[7]=0.05;
    probsTruckNumber[8]=0.05;
    probsTruckNumber[9]=0.05;
    probsTruckNumber[10]=0.05;
    probsTruckNumber[11]=0.05;
    probsTruckNumber[12]=0.05;
    probsTruckNumber[13]=0.05;
    probsTruckNumber[14]=0.05;
    probsTruckNumber[15]=0.05;
    probsTruckNumber[16]=0.05;
    probsTruckNumber[17]=0.05;
    probsTruckNumber[18]=0.05;
    probsTruckNumber[19]=0.05;

    current->truckNumber=draw_multinom(&rCopy, 20, probsTruckNumber);
    }
    
    //return 0; //no change in location
    printf("ex-location: %d . new location: %d\n",currLoc,newLoc);
    printf("Pop ex-loc: %d . Pop new loc: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[newLoc][currGroup]);
}

void outsideVisit(gsl_rng **r, sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20]){ //loc=0,1,2
    sIndiv *current;
    current=*pTargetCopy;
    int currLoc, currGroup, tempLoc;
    gsl_rng *rCopy;
    rCopy=*r;
    int pZone=current->zone;
    int pArea=current->locArea;
    int pUnit=current->locUnit;
    int pCell=current->locCell;
    
/*    if(current->location==1){//prison 1 MIN *from 9010parameters.xlsx
        
        (*pMinCellArray)[pZone][pArea][pUnit][pCell]--;
        
    }else if(current->location==2){//prison 2 MED
        
        (*pMedCellArray)[pZone][pArea][pUnit][pCell]--;
        
    }else if(current->location==3){//prison 3 MAX
        
        (*pMaxCellArray)[pZone][pArea][pUnit][pCell]--;
        
    } */
    
    double probsCourtZone[38];
    
    current->court=1;
    probsCourtZone[0]=0.026;
    probsCourtZone[1]=0.026;
    probsCourtZone[2]=0.026;
    probsCourtZone[3]=0.026;
    probsCourtZone[4]=0.026;
    probsCourtZone[5]=0.026;
    probsCourtZone[6]=0.026;
    probsCourtZone[7]=0.026;
    probsCourtZone[8]=0.026;
    probsCourtZone[9]=0.026;
    probsCourtZone[10]=0.026;
    probsCourtZone[11]=0.026;
    probsCourtZone[12]=0.026;
    probsCourtZone[13]=0.026;
    probsCourtZone[14]=0.026;
    probsCourtZone[15]=0.026;
    probsCourtZone[16]=0.026;
    probsCourtZone[17]=0.026;
    probsCourtZone[18]=0.026;
    probsCourtZone[19]=0.026;
    probsCourtZone[20]=0.026;
    probsCourtZone[21]=0.026;
    probsCourtZone[22]=0.026;
    probsCourtZone[23]=0.026;
    probsCourtZone[24]=0.026;
    probsCourtZone[25]=0.026;
    probsCourtZone[26]=0.026;
    probsCourtZone[27]=0.026;
    probsCourtZone[28]=0.026;
    probsCourtZone[29]=0.026;
    probsCourtZone[30]=0.026;
    probsCourtZone[31]=0.026;
    probsCourtZone[32]=0.026;
    probsCourtZone[33]=0.026;
    probsCourtZone[34]=0.026;
    probsCourtZone[35]=0.026;
    probsCourtZone[36]=0.026;
    probsCourtZone[37]=0.026;
    current->courtNumber=draw_multinom(&rCopy, 38, probsCourtZone);
    
    //set units
    current->courtCell=generateAge(1, 8);
}

void moveHospital(sIndiv **pTargetCopy, int hospitalType){ //1-field hospital 2-community hospital
    sIndiv *current;
    current=*pTargetCopy;
    
    current->hospitalised=hospitalType;
    if(current->isolated==1){
        current->isolated=0;
    }
}

void returnPrison(sIndiv **pTargetCopy, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20]){ //1-field hospital 2-community hospital
    sIndiv *current;
    current=*pTargetCopy;
    
    int pZone=current->zone;
    int pArea=current->locArea;
    int pUnit=current->locUnit;
    int pCell=current->locCell;
    
    if(current->hospitalised==1){
        current->hospitalised=0;
    }else if(current->court==1){
        current->court=0;
        current->courtNumber=99;
    }else if(current->moving==1){
        current->moving=0;
        current->truckNumber=99;
    }
}

void test(sIndiv **pTargetCopy, int currDay){ //loc=0,1,2
    sIndiv *target;
    target=*pTargetCopy;
    
    target->tested=1;
    target->testDate=currDay;
    target->testDaysBeforeResult=5;
}

void testResult(sIndiv **pTargetCopy, int currDay){ //loc=0,1,2
    sIndiv *target;
    target=*pTargetCopy;
    
    if(target->COVID==1){
        target->detected=1;
    }
    target->tested=0;
    target->testDate=0;
    target->testDaysBeforeResult=0;
}

int rapidTest(gsl_rng **r, sIndiv **pTargetCopy, int currDay){ //loc=0,1,2
    sIndiv *target;
    target=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    double probTest[2];
    int binTest;
    int rN=0;
    
    probTest[0]=1-0.50; //not detected
    probTest[1]=0.50; //detected
    binTest=draw_multinom(&rCopy, 2, probTest);
    
    if(binTest==1){ //test success
        if(target->COVID==1){
            target->detected=1;
            target->testResultTrue=1;
            target->TP=1;
            //TP
        }else{
            target->detected=1;
            target->testResultTrue=0;
            target->FP=1;
            //FP
        }
        //target->tested=1;
        target->testDate=currDay;
        //target->testDaysBeforeResult=0;
        rN=1; //test was positive
    }else{
        if(target->COVID==1){
            target->detected=0;
            target->testResultTrue=0;
            target->FN=1;
            //FNffff
        }else{
            target->detected=0;
            target->testResultTrue=1;
            target->TN=1;
            //TN
        }
        target->detected=0;
        target->testResultTrue=1;
        target->tested=1;
        target->testDate=currDay;
        //target->testDaysBeforeResult=0;
    }
    
    return rN;
}

void thermalTest(sIndiv **pTargetCopy){ //loc=0,1,2
    sIndiv *target;
    target=*pTargetCopy;
    
    if(target->COVID==1){
        if(target->severity>=3&&target->severity<=5){ //only detects those with temperature
            target->detected=1;
        }
    }
}

int infect(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*zoneMinCIArray)[30][7], int (*zoneMedCIArray)[30][7], int (*zoneMaxCIArray)[30][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int currDay){
    sIndiv *target, *tempCurrent;
    target=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rStatus, binInf;
    //int zone=target->zone;
    int targetLocation=target->location, aCtr, iTarget, currLoc, currGroup, newGroup=99, eFlag=0;
    //int totalAtRisk=(*pLocArray2)[target->location][3]+(*pLocArray2)[target->location][5]+(*pLocArray2)[target->location][9]+(*pLocArray2)[target->location][11]; //add all at risk except those with immunity
    //sIndiv *refInfectList[totalAtRisk];
    //double probInfectList[totalAtRisk];
    double eventDraw[2]; //proability of recipient to get infected
    double probInf[2]; //probability of source infecting recipient using per event probability
    int success, eventDecision, susCtr;
    int targetArrayLoc;
    int contactIMax=0, contactIMed=0, contactIMin=0, contactPSMax=0, contactPSMed=0, contactPSMin=0, contactHSMax=0, contactHSMed=0, contactHSMin=0, contactEVMax=0, contactEVMed=0, contactEVMin=0, contactFVMax=0, contactFVMed=0, contactFVMin=0;
    int contactCapacityI=0, contactCapacityPS=0, contactCapacityHS=0, contactCapacityEV=0, contactCapacityFV=0;
    int targetAge, inmateFlag=0;
    int infectedCountI=0, infectedCountPS=0, infectedCountHS=0, infectedCountEV=0, infectedCountFV=0;
    double probsRisk[2];
    double ppeEffect;
    
    printf("going through list of individuals\n");
    success=0;
    //probInf[0]=0.0; probInf[1]=0.0;
    eventDraw[0]=0; eventDraw[1]=0;
    //aCtr=0;
    //Go through list of individuals
    
    if(target->indivType==0){ //Inmate
        contactIMax=generateAge(175, 388)+generateAge(30, 76);//gsl_ran_flat(*r, 36, 81);//Interaction within Units/Pods only
        contactIMed=generateAge(107, 238)+generateAge(14, 35);//gsl_ran_flat(*r, 107, 130);//Interaction within Units/Pods
        contactIMin=generateAge(88, 195)+generateAge(15, 38);//gsl_ran_flat(*r, 87, 107);//Interaction within Units/Pods
        contactPSMax=2;
        contactPSMed=1;
        contactPSMin=2;
        probsRisk[0]=1-0.17; //not a patient
        probsRisk[1]=0.17;
        contactHSMax=draw_multinom(&rCopy, 2, probsRisk);
        probsRisk[0]=1-0.07; //not a patient
        probsRisk[1]=0.07;
        contactHSMed=draw_multinom(&rCopy, 2, probsRisk);
        probsRisk[0]=1-0.09; //not a patient
        probsRisk[1]=0.09;
        contactHSMin=draw_multinom(&rCopy, 2, probsRisk);
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=1;
        }else if(targetLocation==2){
            targetArrayLoc=6;
        }else if(targetLocation==3){
            targetArrayLoc=11;
        }
        inmateFlag=1;
    }else if(target->indivType==1){//ps
        contactIMax=5;
        contactIMed=5;
        contactIMin=5;
        contactPSMax=39;
        contactPSMed=14;
        contactPSMin=12;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=2;
        }else if(targetLocation==2){
            targetArrayLoc=7;
        }else if(targetLocation==3){
            targetArrayLoc=12;
        }
    }else if(target->indivType==2){ //hs
        contactIMax=8; //8
        contactIMed=6; //6
        contactIMin=7; //7
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=9; //9
        contactHSMed=7; //7
        contactHSMin=5; //5
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=3;
        }else if(targetLocation==2){
            targetArrayLoc=8;
        }else if(targetLocation==3){
            targetArrayLoc=13;
        }
    }else if(target->indivType==3){ //ev
        contactIMax=1;
        contactIMed=1;
        contactIMin=1;
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=4;
        }else if(targetLocation==2){
            targetArrayLoc=9;
        }else if(targetLocation==3){
            targetArrayLoc=14;
        }
    }else if(target->indivType==4){ //fv
        contactIMax=1;
        contactIMed=1;
        contactIMin=1;
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=5;
        }else if(targetLocation==2){
            targetArrayLoc=10;
        }else if(targetLocation==3){
            targetArrayLoc=15;
        }
    }
    
    if(targetLocation==1){
        contactCapacityI=contactIMin;
        contactCapacityPS=contactPSMin;
        contactCapacityHS=contactHSMin;
        contactCapacityEV=contactEVMin;
        contactCapacityFV=contactFVMin;
    }else if(targetLocation==2){
        contactCapacityI=contactIMed;
        contactCapacityPS=contactPSMed;
        contactCapacityHS=contactHSMed;
        contactCapacityEV=contactEVMed;
        contactCapacityFV=contactFVMed;
    }else if(targetLocation==3){
        contactCapacityI=contactIMax;
        contactCapacityPS=contactPSMax;
        contactCapacityHS=contactHSMax;
        contactCapacityEV=contactEVMax;
        contactCapacityFV=contactFVMax;
    }
    
    tempCurrent=*pHeadCopy; //Target points to the individual at the beginning of the list
    //printf("total at-risk agents in prison %d: %d\n", targetLocation, totalAtRisk);
    
    while(tempCurrent!=NULL){
        //add if injecting
        //FIX location matching to infect other populations
        //
        if((targetLocation==tempCurrent->location)&&(tempCurrent!=target)){//susceptible if same location and not same indiv



            if(tempCurrent->zone==target->zone){  //tempCurrent->zone==target->zone
                if(tempCurrent->hospitalised==0&&tempCurrent->isolated==0&&tempCurrent->court==0){ //not in any hospital
                    //if(tempCurrent->COVID==0){ //not yet been infected and not currently infected
                    if(tempCurrent->indivType==0){
                        if(infectedCountI<contactCapacityI){
                            
                            if(target->indivType==0){
                                if(tempCurrent->locArea==target->locArea){
                                    //herd immunity
                                    if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                        probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[0]=1-probInf[1];
                                    }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                        ppeEffect=0.67;
                                        if(target->severity==4||target->severity==5){ //mod or severe
                                            probInf[1]=0.01*2;//0.001;//0.05;//target = probability of getting infected with COVID-19
                                            probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                        }else{ //mild and asymptomatic
                                            probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//0.05;//probability of getting infected with COVID-19
                                            probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                        }
                                        probInf[0]=1-probInf[1];
                                    }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                        probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[0]=1-probInf[1];
                                    }
                                    
                                    printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                                    
                                    binInf=draw_multinom(&rCopy, 2, probInf);
                                    //printf("Infection commencing\n");
                                    printf("Infect Y or N: %d\n", binInf);
                                    
                                    if(binInf==1){ //Target is to be infected
                                        printf("Infection commencing\n");
                                        
                                        tempCurrent->COVID=1; // infected!!
                                        tempCurrent->severity=0;
                                        currLoc=tempCurrent->location;
                                        tempCurrent->placeInfected=currLoc;
                                        targetAge=tempCurrent->age;
                                        
                                        if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                            (*newInfectedInmatesMin)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=1;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMinCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMinCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMinCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMinCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMinCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMinCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMinCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=2;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=3;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=4;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=5;
                                            }
                                        }else if(currLoc==2){
                                            (*newInfectedInmatesMed)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=6;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMedCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMedCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMedCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMedCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMedCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMedCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMedCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=7;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=8;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=9;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=10;
                                            }
                                        }else if(currLoc==3){
                                            (*newInfectedInmatesMax)++;
                                            
                                            if(tempCurrent->indivType==0){
                                                currLoc=11;
                                                
                                                if(targetAge<=19){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                                    (*totalAge0)++;
                                                }else if(targetAge>=20&&targetAge<=44){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                                    (*totalAge1)++;
                                                }else if(targetAge>=45&&targetAge<=54){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                                    (*totalAge2)++;
                                                }else if(targetAge>=55&&targetAge<=64){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                                    (*totalAge3)++;
                                                }else if(targetAge>=65&&targetAge<=74){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                                    (*totalAge4)++;
                                                }else if(targetAge>=75&&targetAge<=84){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                                    (*totalAge5)++;
                                                }else if(targetAge>=85){
                                                    (*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                                    (*totalAge6)++;
                                                }
                                                
                                            }else if(tempCurrent->indivType==1){
                                                currLoc=12;
                                            }else if(tempCurrent->indivType==2){
                                                currLoc=13;
                                            }else if(tempCurrent->indivType==3){
                                                currLoc=14;
                                            }else if(tempCurrent->indivType==4){
                                                currLoc=15;
                                            }
                                        }
                                        
                                        if(tempCurrent->age<=19){ //means infected in a prison
                                            (*exp19)++;
                                        }else if(tempCurrent->age>=20&&tempCurrent->age<=44){
                                            (*exp20to44)++;
                                        }else if(tempCurrent->age>=45&&tempCurrent->age<=54){
                                            (*exp45to54)++;
                                        }else if(tempCurrent->age>=55&&tempCurrent->age<=64){
                                            (*exp55to64)++;
                                        }else if(tempCurrent->age>=65&&tempCurrent->age<=74){
                                            (*exp65to74)++;
                                        }else if(tempCurrent->age>=75&&tempCurrent->age<=84){
                                            (*exp75to84)++;
                                        }else if(tempCurrent->age>=85){
                                            (*exp85)++;
                                        }
                                        
                                        tempCurrent->timeOfInfection=currDay; //Record date of infection
                                        ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                        currGroup=tempCurrent->group;
                                        
                                        (*pLocArray2)[currLoc][currGroup]--;
                                        
                                        //system("pause");
                                        switch(currGroup){
                                            case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                                break;
                                            case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                                break;
                                            case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                                break;
                                            case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=0;
                                                break;
                                            case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                                break;
                                            case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=1;
                                                break;
                                            case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                                break;
                                            case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                                break;
                                            case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                                break;
                                            case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=6;
                                                break;
                                            case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                                printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                                break;
                                            case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                                newGroup=7;
                                                break;
                                                
                                        }
                                        
                                        //if(currGroup!=newGroup){
                                        tempCurrent->group=newGroup;
                                        //}
                                        
                                        printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                        success++;
                                        
                                        (*pLocArray2)[currLoc][newGroup]++;
                                        printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                        
                                        (*newInfectedInmates)++;
                                    }
                                    
                                    infectedCountI++;
                                }
                            }else{
                                if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                    probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[0]=1-probInf[1];
                                }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                    ppeEffect=0.67;
                                    if(target->severity==4||target->severity==5){
                                        probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                    }else{
                                        probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                        probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                    }
                                    probInf[0]=1-probInf[1];
                                }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                    probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[0]=1-probInf[1];
                                }
                                
                                printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                                
                                binInf=draw_multinom(&rCopy, 2, probInf);
                                //printf("Infection commencing\n");
                                printf("Infect Y or N: %d\n", binInf);
                                
                                if(binInf==1){ //Target is to be infected
                                    printf("Infection commencing\n");
                                    
                                    tempCurrent->COVID=1; // infected!!
                                    tempCurrent->severity=0;
                                    currLoc=tempCurrent->location;
                                    tempCurrent->placeInfected=currLoc;
                                    targetAge=tempCurrent->age;
                                    
                                    if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                        (*newInfectedInmatesMin)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=1;
                                            
                                            if(targetAge<=19){
                                                (*zoneMinCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMinCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMinCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMinCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMinCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMinCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMinCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=2;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=3;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=4;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=5;
                                        }
                                    }else if(currLoc==2){
                                        (*newInfectedInmatesMed)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=6;
                                            
                                            if(targetAge<=19){
                                                (*zoneMedCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMedCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMedCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMedCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMedCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMedCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMedCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=7;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=8;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=9;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=10;
                                        }
                                    }else if(currLoc==3){
                                        (*newInfectedInmatesMax)++;
                                        
                                        if(tempCurrent->indivType==0){
                                            currLoc=11;
                                            
                                            if(targetAge<=19){
                                                (*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                                (*totalAge0)++;
                                            }else if(targetAge>=20&&targetAge<=44){
                                                (*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                                (*totalAge1)++;
                                            }else if(targetAge>=45&&targetAge<=54){
                                                (*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                                (*totalAge2)++;
                                            }else if(targetAge>=55&&targetAge<=64){
                                                (*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                                (*totalAge3)++;
                                            }else if(targetAge>=65&&targetAge<=74){
                                                (*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                                (*totalAge4)++;
                                            }else if(targetAge>=75&&targetAge<=84){
                                                (*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                                (*totalAge5)++;
                                            }else if(targetAge>=85){
                                                (*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                                (*totalAge6)++;
                                            }
                                            
                                        }else if(tempCurrent->indivType==1){
                                            currLoc=12;
                                        }else if(tempCurrent->indivType==2){
                                            currLoc=13;
                                        }else if(tempCurrent->indivType==3){
                                            currLoc=14;
                                        }else if(tempCurrent->indivType==4){
                                            currLoc=15;
                                        }
                                    }
                                    
                                    tempCurrent->timeOfInfection=currDay; //Record date of infection
                                    ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                    currGroup=tempCurrent->group;
                                    
                                    (*pLocArray2)[currLoc][currGroup]--;
                                    
                                    //system("pause");
                                    switch(currGroup){
                                        case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                            break;
                                        case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                            break;
                                        case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                            break;
                                        case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=0;
                                            break;
                                        case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                            break;
                                        case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=1;
                                            break;
                                        case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                            break;
                                        case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                            break;
                                        case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                            break;
                                        case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=6;
                                            break;
                                        case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                            printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                            break;
                                        case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                            newGroup=7;
                                            break;
                                            
                                    }
                                    
                                    //if(currGroup!=newGroup){
                                    tempCurrent->group=newGroup;
                                    //}
                                    
                                    printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                    success++;
                                    
                                    (*pLocArray2)[currLoc][newGroup]++;
                                    printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                    
                                    (*newInfectedInmates)++;
                                }
                                
                                infectedCountI++;
                            }
                                

                        }
                    }else if(tempCurrent->indivType==1){
                        if(infectedCountPS<contactCapacityPS){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.67;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedPSMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedPSMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedPSMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                
                                (*newInfectedPS)++;
                            }
                            
                            infectedCountPS++;
                        }
                    }else if(tempCurrent->indivType==2){
                        if(infectedCountHS<contactCapacityHS){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.67;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedHSMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedHSMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedHSMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedHS)++;
                            }
                            
                            infectedCountHS++;
                        }
                    }else if(tempCurrent->indivType==3){
                        if(infectedCountEV<contactCapacityEV){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.67;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedEVMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedEVMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedEVMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedEV)++;
                            }
                            
                            infectedCountEV++;
                        }
                    }else if(tempCurrent->indivType==4){
                        if(infectedCountFV<contactCapacityFV){
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.67;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=(gsl_ran_flat(*r, 0.02, 0.05))*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedFVMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedFVMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedFVMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                (*newInfectedFV)++;
                            }
                            
                            infectedCountFV++;
                        }
                    }
                    
                    //}//covid
                }//hospitalised
            }//same zone
        }//same location
        tempCurrent=tempCurrent->nextIndiv;
    }
    printf("Number of new infections for this instance: %d", success);
    return success;
}

int infectOutside(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*courtCIArray)[38][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int currDay){
    sIndiv *target, *tempCurrent;
    target=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rStatus, binInf;
    //int zone=target->zone;
    int targetLocation=target->location, aCtr, iTarget, currLoc, currGroup, newGroup=99, eFlag=0;
    //int totalAtRisk=(*pLocArray2)[target->location][3]+(*pLocArray2)[target->location][5]+(*pLocArray2)[target->location][9]+(*pLocArray2)[target->location][11]; //add all at risk except those with immunity
    //sIndiv *refInfectList[totalAtRisk];
    //double probInfectList[totalAtRisk];
    double eventDraw[2]; //proability of recipient to get infected
    double probInf[2]; //probability of source infecting recipient using per event probability
    int success, eventDecision, susCtr;
    int targetArrayLoc;
    int contactIMax=0, contactIMed=0, contactIMin=0, contactPSMax=0, contactPSMed=0, contactPSMin=0, contactHSMax=0, contactHSMed=0, contactHSMin=0, contactEVMax=0, contactEVMed=0, contactEVMin=0, contactFVMax=0, contactFVMed=0, contactFVMin=0, contactMax=0;
    int contactCapacityI=0, contactCapacityPS=0, contactCapacityHS=0, contactCapacityEV=0, contactCapacityFV=0;
    int targetAge, tCourtNumber;
    int infectedCountI=0, infectedCountPS=0, infectedCountHS=0, infectedCountEV=0, infectedCountFV=0;
    printf("going through list of individuals\n");
    success=0;
    //probInf[0]=0.0; probInf[1]=0.0;
    eventDraw[0]=0; eventDraw[1]=0;
    double ppeEffect;
    //aCtr=0;
    //Go through list of individuals
    
    /*
    if(target->indivType==0){ //Inmate
        contactC1=generateAge(36, 81);//gsl_ran_flat(*r, 36, 81);//Interaction within Units/Pods only
        contactC2=generateAge(35, 79);//gsl_ran_flat(*r, 107, 130);//Interaction within Units/Pods
        contactIMin=generateAge(26, 59);//gsl_ran_flat(*r, 87, 107);//Interaction within Units/Pods
        contactPSMax=0;
        contactPSMed=0;
        contactPSMin=0;
        contactHSMax=0;
        contactHSMed=0;
        contactHSMin=0;
        contactEVMax=0;
        contactEVMed=0;
        contactEVMin=0;
        contactFVMax=0;
        contactFVMed=0;
        contactFVMin=0;
        if(targetLocation==1){
            targetArrayLoc=1;
        }else if(targetLocation==2){
            targetArrayLoc=6;
        }else if(targetLocation==3){
            targetArrayLoc=11;
        }
    }
    */
    tCourtNumber=target->courtNumber;
    switch(tCourtNumber){
        case 0:
            contactMax=4;
            break;
        case 1:
            contactMax=3;
            break;
        case 2:
            contactMax=9;
            break;
        case 3:
            contactMax=3;
            break;
        case 4:
            contactMax=3;
            break;
        case 5:
            contactMax=4;
            break;
        case 6:
            contactMax=6;
            break;
        case 7:
            contactMax=3;
            break;
        case 8:
            contactMax=3;
            break;
        case 9:
            contactMax=4;
            break;
        case 10:
            contactMax=2;
            break;
        case 11:
            contactMax=4;
            break;
        case 12:
            contactMax=2;
            break;
        case 13:
            contactMax=6;
            break;
        case 14:
            contactMax=6;
            break;
        case 15:
            contactMax=2;
            break;
        case 16:
            contactMax=2;
            break;
        case 17:
            contactMax=5;
            break;
        case 18:
            contactMax=4;
            break;
        case 19:
            contactMax=3;
            break;
        case 20:
            contactMax=2;
            break;
        case 21:
            contactMax=2;
            break;
        case 22:
            contactMax=6;
            break;
        case 23:
            contactMax=4;
            break;
        case 24:
            contactMax=4;
            break;
        case 25:
            contactMax=6;
            break;
        case 26:
            contactMax=5;
            break;
        case 27:
            contactMax=3;
            break;
        case 28:
            contactMax=3;
            break;
        case 29:
            contactMax=4;
            break;
        case 30:
            contactMax=4;
            break;
        case 31:
            contactMax=4;
            break;
        case 32:
            contactMax=5;
            break;
        case 33:
            contactMax=4;
            break;
        case 34:
            contactMax=6;
            break;
        case 35:
            contactMax=2;
            break;
        case 36:
            contactMax=5;
            break;
        case 37:
            contactMax=2;
            break;
    }

    tempCurrent=*pHeadCopy; //Target points to the individual at the beginning of the list
    //printf("total at-risk agents in prison %d: %d\n", targetLocation, totalAtRisk);
    
    while(tempCurrent!=NULL){
        //add if injecting
        //FIX location matching to infect other populations
        //
        if(tempCurrent->court==1){//susceptible if also in court
            if(tempCurrent->courtNumber==target->courtNumber&&tempCurrent->courtCell==target->courtCell){
                //if(tempCurrent->hospitalised==0&&tempCurrent->isolated==0&&tempCurrent->court==0){ //not in any hospital
                    //if(tempCurrent->COVID==0){ //not yet been infected and not currently infected
                    if(tempCurrent->indivType==0){
                        if(infectedCountI<contactMax){ //infectedCountI<0contactCapacityI
                            
                            //herd immunity
                            if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                                ppeEffect=0.67;
                                if(target->severity==4||target->severity==5){
                                    probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }else{
                                    probInf[1]=gsl_ran_flat(*r, 0.02, 0.05);//0.05;//probability of getting infected with COVID-19
                                    probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                                }
                                probInf[0]=1-probInf[1];
                            }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                                probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                                probInf[0]=1-probInf[1];
                            }
                            
                            printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                            
                            binInf=draw_multinom(&rCopy, 2, probInf);
                            //printf("Infection commencing\n");
                            printf("Infect Y or N: %d\n", binInf);
                            
                            if(binInf==1){ //Target is to be infected
                                printf("Infection commencing\n");
                                
                                tempCurrent->COVID=1; // infected!!
                                tempCurrent->severity=0;
                                currLoc=tempCurrent->location;
                                tempCurrent->placeInfected=currLoc;
                                targetAge=tempCurrent->age;
                                
                                if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                    (*newInfectedInmatesMin)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=1;
                                        
                                        if(targetAge<=19){
                                            //(*zoneMinCIArray)[tempCurrent->zone][0]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][0]++;
                                            (*totalAge0)++;
                                        }else if(targetAge>=20&&targetAge<=44){
                                            //(*zoneMinCIArray)[tempCurrent->zone][1]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][1]++;
                                            (*totalAge1)++;
                                        }else if(targetAge>=45&&targetAge<=54){
                                            //(*zoneMinCIArray)[tempCurrent->zone][2]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][2]++;
                                            (*totalAge2)++;
                                        }else if(targetAge>=55&&targetAge<=64){
                                            //(*zoneMinCIArray)[tempCurrent->zone][3]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][3]++;
                                            (*totalAge3)++;
                                        }else if(targetAge>=65&&targetAge<=74){
                                            //(*zoneMinCIArray)[tempCurrent->zone][4]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][4]++;
                                            (*totalAge4)++;
                                        }else if(targetAge>=75&&targetAge<=84){
                                            //(*zoneMinCIArray)[tempCurrent->zone][5]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][5]++;
                                            (*totalAge5)++;
                                        }else if(targetAge>=85){
                                            //(*zoneMinCIArray)[tempCurrent->zone][6]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][6]++;
                                            (*totalAge6)++;
                                        }
                                        
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=2;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=3;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=4;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=5;
                                    }
                                }else if(currLoc==2){
                                    (*newInfectedInmatesMed)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=6;
                                        
                                        if(targetAge<=19){
                                            //(*zoneMedCIArray)[tempCurrent->zone][0]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][0]++;
                                            (*totalAge0)++;
                                        }else if(targetAge>=20&&targetAge<=44){
                                            //(*zoneMedCIArray)[tempCurrent->zone][1]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][1]++;
                                            (*totalAge1)++;
                                        }else if(targetAge>=45&&targetAge<=54){
                                            //(*zoneMedCIArray)[tempCurrent->zone][2]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][2]++;
                                            (*totalAge2)++;
                                        }else if(targetAge>=55&&targetAge<=64){
                                            //(*zoneMedCIArray)[tempCurrent->zone][3]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][3]++;
                                            (*totalAge3)++;
                                        }else if(targetAge>=65&&targetAge<=74){
                                            //(*zoneMedCIArray)[tempCurrent->zone][4]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][4]++;
                                            (*totalAge4)++;
                                        }else if(targetAge>=75&&targetAge<=84){
                                            //(*zoneMedCIArray)[tempCurrent->zone][5]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][5]++;
                                            (*totalAge5)++;
                                        }else if(targetAge>=85){
                                            //(*zoneMedCIArray)[tempCurrent->zone][6]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][6]++;
                                            (*totalAge6)++;
                                        }
                                        
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=7;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=8;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=9;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=10;
                                    }
                                }else if(currLoc==3){
                                    (*newInfectedInmatesMax)++;
                                    
                                    if(tempCurrent->indivType==0){
                                        currLoc=11;
                                        
                                        if(targetAge<=19){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][0]++;
                                            (*totalAge0)++;
                                        }else if(targetAge>=20&&targetAge<=44){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][1]++;
                                            (*totalAge1)++;
                                        }else if(targetAge>=45&&targetAge<=54){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][2]++;
                                            (*totalAge2)++;
                                        }else if(targetAge>=55&&targetAge<=64){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][3]++;
                                            (*totalAge3)++;
                                        }else if(targetAge>=65&&targetAge<=74){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][4]++;
                                            (*totalAge4)++;
                                        }else if(targetAge>=75&&targetAge<=84){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][5]++;
                                            (*totalAge5)++;
                                        }else if(targetAge>=85){
                                            //(*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                            (*courtCIArray)[tempCurrent->courtNumber][6]++;
                                            (*totalAge6)++;
                                        }
                                        
                                    }else if(tempCurrent->indivType==1){
                                        currLoc=12;
                                    }else if(tempCurrent->indivType==2){
                                        currLoc=13;
                                    }else if(tempCurrent->indivType==3){
                                        currLoc=14;
                                    }else if(tempCurrent->indivType==4){
                                        currLoc=15;
                                    }
                                }
                                
                                tempCurrent->timeOfInfection=currDay; //Record date of infection
                                ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                                currGroup=tempCurrent->group;
                                
                                (*pLocArray2)[currLoc][currGroup]--;
                                
                                //system("pause");
                                switch(currGroup){
                                    case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                        break;
                                    case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                        break;
                                    case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                        break;
                                    case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=0;
                                        break;
                                    case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                        break;
                                    case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=1;
                                        break;
                                    case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                        break;
                                    case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                        break;
                                    case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                        break;
                                    case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=6;
                                        break;
                                    case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                        printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                        break;
                                    case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                        newGroup=7;
                                        break;
                                        
                                }
                                
                                //if(currGroup!=newGroup){
                                tempCurrent->group=newGroup;
                                //}
                                
                                printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                                success++;
                                
                                (*pLocArray2)[currLoc][newGroup]++;
                                printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                                
                                (*newInfectedInmates)++;
                            }
                            
                            infectedCountI++;
                        }
                    }
                    
                    //}//covid
                //}//hospitalised
            }//same zone
        }//same location
        tempCurrent=tempCurrent->nextIndiv;
    }
    printf("Number of new infections for this instance: %d", success);
    return success;
}

int infectInTransit(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*truckCIArray)[20][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int currDay, int *pMoveFlag){
    sIndiv *target, *tempCurrent;
    target=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int rStatus, binInf;
    //int zone=target->zone;
    int targetLocation=target->location, aCtr, iTarget, currLoc, currGroup, newGroup=99, eFlag=0;
    //int totalAtRisk=(*pLocArray2)[target->location][3]+(*pLocArray2)[target->location][5]+(*pLocArray2)[target->location][9]+(*pLocArray2)[target->location][11]; //add all at risk except those with immunity
    //sIndiv *refInfectList[totalAtRisk];
    //double probInfectList[totalAtRisk];
    double eventDraw[2]; //proability of recipient to get infected
    double probInf[2]; //probability of source infecting recipient using per event probability
    int success, eventDecision, susCtr;
    int targetArrayLoc;
    int contactIMax=0, contactIMed=0, contactIMin=0, contactPSMax=0, contactPSMed=0, contactPSMin=0, contactHSMax=0, contactHSMed=0, contactHSMin=0, contactEVMax=0, contactEVMed=0, contactEVMin=0, contactFVMax=0, contactFVMed=0, contactFVMin=0, contactMax=0;
    int contactCapacityI=0, contactCapacityPS=0, contactCapacityHS=0, contactCapacityEV=0, contactCapacityFV=0;
    int targetAge;
    int infectedCountI=0, infectedCountPS=0, infectedCountHS=0, infectedCountEV=0, infectedCountFV=0;
    printf("going through list of individuals\n");
    success=0;
    //probInf[0]=0.0; probInf[1]=0.0;
    eventDraw[0]=0; eventDraw[1]=0;
    double ppeEffect;
    //aCtr=0;
    //Go through list of individuals
    
    ppeEffect=0.67;
    
    
    /*
     if(target->indivType==0){ //Inmate
     contactC1=generateAge(36, 81);//gsl_ran_flat(*r, 36, 81);//Interaction within Units/Pods only
     contactC2=generateAge(35, 79);//gsl_ran_flat(*r, 107, 130);//Interaction within Units/Pods
     contactIMin=generateAge(26, 59);//gsl_ran_flat(*r, 87, 107);//Interaction within Units/Pods
     contactPSMax=0;
     contactPSMed=0;
     contactPSMin=0;
     contactHSMax=0;
     contactHSMed=0;
     contactHSMin=0;
     contactEVMax=0;
     contactEVMed=0;
     contactEVMin=0;
     contactFVMax=0;
     contactFVMed=0;
     contactFVMin=0;
     if(targetLocation==1){
     targetArrayLoc=1;
     }else if(targetLocation==2){
     targetArrayLoc=6;
     }else if(targetLocation==3){
     targetArrayLoc=11;
     }
     }
     
     
     if(targetLocation==1){
     contactCapacityI=contactIMin;
     contactCapacityPS=contactPSMin;
     contactCapacityHS=contactHSMin;
     contactCapacityEV=contactEVMin;
     contactCapacityFV=contactFVMin;
     }else if(targetLocation==2){
     contactCapacityI=contactIMed;
     contactCapacityPS=contactPSMed;
     contactCapacityHS=contactHSMed;
     contactCapacityEV=contactEVMed;
     contactCapacityFV=contactFVMed;
     }else if(targetLocation==3){
     contactCapacityI=contactIMax;
     contactCapacityPS=contactPSMax;
     contactCapacityHS=contactHSMax;
     contactCapacityEV=contactEVMax;
     contactCapacityFV=contactFVMax;
     }
     */
    contactMax=16;
    
    tempCurrent=*pHeadCopy; //Target points to the individual at the beginning of the list
    //printf("total at-risk agents in prison %d: %d\n", targetLocation, totalAtRisk);
    
    while(tempCurrent!=NULL){
        //add if injecting
        //FIX location matching to infect other populations
        //
        if(tempCurrent->moving==1){//susceptible if also in court
            if(tempCurrent->truckNumber==target->truckNumber){
                //if(tempCurrent->hospitalised==0&&tempCurrent->isolated==0&&tempCurrent->court==0){ //not in any hospital
                //if(tempCurrent->COVID==0){ //not yet been infected and not currently infected
                if(tempCurrent->indivType==0){
                    if(infectedCountI<contactMax){ //infectedCountI<0contactCapacityI
                        
                        //herd immunity
                        if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==1){ //immunity
                            probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                            probInf[0]=1-probInf[1];
                        }else if(tempCurrent->COVID==0&&tempCurrent->COVIDAb==0){ //immunity
                            //ppeEffect=0.67;
                            if(target->severity==4||target->severity==5){
                                probInf[1]=0.01*2;//0.05;//probability of getting infected with COVID-19
                                probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                            }else{
                                probInf[1]=gsl_ran_flat(*r, 0.02, 0.05);//0.05;//probability of getting infected with COVID-19
                                probInf[1]=probInf[1]-(probInf[1]*ppeEffect); //assuming face mask at all times
                            }
                            probInf[0]=1-probInf[1];
                        }else if(tempCurrent->COVID==1&&tempCurrent->COVIDAb==0){
                            probInf[1]=0.0*2;//0.05;//probability of getting infected with COVID-19
                            probInf[0]=1-probInf[1];
                        }
                        
                        printf("Infect %d [ID: %d](%p), with probability %f?\n", aCtr, tempCurrent->ID, tempCurrent, probInf[1]);
                        
                        binInf=draw_multinom(&rCopy, 2, probInf);
                        //printf("Infection commencing\n");
                        printf("Infect Y or N: %d\n", binInf);
                        
                        if(binInf==1){ //Target is to be infected
                            printf("Infection commencing\n");
                            
                            tempCurrent->COVID=1; // infected!!
                            tempCurrent->severity=0;
                            currLoc=tempCurrent->location;
                            tempCurrent->placeInfected=currLoc;
                            targetAge=tempCurrent->age;
                            
                            if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==11
                                (*newInfectedInmatesMin)++;
                                
                                if(tempCurrent->indivType==0){
                                    currLoc=1;
                                    
                                    if(targetAge<=19){
                                        //(*zoneMinCIArray)[tempCurrent->zone][0]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][0]++;
                                        (*totalAge0)++;
                                    }else if(targetAge>=20&&targetAge<=44){
                                        //(*zoneMinCIArray)[tempCurrent->zone][1]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][1]++;
                                        (*totalAge1)++;
                                    }else if(targetAge>=45&&targetAge<=54){
                                        //(*zoneMinCIArray)[tempCurrent->zone][2]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][2]++;
                                        (*totalAge2)++;
                                    }else if(targetAge>=55&&targetAge<=64){
                                        //(*zoneMinCIArray)[tempCurrent->zone][3]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][3]++;
                                        (*totalAge3)++;
                                    }else if(targetAge>=65&&targetAge<=74){
                                        //(*zoneMinCIArray)[tempCurrent->zone][4]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][4]++;
                                        (*totalAge4)++;
                                    }else if(targetAge>=75&&targetAge<=84){
                                        //(*zoneMinCIArray)[tempCurrent->zone][5]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][5]++;
                                        (*totalAge5)++;
                                    }else if(targetAge>=85){
                                        //(*zoneMinCIArray)[tempCurrent->zone][6]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][6]++;
                                        (*totalAge6)++;
                                    }
                                    
                                }else if(tempCurrent->indivType==1){
                                    currLoc=2;
                                }else if(tempCurrent->indivType==2){
                                    currLoc=3;
                                }else if(tempCurrent->indivType==3){
                                    currLoc=4;
                                }else if(tempCurrent->indivType==4){
                                    currLoc=5;
                                }
                            }else if(currLoc==2){
                                (*newInfectedInmatesMed)++;
                                
                                if(tempCurrent->indivType==0){
                                    currLoc=6;
                                    
                                    if(targetAge<=19){
                                        //(*zoneMedCIArray)[tempCurrent->zone][0]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][0]++;
                                        (*totalAge0)++;
                                    }else if(targetAge>=20&&targetAge<=44){
                                        //(*zoneMedCIArray)[tempCurrent->zone][1]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][1]++;
                                        (*totalAge1)++;
                                    }else if(targetAge>=45&&targetAge<=54){
                                        //(*zoneMedCIArray)[tempCurrent->zone][2]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][2]++;
                                        (*totalAge2)++;
                                    }else if(targetAge>=55&&targetAge<=64){
                                        //(*zoneMedCIArray)[tempCurrent->zone][3]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][3]++;
                                        (*totalAge3)++;
                                    }else if(targetAge>=65&&targetAge<=74){
                                        //(*zoneMedCIArray)[tempCurrent->zone][4]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][4]++;
                                        (*totalAge4)++;
                                    }else if(targetAge>=75&&targetAge<=84){
                                        //(*zoneMedCIArray)[tempCurrent->zone][5]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][5]++;
                                        (*totalAge5)++;
                                    }else if(targetAge>=85){
                                        //(*zoneMedCIArray)[tempCurrent->zone][6]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][6]++;
                                        (*totalAge6)++;
                                    }
                                    
                                }else if(tempCurrent->indivType==1){
                                    currLoc=7;
                                }else if(tempCurrent->indivType==2){
                                    currLoc=8;
                                }else if(tempCurrent->indivType==3){
                                    currLoc=9;
                                }else if(tempCurrent->indivType==4){
                                    currLoc=10;
                                }
                            }else if(currLoc==3){
                                (*newInfectedInmatesMax)++;
                                
                                if(tempCurrent->indivType==0){
                                    currLoc=11;
                                    
                                    if(targetAge<=19){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][0]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][0]++;
                                        (*totalAge0)++;
                                    }else if(targetAge>=20&&targetAge<=44){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][1]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][1]++;
                                        (*totalAge1)++;
                                    }else if(targetAge>=45&&targetAge<=54){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][2]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][2]++;
                                        (*totalAge2)++;
                                    }else if(targetAge>=55&&targetAge<=64){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][3]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][3]++;
                                        (*totalAge3)++;
                                    }else if(targetAge>=65&&targetAge<=74){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][4]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][4]++;
                                        (*totalAge4)++;
                                    }else if(targetAge>=75&&targetAge<=84){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][5]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][5]++;
                                        (*totalAge5)++;
                                    }else if(targetAge>=85){
                                        //(*zoneMaxCIArray)[tempCurrent->zone][6]++;
                                        (*truckCIArray)[tempCurrent->truckNumber][6]++;
                                        (*totalAge6)++;
                                    }
                                    
                                }else if(tempCurrent->indivType==1){
                                    currLoc=12;
                                }else if(tempCurrent->indivType==2){
                                    currLoc=13;
                                }else if(tempCurrent->indivType==3){
                                    currLoc=14;
                                }else if(tempCurrent->indivType==4){
                                    currLoc=15;
                                }
                            }
                            
                            tempCurrent->timeOfInfection=currDay; //Record date of infection
                            ////////UPDATE SUBGROUP WHEN INFECTED AND UPDATE BUCKETS
                            currGroup=tempCurrent->group;
                            
                            (*pLocArray2)[currLoc][currGroup]--;
                            
                            //system("pause");
                            switch(currGroup){
                                case 0: //[currLoc][0] IDU+; HCV+; ATSI;  already infected
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 0");
                                    break;
                                case 1: //[currLoc][1] IDU+; HCV+; NON-ATSI already infected
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 1");
                                    break;
                                case 2: //[currLoc][2] IDU+; HCV-; ATSI; PREV. EXPOSED not yet infected
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 2");
                                    break;
                                case 3: //[currLoc][3] IDU+; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                    newGroup=0;
                                    break;
                                case 4: //[currLoc][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 4");
                                    break;
                                case 5: //[currLoc][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                    newGroup=1;
                                    break;
                                case 6: //[currLoc][6] IDU-; HCV+; ATSI
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 6");
                                    break;
                                case 7: //[currLoc][7] IDU-; HCV+; NON-ATSI
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 7");
                                    break;
                                case 8: //[currLoc][8] IDU-; HCV-; ATSI; PREV. EXPOSED not yet infected
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 8");
                                    break;
                                case 9: //[currLoc][9] IDU-; HCV-; ATSI; SUSCEPTIBLE not yet infected
                                    newGroup=6;
                                    break;
                                case 10://[currLoc][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED not yet infected
                                    printf("ERROR: NOT SUPPOSED TO HAPPEN 10");
                                    break;
                                case 11://[currLoc][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE not yet infected
                                    newGroup=7;
                                    break;
                                    
                            }
                            
                            //if(currGroup!=newGroup){
                            tempCurrent->group=newGroup;
                            //}
                            
                            printf("Inf Pop ex-group: %d . Inf Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
                            success++;
                            
                            (*pLocArray2)[currLoc][newGroup]++;
                            printf("Infected ID %d on day %d!!!\n", tempCurrent->ID, tempCurrent->timeOfInfection);
                            
                            (*newInfectedInmates)++;
                        }
                        
                        infectedCountI++;
                    }
                }
                
                //}//covid
                //}//hospitalised
            }//same zone
        }//same location
        tempCurrent=tempCurrent->nextIndiv;
    }
    printf("Number of new infections for this instance: %d", success);
    return success;
}

void NtrClear(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]){
    sIndiv *target;
    target=*pTargetCopy;
    int currLoc, newGroup=0, currGroup, currRealLoc;
    
    currLoc=target->location;
    
    if(currLoc==1){ //objLoc==1||objLoc==6||objLoc==1
        if(target->indivType==0){
            currRealLoc=1;
        }else if(target->indivType==1){
            currRealLoc=2;
        }else if(target->indivType==2){
            currRealLoc=3;
        }else if(target->indivType==3){
            currRealLoc=4;
        }else if(target->indivType==4){
            currRealLoc=5;
        }
    }else if(currLoc==2){
        if(target->indivType==0){
            currRealLoc=6;
        }else if(target->indivType==1){
            currRealLoc=7;
        }else if(target->indivType==2){
            currRealLoc=8;
        }else if(target->indivType==3){
            currRealLoc=9;
        }else if(target->indivType==4){
            currRealLoc=10;
        }
    }else if(currLoc==3){
        if(target->indivType==0){
            currRealLoc=11;
        }else if(target->indivType==1){
            currRealLoc=12;
        }else if(target->indivType==2){
            currRealLoc=13;
        }else if(target->indivType==3){
            currRealLoc=14;
        }else if(target->indivType==4){
            currRealLoc=15;
        }
    }
    
    currGroup=target->group;
    
    (*pLocArray2)[currRealLoc][currGroup]--;
    
    switch(currGroup){
        case 0:
            newGroup=2;
            break;
        case 1:
            newGroup=4;
            break;
        case 2:
            printf("CLEAR COVID ERROR");
            break;
        case 3:
            printf("CLEAR COVID ERROR");
            break;
        case 4:
            printf("CLEAR COVID ERROR");
            break;
        case 5:
            printf("CLEAR COVID ERROR");
            break;
        case 6:
            newGroup=8;
            break;
        case 7:
            newGroup=10;
            break;
        case 8:
            printf("CLEAR COVID ERROR");
            break;
        case 9:
            printf("CLEAR COVID ERROR");
            break;
        case 10:
            printf("CLEAR COVID ERROR");
            break;
        case 11:
            printf("CLEAR COVID ERROR");
    }
    
    target->COVID=0;
    target->COVIDAb=1;
    target->severity=5; //clearer
    target->timeOfInfection=0;
    
    target->group=newGroup;
    (*pLocArray2)[currRealLoc][newGroup]++;
    printf("Clr ex-group: %d . Clr new group: %d\n",currGroup,newGroup);
    printf("Clr Pop ex-group: %d . Clr Pop new group: %d\n",(*pLocArray2)[currRealLoc][currGroup],(*pLocArray2)[currRealLoc][newGroup]);
    
    /*
     if(){//condition for number of days or severity?
     currLoc=target->location;
     currGroup=target->group;
     switch(currGroup){
     case 0:
     newGroup=2;
     break;
     case 1:
     newGroup=4;
     break;
     case 2:
     printf("CLEAR COVID ERROR");
     break;
     case 3:
     printf("CLEAR COVID ERROR");
     break;
     case 4:
     printf("CLEAR COVID ERROR");
     break;
     case 5:
     printf("CLEAR COVID ERROR");
     break;
     case 6:
     newGroup=8;
     break;
     case 7:
     newGroup=10;
     break;
     case 8:
     printf("CLEAR COVID ERROR");
     break;
     case 9:
     printf("CLEAR COVID ERROR");
     break;
     case 10:
     printf("CLEAR COVID ERROR");
     break;
     case 11:
     printf("CLEAR COVID ERROR");
     }
     }
     */
    //UPDATE prevention strategy
}

void progress(gsl_rng **r, sIndiv **pTargetCopy, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int *pre19, int *pre20to44, int *pre45to54, int *pre55to64, int *pre65to74, int *pre75to84, int *pre85, int *asy19, int *asy20to44, int *asy45to54, int *asy55to64, int *asy65to74, int *asy75to84, int *asy85, int *mil19, int *mil20to44, int *mil45to54, int *mil55to64, int *mil65to74, int *mil75to84, int *mil85, int *mod19, int *mod20to44, int *mod45to54, int *mod55to64, int *mod65to74, int *mod75to84, int *mod85, int *sev19, int *sev20to44, int *sev45to54, int *sev55to64, int *sev65to74, int *sev75to84, int *sev85, int *cle19, int *cle20to44, int *cle45to54, int *cle55to64, int *cle65to74, int *cle75to84, int *cle85, int currDay){
    
    sIndiv *current;
    current=*pTargetCopy;
    gsl_rng *rCopy;
    rCopy=*r;
    int infTime, age, deltaTime;
    int binEvent;
    double probEvent[2];
    probEvent[1]=0.0;
    probEvent[0]=0.0;
    
    infTime=current->timeOfInfection;
    age=current->age;
    
    deltaTime=currDay-infTime;
    
    //if(target->severity<4){
    //    target->severity++;
    //}
    
    if(current->severity==0){ //exposed -> pre-clinical
        //infTime=infTime+2;
        //if(deltaTime>=2&&deltaTime<=7){ //delta days and severity condition?
        //probEvent[1]=1.0;//gsl_ran_flat(*r, 0.23, 0.43);
        //}else if(deltaTime>7){
        //    probEvent[1]=1.0;
        //}
        
        //probEvent[0]=1-probEvent[1];
        binEvent=1;//draw_multinom(&rCopy, 2, probEvent);
        
        if(binEvent==1){//yes
            current->severity++;
            if(current->indivType==0){
                if(current->age<=19){ //means infected in a prison
                    (*pre19)++;
                }else if(current->age>=20&&current->age<=44){
                    (*pre20to44)++;
                }else if(current->age>=45&&current->age<=54){
                    (*pre45to54)++;
                }else if(current->age>=55&&current->age<=64){
                    (*pre55to64)++;
                }else if(current->age>=65&&current->age<=74){
                    (*pre65to74)++;
                }else if(current->age>=75&&current->age<=84){
                    (*pre75to84)++;
                }else if(current->age>=85){
                    (*pre85)++;
                }
            }
        }
    }else if(current->severity==1){ //pre-clinical -> asymptomatic or mild
        //infTime=infTime+2;
        //currDay>=infTime
        //if(deltaTime>=4&&deltaTime<=7){ //delta days and severity condition?
        if(age>=0&&age<=19){
            probEvent[1]=1-exp(-0.75);//gsl_ran_gaussian(*r, 0.058906);//0.75;//gsl_ran_flat(*r, 0.625, 0.85);
        }else if(age>=20&&age<=44){
            probEvent[1]=1-exp(-0.67);//gsl_ran_gaussian(*r, 0.054375);//0.67;//gsl_ran_flat(*r, 0.56, 0.77);
        }else if(age>=45&&age<=54){
            probEvent[1]=1-exp(-0.55);//gsl_ran_gaussian(*r, 0.058594);//0.55;//gsl_ran_flat(*r, 0.44, 0.657);
        }else if(age>=55&&age<=64){
            probEvent[1]=1-exp(-0.44);//gsl_ran_gaussian(*r, 0.073750);//0.44;//gsl_ran_flat(*r, 0.32, 0.57);
        }else if(age>=65&&age<=74){
            probEvent[1]=1-exp(-0.34);//gsl_ran_gaussian(*r, 0.065781);//0.34;//gsl_ran_flat(*r, 0.21, 0.47);
        }else if(age>=75&&age<=84){
            probEvent[1]=1-exp(-0.31);//gsl_ran_gaussian(*r, 0.063750);//0.31;//gsl_ran_flat(*r, 0.18, 0.43);
        }else{
            probEvent[1]=1-exp(-0.31);//gsl_ran_gaussian(*r, 0.063750);//0.31;//gsl_ran_flat(*r, 0.18, 0.43);
        }
        //}else if(deltaTime>7){
        //probEvent[1]=1.0;
        //}
        
        probEvent[0]=1-probEvent[1];
        binEvent=draw_multinom(&rCopy, 2, probEvent);
        
        if(binEvent==1){//yes == asymptomatic
            current->severity=2;
            if(current->indivType==0){
                if(current->age<=19){ //means infected in a prison
                    (*asy19)++;
                }else if(current->age>=20&&current->age<=44){
                    (*asy20to44)++;
                }else if(current->age>=45&&current->age<=54){
                    (*asy45to54)++;
                }else if(current->age>=55&&current->age<=64){
                    (*asy55to64)++;
                }else if(current->age>=65&&current->age<=74){
                    (*asy65to74)++;
                }else if(current->age>=75&&current->age<=84){
                    (*asy75to84)++;
                }else if(current->age>=85){
                    (*asy85)++;
                }
            }
        }else{ //mild
            current->severity=3;
            if(current->indivType==0){
                if(current->age<=19){ //means infected in a prison
                    (*mil19)++;
                }else if(current->age>=20&&current->age<=44){
                    (*mil20to44)++;
                }else if(current->age>=45&&current->age<=54){
                    (*mil45to54)++;
                }else if(current->age>=55&&current->age<=64){
                    (*mil55to64)++;
                }else if(current->age>=65&&current->age<=74){
                    (*mil65to74)++;
                }else if(current->age>=75&&current->age<=84){
                    (*mil75to84)++;
                }else if(current->age>=85){
                    (*mil85)++;
                }
            }
        }
    }else if(current->severity==3){ //mild to moderate
        //infTime=infTime+5;
        //if(deltaTime>=8){ //delta days and severity condition?
        //probEvent[1]=gsl_ran_flat(*r, 0.15, 0.21);
        //}
        
        //probEvent[0]=1-probEvent[1];
        binEvent=1;//draw_multinom(&rCopy, 2, probEvent);
        
        if(binEvent==1){//yes
            current->severity=4;
            if(current->indivType==0){
                if(current->age<=19){ //means infected in a prison
                    (*mod19)++;
                }else if(current->age>=20&&current->age<=44){
                    (*mod20to44)++;
                }else if(current->age>=45&&current->age<=54){
                    (*mod45to54)++;
                }else if(current->age>=55&&current->age<=64){
                    (*mod55to64)++;
                }else if(current->age>=65&&current->age<=74){
                    (*mod65to74)++;
                }else if(current->age>=75&&current->age<=84){
                    (*mod75to84)++;
                }else if(current->age>=85){
                    (*mod85)++;
                }
            }
        }
        
    }else if(current->severity==4){ //moderate to severe
        //infTime=infTime+5;
        //if(deltaTime>=13){ //delta days and severity condition?
        /*if(age>=0&&age<=19){
         probEvent[1]=gsl_ran_flat(*r, 0.00012, 0.00042);
         }else if(age>=20&&age<=44){
         probEvent[1]=gsl_ran_flat(*r, 0.01731, 0.05937);
         }else if(age>=45&&age<=54){
         probEvent[1]=gsl_ran_flat(*r, 0.03695, 0.1269);
         }else if(age>=55&&age<=64){
         probEvent[1]=gsl_ran_flat(*r, 0.05935, 0.2035);
         }else if(age>=65&&age<=74){
         probEvent[1]=gsl_ran_flat(*r, 0.0844, 0.289);
         }else if(age>=75&&age<=84){
         probEvent[1]=gsl_ran_flat(*r, 0.10435, 0.357);
         }else{
         probEvent[1]=gsl_ran_flat(*r, 0.11, 0.376);
         }*/
        //}
        
        //probEvent[0]=1-probEvent[1];
        binEvent=1;//draw_multinom(&rCopy, 2, probEvent);
        
        if(binEvent==1){//yes
            current->severity=5;
            if(current->indivType==0){
                if(current->age<=19){ //means infected in a prison
                    (*sev19)++;
                }else if(current->age>=20&&current->age<=44){
                    (*sev20to44)++;
                }else if(current->age>=45&&current->age<=54){
                    (*sev45to54)++;
                }else if(current->age>=55&&current->age<=64){
                    (*sev55to64)++;
                }else if(current->age>=65&&current->age<=74){
                    (*sev65to74)++;
                }else if(current->age>=75&&current->age<=84){
                    (*sev75to84)++;
                }else if(current->age>=85){
                    (*sev85)++;
                }
            }
        }
    }
}

void progressHCV(sIndiv **pTargetCopy){
    sIndiv *target;
    target=*pTargetCopy;
    
    if(target->metavir<3){
        target->metavir++;
    }//else
        //printf("HCV liver fibrosis at max!\n");
}

int deathHCV(sIndiv **pTargetCopy){ //NOT IN USE
    sIndiv *target;
    int sRand, clrThresh;
    target=*pTargetCopy;
    
    clrThresh=50; //uniform distribution
    
    do
        sRand=generateRand();
    while(sRand<1||sRand>100);
    
    if(sRand<=clrThresh){
        return 1; //remove Indiv
        //remove in the main program
    }else{
        return 0;
    }
}

//COVID strategies

void startOST(sIndiv **pTargetCopy, int currDay){
    sIndiv *target;
    target=*pTargetCopy;
    
    if(currDay<=4380){
        target->OST=1;
        target->timeStartOST=currDay;
    }else if(currDay>4380){
        target->OST=2;
        target->timeStartOST=currDay;
    }
}

void stopOST(sIndiv **pTargetCopy){
    sIndiv *target;
    target=*pTargetCopy;
    
    target->OST=0;
    target->timeStartOST=0;
}

void startDAA(sIndiv **pTargetCopy, int currDay){
    sIndiv *target;
    target=*pTargetCopy;
    
    target->DAA=1;
    target->timeStartDAA=currDay;
}

void stopDAA(sIndiv **pTargetCopy){
    sIndiv *target;
    target=*pTargetCopy;
    
    target->DAA=0;
    target->timeStartDAA=0;
}

void isolate(sIndiv **pTargetCopy){ //1-field hospital 2-community hospital
    sIndiv *current;
    current=*pTargetCopy;
    
    current->isolated=1; //0 not isolated, 1 for 1 out isolation, 2 for 2 out isolation
}

void endIsolate(sIndiv **pTargetCopy){ //1-field hospital 2-community hospital
    sIndiv *current;
    current=*pTargetCopy;
    
    current->isolated=0;
}

void clearDAA(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]){
    sIndiv *target;
    target=*pTargetCopy;
    int currLoc, newGroup=0, currGroup;
    
    if(target->metavir<5){
    	target->metavir=5; //clear if less than equal to selected threshold
    	target->placeInfected=0; //NA
    	target->timeOfInfection=0;
    	if(target->DAA==1){
        	target->DAA=0; //set DAA to 0
        	target->timeStartDAA=0; //set to 0
    	}
    
    	currLoc=target->location;
    	currGroup=target->group;
    
    	switch(currGroup){
        	case 0:
            	newGroup=2;
            	break;
        	case 1:
            	newGroup=4;
            	break;
        	case 2:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 3:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 4:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 5:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 6:
            	newGroup=8;
            	break;
        	case 7:
            	newGroup=10;
            break;
        	case 8:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 9:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 10:
            	printf("CLEAR HCV ERROR");
            	break;
        	case 11:
        	    printf("CLEAR HCV ERROR");
    	}
    	target->group=newGroup;
    	//update pLocArray; move to HCV-, prev. exposed
    	(*pLocArray2)[currLoc][currGroup]--;
    	(*pLocArray2)[currLoc][newGroup]++;
    	printf("Clr DAA ex-group: %d . Clr new group: %d\n",currGroup,newGroup);
    	printf("Clr DAA Pop ex-group: %d . Clr Pop new group: %d\n",(*pLocArray2)[currLoc][currGroup],(*pLocArray2)[currLoc][newGroup]);
	}
}
