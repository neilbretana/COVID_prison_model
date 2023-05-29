//
//  main.c
//  HCV Model
//
//  Created by Neil Bretana on 15/08/2014.
//  Copyright (c) 2014 Neil Bretana. All rights reserved.
// Quarantine upon entry for 14 days only with PCR tests at day 1 and day 12
// if positive, isolate until clear
// PPE applied standard masks for all 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "individual.h"

#define TCOLCTR 2   //Array of locations containing:
                    //[row][0] IDU+;
                    //[row][1] IDU-;
#define COLCTR 12  //Array of locations containing:
//[row][0] IDU+; HCV+; ATSI                                             // Co-Morbidity+; COVID+; ATSI
//[row][1] IDU+; HCV+; NON-ATSI                                         // Co-Morbidity+; COVID+; Non-ATSI
//[row][2] IDU+; HCV-; ATSI; PREV. EXPOSED                              // Co-Morbidity+; COVID-; ATSI; Antibody+
//[row][3] IDU+; HCV-; ATSI; SUSCEPTIBLE                                // Co-Morbidity+; COVID-; ATSI; Antibody-
//[row][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED                          // Co-Morbidity+; COVID-; Non-ATSI; Antibody+
//[row][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE                            // Co-Morbidity+; COVID-; Non-ATSI; Antibody-
//[row][6] IDU-; HCV+; ATSI                                             // Co-Morbidity-; COVID+; ATSI;
//[row][7] IDU-; HCV+; NON-ATSI                                         // Co-Morbidity-; COVID+; Non-ATSI;
//[row][8] IDU-; HCV-; ATSI; PREV. EXPOSED                              // Co-Morbidity-; COVID-; ATSI; Antibody+
//[row][9] IDU-; HCV-; ATSI; SUSCEPTIBLE                                // Co-Morbidity-; COVID-; ATSI; Antibody-
//[row][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED                         // Co-Morbidity-; COVID-; Non-ATSI; Antibody+
//[row][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE                           // Co-Morbidity-; COVID-; Non-ATSI; Antibody-
//#define ROWPRIS 4 //# prisons +1 community
#define ROWPRIS 16 //community 0, prison 1: 1 inmates, 2 prison staff, 3 healthcare staff, 4 essential visitors 5 family visitors. prison 2: 6, 7, 8, 9, 10. prison 3: 11, 12, 13, 14, 15.
#define NEVENTS 27 //number of events
#define RGROUPS 13 //Risk groups: 0: non-injecting; 1: injecting less than daily; opioid; not sharing 2: injecting less than daily; opioid; sharing less daily 3: injecting less than daily; opioid; sharing daily more 4: injecting less than daily; non-opioid; no sharing 5: injecting less than daily; non-opioid; sharing less daily 6: injecting less than daily; non-opioid; sharing daily or more 7: ...

const gsl_rng_type *T;

int generateRand();
int draw_multinom(gsl_rng **r, int nEvents, double probsInput[]);
int distributePop(int input);
void traverse (sIndiv **pHeadCopy, sIndiv **pTailCopy);
void staffUnisolate (sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay);
int countHCVpris (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countHCVCom (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countOpd (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countOpdNotOST (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countEverIDU (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count0(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count1 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count2(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count3(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count4 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count5(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count6(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count7 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count8(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count9(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count10 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count11(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count12(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countR(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countHCVantibody(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countHCVRNA(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countMovingPop(sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*nMovingPop)[20]);
int countAveStay(sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay);
void countZonePop(sIndiv **pHeadCopy, int (*pZoneArray)[30][3], int (*pZoneCovidArray)[30][3]);
void countSevereInmates(sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pAgeSeverityArray)[7][8]);
void checkIsolate (sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *povercap);
void checkQuarantine (sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *povercap);
void PCRinmate(gsl_rng **r, sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18], int (*pLockdownMinArea)[27][2], int (*pLockdownMedArea)[11][2], int (*pLockdownMaxArea)[18][4], int (*pLockdownMinUnit)[27][2][6], int (*pLockdownMedUnit)[11][2][4], int (*pLockdownMaxUnit)[18][4][5], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *pMoveFlag);
void unlock(gsl_rng **r, sIndiv **pHeadCopy, sIndiv **pTailCopy, int currDay, int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18], int (*pLockdownMinArea)[27][2], int (*pLockdownMedArea)[11][2], int (*pLockdownMaxArea)[18][4], int (*pLockdownMinUnit)[27][2][6], int (*pLockdownMedUnit)[11][2][4], int (*pLockdownMaxUnit)[18][4][5], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *pMoveFlag);
void newIndiv(int *idGlobal, int *HCVentry, int *HCVentryAb, int currDay, sIndiv **pHeadCopy, sIndiv **pTailCopy, gsl_rng **r, int nEvents, int (*pLocArray2)[ROWPRIS][COLCTR], int prison, int category, int seedPop, int seedZone, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int *povercap);
double probMoveCom(gsl_rng **r, sIndiv **pTarget, int currDay, int *pMoveFlag); //14
double probMoveP1(gsl_rng **r, sIndiv **pTarget, int currDay, int *pMoveFlag); //15
double probMoveP2(gsl_rng **r, sIndiv **pTarget); //16
double probMoveP3(gsl_rng **r, sIndiv **pTarget); //17
double probOutsideVisit(gsl_rng **r, sIndiv **pTarget, int currDay, int *pMoveFlag);
double probMoveHospitalF(gsl_rng **r, sIndiv **pTarget, int nMinHospitalF, int nMedHospitalF, int nMaxHospitalF);
double probMoveHospitalC(gsl_rng **r, sIndiv **pTarget, int nMinHospitalF, int nMedHospitalF, int nMaxHospitalF);
double probReturnPrison(gsl_rng **r, sIndiv **pTarget, int nMinHospitalF, int nMedHospitalF, int nMaxHospitalF);
double probTest(gsl_rng **r, sIndiv **pTarget, int currDay, int nInmateTests, int nPStaffTests, int nHStaffTests, int nInmatesTested, int nPStaffTested, int nHStaffTested, int nEVisitorsTested, int nFVisitorsTested);
double probTestResult(gsl_rng **r, sIndiv **pTarget, int currDay);
double probRapidTest(gsl_rng **r, sIndiv **pTarget, int currDay, int nTests);
double probThermalTest(gsl_rng **r, sIndiv **pTarget, int currDay, int nTests);
double probInfectIndiv(sIndiv **pTarget, int currDay, int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18], int (*pLockdownMinArea)[27][2], int (*pLockdownMedArea)[11][2], int (*pLockdownMaxArea)[18][4], int (*pLockdownMinUnit)[27][2][6], int (*pLockdownMedUnit)[11][2][4], int (*pLockdownMaxUnit)[18][4][5], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20]); //18
double probNtrClear(gsl_rng **r, sIndiv **pTarget, int currDay); //19
double probProgress(gsl_rng **r, sIndiv **pTarget, int currDay); //20
double probCOVIDdeath(gsl_rng **r, sIndiv **pTarget, int currDay); //21
double probNtrDeath(gsl_rng **r, sIndiv **pTarget); //0
double probIsolate(gsl_rng **r, sIndiv **pTarget, int currDay); //26
double probEndIsolate(gsl_rng **r, sIndiv **pTarget, int currDay);
double probQuarantine(gsl_rng **r, sIndiv **pTarget, int currDay);
double probEndQuarantine(gsl_rng **r, sIndiv **pTarget, int currDay);
double probCohorting(gsl_rng **r, sIndiv **pTarget, int currDay);
double probEndCohorting(gsl_rng **r, sIndiv **pTarget, int currDay);
double probPPE(gsl_rng **r, sIndiv **pTarget, int currDay);
double probEndPPE(gsl_rng **r, sIndiv **pTarget, int currDay);
double probGroupVulnerable(gsl_rng **r, sIndiv **pTarget, int currDay);
double probEndGroupVulnerable(gsl_rng **r, sIndiv **pTarget, int currDay);

void age(sIndiv **pTargetCopy);
void removeIndiv(sIndiv **pTarget, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20]);
void changeRisk(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int newRisk);
void moveLocation(gsl_rng **r, sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int newLoc, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20], int (*pIsolateMinZone)[27], int (*pIsolateMedZone)[11], int (*pIsolateMaxZone)[18], int (*pIsolateMinArea)[27][2], int (*pIsolateMedArea)[11][2], int (*pIsolateMaxArea)[18][4], int (*pIsolateMinUnit)[27][2][6], int (*pIsolateMedUnit)[11][2][4], int (*pIsolateMaxUnit)[18][4][5], int (*pIsolateMinCell)[27][2][6][13], int (*pIsolateMedCell)[11][2][4][19], int (*pIsolateMaxCell)[18][4][5][20], int (*pLockdownMinZone)[27], int (*pLockdownMedZone)[11], int (*pLockdownMaxZone)[18]);
void outsideVisit(gsl_rng **r, sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20]);
void moveHospital(sIndiv **pTargetCopy, int hospitalType);
void returnPrison(sIndiv **pTargetCopy, int (*pMinCellArray)[27][2][6][13], int (*pMedCellArray)[11][2][4][19], int (*pMaxCellArray)[18][4][5][20]);
void test(sIndiv **pTargetCopy, int currDay);
void testResult(sIndiv **pTargetCopy, int currDay);
int rapidTest(gsl_rng **r, sIndiv **pTargetCopy, int currDay);
int rapidTestPStaff(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pTailCopy, int currDay, int *PSTP, int *PSTN, int *PSFP, int *PSFN);
int rapidTestHStaff(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pTailCopy, int currDay, int *HSTP, int *HSTN, int *HSFP, int *HSFN);
void thermalTest(sIndiv **pTargetCopy);
int infect(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*zoneMinCIArray)[30][7], int (*zoneMedCIArray)[30][7], int (*zoneMaxCIArray)[30][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax, int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int currDay);
int infectOutside(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*courtCIArray)[38][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax, int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6,  int currDay);
int infectInTransit(gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*truckCIArray)[20][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax, int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6,  int currDay, int *pMoveFlag);
int infectRedZone (gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*zoneMinCIArray)[30][7], int (*zoneMedCIArray)[30][7], int (*zoneMaxCIArray)[30][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int currDay);
int StaffInfectRedZone (gsl_rng **r, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int (*zoneMinCIArray)[30][7], int (*zoneMedCIArray)[30][7], int (*zoneMaxCIArray)[30][7], int *newInfectedInmates, int *newInfectedPS, int *newInfectedHS, int *newInfectedEV, int *newInfectedFV, int *newInfectedInmatesMin, int *newInfectedPSMin, int *newInfectedHSMin, int *newInfectedEVMin, int *newInfectedFVMin, int *newInfectedInmatesMed, int *newInfectedPSMed, int *newInfectedHSMed, int *newInfectedEVMed, int *newInfectedFVMed, int *newInfectedInmatesMax, int *newInfectedPSMax, int *newInfectedHSMax, int *newInfectedEVMax, int *newInfectedFVMax,  int *totalAge0, int *totalAge1, int *totalAge2, int *totalAge3, int *totalAge4, int *totalAge5, int *totalAge6, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int currDay);
void clearHCV(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]);
void NtrClear(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]);
void progress(gsl_rng **r, sIndiv **pTargetCopy, int *exp19, int *exp20to44, int *exp45to54, int *exp55to64, int *exp65to74, int *exp75to84, int *exp85, int *pre19, int *pre20to44, int *pre45to54, int *pre55to64, int *pre65to74, int *pre75to84, int *pre85, int *asy19, int *asy20to44, int *asy45to54, int *asy55to64, int *asy65to74, int *asy75to84, int *asy85, int *mil19, int *mil20to44, int *mil45to54, int *mil55to64, int *mil65to74, int *mil75to84, int *mil85, int *mod19, int *mod20to44, int *mod45to54, int *mod55to64, int *mod65to74, int *mod75to84, int *mod85, int *sev19, int *sev20to44, int *sev45to54, int *sev55to64, int *sev65to74, int *sev75to84, int *sev85, int *cle19, int *cle20to44, int *cle45to54, int *cle55to64, int *cle65to74, int *cle75to84, int *cle85, int currDay);
//int deathHCV(sIndiv **pTargetCopy);
void startOST(sIndiv **pTargetCopy, int currDay);
void stopOST(sIndiv **pTargetCopy);
void startDAA(sIndiv **pTargetCopy, int currDay);
void stopDAA(sIndiv **pTargetCopy);
void clearDAA(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]);

int main(int argc, const char * argv[])
{
    //int newInfected, newCleared; //new cases per week/month
    FILE *fp, *fw, *fr, *fm;
    int idG=0, totalInfectedCommunity=0, daysSim, i, j, k, l, normCount, nShare, totalIndivDAA=0, totalIndivOST=0, totalClrN=0, totalClrD=0, fi, ni, newCases=0, success, currDay=1, totalPrisonPop=0, dead=0, deadHCV=0, deadHCVCom=0, deadHCVPris=0, prisonHCV=0, communityHCV=0, prisonOpd=0, prisonOpdNotOST=0, prisonEverIDU=0, releasedHCVAb=0, releasedHCVRNA=0, releasedHCVCom=0, releasedHCVPris=0, released=0, eIndex=999, eInf=999, deadFlag=0, eFlag, mFlag=0, dFlag=0, cnFlag=0, cdFlag=0, cEvents, cIDUp=0, cIDUn=0, NSWpop=6774000, inPop, HCVe, HCVeAb, exitAb=0, exitRNA=0; //total prison population, number of days to simulate; NSW pop=2005
    int E0=0, E1=0, E2=0, E3=0, E4=0, E5=0, E6=0, E7=0, E8=0, E9=0, E10=0, E11=0, E12=0, E13=0, E14=0, E15=0, E16=0, E17=0, E18=0, E19=0, E20=0, E21=0, E22=0, E23=0, E24=0, E25=0, E26=0, R0=0, R1=0, R2=0, R3=0, R4=0, R5=0, R6=0, R7=0, R8=0, R9=0, R10=0, R11=0, R12=0;//Event counters
    int nOST=0, nDAA=0, nReinfected, inIncrease, nDailyOST=0, nDaily=0, nOSTDailyCap=0, nDailyCap=0, dailyInPop=0, hcvAntibody=0, hcvRNA=0, aveLengthStay=0, releasedPSCov=0, releasedHSCov=0, releasedEVCov=0, releasedFVCov=0;
    int nMinHospitalF=0, nMedHospitalF=0, nMaxHospitalF=0, nMinHospitalC=0, nMedHospitalC=0, nMaxHospitalC=0;
    int nInmateTests, nPStaffTests, nHStaffTests, nRTests, nInmatesTested=0, nPStaffTested=0, nHStaffTested=0, nEVisitorsTested=0, nFVisitorsTested=0;//per day
    int nInmatesDetected=0, nInmatesUndetected=0, nPStaffDetected=0, nPStaffUndetected=0, nHStaffDetected=0, nHStaffUndetected=0, nEVisitorsDetected=0, nEVisitorsUndetected=0, nFVisitorsDetected=0, nFVisitorsUndetected=0; //per day
    int releasedI=0, releasedICov=0, releasedPS=0, releasedHS=0, releasedEV=0, releasedFV=0, infectedOutside=0, prevStatus=0, newStatus=0, nInmatesThermal=0, nPStaffThermal=0, nHStaffThermal=0, nEVisitorsThermal=0, nFVisitorsThermal=0, nInmatesRapidTested=0, nPStaffRapidTested=0, nHStaffRapidTested=0, nEVisitorsRapidTested=0, nFVisitorsRapidTested=0, inCtr; //per day
    int newInfectedInmates=0, newInfectedPS=0, newInfectedHS=0, newInfectedEV=0, newInfectedFV=0;
    int newInfectedInmatesMin=0, newInfectedPSMin=0, newInfectedHSMin=0, newInfectedEVMin=0, newInfectedFVMin=0;
    int newInfectedInmatesMed=0, newInfectedPSMed=0, newInfectedHSMed=0, newInfectedEVMed=0, newInfectedFVMed=0;
    int newInfectedInmatesMax=0, newInfectedPSMax=0, newInfectedHSMax=0, newInfectedEVMax=0, newInfectedFVMax=0;
    int i19=0, i20to44=0, i45to54=0, i55to64=0, i65to74=0, i75to84=0, i85=0, ps19=0, ps20to44=0, ps45to54=0, ps55to64=0, ps65to74=0, ps75to84=0, ps85=0, hs19=0, hs20to44=0, hs45to54=0, hs55to64=0, hs65to74=0, hs75to84=0, hs85=0, ev19=0, ev20to44=0, ev45to54=0, ev55to64=0, ev65to74=0, ev75to84=0, ev85=0, fv19=0, fv20to44=0, fv45to54=0, fv55to64=0, fv65to74=0, fv75to84=0, fv85=0;
    int inmateDeathCOVID=0, PSDeathCOVID=0, HSDeathCOVID=0, EVDeathCOVID=0, FVDeathCOVID=0, totalDeathCOVID=0;
    int finalSevExp=0, finalSevPre=0, finalSevAsy=0, finalSevMil=0, finalSevMod=0, finalSevSev=0, finalSevNon=0;
    int totalTestPS=0, totalTestHS=0, totalPSTP=0, totalHSTP=0, totalPSTN=0, totalHSTN=0, totalPSFP=0, totalHSFP=0, totalPSFN=0, totalHSFN=0, rapidResult=0;
    //add 0:exposed -> 1:pre-clinical -> 2:asymptomatic -> 3:mild -> 4:mod -> 5:severe, 6:clearer, 7: no infection
	gsl_rng *r;
    //int locArray[ROWPRIS][COLCTR]; //Array of locations containing: //[row][0] IDU+ //[row][1] IDU-
    //int (*pLocArray)[ROWPRIS][COLCTR]=&locArray;
    //double probLocMat[ROWPRIS][ROWPRIS]; //[row][col] for movement between locations
    //double probRiskMat[RGROUPS][RGROUPS]; //[row][col] for movement between locations
    double riskNormDenominator, movesNormDenominator, normDenominator, probEvents[NEVENTS], probInf[2];
    //arrays for events
	double eventProbs[27], eventDraw[2]; 
	int eventDecision[27], eventOrder[27], outArray[2];
	int eCtr, iFlag, iCrawl, iCount; 
    int tempLocArray[ROWPRIS][TCOLCTR]; //Array of locations containing: //[row][0] IDU+ //[row][1] IDU-
    int (*pTempLocArray)[ROWPRIS][TCOLCTR]=&tempLocArray;
    
    int locArray[ROWPRIS][COLCTR]; //Array of locations containing: [row][0] IDU+; HCV+; ATSI, [row][1] IDU+; HCV+; NON-ATSI, [row][2] IDU+; HCV-; ATSI; PREV. EXPOSED, [row][3] IDU+; HCV-; ATSI; SUSCEPTIBLE, [row][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED, [row][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE, [row][6] IDU-; HCV+; ATSI, [row][7] IDU-; HCV+; NON-ATSI, [row][8] IDU-; HCV-; ATSI; PREV. EXPOSED, [row][9] IDU-; HCV-; ATSI; SUSCEPTIBLE, [row][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED, [row][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE
    int (*pLocArray)[ROWPRIS][COLCTR]=&locArray;
    int exp19=0, exp20to44=0, exp45to54=0, exp55to64=0, exp65to74=0, exp75to84=0, exp85=0;
    int pre19=0, pre20to44=0, pre45to54=0, pre55to64=0, pre65to74=0, pre75to84=0, pre85=0;
    int asy19=0, asy20to44=0, asy45to54=0, asy55to64=0, asy65to74=0, asy75to84=0, asy85=0;
    int mil19=0, mil20to44=0, mil45to54=0, mil55to64=0, mil65to74=0, mil75to84=0, mil85=0;
    int mod19=0, mod20to44=0, mod45to54=0, mod55to64=0, mod65to74=0, mod75to84=0, mod85=0;
    int sev19=0, sev20to44=0, sev45to54=0, sev55to64=0, sev65to74=0, sev75to84=0, sev85=0;
    int cle19=0, cle20to44=0, cle45to54=0, cle55to64=0, cle65to74=0, cle75to84=0, cle85=0;
 
    printf( "Creating indivArray\n");
    
    int indivArray[30][5]; //[correctional centres][indivTypes]
    int (*pIndivArray)[30][5]=&indivArray;
    
    int indivCtrR;
    int indivCtrC;
    for(indivCtrR=0; indivCtrR<30; indivCtrR++){
        for(indivCtrC=0; indivCtrC<5; indivCtrC++){
            indivArray[indivCtrR][indivCtrC]=0;
        }
    }
 
    /*
    int isolateMaxZone[30], isolateMedZone[30], isolateMinZone[30], isolateMinArea[4], isolateMedArea[4], isolateMaxArea[4], isolateMinUnit[6], isolateMedUnit[6], isolateMaxUnit[6], isolateMinCell[20], isolateMedCell[20], isolateMaxCell[20], timerMax[30], timerMed[30], timerMin[30];
    
    int isolateMaxArray[30][4][6][20]; //zone, area, unit, cell
    int isolateMedArray[30][4][6][20];
    int isolateMinArray[30][4][6][20];
    
    int iZoneCtr, iAreaCtr, iUnitCtr, iCellCtr, iMCtr;
    //int iMCtr;
    for(iZoneCtr=0; iZoneCtr<30; iZoneCtr++){
        for(iAreaCtr=0; iAreaCtr<4; iAreaCtr++){
            for(iUnitCtr=0; iUnitCtr<6; iUnitCtr++){
                for(iCellCtr=0; iCellCtr<20; iCellCtr++){
                    //printf( "%d, %d, %d, %d\n", iZoneCtr, iAreaCtr, iUnitCtr, iCellCtr);
                    isolateMaxArray[iZoneCtr][iAreaCtr][iUnitCtr][iCellCtr]=0;
                    isolateMedArray[iZoneCtr][iAreaCtr][iUnitCtr][iCellCtr]=0;
                    isolateMinArray[iZoneCtr][iAreaCtr][iUnitCtr][iCellCtr]=0;
                }
            }
        }
    }
     */
    
    int MinCellArray[27][2][6][13]; //[correctional centres][indivTypes]
    int (*pMinCellArray)[27][2][6][13]=&MinCellArray;

    int isolateMinZone[27]; //[correctional centres][indivTypes]
    int isolateMinArea[27][2];
    int isolateMinUnit[27][2][6];
    int isolateMinCell[27][2][6][13];
    
    int (*pIsolateMinZone)[27]=&isolateMinZone;
    int (*pIsolateMinArea)[27][2]=&isolateMinArea;
    int (*pIsolateMinUnit)[27][2][6]=&isolateMinUnit;
    int (*pIsolateMinCell)[27][2][6][13]=&isolateMinCell;

    int lockdownMinZone[27]; //to check if a prison is in lockdown
    int lockdownMinArea[27][2];
    int lockdownMinUnit[27][2][6];
    
    int (*pLockdownMinZone)[27]=&lockdownMinZone;
    int (*pLockdownMinArea)[27][2]=&lockdownMinArea;
    int (*pLockdownMinUnit)[27][2][6]=&lockdownMinUnit;
    
    int mcCtr1, mcCtr2, mcCtr3, mcCtr4;
    for(mcCtr1=0; mcCtr1<27; mcCtr1++){ //correctional centres
        isolateMinZone[mcCtr1]=0; //isolate prison or not?
        lockdownMinZone[mcCtr1]=0;
        for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
            isolateMinArea[mcCtr1][mcCtr2]=0;
            lockdownMinArea[mcCtr1][mcCtr2]=0;
            for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                isolateMinUnit[mcCtr1][mcCtr2][mcCtr3]=0;
                lockdownMinUnit[mcCtr1][mcCtr2][mcCtr3]=0;
                for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                    isolateMinCell[mcCtr1][mcCtr2][mcCtr3][mcCtr4]=0; //by default 2-out cells
                    MinCellArray[mcCtr1][mcCtr2][mcCtr3][mcCtr4]=0;
                }
            }
        }
    }
    
    int MedCellArray[11][2][4][19]; //[correctional centres][indivTypes]
    int (*pMedCellArray)[11][2][4][19]=&MedCellArray;

    int isolateMedZone[11]; //[correctional centres][indivTypes]
    int isolateMedArea[11][2];
    int isolateMedUnit[11][2][4];
    int isolateMedCell[11][2][4][19];

    int lockdownMedZone[11]; //[correctional centres][indivTypes]
    int lockdownMedArea[11][2];
    int lockdownMedUnit[11][2][4];
    
    int (*pIsolateMedZone)[11]=&isolateMedZone;
    int (*pIsolateMedArea)[11][2]=&isolateMedArea;
    int (*pIsolateMedUnit)[11][2][4]=&isolateMedUnit;
    int (*pIsolateMedCell)[11][2][4][19]=&isolateMedCell;

    int (*pLockdownMedZone)[11]=&lockdownMedZone;
    int (*pLockdownMedArea)[11][2]=&lockdownMedArea;
    int (*pLockdownMedUnit)[11][2][4]=&lockdownMedUnit;
    
    for(mcCtr1=0; mcCtr1<11; mcCtr1++){ //correctional centres
        isolateMedZone[mcCtr1]=0; //isolate prison or not?
        lockdownMedZone[mcCtr1]=0;
        for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
            isolateMedArea[mcCtr1][mcCtr2]=0;
            lockdownMedArea[mcCtr1][mcCtr2]=0;
            for(mcCtr3=0; mcCtr3<4; mcCtr3++){ //min med max
                isolateMedUnit[mcCtr1][mcCtr2][mcCtr3]=0;
                lockdownMedUnit[mcCtr1][mcCtr2][mcCtr3]=0;
                for(mcCtr4=0; mcCtr4<19; mcCtr4++){ //min med max
                    isolateMedCell[mcCtr1][mcCtr2][mcCtr3][mcCtr4]=0;//by default 2-out cells
                    MedCellArray[mcCtr1][mcCtr2][mcCtr3][mcCtr4]=0;
                }
            }
        }
    }
    
    int MaxCellArray[18][4][5][20]; //[correctional centres][indivTypes]
    int (*pMaxCellArray)[18][4][5][20]=&MaxCellArray;

    int isolateMaxZone[18]; //[correctional centres][indivTypes]
    int isolateMaxArea[18][4];
    int isolateMaxUnit[18][4][5];
    int isolateMaxCell[18][4][5][20];

    int lockdownMaxZone[18]; //[correctional centres][indivTypes]
    int lockdownMaxArea[18][4];
    int lockdownMaxUnit[18][4][5];
    
    int (*pIsolateMaxZone)[18]=&isolateMaxZone;
    int (*pIsolateMaxArea)[18][4]=&isolateMaxArea;
    int (*pIsolateMaxUnit)[18][4][5]=&isolateMaxUnit;
    int (*pIsolateMaxCell)[18][4][5][20]=&isolateMaxCell;

    int (*pLockdownMaxZone)[18]=&lockdownMaxZone;
    int (*pLockdownMaxArea)[18][4]=&lockdownMaxArea;
    int (*pLockdownMaxUnit)[18][4][5]=&lockdownMaxUnit;
    
    int overcap=0;
    int *povercap=&overcap;
    
    int infectedRedZone=0;
    
    for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
        isolateMaxZone[mcCtr1]=0; //isolate yes, or no?
        lockdownMaxZone[mcCtr1]=0;
        for(mcCtr2=0; mcCtr2<4; mcCtr2++){ //min med max
            isolateMaxArea[mcCtr1][mcCtr2]=0;
            lockdownMaxArea[mcCtr1][mcCtr2]=0;
            for(mcCtr3=0; mcCtr3<5; mcCtr3++){ //min med max
                isolateMaxUnit[mcCtr1][mcCtr2][mcCtr3]=0;
                lockdownMaxUnit[mcCtr1][mcCtr2][mcCtr3]=0;
                for(mcCtr4=0; mcCtr4<20; mcCtr4++){ //min med max
                    isolateMaxCell[mcCtr1][mcCtr2][mcCtr3][mcCtr4]=0; //by default 2-out cells
                    MaxCellArray[mcCtr1][mcCtr2][mcCtr3][mcCtr4]=0;
                }
            }
        }
    }
    
    printf( "Creating zoneArray\n");
    
    int zoneArray[30][3];
    int zoneCovidArray[30][3];
    //int (*pZoneArray)[30][3]=&zoneArray;
    
    int zaCtrC, zaCtrR;
    for(zaCtrR=0; zaCtrR<30; zaCtrR++){ //correctional centres
        for(zaCtrC=0; zaCtrC<3; zaCtrC++){ //min med max
            zoneArray[zaCtrR][zaCtrC]=0;
            zoneCovidArray[zaCtrR][zaCtrC]=0;
        }
    }

    printf( "Creating zoneCIArray\n");

    int zoneMinDIArray[30][7]; //daily new infections
    int zoneMedDIArray[30][7];
    int zoneMaxDIArray[30][7];
    int zoneMinDDArray[30][7]; //daily deaths
    int zoneMedDDArray[30][7];
    int zoneMaxDDArray[30][7];
    int zoneMinCIArray[30][7]; //cumulative infections
    int zoneMedCIArray[30][7];
    int zoneMaxCIArray[30][7];
    int zoneMinCDArray[30][7]; //cumulative deaths
    int zoneMedCDArray[30][7];
    int zoneMaxCDArray[30][7];
    
    int zoneMaxPrevCI[30][7];
    int zoneMedPrevCI[30][7];
    int zoneMinPrevCI[30][7];
    //int (*pZoneArray)[30][3]=&zoneArray;
    
    int zInfCtrC, zInfCtrR;
    for(zInfCtrR=0; zInfCtrR<30; zInfCtrR++){ //correctional centres
        for(zInfCtrC=0; zInfCtrC<7; zInfCtrC++){ //min med max
            zoneMinDIArray[zInfCtrR][zInfCtrC]=0; //daily new
            zoneMedDIArray[zInfCtrR][zInfCtrC]=0;
            zoneMaxDIArray[zInfCtrR][zInfCtrC]=0;
            zoneMinDDArray[zInfCtrR][zInfCtrC]=0; //daily deaths
            zoneMedDDArray[zInfCtrR][zInfCtrC]=0;
            zoneMaxDDArray[zInfCtrR][zInfCtrC]=0;
            zoneMinCIArray[zInfCtrR][zInfCtrC]=0;
            zoneMedCIArray[zInfCtrR][zInfCtrC]=0;
            zoneMaxCIArray[zInfCtrR][zInfCtrC]=0;
            zoneMinCDArray[zInfCtrR][zInfCtrC]=0;
            zoneMedCDArray[zInfCtrR][zInfCtrC]=0;
            zoneMaxCDArray[zInfCtrR][zInfCtrC]=0;
            
            zoneMaxPrevCI[zInfCtrR][zInfCtrC]=0;
            zoneMedPrevCI[zInfCtrR][zInfCtrC]=0;
            zoneMinPrevCI[zInfCtrR][zInfCtrC]=0;
        }
    }

    int totalAge0=0, totalAge1=0, totalAge2=0, totalAge3=0, totalAge4=0, totalAge5=0, totalAge6=0;
    int totalAD0=0, totalAD1=0, totalAD2=0, totalAD3=0, totalAD4=0, totalAD5=0, totalAD6=0;

    printf( "Creating courtCIArray\n");
    
    int courtCIArray[38][7];
    int courtCDArray[38][7];
    int courtDIArray[38][7]; //daily infect
    int courtDDArray[38][7]; //daily death
    int nCourtPop[38];
    
    int courtPrevCI[38][7];
    
    int moveFlag=0;
    
    int cInfCtrC, cInfCtrR;
    for(cInfCtrR=0; cInfCtrR<38; cInfCtrR++){ //correctional centres
        nCourtPop[cInfCtrR]=0;
        printf( "Creating nCourtPop=%d\n", cInfCtrC);
        for(cInfCtrC=0; cInfCtrC<7; cInfCtrC++){ //min med max
            printf( "Creating nCourtPop=%d\n", cInfCtrC);
            courtCIArray[cInfCtrR][cInfCtrC]=0;
            courtCDArray[cInfCtrR][cInfCtrC]=0;
            courtDIArray[cInfCtrR][cInfCtrC]=0;
            courtDDArray[cInfCtrR][cInfCtrC]=0;
            
            courtPrevCI[cInfCtrR][cInfCtrC]=0;
        }
    }

    printf( "Creating truckCIArray\n");
    
    int truckCIArray[20][7];
    int truckDIArray[20][7];
    int nMovingPop[20];
    
    int truckPrevCI[20][7];
    
    int tInfCtrC, tInfCtrR;
    for(tInfCtrR=0; tInfCtrR<20; tInfCtrR++){ //correctional centres
        nMovingPop[tInfCtrR]=0;
        for(tInfCtrC=0; tInfCtrC<7; tInfCtrC++){ //min med max
            truckCIArray[tInfCtrR][tInfCtrC]=0;
            truckDIArray[tInfCtrR][tInfCtrC]=0;
            
            truckPrevCI[tInfCtrR][tInfCtrC]=0;
        }
    }
    
    //prepare population age breakdown
    int ageGroups=7;
    int iAgeArray[ageGroups][3]; //Array containing: [age group][0] with all, [1] with COVID, [2] with detected COVID
    //int (*pAgeArray)[ROWPRIS][COLCTR]=&ageArray;
    int psAgeArray[ageGroups][3];
    int hsAgeArray[ageGroups][3];
    int evAgeArray[ageGroups][3];
    int fvAgeArray[ageGroups][3];

    printf( "Creating age arrays\n");
    
    int agCtrC, agCtrR;
    for(agCtrR=0; agCtrR<ageGroups; agCtrR++){ //age groups
        for(agCtrC=0; agCtrC<3; agCtrC++){ //populations
            iAgeArray[agCtrR][agCtrC]=0;
            psAgeArray[agCtrR][agCtrC]=0;
            hsAgeArray[agCtrR][agCtrC]=0;
            evAgeArray[agCtrR][agCtrC]=0;
            fvAgeArray[agCtrR][agCtrC]=0;
        }
    }
    
    printf( "Creating severityArray\n");
    
    //age breakdown according to severity
    int ctrSeverity=8;  //add 0:exposed -> 1:pre-clinical -> 2:asymptomatic -> 3:mild -> 4:mod -> 5:severe, 6:clearer, 7: no infection
    int iAgeSeverityArray[ageGroups][ctrSeverity]; //7 rows by 8 columns

    int asCtrC, asCtrR;
    for(asCtrR=0; asCtrR<ageGroups; asCtrR++){ //age groups
        for(asCtrC=0; asCtrC<ctrSeverity; asCtrC++){ //populations
            iAgeSeverityArray[asCtrR][asCtrC]=0;
        }
    }
    
    int tAge=0, tCov=0, tDetCov=0, tAgeGroup=0, tAgePop=0;
    double fd;
    float rProp, lb, ub, inPopRate = 56.3;//40.8127854;//2.68;// inPopRateMin=2.67, inPopRateMax=3.236; //For population
    char option1[2]={"b"}, option2[2]={"d"}; //basic or detailed input

    printf( "Generating r\n");
    
    //Create a generator chosen by the environment variable gsl_rng_type *
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, (unsigned long)time(NULL));
    
    //printf("Days:");
    //scanf("%d", &daysSim);
    
    srand((unsigned int)time(NULL));
    
    //Set linked list parameters
    sIndiv *nHead, *nTail, *target;
    nHead=NULL;
    nTail=NULL;
    
    //printf("int: %f\n", gsl_ran_exponential(r, 0.0004593));
    //printf("int: %f\n", gsl_ran_beta(r, 1.28527, 2797.038));
    //printf("int: %f\n", gsl_ran_exponential(r, 0.0004593));
    //printf("int: %f\n", gsl_ran_beta(r, 1.28527, 2797.038));
    //printf("int: %f\n", gsl_ran_exponential(r, 0.0004593));
    //printf("int: %f\n", gsl_ran_beta(r, 1.28527, 2797.038));
    //exit(0);
    
    for(i=0;i<ROWPRIS;i++){ //read through prisons
        for(j=0;j<COLCTR;j++){
            locArray[i][j]=0;
        }
    }
            printf( "starting\n");
    if ( argc != 4) // There should be one option
    {
        printf( "Missing arguments\n");
        printf( "Input file name, Output file name, Run option\n");
        exit(0);
    }
    
    printf( "Input read: done");
    //printf( "Input: %s\n", argv[1] );
    //printf( "Output: %s\n", argv[2] );
    //printf( "Option: %s\n", argv[3] ); //basic or detailed; use d
    //exit(0);
     
    //use files
    //fp=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/input.txt", "r");
    fp=fopen(argv[1], "r");
    //fw=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/output.xls", "w");
    fw=fopen(argv[2], "w");
    
    //argv[3]="d";

    if(fw==NULL){
        printf("Cannot open file!\n");
        exit(1);
    }
    //fprintf(fw, "TIME\tC IDU+HCV+ATSI+\tC IDU+HCV+ATSI-\tC IDU+HCV-ATSI+PE\tC IDU+HCV-ATSI+S\tC IDU+HCV-ATSI-PE\tC IDU+HCV-ATSI-S\tC IDU-HCV+ATSI+\tC IDU-HCV+ATSI-\tC IDU-HCV-ATSI+PE\tC IDU-HCV-ATSI+S\tC IDU-HCV-ATSI-PE\tC IDU-HCV-ATSI-S\tP1 IDU+HCV+ATSI+\tP1 IDU+HCV+ATSI-\tP1 IDU+HCV-ATSI+PE\tP1 IDU+HCV-ATSI+S\tP1 IDU+HCV-ATSI-PE\tP1 IDU+HCV-ATSI-S\tP1 IDU-HCV+ATSI+\tP1 IDU-HCV+ATSI-\tP1 IDU-HCV-ATSI+PE\tP1 IDU-HCV-ATSI+S\tP1 IDU-HCV-ATSI-PE\tP1 IDU-HCV-ATSI-S\tP2 IDU+HCV+ATSI+\tP2 IDU+HCV+ATSI-\tP2 IDU+HCV-ATSI+PE\tP2 IDU+HCV-ATSI+S\tP2 IDU+HCV-ATSI-PE\tP2 IDU+HCV-ATSI-S\tP2 IDU-HCV+ATSI+\tP2 IDU-HCV+ATSI-\tP2 IDU-HCV-ATSI+PE\tP2 IDU-HCV-ATSI+S\tP2 IDU-HCV-ATSI-PE\tP2 IDU-HCV-ATSI-S\tP3 IDU+HCV+ATSI+\tP3 IDU+HCV+ATSI-\tP3 IDU+HCV-ATSI+PE\tP3 IDU+HCV-ATSI+S\tP3 IDU+HCV-ATSI-PE\tP3 IDU+HCV-ATSI-S\tP3 IDU-HCV+ATSI+\tP3 IDU-HCV+ATSI-\tP3 IDU-HCV-ATSI+PE\tP3 IDU-HCV-ATSI+S\tP3 IDU-HCV-ATSI-PE\tP3 IDU-HCV-ATSI-S\tNEW CASES\tPRISON HCV\tCOMMUNITY HCV\tRELEASED\tRELEASED HCV Ab\tRELEASED HCV RNA\tRELEASED HCV COM\tRELEASED HCV PRIS\tCLEAR HCV\tDEATHS\tHCV DEATHS\tHCV COM DEATH\tHCV PRISON DEATH\tTOTAL PRISON POP\tE0\tE1\tE2\tE3\tE4\tE5\tE6\tE7\tE8\tE9\tE10\tE11\tE12\tE13\tE14\tE15\tE16\tE17\tE18\tE19\tE20\tE21\tE22\tE23\tE24\tE25\tE26\tR0\tR1\tR2\tR3\tR4\tR5\tR6\tR7\tR8\tR9\tR10\tR11\tR12\teverIDU\tOPD\tOpdNotOST\ttotalIndivOST\ttotalIndivDAA\ttotalClrD\tnReinfected\tDAAremaining\tOSTremaining\tHCVentryRNA\tHCVentryAb\tNewEntry\thcvAntibody\thcvRNA\tAveStay\n");

    fprintf(fw, "TIME\tC CM+COVID+ATSI+\tC CM+COVID+ATSI-\tC CM+COVID-ATSI+Ab+\tC CM+COVID-ATSI+Ab-\tC CM+COVID-ATSI-Ab+\tC CM+COVID-ATSI-Ab-\tC CM-COVID+ATSI+\tC CM-COVID+ATSI-\tC CM-COVID-ATSI+Ab+\tC CM-COVID-ATSI+Ab-\tC CM-COVID-ATSI-Ab+\tC CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P1I CM+COVID+ATSI+\tP1I  CM+COVID+ATSI-\tP1I  CM+COVID-ATSI+Ab+\tP1I  CM+COVID-ATSI+Ab-\tP1I  CM+COVID-ATSI-Ab+\tP1I  CM+COVID-ATSI-Ab-\tP1I  CM-COVID+ATSI+\tP1I  CM-COVID+ATSI-\tP1I  CM-COVID-ATSI+Ab+\tP1I  CM-COVID-ATSI+Ab-\tP1I  CM-COVID-ATSI-Ab+\tP1I  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P1PS CM+COVID+ATSI+\tP1PS  CM+COVID+ATSI-\tP1PS  CM+COVID-ATSI+Ab+\tP1PS  CM+COVID-ATSI+Ab-\tP1PS  CM+COVID-ATSI-Ab+\tP1PS  CM+COVID-ATSI-Ab-\tP1PS  CM-COVID+ATSI+\tP1PS  CM-COVID+ATSI-\tP1PS  CM-COVID-ATSI+Ab+\tP1PS  CM-COVID-ATSI+Ab-\tP1PS  CM-COVID-ATSI-Ab+\tP1PS  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P1HS CM+COVID+ATSI+\tP1HS  CM+COVID+ATSI-\tP1HS  CM+COVID-ATSI+Ab+\tP1HS  CM+COVID-ATSI+Ab-\tP1HS  CM+COVID-ATSI-Ab+\tP1HS  CM+COVID-ATSI-Ab-\tP1HS  CM-COVID+ATSI+\tP1HS  CM-COVID+ATSI-\tP1HS  CM-COVID-ATSI+Ab+\tP1HS  CM-COVID-ATSI+Ab-\tP1HS  CM-COVID-ATSI-Ab+\tP1HS  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P1EV CM+COVID+ATSI+\tP1EV  CM+COVID+ATSI-\tP1EV  CM+COVID-ATSI+Ab+\tP1EV  CM+COVID-ATSI+Ab-\tP1EV  CM+COVID-ATSI-Ab+\tP1EV  CM+COVID-ATSI-Ab-\tP1EV  CM-COVID+ATSI+\tP1EV  CM-COVID+ATSI-\tP1EV  CM-COVID-ATSI+Ab+\tP1EV  CM-COVID-ATSI+Ab-\tP1EV  CM-COVID-ATSI-Ab+\tP1EV  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P1FV CM+COVID+ATSI+\tP1FV  CM+COVID+ATSI-\tP1FV  CM+COVID-ATSI+Ab+\tP1FV  CM+COVID-ATSI+Ab-\tP1FV  CM+COVID-ATSI-Ab+\tP1FV  CM+COVID-ATSI-Ab-\tP1FV  CM-COVID+ATSI+\tP1FV  CM-COVID+ATSI-\tP1FV  CM-COVID-ATSI+Ab+\tP1FV  CM-COVID-ATSI+Ab-\tP1FV  CM-COVID-ATSI-Ab+\tP1FV  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P2I CM+COVID+ATSI+\tP2I  CM+COVID+ATSI-\tP2I  CM+COVID-ATSI+Ab+\tP2I  CM+COVID-ATSI+Ab-\tP2I  CM+COVID-ATSI-Ab+\tP2I  CM+COVID-ATSI-Ab-\tP2I  CM-COVID+ATSI+\tP2I  CM-COVID+ATSI-\tP2I  CM-COVID-ATSI+Ab+\tP2I  CM-COVID-ATSI+Ab-\tP2I  CM-COVID-ATSI-Ab+\tP2I  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P2PS CM+COVID+ATSI+\tP2PS  CM+COVID+ATSI-\tP2PS  CM+COVID-ATSI+Ab+\tP2PS  CM+COVID-ATSI+Ab-\tP2PS  CM+COVID-ATSI-Ab+\tP2PS  CM+COVID-ATSI-Ab-\tP2PS  CM-COVID+ATSI+\tP2PS  CM-COVID+ATSI-\tP2PS  CM-COVID-ATSI+Ab+\tP2PS  CM-COVID-ATSI+Ab-\tP2PS  CM-COVID-ATSI-Ab+\tP2PS  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P2HS CM+COVID+ATSI+\tP2HS  CM+COVID+ATSI-\tP2HS  CM+COVID-ATSI+Ab+\tP2HS  CM+COVID-ATSI+Ab-\tP2HS  CM+COVID-ATSI-Ab+\tP2HS  CM+COVID-ATSI-Ab-\tP2HS  CM-COVID+ATSI+\tP2HS  CM-COVID+ATSI-\tP2HS  CM-COVID-ATSI+Ab+\tP2HS  CM-COVID-ATSI+Ab-\tP2HS  CM-COVID-ATSI-Ab+\tP2HS  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P2EV CM+COVID+ATSI+\tP2EV  CM+COVID+ATSI-\tP2EV  CM+COVID-ATSI+Ab+\tP2EV  CM+COVID-ATSI+Ab-\tP2EV  CM+COVID-ATSI-Ab+\tP2EV  CM+COVID-ATSI-Ab-\tP2EV  CM-COVID+ATSI+\tP2EV  CM-COVID+ATSI-\tP2EV  CM-COVID-ATSI+Ab+\tP2EV  CM-COVID-ATSI+Ab-\tP2EV  CM-COVID-ATSI-Ab+\tP2EV  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P2FV CM+COVID+ATSI+\tP2FV  CM+COVID+ATSI-\tP2FV  CM+COVID-ATSI+Ab+\tP2FV  CM+COVID-ATSI+Ab-\tP2FV  CM+COVID-ATSI-Ab+\tP2FV  CM+COVID-ATSI-Ab-\tP2FV  CM-COVID+ATSI+\tP2FV  CM-COVID+ATSI-\tP2FV  CM-COVID-ATSI+Ab+\tP2FV  CM-COVID-ATSI+Ab-\tP2FV  CM-COVID-ATSI-Ab+\tP2FV  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P3I CM+COVID+ATSI+\tP3I  CM+COVID+ATSI-\tP3I  CM+COVID-ATSI+Ab+\tP3I  CM+COVID-ATSI+Ab-\tP3I  CM+COVID-ATSI-Ab+\tP3I  CM+COVID-ATSI-Ab-\tP3I  CM-COVID+ATSI+\tP3I  CM-COVID+ATSI-\tP3I  CM-COVID-ATSI+Ab+\tP3I  CM-COVID-ATSI+Ab-\tP3I  CM-COVID-ATSI-Ab+\tP3I  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P3PS CM+COVID+ATSI+\tP3PS  CM+COVID+ATSI-\tP3PS  CM+COVID-ATSI+Ab+\tP3PS  CM+COVID-ATSI+Ab-\tP3PS  CM+COVID-ATSI-Ab+\tP3PS  CM+COVID-ATSI-Ab-\tP3PS  CM-COVID+ATSI+\tP3PS  CM-COVID+ATSI-\tP3PS  CM-COVID-ATSI+Ab+\tP3PS  CM-COVID-ATSI+Ab-\tP3PS  CM-COVID-ATSI-Ab+\tP3PS  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P3HS CM+COVID+ATSI+\tP3HS  CM+COVID+ATSI-\tP3HS  CM+COVID-ATSI+Ab+\tP3HS  CM+COVID-ATSI+Ab-\tP3HS  CM+COVID-ATSI-Ab+\tP3HS  CM+COVID-ATSI-Ab-\tP3HS  CM-COVID+ATSI+\tP3HS  CM-COVID+ATSI-\tP3HS  CM-COVID-ATSI+Ab+\tP3HS  CM-COVID-ATSI+Ab-\tP3HS  CM-COVID-ATSI-Ab+\tP3HS  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P3EV CM+COVID+ATSI+\tP3EV  CM+COVID+ATSI-\tP3EV  CM+COVID-ATSI+Ab+\tP3EV  CM+COVID-ATSI+Ab-\tP3EV  CM+COVID-ATSI-Ab+\tP3EV  CM+COVID-ATSI-Ab-\tP3EV  CM-COVID+ATSI+\tP3EV  CM-COVID+ATSI-\tP3EV  CM-COVID-ATSI+Ab+\tP3EV  CM-COVID-ATSI+Ab-\tP3EV  CM-COVID-ATSI-Ab+\tP3EV  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "P3FV CM+COVID+ATSI+\tP3FV  CM+COVID+ATSI-\tP3FV  CM+COVID-ATSI+Ab+\tP3FV  CM+COVID-ATSI+Ab-\tP3FV  CM+COVID-ATSI-Ab+\tP3FV  CM+COVID-ATSI-Ab-\tP3FV  CM-COVID+ATSI+\tP3FV  CM-COVID+ATSI-\tP3FV  CM-COVID-ATSI+Ab+\tP3FV  CM-COVID-ATSI+Ab-\tP3FV  CM-COVID-ATSI-Ab+\tP3FV  CM-COVID-ATSI-Ab-\t");
    fprintf(fw, "NEW TOTAL CASES\tNEW CASES INMATES\tNEW CASES PS\tNEW CASES HS\tNEW CASES EV\tNEW CASES FV\t");
    fprintf(fw, "NEW CASES INMATES MIN\tNEW CASES PS MIN\tNEW CASES HS MIN\tNEW CASES EV MIN\tNEW CASES FV MIN\t");
    fprintf(fw, "NEW CASES INMATES MED\tNEW CASES PS MED\tNEW CASES HS MED\tNEW CASES EV MED\tNEW CASES FV MED\t");
    fprintf(fw, "NEW CASES INMATES MAX\tNEW CASES PS MAX\tNEW CASES HS MAX\tNEW CASES EV MAX\tNEW CASES FV MAX\t");
    fprintf(fw, "InfectedOutVisits\treleasedI\treleasedICov\treleasedPS\treleasedPSCov\treleasedHS\treleasedHSCov\t releasedEV\treleasedEVCov\treleasedFV\treleasedFVCov\t");
    fprintf(fw, "nInmatesTested\tnPStaffTested\tnHStaffTested\tnEVisitorsTested\tnFVisitorsTested\tnInmatesThermal\tnPStaffThermal\tnHStaffThermal\tnEVisitorsThermal\tnFVisitorsThermal\tnInmatesDetected\t nPStaffDetected\tnHStaffDetected\tnEVisitorsDetected\tnFVisitorsDetected\tnInmatesRapidTested\tnPStaffRapidTested\t nHStaffRapidTested\tnEVisitorsRapidTested\tnFVisitorsRapidTested\t");
    fprintf(fw, "nMinHospitalF\tnMedHospitalF\tnMaxHospitalF\tnMinHospitalC\tnMedHospitalC\tnMaxHospitalC\t");
    //int ageArray[ageGroups][2]; //Array containing: [age group][0] with all, [1] with COVID, [2] with detected
    fprintf(fw, "i19\t, i20to44\t, i45to54\t, i55to64\t, i65to74\t, i75to84\t, i85\t, i19c\t, i20to44c\t, i45to54c\t, i55to64c\t, i65to74c\t, i75to84c\t, i85c\t, i19d\t, i20to44d\t, i45to54d\t, i55to64d\t, i65to74d\t, i75to84d\t, i85d\t");
    fprintf(fw, "ps19\t, ps20to44\t, ps45to54\t, ps55to64\t, ps65to74\t, ps75to84\t, ps85\t, ps19c\t, ps20to44c\t, ps45to54c\t, ps55to64c\t, ps65to74c\t, ps75to84c\t, ps85c\t, ps19d\t, ps20to44d\t, ps45to54d\t, ps55to64d\t, ps65to74d\t, ps75to84d\t, ps85d\t");
    fprintf(fw, "hs19\t, hs20to44\t, hs45to54\t, hs55to64\t, hs65to74\t, hs75to84\t, hs85\t, hs19c\t, hs20to44c\t, hs45to54c\t, hs55to64c\t, hs65to74c\t, hs75to84c\t, hs85c\t, hs19d\t, hs20to44d\t, hs45to54d\t, hs55to64d\t, hs65to74d\t, hs75to84d\t, hs85d\t");
    fprintf(fw, "ev19\t, ev20to44\t, ev45to54\t, ev55to64\t, ev65to74\t, ev75to84\t, ev85\t, ev19c\t, ev20to44c\t, ev45to54c\t, ev55to64c\t, ev65to74c\t, ev75to84c\t, ev85c\t, ev19d\t, ev20to44d\t, ev45to54d\t, ev55to64d\t, ev65to74d\t, ev75to84d\t, ev85d\t");
    fprintf(fw, "fv19\t, fv20to44\t, fv45to54\t, fv55to64\t, fv65to74\t, fv75to84\t, fv85\t, fv19c\t, fv20to44c\t, fv45to54c\t, fv55to64c\t, fv65to74c\t, fv75to84c\t, fv85c\t, fv19d\t, fv20to44d\t, fv45to54d\t, fv55to64d\t, fv65to74d\t, fv75to84d\t, fv85d\t");
    fprintf(fw, "COMMUNITY HCV\tRELEASED\tRELEASED HCV Ab\tRELEASED HCV RNA\tRELEASED HCV COM\tRELEASED HCV PRIS\tCLEAR HCV\tDEATHS\tHCV DEATHS\tHCV COM DEATH\tHCV PRISON DEATH\tTOTAL PRISON POP\tE0\tE1\tE2\tE3\tE4\tE5\tE6\tE7\tE8\tE9\tE10\tE11\tE12\tE13\tE14\tE15\tE16\tE17\tE18\tE19\tE20\tE21\tE22\tE23\tE24\tE25\tE26\tR0\tR1\tR2\tR3\tR4\tR5\tR6\tR7\tR8\tR9\tR10\tR11\tR12\teverIDU\tOPD\tOpdNotOST\ttotalIndivOST\ttotalIndivDAA\ttotalClrD\tnReinfected\tDAAremaining\tOSTremaining\tHCVentryRNA\tHCVentryAb\tNewEntry\thcvAntibody\thcvRNA\tAveStay\t");
      fprintf(fw, "zone100\t, zone101\t,zone102\t,zone103\t,zone104\t,zone105\t,zone106\t,zone107\t,zone108\t,zone109\t,zone110\t, zone111\t, zone112\t,zone113\t,zone114\t,zone115\t,zone116\t,zone117\t,zone118\t,zone119\t,zone120\t,zone121\t, zone122\t,zone123\t,zone124\t,zone125\t, zone126\t,zone127\t,zone128\t,zone129\t");
        fprintf(fw, "zone200\t, zone201\t,zone202\t,zone203\t,zone204\t,zone205\t,zone206\t,zone207\t,zone208\t,zone209\t,zone210\t, zone211\t, zone212\t,zone213\t,zone214\t,zone215\t,zone216\t,zone217\t,zone218\t,zone219\t,zone220\t,zone221\t, zone222\t,zone223\t,zone224\t,zone225\t, zone226\t,zone227\t,zone228\t,zone229\t");
        fprintf(fw, "zone300\t, zone301\t,zone302\t,zone303\t,zone304\t,zone305\t,zone306\t,zone307\t,zone308\t,zone309\t,zone310\t, zone311\t, zone312\t,zone313\t,zone314\t,zone315\t,zone316\t,zone317\t,zone318\t,zone319\t,zone320\t,zone321\t, zone322\t,zone323\t,zone324\t,zone325\t, zone326\t,zone327\t,zone328\t,zone329\t");
    fprintf(fw, "zc100\tzc101\tzc102\tzc103\tzc104\tzc105\tzc106\tzc107\tzc108\tzc109\tzc110\t zc111\t zc112\tzc113\tzc114\tzc115\tzc116\tzc117\tzc118\tzc119\tzc120\tzc121\t zc122\tzc123\tzc124\tzc125\tzc126\tzc127\tzc128\tzc129\t");
    fprintf(fw, "zc200\tzc201\tzc202\tzc203\tzc204\tzc205\tzc206\tzc207\tzc208\tzc209\tzc210\t zc211\t zc212\tzc213\tzc214\tzc215\tzc216\tzc217\tzc218\tzc219\tzc220\tzc221\t zc222\tzc223\tzc224\tzc225\tzc226\tzc227\tzc228\tzc229\t");
    fprintf(fw, "zc300\tzc301\tzc302\tzc303\tzc304\tzc305\tzc306\tzc307\tzc308\tzc309\tzc310\tzc311\t zc312\tzc313\tzc314\tzc315\tzc316\tzc317\tzc318\tzc319\tzc320\tzc321\tzc322\tzc323\tzc324\tzc325\tzc326\tzc327\tzc328\tzc329\t");
    fprintf(fw, "i19s0\ti19s1\ti19s2\ti19s3\ti19s4\ti19s5\t");
    fprintf(fw, "i20to44s0\ti20to44s1\tii20to44s2\ti20to44s3\ti20to44s4\ti20to44s5\t");
    fprintf(fw, "i45to54s0\ti45to54s1\ti45to54s2\ti45to54s3\ti45to54s4\ti45to54s5\t");
    fprintf(fw, "i55to64s0\ti55to64s1\ti55to64s2\ti55to64s3\ti55to64s4\ti55to64s5\t");
    fprintf(fw, "i65to74s0\ti65to74s1\ti65to74s2\ti65to74s3\ti65to74s4\ti65to74s5\t");
    fprintf(fw, "i75to84s0\ti75to84s1\ti75to84s2\ti75to84s3\ti75to84s4\ti75to84s5\t");
    fprintf(fw, "i85s0\ti85s1\ti85s2\ti85s3\ti85s4\ti85s5\t");
    fprintf(fw, "z100_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z101_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z102_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z103_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z104_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z105_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z106_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z107_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z108_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z109_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z110_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z111_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z112_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z113_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z114_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z115_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z116_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z117_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z118_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z119_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z120_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z121_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z122_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z123_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z124_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z125_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z126_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z127_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z128_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z129_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t");
    fprintf(fw, "z200_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z201_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z202_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z203_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z204_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z205_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z206_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z207_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z208_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z209_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z210_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z211_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z212_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z213_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z214_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z215_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z216_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z217_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z218_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z219_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z220_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z221_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z222_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z223_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z224_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z225_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z226_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z227_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z228_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z229_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t");
    fprintf(fw, "z300_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z301_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z302_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z303_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z304_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z305_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z306_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z307_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z308_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z309_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z310_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z311_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z312_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z313_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z314_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z315_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z316_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z317_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z318_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z319_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z320_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z321_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z322_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z323_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z324_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z325_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z326_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z327_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z328_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, z329_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t");
    fprintf(fw, "zd100_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd101_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd102_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd103_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd104_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t,  zd105_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd106_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd107_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd108_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd109_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd110_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd111_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd112_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd113_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd114_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd115_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd116_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd117_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd118_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd119_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd120_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd121_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd122_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd123_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd124_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd125_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd126_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd127_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd128_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd129_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t");
    fprintf(fw, "zd200_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd201_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd202_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd203_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd204_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t,  zd205_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd206_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd207_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd208_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd209_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd210_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd211_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd212_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd213_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd214_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd215_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd216_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd217_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd218_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd219_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd220_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd221_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd222_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd223_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd224_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd225_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd226_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd227_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd228_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd229_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t");
    fprintf(fw, "zd300_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd301_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd302_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd303_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd304_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t,  zd305_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd306_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd307_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd308_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd309_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd310_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd311_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd312_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd313_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd314_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd315_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd316_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd317_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd318_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd319_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd320_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd321_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd322_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd323_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd324_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd325_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd326_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd327_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd328_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t, zd329_19\t, 20to44\t, 45to54\t, 55to64\t, 65to74\t, 75to84\t, 85\t");
    fprintf(fw, "inf 19\tinf 20 to 44\tinf 45 to 54\tinf 55 to 64\tinf 65to74\tinf 75to84\tinf 85\t");
    fprintf(fw, "dth 19\tdth 20 to 44\tdth 45 to 54\tdth 55 to 64\tdth 65to74\tdth 75to84\tdth 85\t");
    
    fprintf(fw, "court1\tcourt2\tcourt3\tcourt4\tcourt5\tcourt6\tcourt7\tcourt8\tcourt9\tcourt10\tcourt11\tcourt12\tcourt13\tcourt14\tcourt15\tcourt16\tcourt17\tcourt18\tcourt19\tcourt20\tcourt21\tcourt22\tcourt23\tcourt24\tcourt25\tcourt26\tcourt27\tcourt28\tcourt29\tcourt30\tcourt31\tcourt32\tcourt33\tcourt34\tcourt35\tcourt36\tcourt37\tcourt38\t");
    fprintf(fw, "ci1\tci2\tci3\tci4\tci5\tci6\tci7\tci8\tci9\tci10\tci11\tci12\tci13\tci14\tci15\tci16\tci17\tci18\tci19\tci20\tci21\tci22\tci23\tci24\tci25\tci26\tci27\tci28\tci29\tci30\tci31\tci32\tci33\tci34\tci35\tci36\tci37\tci38\t");
    fprintf(fw, "truck1\ttruck2\ttruck3\ttruck4\ttruck5\ttruck6\ttruck7\ttruck8\ttruck9\ttruck10\ttruck11\ttruck12\ttruck13\ttruck14\ttruck15\ttruck16\ttruck17\ttruck18\ttruck19\ttruck20\t");
    fprintf(fw, "ti1\tti2\tti3\tti4\tti5\tti6\tti7\tti8\tti9\tti10\tti11\tti12\tti13\tti14\tti15\tti16\tti17\tti18\tti19\tti20\t");
    //print exposed19 \texposed20to44 \texposed45to54 \texposed55to64 \texposed65to74 \texposed75-84 \texposed85
    fprintf(fw, "exposed19\texposed20to44\texposed45to54\texposed55to64\texposed65to74\texposed75to84\texposed85\t");
    
    //print asymptomatic19 \tasymptomatic20to44 \tasymptomatic45to54 \tasymptomatic55to64 \tasymptomatic65to74 \tasymptomatic75-84 \tasymptomatic85
    fprintf(fw, "preclinical19\tpreclinical20to44\tpreclinical45to54\tpreclinical55to64\tpreclinical65to74\tpreclinical75-84\tpreclinical85\t");
    
    fprintf(fw, "asymptomatic19\tasymptomatic20to44\tasymptomatic45to54\tasymptomatic55to64\tasymptomatic65to74\tasymptomatic75-84\tasymptomatic85\t");

    //print mild19 \tmild20to44 \tmild45to54 \tmild55to64 \tmild65to74 \tmild75-84 \tmild85
    fprintf(fw, "mild19\tmild20to44\tmild45to54\tmild55to64\tmild65to74\tmild75-84\tmild85\t");
    
    //print moderate19 \tmoderate20to44 \tmoderate45to54 \tmoderate55to64 \tmoderate65to74 \tmoderate75-84 \tmoderate85
    fprintf(fw, "moderate19\tmoderate20to44\tmoderate45to54\tmoderate55to64\tmoderate65to74\tmoderate75to84\tmoderate85\t");
    
    //print severe19 \tsevere20to44 \tsevere45to54 \tsevere55to64 \tsevere65to74 \tsevere75-84 \tsevere85
    fprintf(fw, "severe19\tsevere20to44\tsevere45to54\tsevere55to64\tsevere65to74\tsevere75to84\tsevere85\t");
    
    //print cleared19 \tcleared20to44 \tcleared45to54 \tcleared55to64 \tcleared65to74 \tcleared75-84 \tcleared85
    fprintf(fw, "cleared19\tcleared20to44\tcleared45to54\tcleared55to64\tcleared65to74\tcleared75to84\tcleared85\t");
    
    //int exp19=0, exp20to44=0, exp45to54=0, exp55to64=0, exp65to74=0, exp75to84=0, exp85=0;
    //int asy19=0, asy20to44=0, asy45to54=0, asy55to64=0, asy65to74=0, asy75to84=0, asy85=0;
    //int mil19=0, mil20to44=0, mil45to54=0, mil55to64=0, mil65to74=0, mil75to84=0, mil85=0;
    //int mod19=0, mod20to44=0, mod45to54=0, mod55to64=0, mod65to74=0, mod75to84=0, mod85=0;
    //int sev19=0, sev20to44=0, sev45to54=0, sev55to64=0, sev65to74=0, sev75to84=0, sev85=0;
    //int cle19=0, cle20to44=0, cle45to54=0, cle55to64=0, cle65to74=0, cle75to84=0, cle85=0;
    
    fprintf(fw, "z100_total\tz101_total\tz102_total\tz103_total\tz104_total\tz105_total\tz106_total\tz107_total\tz108_total\tz109_total\tz110_total\tz111_total\tz112_total\tz113_total\tz114_total\tz115_total\tz116_total\tz117_total\tz118_total\tz119_total\tz120_total\tz121_total\tz122_total\tz123_total\tz124_total\tz125_total\tz126_total\tz127_total\tz128_total\tz129_total\t");
    
    fprintf(fw, "z200_total\tz201_total\tz202_total\tz203_total\tz204_total\tz205_total\tz206_total\tz207_total\tz208_total\tz209_total\tz210_total\tz211_total\tz212_total\tz213_total\tz214_total\tz215_total\tz216_total\tz217_total\tz218_total\tz219_total\tz220_total\tz221_total\tz222_total\tz223_total\tz224_total\tz225_total\tz226_total\tz227_total\tz228_total\tz229_total\t");
    
    fprintf(fw, "z300_total\tz301_total\tz302_total\tz303_total\tz304_total\tz305_total\tz306_total\tz307_total\tz308_total\tz309_total\tz310_total\tz311_total\tz312_total\tz313_total\tz314_total\tz315_total\tz316_total\tz317_total\tz318_total\tz319_total\tz320_total\tz321_total\tz322_total\tz323_total\tz324_total\tz325_total\tz326_total\tz327_total\tz328_total\tz329_total\t");
    
        //int finalSevExp=0, finalSevPre=0, finalSecAsy=0, finalSecMil=0, finalSevMod=0, finalSevSev=0, finalSevNon=0;
    fprintf(fw, "finalExposed\tfinalPreClinical\tfinalAsymptomatic\tfinalMild\tfinalModerate\tfinalSevere\tfinalNonInfected\t");
    
    fprintf(fw, "z100_dailyInf\tz101_dailyInf\tz102_dailyInf\tz103_dailyInf\tz104_dailyInf\tz105_dailyInf\tz106_dailyInf\tz107_dailyInf\tz108_dailyInf\tz109_dailyInf\tz110_dailyInf\tz111_dailyInf\tz112_dailyInf\tz113_dailyInf\tz114_dailyInf\tz115_dailyInf\tz116_dailyInf\tz117_dailyInf\tz118_dailyInf\tz119_dailyInf\tz120_dailyInf\tz121_dailyInf\tz122_dailyInf\tz123_dailyInf\tz124_dailyInf\tz125_dailyInf\tz126_dailyInf\tz127_dailyInf\tz128_dailyInf\tz129_dailyInf\t");
    
    fprintf(fw, "z200_dailyInf\tz201_dailyInf\tz202_dailyInf\tz203_dailyInf\tz204_dailyInf\tz205_dailyInf\tz206_dailyInf\tz207_dailyInf\tz208_dailyInf\tz209_dailyInf\tz210_dailyInf\tz211_dailyInf\tz212_dailyInf\tz213_dailyInf\tz214_dailyInf\tz215_dailyInf\tz216_dailyInf\tz217_dailyInf\tz218_dailyInf\tz219_dailyInf\tz220_dailyInf\tz221_dailyInf\tz222_dailyInf\tz223_dailyInf\tz224_dailyInf\tz225_dailyInf\tz226_dailyInf\tz227_dailyInf\tz228_dailyInf\tz229_dailyInf\t");
    
    fprintf(fw, "z300_dailyInf\tz301_dailyInf\tz302_dailyInf\tz303_dailyInf\tz304_dailyInf\tz305_dailyInf\tz306_dailyInf\tz307_dailyInf\tz308_dailyInf\tz309_dailyInf\tz310_dailyInf\tz311_dailyInf\tz312_dailyInf\tz313_dailyInf\tz314_dailyInf\tz315_dailyInf\tz316_dailyInf\tz317_dailyInf\tz318_dailyInf\tz319_dailyInf\tz320_dailyInf\tz321_dailyInf\tz322_dailyInf\tz323_dailyInf\tz324_dailyInf\tz325_dailyInf\tz326_dailyInf\tz327_dailyInf\tz328_dailyInf\tz329_dailyInf\t");
    
    fprintf(fw, "cDI1\tcDI2\tcDI3\tcDI4\tcDI5\tcDI6\tcDI7\tcDI8\tcDI9\tcDI10\tcDI11\tcDI12\tcDI13\tcDI14\tcDI15\tcDI16\tcDI17\tcDI18\tcDI19\tcDI20\tcDI21\tcDI22\tcDI23\tcDI24\tcDI25\tcDI26\tcDI27\tcDI28\tcDI29\tcDI30\tcDI31\tcDI32\tcDI33\tcDI34\tcDI35\tcDI36\tcDI37\tcDI38\t");
    
    fprintf(fw, "tDI1\ttDI2\ttDI3\ttDI4\ttDI5\ttDI6\ttDI7\ttDI8\ttDI9\ttDI10\ttDI11\ttDI12\ttDI13\ttDI14\ttDI15\ttDI16\ttDI17\ttDI18\ttDI19\ttDI20\t");
   
    fprintf(fw, "TotalRTestPS\tTotalRTestHS\tTotalRTTPPS\tTotalRTTPHS\tTotalRTTNPS\tTotalRTTNHS\tTotalRTFPPS\tTotalRTFPHS\tTotalRTFNPS\tTotalRTFNHS\t");
    
    for(mcCtr1=0; mcCtr1<27; mcCtr1++){ //correctional centres
        for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
            for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                    fprintf(fw, "%d\t", mcCtr4);
                }
            }
        }
    }

    for(mcCtr1=0; mcCtr1<11; mcCtr1++){ //correctional centres
        for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
            for(mcCtr3=0; mcCtr3<4; mcCtr3++){ //min med max
                for(mcCtr4=0; mcCtr4<19; mcCtr4++){ //min med max
                    fprintf(fw, "%d\t", mcCtr4);
                }
            }
        }
    }

    for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
        for(mcCtr2=0; mcCtr2<4; mcCtr2++){ //min med max
            for(mcCtr3=0; mcCtr3<5; mcCtr3++){ //min med max
                for(mcCtr4=0; mcCtr4<20; mcCtr4++){ //min med max
                    fprintf(fw, "%d\t", mcCtr4);
                }
            }
        }
    }
    
    fprintf(fw, "overcap\t");
    fprintf(fw, "infectedRZ\t");
    
    for(mcCtr1=0; mcCtr1<27; mcCtr1++){ //correctional centres
        //for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
          //  for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
            //    for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                    fprintf(fw, "%d\t", mcCtr1);
              //  }
            //}
        //}
    }

    for(mcCtr1=0; mcCtr1<11; mcCtr1++){ //correctional centres
        //for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
           //for(mcCtr3=0; mcCtr3<4; mcCtr3++){ //min med max
                //for(mcCtr4=0; mcCtr4<19; mcCtr4++){ //min med max
                    fprintf(fw, "%d\t", mcCtr1);
                //}
            //}
        //}
    }

    for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
        //for(mcCtr2=0; mcCtr2<4; mcCtr2++){ //min med max
            //for(mcCtr3=0; mcCtr3<5; mcCtr3++){ //min med max
                //for(mcCtr4=0; mcCtr4<20; mcCtr4++){ //min med max
                    fprintf(fw, "%d\t", mcCtr1);
                //}
            //}
        //}
    }
    
    fprintf(fw, "\n");
    
//    totalTestPS=0, totalTestHS=0, totalTPPS=0, totalTPHS=0, totalTNPS=0, totalTNHS=0, totalFPPS=0, totalFPHS=0, totalFNPS=0, totalFNHS=0;
    
    if(fp==NULL){
        printf("Cannot open file\n");
        exit(0);
    }else{
        //check if correct file format
        
        //if file correct
        fscanf(fp, "%d,", &daysSim);
        printf("Days: %d\n", daysSim);
        
        //option a
        //    for(i=0;i<4;i++){ //read through all 4 populations
        //        fscanf(fp, "%d,", &fi); //read input
        //        printf("No. of people in location %d group %d: %d\n", i, j, fi); //print input
        //    }
        
        if(strcmp(argv[3],option1)==0){ //basic populations USE COLCTR to read
            //fp=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/inputII.txt", "r");
            for(i=0;i<ROWPRIS;i++){ //read through prison locations 16 community 0, prison 1: 1 inmates, 2 prison staff, 3 healthcare staff, 4 essential visitors 5 family visitors. prison 2: 6, 7, 8, 9, 10. prison 3: 11, 12, 13, 14, 15.
                //for(j=0;j<TCOLCTR;j++){ //0 and 1, if 0 IDU+, 1 IDU-
                    fscanf(fp, "%d,", &fi);
                if(i!=0){
                    printf("No. of people in location %d : %d\n", i, fi); //
                    for(inCtr=0;inCtr<fi;inCtr++){
                        newIndiv(&idG, &HCVe, &HCVeAb, 0, &nHead, &nTail, &r, ROWPRIS, pLocArray, i, 99, 1, 99, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap); //seedPop is true
                    }
                }
                //}
            }
        }else if(strcmp(argv[3],option2)==0){ //detailed populations USE COLCTR2 to read
            //fp=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/input.txt", "r");
            
            for(i=0;i<ROWPRIS;i++){ //read through prisons
                for(j=0;j<COLCTR;j++){ //subpopulation array index
                    fscanf(fp, "%d,", &fi);
                    if(i!=0){ //prisons 1-3 (1 (min): 1 (inmate), 2 (prison staff), 3 (health staff), 4 (visitor); 2: 5, 6, 7, 8; 3: 9, 10, 11, 12)
                        printf("No. of people in location %d group %d: %d\n", i, j, fi); //
                        for(l=0; l<fi; l++){
                            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, i, j, 1, 99, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap); //creating individual agents
                        }
                    }else{ //community
                        locArray[i][j]=fi;
                        printf("No. of people in location %d group %d: %d\n", i, j, locArray[i][j]);
                    }
                }
            }

        }
        
        printf( "Initialise populations: done \n");
        fclose(fp);
        
        //else if incorrect
            //exit
    }
    
    int pK;
    //enter prison staff
    for(k=0;k<27;k++){
        for(pK=0; pK<127; pK++){
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 2, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
        }
    }
    for(k=0;k<11;k++){
        for(pK=0; pK<120; pK++){
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 7, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
        }
    }
    for(k=0;k<18;k++){
        for(pK=0; pK<120; pK++){
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 12, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
        }
    }
    
    int sK;
    //enter hs staff //about 3 healthcare staff per prison in Min
    for(k=0;k<27;k++){ //
        for(sK=0; sK<3; sK++){ //3
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 3, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
        }
    }
    for(k=0;k<11;k++){ //about 6 per Med prison
        for(sK=0; sK<6; sK++){ //6
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 8, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
        }
    }
    /*for(k=0;k<18;k++){ //about 7 per Max prison
        for(sK=0; sK<7; sK++){
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 13, 99, 1, k);//new prisoner to be moved from com to a random prison based on **
        }
    }*/
    for(k=0;k<18;k++){ //about 7 per Max prison
        if(k==18){
            for(sK=0; sK<7; sK++){ //7
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 13, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }else{
            for(sK=0; sK<7; sK++){ //7
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 13, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
    }
    //newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 13, 99, 0, 1);

    
    //simulation starts here
    //loop depending on number of days being simulated
    while(currDay!=daysSim+1){
        //printf("DAY %d\n", currDay);
        
        HCVe=0;
        HCVeAb=0;
        totalClrD=0;
        //totalIndivDAA=0;
        //totalIndivOST=0;
        nDaily=0;
        nDailyOST=0;
        //Move non-individual subpopulations from community to prisons
        //i=0; //community
        //function creates new individuals from the community, choose from 0-23
        //lb=0;//Min influx from DCS
        //ub=23;//Max influx from DCS
        //lb=1;//Min Max from influx frquencies Influx0515 file MONTHLY though
        //ub=703;
        //inPop=gsl_ran_flat(r, lb, ub); // uniform distribution between min and max (min= 4.793333e-06; max=6.74e-06) <- decimal percent
        //every 365th day
        
        //inPopRate=gsl_ran_flat(r, 2.1197, 3.236);
        
        /*if(currDay%365==0){
            if(currDay<3286){
                inPopRate = inPopRate + (inPopRate*0.01005349); //(inPopRate*0.025026136); //0.0409 //0.00419//0.0169
                //inPopRateMin= inPopRate + (inPopRateMin*0.0419);//(inPopRate*0.0419); //0.0169
                //inPopRateMax= inPopRate + (inPopRateMax*0.0419);
            }else if(currDay>3285&&currDay<4381){
                inPopRate = inPopRate + (inPopRate*0.20030862);
            }else if(currDay>4380){
                inPopRate = inPopRate + (inPopRate*0.025026136); // averaged 06-16
            }
        }
         */
        
        //if(currDay==1){
            //inPop=1;
            inPop=gsl_ran_flat(r, 56, 62);;
        //}else{
        //    inPop=0;
        //}
        
        //else{
         //   inPopRate = inPopRate + (inPopRate*0.025026136);
         //   inPop=gsl_ran_exponential(r, inPopRate);
        //}

        //inPop=gsl_ran_exponential(r, inPopRate);//inPop=gsl_ran_flat(r, inPopRateMin, inPopRateMax);;//inPop=gsl_ran_exponential(r, inPopRate); //3.0); //2.67); //See influxtrnds.xls 2.67 is the mean influx in days based on data from DCS
        //inIncrease=gsl_ran_flat(r, 0.0, 0.027);
        //inPop=inPop+(inPop*inIncrease); //accomodate increase in new incarcerations
        //printf("rProp %f\n", rProp);
        //inPop=NSWpop*rProp; // get proportion from NSW population
        //printf("inPop %d\n", inPop);
        
        //set daily new cases per prison counter
        for(zInfCtrR=0; zInfCtrR<30; zInfCtrR++){ //correctional centres
            for(zInfCtrC=0; zInfCtrC<7; zInfCtrC++){ //min med max
                zoneMaxPrevCI[zInfCtrR][zInfCtrC]=zoneMaxCIArray[zInfCtrR][zInfCtrC];
                zoneMedPrevCI[zInfCtrR][zInfCtrC]=zoneMedCIArray[zInfCtrR][zInfCtrC];
                zoneMinPrevCI[zInfCtrR][zInfCtrC]=zoneMinCIArray[zInfCtrR][zInfCtrC];
            }
        }
        
        for(cInfCtrR=0; cInfCtrR<38; cInfCtrR++){ //correctional centres
            for(cInfCtrC=0; cInfCtrC<7; cInfCtrC++){ //min med max
                courtPrevCI[cInfCtrR][cInfCtrC]=courtCIArray[cInfCtrR][cInfCtrC];
            }
        }

        for(tInfCtrR=0; tInfCtrR<20; tInfCtrR++){ //correctional centres
            for(tInfCtrC=0; tInfCtrC<7; tInfCtrC++){ //min med max
                truckPrevCI[tInfCtrR][tInfCtrC]=truckCIArray[tInfCtrR][tInfCtrC];
            }
        }
        
        if(currDay==1){ //put in new infections in the first 7 days
            for(k=0;k<inPop-1;k++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 99, 99, 1, 99, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
            newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 99, 99, 0, 99, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap); //set as covid-19 starter
        }else{
            for(k=0;k<inPop;k++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 99, 99, 1, 99, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
        
        //enter essential visitors
        for(k=0;k<27;k++){ //
            for(sK=0; sK<1; sK++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 4, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
        for(k=0;k<11;k++){ //about 6 per Med prison
            for(sK=0; sK<1; sK++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 9, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
        for(k=0;k<18;k++){ //about 7 per Max prison
            for(sK=0; sK<2; sK++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 14, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
        
        //enter family visitor
        for(k=0;k<27;k++){ //
            for(sK=0; sK<13; sK++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 5, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
        for(k=0;k<11;k++){ //about 6 per Med prison
            for(sK=0; sK<17; sK++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 10, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }
        for(k=0;k<18;k++){ //about 7 per Max prison
            for(sK=0; sK<34; sK++){
                newIndiv(&idG, &HCVe, &HCVeAb, currDay, &nHead, &nTail, &r, ROWPRIS, pLocArray, 15, 99, 1, k, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);//new prisoner to be moved from com to a random prison based on **
            }
        }

        dailyInPop=k;

        //unisolate recovered staff every morning
        //staffUnisolate(&nHead, &nTail, currDay);
        //checkQuarantine(&nHead, &nTail, currDay, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);
        //checkIsolate(&nHead, &nTail, currDay, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &overcap);
        
        //rapidTest staff every morning before we begin
        /*if(currDay%2==1){
            totalTestPS=rapidTestPStaff(&r, &nHead, &nTail, currDay, &totalPSTP, &totalPSTN, &totalPSFP, &totalPSFN);
            totalTestHS=rapidTestHStaff(&r, &nHead, &nTail, currDay, &totalHSTP, &totalHSTN, &totalHSFP, &totalHSFN);
  //      }
        */
        
        //PCR test prison community (non-isolated/non-quarantine) & check results and isolate
        //for lockdown and stop movement
        //PCRinmate
        PCRinmate(&r, &nHead, &nTail, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &moveFlag);
        //checkLocks
        unlock(&r, &nHead, &nTail, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, &moveFlag);
        
/*        if(currDay>=1&&currDay<=4380){ //asume OST is captured in current parameters from HITS-p
            if(currDay%365==1){
                nOST=0;
                nOSTDailyCap=0;
            }
        }else if(currDay>4380){ //From 2018 onwards //Here is where you set the OST SCENARIO
            if(currDay%365==1){
                nOST=0;
                nOSTDailyCap=0;
            }
        }
*/
 
        target=nHead; //Target points to the individual at the beginning of the list
        while(target!=NULL){
            //age if 365
            printf("Modifying individual \n");
            
			if(currDay%365==0){
                age(&target);
            }
            
            //Fill up probability array of events
            if(target->indivType==0){ //if inmate
            probEvents[0]=probMoveCom(&r, &target, currDay, &moveFlag);//*rate of release to community
            probEvents[1]=probMoveP1(&r, &target, currDay, &moveFlag);//*rate of movement to Prison 1
            probEvents[2]=0.0;//probMoveP2(&r, &target);//*rate of movement to Prison 2
            probEvents[3]=0.0;//probMoveP3(&r, &target);//*rate of movement to Prison 3 than daily
            probEvents[4]=probOutsideVisit(&r, &target, currDay, &moveFlag); //and probability of getting infected
            probEvents[5]=0.0;//probMoveHospitalF;
            probEvents[6]=0.0;//probMoveHospitalC;
            probEvents[7]=probReturnPrison(&r, &target, nMinHospitalF, nMedHospitalF, nMaxHospitalF);
            probEvents[8]=0.0;//probTest(&r, &target, currDay, 999, 999, 999, nInmatesTested, nPStaffTested, nHStaffTested, nEVisitorsTested, nFVisitorsTested);
            probEvents[9]=0.0;//probTestResult(&r, &target, currDay);
            probEvents[10]=0.0;//probThermalTest;
            probEvents[11]=0.0;//probRapidTest(&r, &target, currDay);
            probEvents[12]=probInfectIndiv(&target, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell);
            probEvents[13]=probProgress(&r, &target, currDay);
            probEvents[14]=probNtrClear(&r, &target, currDay);
            probEvents[15]=probCOVIDdeath(&r, &target, currDay);
            probEvents[16]=0.0;///probNtrDeath(&r, &target);
            probEvents[17]=0.0;///probIsolate
            probEvents[18]=0.0;///probEndIsolate
            probEvents[19]=0.0;///probQuarantine
            probEvents[20]=0.0;///probEndQuarantine
            probEvents[21]=0.0;//probCohorting
            probEvents[22]=0.0;//probEndCohorting
            probEvents[23]=0.0;//probdPPE
            probEvents[24]=0.0;//probEndPPE
            probEvents[25]=0.0;//
            probEvents[26]=0.0;//
            }else if(target->indivType==1){ //if prison staff
                probEvents[0]=0.0;//probMoveCom(&r, &target, currDay);//*rate of release to community
                probEvents[1]=0.0;//probMoveP1(&r, &target);//*rate of movement to Prison 1
                probEvents[2]=0.0;//probMoveP2(&r, &target);//*rate of movement to Prison 2
                probEvents[3]=0.0;//probMoveP3(&r, &target);//*rate of movement to Prison 3 than daily
                probEvents[4]=0.0;//probOutsideVisit; //and probability of getting infected
                probEvents[5]=0.0;//probMoveHospitalF;
                probEvents[6]=0.0;//probMoveHospitalC;
                probEvents[7]=0.0;// probReturnPrison;
                probEvents[8]=0.0;////probTest;
                probEvents[9]=0.0;//probTestResult;
                probEvents[10]=0.0;//probThermalTest;
                probEvents[11]=0.0;//probRapidTest; //rapid test everyone at the beginning of the day
                probEvents[12]=probInfectIndiv(&target, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell);
                probEvents[13]=probProgress(&r, &target, currDay);
                probEvents[14]=probNtrClear(&r, &target, currDay);
                probEvents[15]=probCOVIDdeath(&r, &target, currDay);
                probEvents[16]=0.0;///probNtrDeath(&r, &target);
                probEvents[17]=0.0;//probIsolate(&r, &target, currDay);
                probEvents[18]=0.0;//probEndIsolate(&r, &target, currDay);
                probEvents[19]=0.0;///probQuarantine
                probEvents[20]=0.0;///probEndQuarantine
                probEvents[21]=0.0;//probCohorting
                probEvents[22]=0.0;//probEndCohorting
                probEvents[23]=0.0;//probdPPE
                probEvents[24]=0.0;//probEndPPE
                probEvents[25]=0.0;//
                probEvents[26]=0.0;//
            }else if(target->indivType==2){ //if healthcare staff
                probEvents[0]=0.0;//probMoveCom(&r, &target, currDay);//*rate of release to community
                probEvents[1]=0.0;//probMoveP1(&r, &target);//*rate of movement to Prison 1
                probEvents[2]=0.0;//probMoveP2(&r, &target);//*rate of movement to Prison 2
                probEvents[3]=0.0;//probMoveP3(&r, &target);//*rate of movement to Prison 3 than daily
                probEvents[4]=0.0;//probOutsideVisit; //and probability of getting infected
                probEvents[5]=0.0;//probMoveHospitalF;
                probEvents[6]=0.0;//probMoveHospitalC;
                probEvents[7]=0.0;// probReturnPrison;
                probEvents[8]=0.0;////probTest;
                probEvents[9]=0.0;//probTestResult;
                probEvents[10]=0.0;//probThermalTest;
                probEvents[11]=0.0;//probRapidTest;
                probEvents[12]=probInfectIndiv(&target, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell);
                probEvents[13]=probProgress(&r, &target, currDay);
                probEvents[14]=0.0;//probNtrClear(&r, &target, currDay);
                probEvents[15]=probCOVIDdeath(&r, &target, currDay);
                probEvents[16]=0.0;///probNtrDeath(&r, &target);
                probEvents[17]=0.0;//probIsolate(&r, &target, currDay);
                probEvents[18]=0.0;//probEndIsolate(&r, &target, currDay);
                probEvents[19]=0.0;///probQuarantine
                probEvents[20]=0.0;///probEndQuarantine
                probEvents[21]=0.0;//probCohorting
                probEvents[22]=0.0;//probEndCohorting
                probEvents[23]=0.0;//probdPPE
                probEvents[24]=0.0;//probEndPPE
                probEvents[25]=0.0;//
                probEvents[26]=0.0;//
            }else if(target->indivType==3){ //if essential
                probEvents[0]=1.0;//probMoveCom(&r, &target, currDay);//*rate of release to community
                probEvents[1]=0.0;//
                probEvents[2]=0.0;//
                probEvents[3]=0.0;//
                probEvents[4]=0.0;//
                probEvents[5]=0.0;//
                probEvents[6]=0.0;//
                probEvents[7]=0.0;///
                probEvents[8]=0.0;//
                probEvents[9]=0.0;//
                probEvents[10]=0.0;//probThermalTest;
                probEvents[11]=0.0;//probRapidTest;
                probEvents[12]=probInfectIndiv(&target, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell);
                probEvents[13]=probProgress(&r, &target, currDay);
                probEvents[14]=probNtrClear(&r, &target, currDay);
                probEvents[15]=probCOVIDdeath(&r, &target, currDay);
                probEvents[16]=0.0;//
                probEvents[17]=0.0;//
                probEvents[18]=0.0;///
                probEvents[19]=0.0;///
                probEvents[20]=0.0;//
                probEvents[21]=0.0;//
                probEvents[22]=0.0;//
                probEvents[23]=0.0;//probdPPE
                probEvents[24]=0.0;//probEndPPE
                probEvents[25]=0.0;//
                probEvents[26]=0.0;//
            }else if(target->indivType==4){ //if family visitor
                probEvents[0]=1.0;//probMoveCom(&r, &target, currDay);//*rate of release to community
                probEvents[1]=0.0;//
                probEvents[2]=0.0;//
                probEvents[3]=0.0;//
                probEvents[4]=0.0;//
                probEvents[5]=0.0;//
                probEvents[6]=0.0;//
                probEvents[7]=0.0;///
                probEvents[8]=0.0;//
                probEvents[9]=0.0;//
                probEvents[10]=0.0;//probThermalTest;
                probEvents[11]=0.0;//probRapidTest;
                probEvents[12]=probInfectIndiv(&target, currDay, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone, pLockdownMinArea, pLockdownMedArea, pLockdownMaxArea, pLockdownMinUnit, pLockdownMedUnit, pLockdownMaxUnit, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell);
                probEvents[13]=probProgress(&r, &target, currDay);
                probEvents[14]=probNtrClear(&r, &target, currDay);
                probEvents[15]=probCOVIDdeath(&r, &target, currDay);
                probEvents[16]=0.0;//
                probEvents[17]=0.0;//
                probEvents[18]=0.0;///
                probEvents[19]=0.0;///
                probEvents[20]=0.0;//
                probEvents[21]=0.0;//
                probEvents[22]=0.0;//
                probEvents[23]=0.0;//probdPPE
                probEvents[24]=0.0;//probEndPPE
                probEvents[25]=0.0;//
                probEvents[26]=0.0;//
        }
        
            //1. fill up eventDraw 1 or 2; will a particular event occur or not? 
            eventDraw[0]=0; eventDraw[1]=0; outArray[0]=0; outArray[1]=0;
			for(eCtr=0; eCtr<27; eCtr++){
            	eventDecision[eCtr]=0;
            	eventOrder[eCtr]=eCtr;
			}
			printf("event draw array prepared \n");
            
            //2. Fill-up 1/0s for eventDecision, based on eventDraw
            for(eCtr=0;eCtr<27;eCtr++){
            	eventDraw[0]=1-probEvents[eCtr];
            	eventDraw[1]=probEvents[eCtr];
            	//gsl_ran_multinomial(r, 2, 1, eventDraw, outArray);
            	eventDecision[eCtr]=draw_multinom(&r, 2, eventDraw);
            //	printf("event decision %d: %d\n", eCtr, eventDecision[eCtr]);
            	//iCrawl=0;
            	//iFlag=0;
            	//while(iFlag==0){
            	//	if(outArray[iCrawl]==1){
            	//		eventDecision[eCtr]=iCrawl;
            	//		iFlag=1;
				//	}
				//	iCrawl++;
				//}
			}
			printf("event decision array prepared \n");
            
			//3. Randomise eventOrder array
			gsl_ran_shuffle(r, eventOrder, 27, sizeof(int));
			iCount=0;
			eFlag=0; //0 to keep drawing events; 1 if individual is dead 
			
			//for(eCtr=0;eCtr<27;eCtr++){
			//	printf("eventOrder %d: %d\n", eCtr, eventOrder[eCtr]);
			//}
			
			//for(eCtr=0;eCtr<27;eCtr++){
			//	printf("eventDecision %d: %d\n", eCtr, eventDecision[eCtr]);
			//}
			
            while(iCount<27 && eFlag==0){
            	iFlag=0; //flag for event search
            	iCrawl=0; //index of event to initiate
                //printf("Looking for event order %d\n", iCount);
            	while(iFlag==0){
            		if(eventOrder[iCrawl]==iCount){
            			//printf("to perform event %d: y or n?\n", iCrawl);
            			//if(iCrawl>27){
            			//	exit(1);
						//}
						if(eventDecision[iCrawl]==1){
							printf("PERFORM!\n");
            				 switch(iCrawl){
                                 case 0://move to community
                                     printf("case 0: release\n");
                                     moveLocation(&r, &target, pLocArray, 0, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone); //Community
                                     
                                     if(target->indivType==0){
                                         releasedI++;
                                         if(target->COVID==1){
                                             releasedICov++; //inmates
                                         }
                                     }else if(target->indivType==1){
                                         releasedPS++;
                                         if(target->COVID==1){
                                             releasedPSCov++; //inmates
                                         }
                                     }else if(target->indivType==2){
                                         releasedHS++;
                                         if(target->COVID==1){
                                             releasedHSCov++; //inmates
                                         }
                                     }else if(target->indivType==3){
                                         releasedEV++;
                                         if(target->COVID==1){
                                             releasedEVCov++; //inmates
                                         }
                                     }else if(target->indivType==4){
                                         releasedFV++;
                                         if(target->COVID==1){
                                             releasedFVCov++; //inmates
                                         }
                                     }
                                     
                                     removeIndiv(&target, &nHead, &nTail, pLocArray, pMinCellArray, pMedCellArray, pMaxCellArray); //remove from list and release to community
                                     
                                     deadFlag=1; //totalIndivAgents--; eFlag=1; //to move to the next indiv
                                     E0++; //how many times this event is called
                                     eFlag=1;
                                     
                                     printf("released to community\n");
                                     break;
                                 case 1://move to prison 1
                                     printf("case 1\n");
                                     //if(target->location!=1){
                                     moveLocation(&r, &target, pLocArray, 1, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone); //Prison 1
                                     //nMovingPop[target->truckNumber]++;
                                     if (target->COVID==1){
                                         
                                         success=infectInTransit(&r, &target, &nHead, &nTail, pLocArray, &truckCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, currDay, &moveFlag);
                                         if(success>0){ // success indicates number of new cases from this individual
                                             newCases=newCases+success;
                                         }
                                         
                                     }
                                     //}
                                     E1++;
                                     printf("move prison 1\n");
                                     break;
                                 case 2://move to prison 2
                                     printf("case 2\n");
                                     //if(target->location!=2){
                                     moveLocation(&r, &target, pLocArray, 2, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone); //Prison 2
                                     //nMovingPop[target->truckNumber]++;
                                     if (target->COVID==1){
                                         success=infectInTransit(&r, &target, &nHead, &nTail, pLocArray, &truckCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, currDay, &moveFlag);
                                         if(success>0){ // success indicates number of new cases from this individual
                                             newCases=newCases+success;
                                         }
                                     }
                                     //}
                                     E2++;
                                     printf("move prison 2\n");
                                     break;
                                 case 3://move to prison 3
                                     printf("case 3\n");
                                     //if(target->location!=3){
                                     moveLocation(&r, &target, pLocArray, 3, pMinCellArray, pMedCellArray, pMaxCellArray, pIsolateMinZone, pIsolateMedZone, pIsolateMaxZone, pIsolateMinArea, pIsolateMedArea, pIsolateMaxArea, pIsolateMinUnit, pIsolateMedUnit, pIsolateMaxUnit, pIsolateMinCell, pIsolateMedCell, pIsolateMaxCell, pLockdownMinZone, pLockdownMedZone, pLockdownMaxZone); //Prison 3
                                     //nMovingPop[target->truckNumber]++;
                                     if (target->COVID==1){
                                         success=infectInTransit(&r, &target, &nHead, &nTail, pLocArray, &truckCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, currDay, &moveFlag);
                                         if(success>0){ // success indicates number of new cases from this individual
                                             newCases=newCases+success;
                                         }
                                     }
                                     //}
                                     E3++;
                                     printf("move prison 3\n");
                                     break;
                                 case 4://outside visit
                                     printf("case 4\n");
                                     prevStatus=0;
                                     newStatus=0;
                                     prevStatus=target->COVID;
                                     
                                     outsideVisit(&r, &target, pLocArray, pMinCellArray, pMedCellArray, pMaxCellArray);
                                     nCourtPop[target->courtNumber]++;
                                     //chance of infecting
                                     if (target->COVID==1){
                                         
                                         success=infectOutside(&r, &target, &nHead, &nTail, pLocArray, &courtCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, currDay);
                                         if(success>0){ // success indicates number of new cases from this individual
                                             newCases=newCases+success;
                                         }
                                         
                                     }
                                     E4++;
                                     printf("outside visit 4\n");
                                     break;
                                 case 5://Move to field hospital
                                     printf("case 5\n");
                                     
                                     moveHospital(&target, 1); //Field hospital
                                     
                                     if(target->location==1){
                                         nMinHospitalF++; //total in min field hospitals
                                     }else if(target->location==2){
                                         nMedHospitalF++;
                                     }else if(target->location==3){
                                         nMaxHospitalF++;
                                     }

                                     E5++;
                                     printf("move to field hospital 5\n");
                                     break;
                                 case 6://Move to community hospital
                                     printf("case 6\n");
                                     
                                     moveHospital(&target, 2); //Community hospital
                                     
                                     if(target->location==1){
                                         nMinHospitalC++;
                                     }else if(target->location==2){
                                         nMedHospitalC++;
                                     }else if(target->location==3){
                                         nMaxHospitalC++;
                                     }
                                     
                                     E6++;
                                     printf("move to community hospital 6\n");
                                     break;
                                 case 7://Return to Prison cell
                                     printf("case 7\n");
                                     
                                     if(target->hospitalised==1){ //in field hospital
                                         if(target->location==1){
                                             nMinHospitalF--;
                                         }else if(target->location==2){
                                             nMedHospitalF--;
                                         }else if(target->location==3){
                                             nMaxHospitalF--;
                                         }
                                     }else if(target->hospitalised==2){ //in field hospital
                                         if(target->location==1){
                                             nMinHospitalC--;
                                         }else if(target->location==2){
                                             nMedHospitalC--;
                                         }else if(target->location==3){
                                             nMaxHospitalC--;
                                         }
                                     }
                                     
                                     if(target->court==1){
                                         nCourtPop[target->courtNumber]--;
                                     }
                                     
                                     if(target->moving==1){
                                         //nMovingPop[target->truckNumber]--;
                                     }
                                     
                                     returnPrison(&target, pMinCellArray, pMedCellArray, pMaxCellArray); //return to cell
                                     
                                     E7++;
                                     printf("return to prison cell 7\n");
                                     break;
                                 case 8://RNA Test
                                     printf("case 8\n");
                                     
                                     test(&target, currDay); //test
                                     
                                     if(target->indivType==0){
                                         nInmatesTested++;
                                     }else if(target->indivType==1){
                                         nPStaffTested++;
                                     }else if(target->indivType==2){
                                         nHStaffTested++;
                                     }else if(target->indivType==3){
                                         nEVisitorsTested++;
                                     }else if(target->indivType==4){
                                         nFVisitorsTested++;
                                     }
                                     
                                     E8++;
                                     printf("RNA test 5\n");
                                     break;
                                case 9: //RNA test result
                                     printf("case 9\n");
                                            
                                     testResult(&target, currDay);
                                            
                                     if(target->indivType==0){
                                        nInmatesDetected++;
                                     }else if(target->indivType==1){
                                        nPStaffDetected++;
                                     }else if(target->indivType==2){
                                        nHStaffDetected++;
                                     }else if(target->indivType==3){
                                        nEVisitorsDetected++;
                                     }else if(target->indivType==4){
                                        nFVisitorsDetected++;
                                     }
                                            
                                     E9++;
                                     printf("RNA Test result\n");
                                     break;
                                case 10://Thermal Test
                                     printf("case 10\n");
                                            
                                     thermalTest(&target); //test
                                    
                                     if(target->indivType==0){
                                        nInmatesThermal++;
                                     }else if(target->indivType==1){
                                        nPStaffThermal++;
                                     }else if(target->indivType==2){
                                        nHStaffThermal++;
                                     }else if(target->indivType==3){
                                        nEVisitorsThermal++;
                                     }else if(target->indivType==4){
                                        nFVisitorsThermal++;
                                     }
                                            
                                     E10++;
                                     printf("Thermal test\n");
                                     break;
                                case 11: //rapid test
                                     printf("case 11\n");
                                            
                                     /*rapidResult=rapidTest(&r, &target, currDay);
                                            
                                     if(target->indivType==0){
                                        nInmatesRapidTested++;
                                     }else if(target->indivType==1){
                                        nPStaffRapidTested++;
                                     }else if(target->indivType==2){
                                        nHStaffRapidTested++;
                                     }else if(target->indivType==3){
                                        nEVisitorsRapidTested++;
                                     }else if(target->indivType==4){
                                        nFVisitorsRapidTested++;
                                     }
                                            
                                     E11++;
                                            
                                     printf("Rapid test\n");
                                     if(rapidResult==1){//if target is detected == 1, then isolate
                                         if(target->indivType==1){
                                             target->isolated=1;
                                         }
                                         if(target->indivType==2){
                                             target->isolated=1;
                                         }
                                     
                                         //totalTestPS=0,totalTestHS=0, totalTPPS=0, totalTPHS=0, totalTNPS=0, totalTNHS=0, totalFPPS=0, totalFPHS=0, totalFNPS=0, totalFNHS=0;
                                     }
                                     
                                     if(target->TP==1){
                                         if(target->indivType==1){
                                             totalTPPS++;
                                         }
                                         if(target->indivType==2){
                                             totalTPHS++;
                                         }
                                     }else if(target->FP==1){
                                         if(target->indivType==1){
                                             totalFPPS++;
                                         }
                                         if(target->indivType==2){
                                             totalFPHS++;
                                         }
                                     }else if(target->TN==1){
                                         if(target->indivType==1){
                                             totalTNPS++;
                                         }
                                         if(target->indivType==2){
                                             totalTNHS++;
                                         }
                                     }else if(target->FP==1){
                                         if(target->indivType==1){
                                             totalFPPS++;
                                         }
                                         if(target->indivType==2){
                                             totalFPHS++;
                                         }
                                     }
                                     */
                                     E11++;
                                     break;
                                case 12://infect others 0.05
                                     printf("case 12\n");
                                     success=0;
                                     
                                     if(target->COVID==1){ //add severity condition if needed?
                                         if(target->isolated!=0&&target->indivType==0){
                                             success=infectRedZone(&r, &target, &nHead, &nTail, pLocArray, &zoneMinCIArray, &zoneMedCIArray, &zoneMaxCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, &exp19, &exp20to44, &exp45to54, &exp55to64, &exp65to74, &exp75to84, &exp85, currDay);
                                             if(success>0){ // success indicates number of new cases from this individual
                                                 newCases=newCases+success;
                                                 infectedRedZone=infectedRedZone+success;
                                             }
                                         }else if(target->isolated==0){
                                             success=infect(&r, &target, &nHead, &nTail, pLocArray, &zoneMinCIArray, &zoneMedCIArray, &zoneMaxCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, &exp19, &exp20to44, &exp45to54, &exp55to64, &exp65to74, &exp75to84, &exp85, currDay);
                                             if(success>0){ // success indicates number of new cases from this individual
                                                 newCases=newCases+success;
                                             }
                                         }
                                         
                                         
                                         //for staff, infect inmates in red zone
                                         if(target->indivType==1||target->indivType==2){
                                             //staff to infect inmates in red zones
                                             if(target->isolated==0){
                                                 success=StaffInfectRedZone(&r, &target, &nHead, &nTail, pLocArray, &zoneMinCIArray, &zoneMedCIArray, &zoneMaxCIArray, &newInfectedInmates, &newInfectedPS, &newInfectedHS, &newInfectedEV, &newInfectedFV, &newInfectedInmatesMin, &newInfectedPSMin, &newInfectedHSMin, &newInfectedEVMin, &newInfectedFVMin, &newInfectedInmatesMed, &newInfectedPSMed, &newInfectedHSMed, &newInfectedEVMed, &newInfectedFVMed, &newInfectedInmatesMax, &newInfectedPSMax, &newInfectedHSMax, &newInfectedEVMax, &newInfectedFVMax, &totalAge0, &totalAge1, &totalAge2, &totalAge3, &totalAge4, &totalAge5, &totalAge6, &exp19, &exp20to44, &exp45to54, &exp55to64, &exp65to74, &exp75to84, &exp85, currDay);
                                                 if(success>0){ // success indicates number of new cases from this individual
                                                     newCases=newCases+success;
                                                     infectedRedZone=infectedRedZone+success;
                                                 }
                                             }
                                         }
                                     }
                                     
                                     
                                     
                                     //separate infection functions for inmates, hs, ps, ev, and fv
                                     E12++;
                                     printf("Infect others\n");
                                     break;
                                case 13: //DISEASE PROGRESS
                                     printf("case 13\n");
                                     progress(&r, &target, &exp19, &exp20to44, &exp45to54, &exp55to64, &exp65to74, &exp75to84, &exp85, &pre19, &pre20to44, &pre45to54, &pre55to64, &pre65to74, &pre75to84, &pre85, &asy19, &asy20to44, &asy45to54, &asy55to64, &asy65to74, &asy75to84, &asy85, &mil19, &mil20to44, &mil45to54, &mil55to64, &mil65to74, &mil75to84, &mil85, &mod19, &mod20to44, &mod45to54, &mod55to64, &mod65to74, &mod75to84, &mod85, &sev19, &sev20to44, &sev45to54, &sev55to64, &sev65to74, &sev75to84, &sev85, &cle19, &cle20to44, &cle45to54, &cle55to64, &cle65to74, &cle75to84, &cle85, currDay);
                                     
                                     /*if(target->indivType==0){
                                         if(target->age<=19){ //means infected in a prison
                                             if(target->severity==1){ //asy
                                                 asy19++;
                                             }else if(target->severity==2){ //mil
                                                 mil19++;
                                             }else if(target->severity==3){ //mod
                                                 mod19++;
                                             }else if(target->severity==4){ //sev
                                                 sev19++;
                                             }
                                         }else if(target->age>=20&&target->age<=44){
                                             if(target->severity==1){ //asy
                                                 asy20to44++;
                                             }else if(target->severity==2){ //mil
                                                 mil20to44++;
                                             }else if(target->severity==3){ //mod
                                                 mod20to44++;
                                             }else if(target->severity==4){ //sev
                                                 sev20to44++;
                                             }
                                         }else if(target->age>=45&&target->age<=54){
                                             if(target->severity==1){ //asy
                                                 asy45to54++;
                                             }else if(target->severity==2){ //mil
                                                 mil45to54++;
                                             }else if(target->severity==3){ //mod
                                                 mod45to54++;
                                             }else if(target->severity==4){ //sev
                                                 sev45to54++;
                                             }
                                         }else if(target->age>=55&&target->age<=64){
                                             if(target->severity==1){ //asy
                                                 asy55to64++;
                                             }else if(target->severity==2){ //mil
                                                 mil55to64++;
                                             }else if(target->severity==3){ //mod
                                                 mod55to64++;
                                             }else if(target->severity==4){ //sev
                                                 sev55to64++;
                                             }
                                         }else if(target->age>=65&&target->age<=74){
                                             if(target->severity==1){ //asy
                                                 asy65to74++;
                                             }else if(target->severity==2){ //mil
                                                 mil65to74++;
                                             }else if(target->severity==3){ //mod
                                                 mod65to74++;
                                             }else if(target->severity==4){ //sev
                                                 sev65to74++;
                                             }
                                         }else if(target->age>=75&&target->age<=84){
                                             if(target->severity==1){ //asy
                                                 asy75to84++;
                                             }else if(target->severity==2){ //mil
                                                 mil75to84++;
                                             }else if(target->severity==3){ //mod
                                                 mod75to84++;
                                             }else if(target->severity==4){ //sev
                                                 sev75to84++;
                                             }
                                         }else if(target->age>=85){
                                             if(target->severity==1){ //asy
                                                 asy85++;
                                             }else if(target->severity==2){ //mil
                                                 mil85++;
                                             }else if(target->severity==3){ //mod
                                                 mod85++;
                                             }else if(target->severity==4){ //sev
                                                 sev85++;
                                             }
                                         }
                                     }*/
                                     
                                     E13++;
                                     printf("COVID progress\n");
                                     break;
                                case 14: //NATURAL CLEARANCE
                                     printf("case 14\n");
                                     
                                    //int finalSevExp=0, finalSevPre=0, finalSevAsy=0, finalSevMil=0, finalSevMod=0, finalSevSev=0, finalSevNon=0;
                                     if(target->severity==0){
                                         finalSevExp++;
                                     }else if(target->severity==1){
                                         finalSevPre++;
                                     }else if(target->severity==2){
                                         finalSevAsy++;
                                     }else if(target->severity==3){
                                         finalSevMil++;
                                     }else if(target->severity==4){
                                         finalSevMod++;
                                     }else if(target->severity==5){
                                         finalSevSev++;
                                     }else{
                                         finalSevNon++;
                                     }
                                     
                                     NtrClear(&target, pLocArray);
                                     totalClrN++;
                                     
                                     if(target->indivType==0){
                                         if(target->age<=19){ //means infected in a prison
                                             cle19++;
                                         }else if(target->age>=20&&target->age<=44){
                                             cle20to44++;
                                         }else if(target->age>=45&&target->age<=54){
                                             cle45to54++;
                                         }else if(target->age>=55&&target->age<=64){
                                             cle55to64++;
                                         }else if(target->age>=65&&target->age<=74){
                                             cle65to74++;
                                         }else if(target->age>=75&&target->age<=84){
                                             cle75to84++;
                                         }else if(target->age>=85){
                                             cle85++;
                                         }
                                     }
                                     
                                     E14++;
                                     printf("natural clearance\n");
                                     break;
                                case 15: //COVID death
                                     printf("case 15\n");
                                     
                                     //int finalSevExp=0, finalSevPre=0, finalSevAsy=0, finalSevMil=0, finalSevMod=0, finalSevSev=0, finalSevNon=0;
                                     if(target->severity==0){
                                         finalSevExp++;
                                     }else if(target->severity==1){
                                         finalSevPre++;
                                     }else if(target->severity==2){
                                         finalSevAsy++;
                                     }else if(target->severity==3){
                                         finalSevMil++;
                                     }else if(target->severity==4){
                                         finalSevMod++;
                                     }else if(target->severity==5){
                                         finalSevSev++;
                                     }else{
                                         finalSevNon++;
                                     }
                                     
                                     if(target->indivType==0){
                                         inmateDeathCOVID++;
                                         
                                         if(target->location==1){
                                             if(target->age<=19){
                                                 zoneMinCDArray[target->zone][0]++;
                                                 totalAD0++;
                                             }else if(target->age>=20&&target->age<=44){
                                                 zoneMinCDArray[target->zone][1]++;
                                                 totalAD1++;
                                             }else if(target->age>=45&&target->age<=54){
                                                 zoneMinCDArray[target->zone][2]++;
                                                 totalAD2++;
                                             }else if(target->age>=55&&target->age<=64){
                                                 zoneMinCDArray[target->zone][3]++;
                                                 totalAD3++;
                                             }else if(target->age>=65&&target->age<=74){
                                                 zoneMinCDArray[target->zone][4]++;
                                                 totalAD4++;
                                             }else if(target->age>=75&&target->age<=84){
                                                 zoneMinCDArray[target->zone][5]++;
                                                 totalAD5++;
                                             }else if(target->age>=85){
                                                 zoneMinCDArray[target->zone][6]++;
                                                 totalAD6++;
                                             }
                                         }else if(target->location==2){
                                             if(target->age<=19){
                                                 zoneMedCDArray[target->zone][0]++;
                                                 totalAD0++;
                                             }else if(target->age>=20&&target->age<=44){
                                                 zoneMedCDArray[target->zone][1]++;
                                                 totalAD1++;
                                             }else if(target->age>=45&&target->age<=54){
                                                 zoneMedCDArray[target->zone][2]++;
                                                 totalAD2++;
                                             }else if(target->age>=55&&target->age<=64){
                                                 zoneMedCDArray[target->zone][3]++;
                                                 totalAD3++;
                                             }else if(target->age>=65&&target->age<=74){
                                                 zoneMedCDArray[target->zone][4]++;
                                                 totalAD4++;
                                             }else if(target->age>=75&&target->age<=84){
                                                 zoneMedCDArray[target->zone][5]++;
                                                 totalAD5++;
                                             }else if(target->age>=85){
                                                 zoneMedCDArray[target->zone][6]++;
                                                 totalAD6++;
                                             }
                                         }else if(target->location==3){
                                             if(target->age<=19){
                                                 zoneMaxCDArray[target->zone][0]++;
                                                 totalAD0++;
                                             }else if(target->age>=20&&target->age<=44){
                                                 zoneMaxCDArray[target->zone][1]++;
                                                 totalAD1++;
                                             }else if(target->age>=45&&target->age<=54){
                                                 zoneMaxCDArray[target->zone][2]++;
                                                 totalAD2++;
                                             }else if(target->age>=55&&target->age<=64){
                                                 zoneMaxCDArray[target->zone][3]++;
                                                 totalAD3++;
                                             }else if(target->age>=65&&target->age<=74){
                                                 zoneMaxCDArray[target->zone][4]++;
                                                 totalAD4++;
                                             }else if(target->age>=75&&target->age<=84){
                                                 zoneMaxCDArray[target->zone][5]++;
                                                 totalAD5++;
                                             }else if(target->age>=85){
                                                 zoneMaxCDArray[target->zone][6]++;
                                                 totalAD6++;
                                             }
                                         }
                                         
                                     }else if(target->indivType==1){
                                         PSDeathCOVID++;
                                     }else if(target->indivType==2){
                                         HSDeathCOVID++;
                                     }else if(target->indivType==3){
                                         EVDeathCOVID++;
                                     }else if(target->indivType==4){
                                         FVDeathCOVID++;
                                     }
                                     
                                     removeIndiv(&target, &nHead, &nTail, pLocArray, pMinCellArray, pMedCellArray, pMaxCellArray); //remove from list
                                     totalDeathCOVID++;
                    
                                     dead++; deadFlag=1;
                                     E15++; eFlag=1;
                                     printf("COVID death\n");
                                     break;
                                case 16: //natural death
                                     printf("case 16\n");
                                     removeIndiv(&target, &nHead, &nTail, pLocArray, pMinCellArray, pMedCellArray, pMaxCellArray);
                                     dead++; deadFlag=1;
                                     E16++; eFlag=1;
                                     printf("natural death\n");
                                     break;
                                case 17: //Isolate
                                     printf("case 17\n");
                                     E17++;
                                     printf("Isolate\n");
                                     break;
                                case 18: //EndIsolate
                                     printf("case 18\n");
                                     E18++;
                                     printf("End isolation\n");
                                     break;
                                case 19: //Quarantine
                                     printf("case 19\n");
                                     E19++;
                                     printf("Quarantine\n");
                                     break;
                                case 20: //EndQuarantine
                                     printf("case 20\n");
                                     E20++;
                                     printf("End quarantine\n");
                                     break;
                                case 21: //Cohorting
                                     printf("case 21\n");
                                     E21++;
                                     printf("Cohorting\n");
                                     break;
                                case 22: //EndCohorting
                                     printf("case 22\n");
                                     E22++;
                                     printf("End cohorting\n");
                                     break;
                                case 23: //PPE
                                     printf("case 23\n");
                                     E23++;
                                     printf("PPE\n");
                                     break;
                                case 24: //End PPE
                                     printf("case 23\n");
                                     E23++;
                                     printf("PPE\n");
                                     break;
                    			case 25: //Blank
                                    printf("case 25\n");
                        			E25++;
                        			printf("Blank case\n");
                        			break;
                    			case 26: //Blank
                                     printf("case 26\n");
									 E26++;
                                     printf("Blank case\n");
                        			 break;
                             }
						}
						iCount++;
						iFlag=1;
					}
					iCrawl++;
				}
				//exit(1);
			}
            
            //Move to the next individual (if current individual is still alive)
            if(deadFlag==0){
                //probsAge[0]=;//Age 0-19
                //probsAge[1]=;//Age 20-44
                //probsAge[2]=;//Age 45-54
                //probsAge[3]=;//Age 55-64
                //probsAge[4]=;//Age 65-74
                //probsAge[5]=;//Age 75-84
                //probsAge[6]=;//Age >= 85
                //update population age breakdown
                //int ageArray[ageGroups][2]; //Array containing: [age group][0] with all, [1] with COVID, [2] with detected COVID
                tAge=0; tCov=0; tDetCov=0;
                tAge=target->age;
                tCov=target->COVID;
                tDetCov=target->detected;
                
                if(tAge<=19){
                    if(target->indivType==0){//inmate
                        iAgeArray[0][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[0][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[0][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[0][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[0][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[0][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[0][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[0][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[0][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[0][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[0][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[0][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[0][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[0][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[0][2]++;
                            }
                        }
                    }
                }else if(tAge>=20&&tAge<=44){
                    if(target->indivType==0){//inmate
                        iAgeArray[1][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[1][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[1][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[1][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[1][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[1][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[1][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[1][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[1][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[1][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[1][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[1][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[1][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[1][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[1][2]++;
                            }
                        }
                    }
                }else if(tAge>=45&&tAge<=54){
                    if(target->indivType==0){//inmate
                        iAgeArray[2][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[2][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[2][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[2][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[2][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[2][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[2][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[2][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[2][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[2][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[2][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[2][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[2][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[2][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[2][2]++;
                            }
                        }
                    }
                }else if(tAge>=55&&tAge<=64){
                    if(target->indivType==0){//inmate
                        iAgeArray[3][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[3][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[3][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[3][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[3][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[3][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[3][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[3][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[3][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[3][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[3][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[3][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[3][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[3][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[3][2]++;
                            }
                        }
                    }
                }else if(tAge>=65&&tAge<=74){
                    if(target->indivType==0){//inmate
                        iAgeArray[4][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[4][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[4][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[4][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[4][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[4][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[4][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[4][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[4][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[4][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[4][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[4][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[4][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[4][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[4][2]++;
                            }
                        }
                    }
                }else if(tAge>=75&&tAge<=84){
                    if(target->indivType==0){//inmate
                        iAgeArray[5][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[5][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[5][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[5][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[5][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[5][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[5][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[5][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[5][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[5][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[5][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[5][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[5][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[5][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[5][2]++;
                            }
                        }
                    }
                }else if(tAge>=85){
                    if(target->indivType==0){//inmate
                        iAgeArray[6][0]++;
                    }else if(target->indivType==1){//ps staff
                        psAgeArray[6][0]++;
                    }else if(target->indivType==2){//hs staff
                        hsAgeArray[6][0]++;
                    }else if(target->indivType==3){//essential visitors
                        evAgeArray[6][0]++;
                    }else if(target->indivType==4){//family visitor
                        fvAgeArray[6][0]++;
                    }
                    
                    if(tCov==1){ //
                        if(target->indivType==0){//inmate
                            iAgeArray[6][1]++;
                        }else if(target->indivType==1){//ps staff
                            psAgeArray[6][1]++;
                        }else if(target->indivType==2){//hs staff
                            hsAgeArray[6][1]++;
                        }else if(target->indivType==3){//essential visitors
                            evAgeArray[6][1]++;
                        }else if(target->indivType==4){//family visitor
                            fvAgeArray[6][1]++;
                        }
                        
                        if(tDetCov==1){
                            if(target->indivType==0){//inmate
                                iAgeArray[6][2]++;
                            }else if(target->indivType==1){//ps staff
                                psAgeArray[6][2]++;
                            }else if(target->indivType==2){//hs staff
                                hsAgeArray[6][2]++;
                            }else if(target->indivType==3){//essential visitors
                                evAgeArray[6][2]++;
                            }else if(target->indivType==4){//family visitor
                                fvAgeArray[6][2]++;
                            }
                        }
                    }
                }
                
                target=target->nextIndiv;
            }
            deadFlag=0;

        }
        
        //while loop for events ends here
        //print counts and individuals
        fprintf(fw, "%d\t", currDay);
        for(i=0; i<ROWPRIS; i++) //15
        {
            for(j=0; j<COLCTR; j++){ //12
                //printf("LOCATION %d GROUP %d >>> %d \n", i, j, locArray[i][j]);
                fprintf(fw, "%d\t", locArray[i][j]);
                if(i==1&&i==6&&i==11){
                    totalPrisonPop=totalPrisonPop+locArray[i][j];
                }
            }
        }
        
        //count totalPrisPop
        
        
        //count prevalence in prison
        if(nTail!=NULL){
            prisonHCV=countHCVpris(&nHead, &nTail);
        }
        
        //count prevalence in prison obtained in community
        if(nTail!=NULL){
            communityHCV=countHCVCom(&nHead, &nTail);
        }
        
        if(nTail!=NULL){
            prisonOpd=countOpd(&nHead, &nTail);
        }

        if(nTail!=NULL){
            prisonOpdNotOST=countOpdNotOST(&nHead, &nTail);
        }
        
        if(nTail!=NULL){
            prisonEverIDU=countEverIDU(&nHead, &nTail);
        }
        
        if(nTail!=NULL){
            hcvAntibody=countHCVantibody(&nHead, &nTail);
        }
        
        if(nTail!=NULL){
            hcvRNA=countHCVRNA(&nHead, &nTail);
        }

        if(nTail!=NULL){
            aveLengthStay=countAveStay(&nHead, &nTail, currDay);
        }
        
        if(nTail!=NULL){
            countZonePop(&nHead, &zoneArray, &zoneCovidArray);
        }
        
        if(nTail!=NULL){
            countSevereInmates(&nHead, &nTail, &iAgeSeverityArray);
        }
        
        if(nTail!=NULL){
            countMovingPop(&nHead, &nTail, &nMovingPop);
        }
        
        //totalTestPS=0, totalTestHS=0, totalTPPS=0, totalTPHS=0, totalTNPS=0, totalTNHS=0, totalFPPS=0, totalFPHS=0, totalFNPS=0, totalFNHS=0;
        
        //Calculate new daily infections
        for(zInfCtrR=0; zInfCtrR<30; zInfCtrR++){ //correctional centres
            for(zInfCtrC=0; zInfCtrC<7; zInfCtrC++){ //min med max
                zoneMaxDIArray[zInfCtrR][zInfCtrC]=zoneMaxCIArray[zInfCtrR][zInfCtrC]-zoneMaxPrevCI[zInfCtrR][zInfCtrC];
                zoneMedDIArray[zInfCtrR][zInfCtrC]=zoneMedCIArray[zInfCtrR][zInfCtrC]-zoneMedPrevCI[zInfCtrR][zInfCtrC];
                zoneMinDIArray[zInfCtrR][zInfCtrC]=zoneMinCIArray[zInfCtrR][zInfCtrC]-zoneMinPrevCI[zInfCtrR][zInfCtrC];
            }
        }
        
        for(cInfCtrR=0; cInfCtrR<38; cInfCtrR++){ //correctional centres
            for(cInfCtrC=0; cInfCtrC<7; cInfCtrC++){ //min med max
                courtDIArray[cInfCtrR][cInfCtrC]=courtCIArray[cInfCtrR][cInfCtrC]-courtPrevCI[cInfCtrR][cInfCtrC];
            }
        }


        for(tInfCtrR=0; tInfCtrR<20; tInfCtrR++){ //correctional centres
            for(tInfCtrC=0; tInfCtrC<7; tInfCtrC++){ //min med max
                truckDIArray[tInfCtrR][tInfCtrC]=truckCIArray[tInfCtrR][tInfCtrC]-truckPrevCI[tInfCtrR][tInfCtrC];
            }
        }
        
        //count Rs
        /*if(nTail!=NULL){
            R0=count0(&nHead, &nTail);
            R1=count1(&nHead, &nTail);
            R2=count2(&nHead, &nTail);
            R3=count3(&nHead, &nTail);
            R4=count4(&nHead, &nTail);
            R5=count5(&nHead, &nTail);
            R6=count6(&nHead, &nTail);
            R7=count7(&nHead, &nTail);
            R8=count8(&nHead, &nTail);
            R9=count9(&nHead, &nTail);
            R10=count10(&nHead, &nTail);
            R11=count11(&nHead, &nTail);
            R12=count12(&nHead, &nTail);
            nReinfected=countR(&nHead, &nTail);
        }
         */
        
        //fprintf(fw, "NEW CASES\t
        //InfectedOutVisits\t
        //releasedI\t
        //releasedICov\t
        //releasedPS\t
        //releasedPSCov\t
        //releasedHS\t
        //releasedHSCov\t
        //releasedEV\t
        //releasedEVCov\t
        //releasedFV\t
        //releasedFVCov\t");
        
        fprintf(fw, "%d\t", newCases); //print new infections per day
        fprintf(fw, "%d\t", newInfectedInmates); //print new infections per day
        fprintf(fw, "%d\t", newInfectedPS); //print new infections per day
        fprintf(fw, "%d\t", newInfectedHS); //print new infections per day
        fprintf(fw, "%d\t", newInfectedEV); //print new infections per day
        fprintf(fw, "%d\t", newInfectedFV); //print new infections per day
        fprintf(fw, "%d\t", newInfectedInmatesMin); //print new infections per day
        fprintf(fw, "%d\t", newInfectedPSMin); //print new infections per day
        fprintf(fw, "%d\t", newInfectedHSMin); //print new infections per day
        fprintf(fw, "%d\t", newInfectedEVMin); //print new infections per day
        fprintf(fw, "%d\t", newInfectedFVMin); //print new infections per day
        fprintf(fw, "%d\t", newInfectedInmatesMed); //print new infections per day
        fprintf(fw, "%d\t", newInfectedPSMed); //print new infections per day
        fprintf(fw, "%d\t", newInfectedHSMed); //print new infections per day
        fprintf(fw, "%d\t", newInfectedEVMed); //print new infections per day
        fprintf(fw, "%d\t", newInfectedFVMed); //print new infections per day
        fprintf(fw, "%d\t", newInfectedInmatesMax); //print new infections per day
        fprintf(fw, "%d\t", newInfectedPSMax); //print new infections per day
        fprintf(fw, "%d\t", newInfectedHSMax); //print new infections per day
        fprintf(fw, "%d\t", newInfectedEVMax); //print new infections per day
        fprintf(fw, "%d\t", newInfectedFVMax); //print new infections per day
        fprintf(fw, "%d\t", infectedOutside); //numbers infected during outside visit
        fprintf(fw, "%d\t", releasedI);
        fprintf(fw, "%d\t", releasedICov);
        fprintf(fw, "%d\t", releasedPS);
        fprintf(fw, "%d\t", releasedPSCov);
        fprintf(fw, "%d\t", releasedHS);
        fprintf(fw, "%d\t", releasedHSCov);
        fprintf(fw, "%d\t", releasedEV);
        fprintf(fw, "%d\t", releasedEVCov);
        fprintf(fw, "%d\t", releasedFV);
        fprintf(fw, "%d\t", releasedFVCov);
        
        fprintf(fw, "%d\t", nInmatesTested);
        fprintf(fw, "%d\t", nPStaffTested);
        fprintf(fw, "%d\t", nHStaffTested);
        fprintf(fw, "%d\t", nEVisitorsTested);
        fprintf(fw, "%d\t", nFVisitorsTested);
        
        fprintf(fw, "%d\t", nInmatesThermal);
        fprintf(fw, "%d\t", nPStaffThermal);
        fprintf(fw, "%d\t", nHStaffThermal);
        fprintf(fw, "%d\t", nEVisitorsThermal);
        fprintf(fw, "%d\t", nFVisitorsThermal);
        
        fprintf(fw, "%d\t", nInmatesDetected);
        fprintf(fw, "%d\t", nPStaffDetected);
        fprintf(fw, "%d\t", nHStaffDetected);
        fprintf(fw, "%d\t", nEVisitorsDetected);
        fprintf(fw, "%d\t", nFVisitorsDetected);
        
        fprintf(fw, "%d\t", nInmatesRapidTested);
        fprintf(fw, "%d\t", nPStaffRapidTested);
        fprintf(fw, "%d\t", nHStaffRapidTested);
        fprintf(fw, "%d\t", nEVisitorsRapidTested);
        fprintf(fw, "%d\t", nFVisitorsRapidTested);

        fprintf(fw, "%d\t", nMinHospitalF);
        fprintf(fw, "%d\t", nMedHospitalF);
        fprintf(fw, "%d\t", nMaxHospitalF);
        fprintf(fw, "%d\t", nMinHospitalC);
        fprintf(fw, "%d\t", nMedHospitalC);
        fprintf(fw, "%d\t", nMaxHospitalC);
        
        fprintf(fw, "%d\t", iAgeArray[0][0]);
        fprintf(fw, "%d\t", iAgeArray[1][0]);
        fprintf(fw, "%d\t", iAgeArray[2][0]);
        fprintf(fw, "%d\t", iAgeArray[3][0]);
        fprintf(fw, "%d\t", iAgeArray[4][0]);
        fprintf(fw, "%d\t", iAgeArray[5][0]);
        fprintf(fw, "%d\t", iAgeArray[6][0]);
        fprintf(fw, "%d\t", iAgeArray[0][1]);
        fprintf(fw, "%d\t", iAgeArray[1][1]);
        fprintf(fw, "%d\t", iAgeArray[2][1]);
        fprintf(fw, "%d\t", iAgeArray[3][1]);
        fprintf(fw, "%d\t", iAgeArray[4][1]);
        fprintf(fw, "%d\t", iAgeArray[5][1]);
        fprintf(fw, "%d\t", iAgeArray[6][1]);
        fprintf(fw, "%d\t", iAgeArray[0][2]);
        fprintf(fw, "%d\t", iAgeArray[1][2]);
        fprintf(fw, "%d\t", iAgeArray[2][2]);
        fprintf(fw, "%d\t", iAgeArray[3][2]);
        fprintf(fw, "%d\t", iAgeArray[4][2]);
        fprintf(fw, "%d\t", iAgeArray[5][2]);
        fprintf(fw, "%d\t", iAgeArray[6][2]);
        
        fprintf(fw, "%d\t", psAgeArray[0][0]);
        fprintf(fw, "%d\t", psAgeArray[1][0]);
        fprintf(fw, "%d\t", psAgeArray[2][0]);
        fprintf(fw, "%d\t", psAgeArray[3][0]);
        fprintf(fw, "%d\t", psAgeArray[4][0]);
        fprintf(fw, "%d\t", psAgeArray[5][0]);
        fprintf(fw, "%d\t", psAgeArray[6][0]);
        fprintf(fw, "%d\t", psAgeArray[0][1]);
        fprintf(fw, "%d\t", psAgeArray[1][1]);
        fprintf(fw, "%d\t", psAgeArray[2][1]);
        fprintf(fw, "%d\t", psAgeArray[3][1]);
        fprintf(fw, "%d\t", psAgeArray[4][1]);
        fprintf(fw, "%d\t", psAgeArray[5][1]);
        fprintf(fw, "%d\t", psAgeArray[6][1]);
        fprintf(fw, "%d\t", psAgeArray[0][2]);
        fprintf(fw, "%d\t", psAgeArray[1][2]);
        fprintf(fw, "%d\t", psAgeArray[2][2]);
        fprintf(fw, "%d\t", psAgeArray[3][2]);
        fprintf(fw, "%d\t", psAgeArray[4][2]);
        fprintf(fw, "%d\t", psAgeArray[5][2]);
        fprintf(fw, "%d\t", psAgeArray[6][2]);
        
        fprintf(fw, "%d\t", hsAgeArray[0][0]);
        fprintf(fw, "%d\t", hsAgeArray[1][0]);
        fprintf(fw, "%d\t", hsAgeArray[2][0]);
        fprintf(fw, "%d\t", hsAgeArray[3][0]);
        fprintf(fw, "%d\t", hsAgeArray[4][0]);
        fprintf(fw, "%d\t", hsAgeArray[5][0]);
        fprintf(fw, "%d\t", hsAgeArray[6][0]);
        fprintf(fw, "%d\t", hsAgeArray[0][1]);
        fprintf(fw, "%d\t", hsAgeArray[1][1]);
        fprintf(fw, "%d\t", hsAgeArray[2][1]);
        fprintf(fw, "%d\t", hsAgeArray[3][1]);
        fprintf(fw, "%d\t", hsAgeArray[4][1]);
        fprintf(fw, "%d\t", hsAgeArray[5][1]);
        fprintf(fw, "%d\t", hsAgeArray[6][1]);
        fprintf(fw, "%d\t", hsAgeArray[0][2]);
        fprintf(fw, "%d\t", hsAgeArray[1][2]);
        fprintf(fw, "%d\t", hsAgeArray[2][2]);
        fprintf(fw, "%d\t", hsAgeArray[3][2]);
        fprintf(fw, "%d\t", hsAgeArray[4][2]);
        fprintf(fw, "%d\t", hsAgeArray[5][2]);
        fprintf(fw, "%d\t", hsAgeArray[6][2]);

        fprintf(fw, "%d\t", evAgeArray[0][0]);
        fprintf(fw, "%d\t", evAgeArray[1][0]);
        fprintf(fw, "%d\t", evAgeArray[2][0]);
        fprintf(fw, "%d\t", evAgeArray[3][0]);
        fprintf(fw, "%d\t", evAgeArray[4][0]);
        fprintf(fw, "%d\t", evAgeArray[5][0]);
        fprintf(fw, "%d\t", evAgeArray[6][0]);
        fprintf(fw, "%d\t", evAgeArray[0][1]);
        fprintf(fw, "%d\t", evAgeArray[1][1]);
        fprintf(fw, "%d\t", evAgeArray[2][1]);
        fprintf(fw, "%d\t", evAgeArray[3][1]);
        fprintf(fw, "%d\t", evAgeArray[4][1]);
        fprintf(fw, "%d\t", evAgeArray[5][1]);
        fprintf(fw, "%d\t", evAgeArray[6][1]);
        fprintf(fw, "%d\t", evAgeArray[0][2]);
        fprintf(fw, "%d\t", evAgeArray[1][2]);
        fprintf(fw, "%d\t", evAgeArray[2][2]);
        fprintf(fw, "%d\t", evAgeArray[3][2]);
        fprintf(fw, "%d\t", evAgeArray[4][2]);
        fprintf(fw, "%d\t", evAgeArray[5][2]);
        fprintf(fw, "%d\t", evAgeArray[6][2]);

        fprintf(fw, "%d\t", fvAgeArray[0][0]);
        fprintf(fw, "%d\t", fvAgeArray[1][0]);
        fprintf(fw, "%d\t", fvAgeArray[2][0]);
        fprintf(fw, "%d\t", fvAgeArray[3][0]);
        fprintf(fw, "%d\t", fvAgeArray[4][0]);
        fprintf(fw, "%d\t", fvAgeArray[5][0]);
        fprintf(fw, "%d\t", fvAgeArray[6][0]);
        fprintf(fw, "%d\t", fvAgeArray[0][1]);
        fprintf(fw, "%d\t", fvAgeArray[1][1]);
        fprintf(fw, "%d\t", fvAgeArray[2][1]);
        fprintf(fw, "%d\t", fvAgeArray[3][1]);
        fprintf(fw, "%d\t", fvAgeArray[4][1]);
        fprintf(fw, "%d\t", fvAgeArray[5][1]);
        fprintf(fw, "%d\t", fvAgeArray[6][1]);
        fprintf(fw, "%d\t", fvAgeArray[0][2]);
        fprintf(fw, "%d\t", fvAgeArray[1][2]);
        fprintf(fw, "%d\t", fvAgeArray[2][2]);
        fprintf(fw, "%d\t", fvAgeArray[3][2]);
        fprintf(fw, "%d\t", fvAgeArray[4][2]);
        fprintf(fw, "%d\t", fvAgeArray[5][2]);
        fprintf(fw, "%d\t", fvAgeArray[6][2]);
        
        fprintf(fw, "%d\t", communityHCV); //prevalence of those infected in community
        fprintf(fw, "%d\t", released); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCVAb); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCVRNA); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCVCom); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCVPris); //print HCV-infected individuals released per day
            
        fprintf(fw, "%d\t", totalClrN); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", dead); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", deadHCV); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", deadHCVCom); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", deadHCVPris); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", totalPrisonPop);
        fprintf(fw, "%d\t", E0);
        fprintf(fw, "%d\t", E1);
        fprintf(fw, "%d\t", E2);
        fprintf(fw, "%d\t", E3);
        fprintf(fw, "%d\t", E4);
        fprintf(fw, "%d\t", E5);
        fprintf(fw, "%d\t", E6);
        fprintf(fw, "%d\t", E7);
        fprintf(fw, "%d\t", E8);
        fprintf(fw, "%d\t", E9);
        fprintf(fw, "%d\t", E10);
        fprintf(fw, "%d\t", E11);
        fprintf(fw, "%d\t", E12);
        fprintf(fw, "%d\t", E13);
        fprintf(fw, "%d\t", E14);
        fprintf(fw, "%d\t", E15);
        fprintf(fw, "%d\t", E16);
        fprintf(fw, "%d\t", E17);
        fprintf(fw, "%d\t", E18);
        fprintf(fw, "%d\t", E19);
        fprintf(fw, "%d\t", E20);
        fprintf(fw, "%d\t", E21);
        fprintf(fw, "%d\t", E22);
        fprintf(fw, "%d\t", E23);
        fprintf(fw, "%d\t", E24); //How many got DAA
        fprintf(fw, "%d\t", E25);
        fprintf(fw, "%d\t", E26); //total DAA clearance
        fprintf(fw, "%d\t", R0);
        fprintf(fw, "%d\t", R1);
        fprintf(fw, "%d\t", R2);
        fprintf(fw, "%d\t", R3);
        fprintf(fw, "%d\t", R4);
        fprintf(fw, "%d\t", R5);
        fprintf(fw, "%d\t", R6);
        fprintf(fw, "%d\t", R7);
        fprintf(fw, "%d\t", R8);
        fprintf(fw, "%d\t", R9);
        fprintf(fw, "%d\t", R10);
        fprintf(fw, "%d\t", R11);
        fprintf(fw, "%d\t", R12);

        fprintf(fw, "%d\t", prisonEverIDU); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", prisonOpd); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", prisonOpdNotOST); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", totalIndivOST); //total under OST in prison
        fprintf(fw, "%d\t", totalIndivDAA); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", totalClrD); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", nReinfected); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", nDAA); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nOST); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", HCVe); //number of those infected with HCV RNA upon entry
        fprintf(fw, "%d\t", HCVeAb); //number of those infected with HCV Ab upon entry
        fprintf(fw, "%d\t", dailyInPop); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", hcvAntibody); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", hcvRNA); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", aveLengthStay); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneArray[29][2]); //number of those infected with HCV upon entry
 
        /*for(zInfCtrR=0; zInfCtrR<30; zInfCtrR++){ //correctional centres
         totalAge0 = totalAge0+zoneMinCIArray[zInfCtrR][0]+zoneMedCIArray[zInfCtrR][0]+zoneMaxCIArray[zInfCtrR][0];
         totalAge1 = totalAge1+zoneMinCIArray[zInfCtrR][1]+zoneMedCIArray[zInfCtrR][1]+zoneMaxCIArray[zInfCtrR][1];
         totalAge2 = totalAge2+zoneMinCIArray[zInfCtrR][2]+zoneMedCIArray[zInfCtrR][2]+zoneMaxCIArray[zInfCtrR][2];
         totalAge3 = totalAge3+zoneMinCIArray[zInfCtrR][3]+zoneMedCIArray[zInfCtrR][3]+zoneMaxCIArray[zInfCtrR][3];
         totalAge4 = totalAge4+zoneMinCIArray[zInfCtrR][4]+zoneMedCIArray[zInfCtrR][4]+zoneMaxCIArray[zInfCtrR][4];
         totalAge5 = totalAge5+zoneMinCIArray[zInfCtrR][5]+zoneMedCIArray[zInfCtrR][5]+zoneMaxCIArray[zInfCtrR][5];
         totalAge6 = totalAge6+zoneMinCIArray[zInfCtrR][6]+zoneMedCIArray[zInfCtrR][6]+zoneMaxCIArray[zInfCtrR][6];
         totalAD0 = totalAD0+zoneMinCDArray[zInfCtrR][0]+zoneMedCDArray[zInfCtrR][0]+zoneMaxCDArray[zInfCtrR][0];
         totalAD1 = totalAD1+zoneMinCDArray[zInfCtrR][1]+zoneMedCDArray[zInfCtrR][1]+zoneMaxCDArray[zInfCtrR][1];
         totalAD2 = totalAD2+zoneMinCDArray[zInfCtrR][2]+zoneMedCDArray[zInfCtrR][2]+zoneMaxCDArray[zInfCtrR][2];
         totalAD3 = totalAD3+zoneMinCDArray[zInfCtrR][3]+zoneMedCDArray[zInfCtrR][3]+zoneMaxCDArray[zInfCtrR][3];
         totalAD4 = totalAD4+zoneMinCDArray[zInfCtrR][4]+zoneMedCDArray[zInfCtrR][4]+zoneMaxCDArray[zInfCtrR][4];
         totalAD5 = totalAD5+zoneMinCDArray[zInfCtrR][5]+zoneMedCDArray[zInfCtrR][5]+zoneMaxCDArray[zInfCtrR][5];
         totalAD6 = totalAD6+zoneMinCDArray[zInfCtrR][6]+zoneMedCDArray[zInfCtrR][6]+zoneMaxCDArray[zInfCtrR][6];
         }*/
        for(zaCtrR=0; zaCtrR<30; zaCtrR++){ //correctional centres
                zoneCovidArray[zaCtrR][0]+=zoneMinCIArray[zaCtrR][0]+zoneMinCIArray[zaCtrR][1]+zoneMinCIArray[zaCtrR][2]+zoneMinCIArray[zaCtrR][3]+zoneMinCIArray[zaCtrR][4]+zoneMinCIArray[zaCtrR][5]+zoneMinCIArray[zaCtrR][6];
                zoneCovidArray[zaCtrR][1]+=zoneMedCIArray[zaCtrR][0]+zoneMedCIArray[zaCtrR][1]+zoneMedCIArray[zaCtrR][2]+zoneMedCIArray[zaCtrR][3]+zoneMedCIArray[zaCtrR][4]+zoneMedCIArray[zaCtrR][5]+zoneMedCIArray[zaCtrR][6];
                zoneCovidArray[zaCtrR][2]+=zoneMaxCIArray[zaCtrR][0]+zoneMaxCIArray[zaCtrR][1]+zoneMaxCIArray[zaCtrR][2]+zoneMaxCIArray[zaCtrR][3]+zoneMaxCIArray[zaCtrR][4]+zoneMaxCIArray[zaCtrR][5]+zoneMaxCIArray[zaCtrR][6];
        }
        
        fprintf(fw, "%d\t", zoneCovidArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneCovidArray[29][2]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", iAgeSeverityArray[0][0]); //age group 0
        fprintf(fw, "%d\t", iAgeSeverityArray[0][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[0][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[0][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[0][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[0][5]);
        fprintf(fw, "%d\t", iAgeSeverityArray[1][0]); //age group 1
        fprintf(fw, "%d\t", iAgeSeverityArray[1][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[1][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[1][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[1][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[1][5]);
        fprintf(fw, "%d\t", iAgeSeverityArray[2][0]); //age group 1
        fprintf(fw, "%d\t", iAgeSeverityArray[2][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[2][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[2][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[2][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[2][5]);
        fprintf(fw, "%d\t", iAgeSeverityArray[3][0]); //age group 0
        fprintf(fw, "%d\t", iAgeSeverityArray[3][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[3][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[3][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[3][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[3][5]);
        fprintf(fw, "%d\t", iAgeSeverityArray[4][0]); //age group 1
        fprintf(fw, "%d\t", iAgeSeverityArray[4][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[4][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[4][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[4][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[4][5]);
        fprintf(fw, "%d\t", iAgeSeverityArray[5][0]); //age group 1
        fprintf(fw, "%d\t", iAgeSeverityArray[5][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[5][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[5][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[5][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[5][5]);
        fprintf(fw, "%d\t", iAgeSeverityArray[6][0]); //age group 1
        fprintf(fw, "%d\t", iAgeSeverityArray[6][1]);
        fprintf(fw, "%d\t", iAgeSeverityArray[6][2]);
        fprintf(fw, "%d\t", iAgeSeverityArray[6][3]);
        fprintf(fw, "%d\t", iAgeSeverityArray[6][4]);
        fprintf(fw, "%d\t", iAgeSeverityArray[6][5]);
        
        fprintf(fw, "%d\t", zoneMinCIArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[0][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[0][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[0][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[0][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][6]); //number of those infected with HCV upon entry
    
        fprintf(fw, "%d\t", zoneMinCIArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCIArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][6]); //number of those infected with HCV upon entry
 
        fprintf(fw, "%d\t", zoneMedCIArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[0][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[0][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[0][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[0][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCIArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[0][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[0][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[0][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[0][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][6]); //number of those infected with HCV upon entry

        fprintf(fw, "%d\t", zoneMinCDArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[0][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[0][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[0][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[0][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[1][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[1][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[1][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[1][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[2][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[2][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[2][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[2][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[3][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[3][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[3][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[3][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[4][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[4][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[4][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[4][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[5][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[5][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[5][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[5][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[6][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[6][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[6][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[6][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[7][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[7][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[7][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[7][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[8][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[8][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[8][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[8][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[9][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[9][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[9][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[9][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[10][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[10][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[10][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[10][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[11][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[11][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[11][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[11][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[12][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[12][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[12][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[12][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[13][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[13][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[13][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[13][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[14][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[14][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[14][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[14][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[15][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[15][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[15][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[15][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[16][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[16][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[16][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[16][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[17][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[17][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[17][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[17][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[18][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[18][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[18][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[18][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[19][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[19][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[19][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[20][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[20][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[20][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[20][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[21][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[21][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[21][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[21][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[22][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[22][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[22][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[22][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[23][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[23][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[23][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[23][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[24][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[24][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[24][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[24][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[25][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[25][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[25][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[25][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[26][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[26][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[26][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[26][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[27][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[27][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[27][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[27][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[28][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[28][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[28][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[28][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMinCDArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[29][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[29][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[29][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[29][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCDArray[29][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[0][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[0][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[0][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[0][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[1][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[1][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[1][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[1][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[2][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[2][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[2][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[2][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[3][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[3][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[3][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[3][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[4][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[4][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[4][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[4][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[5][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[5][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[5][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[5][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[6][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[6][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[6][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[6][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[7][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[7][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[7][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[7][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[8][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[8][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[8][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[8][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[9][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[9][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[9][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[9][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[10][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[10][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[10][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[10][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[11][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[11][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[11][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[11][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[12][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[12][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[12][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[12][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[13][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[13][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[13][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[13][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[14][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[14][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[14][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[14][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[15][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[15][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[15][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[15][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[16][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[16][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[16][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[16][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[17][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[17][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[17][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[17][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[18][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[18][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[18][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[18][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[19][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[19][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[19][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[20][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[20][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[20][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[20][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[21][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[21][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[21][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[21][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[22][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[22][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[22][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[22][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[23][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[23][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[23][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[23][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[24][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[24][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[24][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[24][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[25][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[25][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[25][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[25][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[26][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[26][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[26][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[26][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[27][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[27][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[27][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[27][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[28][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[28][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[28][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[28][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMedCDArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[29][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[29][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[29][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[29][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCDArray[29][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[0][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[0][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[0][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[0][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[0][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[0][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[0][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[1][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[1][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[1][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[1][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[1][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[1][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[1][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[2][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[2][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[2][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[2][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[2][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[2][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[2][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[3][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[3][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[3][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[3][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[3][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[3][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[3][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[4][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[4][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[4][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[4][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[4][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[4][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[4][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[5][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[5][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[5][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[5][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[5][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[5][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[5][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[6][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[6][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[6][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[6][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[6][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[6][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[6][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[7][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[7][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[7][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[7][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[7][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[7][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[7][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[8][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[8][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[8][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[8][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[8][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[8][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[8][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[9][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[9][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[9][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[9][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[9][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[9][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[9][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[10][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[10][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[10][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[10][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[10][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[10][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[10][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[11][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[11][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[11][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[11][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[11][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[11][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[11][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[12][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[12][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[12][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[12][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[12][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[12][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[12][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[13][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[13][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[13][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[13][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[13][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[13][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[13][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[14][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[14][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[14][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[14][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[14][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[14][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[14][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[15][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[15][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[15][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[15][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[15][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[15][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[15][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[16][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[16][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[16][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[16][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[16][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[16][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[16][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[17][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[17][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[17][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[17][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[17][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[17][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[17][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[18][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[18][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[18][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[18][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[18][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[18][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[18][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[19][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[19][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[19][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[19][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[19][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[19][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[20][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[20][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[20][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[20][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[20][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[20][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[20][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[21][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[21][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[21][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[21][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[21][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[21][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[21][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[22][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[22][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[22][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[22][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[22][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[22][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[22][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[23][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[23][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[23][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[23][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[23][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[23][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[23][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[24][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[24][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[24][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[24][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[24][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[24][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[24][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[25][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[25][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[25][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[25][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[25][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[25][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[25][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[26][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[26][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[26][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[26][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[26][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[26][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[26][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[27][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[27][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[27][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[27][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[27][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[27][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[27][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[28][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[28][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[28][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[28][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[28][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[28][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[28][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCDArray[29][0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[29][1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[29][2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[29][3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[29][4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[29][5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCDArray[29][6]); //number of those infected with HCV upon entry
        
        /*for(zInfCtrR=0; zInfCtrR<30; zInfCtrR++){ //correctional centres
            totalAge0 = totalAge0+zoneMinCIArray[zInfCtrR][0]+zoneMedCIArray[zInfCtrR][0]+zoneMaxCIArray[zInfCtrR][0];
            totalAge1 = totalAge1+zoneMinCIArray[zInfCtrR][1]+zoneMedCIArray[zInfCtrR][1]+zoneMaxCIArray[zInfCtrR][1];
            totalAge2 = totalAge2+zoneMinCIArray[zInfCtrR][2]+zoneMedCIArray[zInfCtrR][2]+zoneMaxCIArray[zInfCtrR][2];
            totalAge3 = totalAge3+zoneMinCIArray[zInfCtrR][3]+zoneMedCIArray[zInfCtrR][3]+zoneMaxCIArray[zInfCtrR][3];
            totalAge4 = totalAge4+zoneMinCIArray[zInfCtrR][4]+zoneMedCIArray[zInfCtrR][4]+zoneMaxCIArray[zInfCtrR][4];
            totalAge5 = totalAge5+zoneMinCIArray[zInfCtrR][5]+zoneMedCIArray[zInfCtrR][5]+zoneMaxCIArray[zInfCtrR][5];
            totalAge6 = totalAge6+zoneMinCIArray[zInfCtrR][6]+zoneMedCIArray[zInfCtrR][6]+zoneMaxCIArray[zInfCtrR][6];
            totalAD0 = totalAD0+zoneMinCDArray[zInfCtrR][0]+zoneMedCDArray[zInfCtrR][0]+zoneMaxCDArray[zInfCtrR][0];
            totalAD1 = totalAD1+zoneMinCDArray[zInfCtrR][1]+zoneMedCDArray[zInfCtrR][1]+zoneMaxCDArray[zInfCtrR][1];
            totalAD2 = totalAD2+zoneMinCDArray[zInfCtrR][2]+zoneMedCDArray[zInfCtrR][2]+zoneMaxCDArray[zInfCtrR][2];
            totalAD3 = totalAD3+zoneMinCDArray[zInfCtrR][3]+zoneMedCDArray[zInfCtrR][3]+zoneMaxCDArray[zInfCtrR][3];
            totalAD4 = totalAD4+zoneMinCDArray[zInfCtrR][4]+zoneMedCDArray[zInfCtrR][4]+zoneMaxCDArray[zInfCtrR][4];
            totalAD5 = totalAD5+zoneMinCDArray[zInfCtrR][5]+zoneMedCDArray[zInfCtrR][5]+zoneMaxCDArray[zInfCtrR][5];
            totalAD6 = totalAD6+zoneMinCDArray[zInfCtrR][6]+zoneMedCDArray[zInfCtrR][6]+zoneMaxCDArray[zInfCtrR][6];
        }*/
        
        fprintf(fw, "%d\t", totalAge0); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAge1); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAge2); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAge3); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAge4); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAge5); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAge6); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", totalAD0); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAD1); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAD2); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAD3); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAD4); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAD5); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalAD6); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", nCourtPop[0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[7]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[8]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[9]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[10]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[11]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[12]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[13]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[14]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[15]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[16]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[17]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[18]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[19]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[20]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[21]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[22]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[23]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[24]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[25]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[26]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[27]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[28]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[29]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[30]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[31]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[32]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[33]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[34]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[35]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[36]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nCourtPop[37]); //number of those infected with HCV upon entry
 
        fprintf(fw, "%d\t", courtCIArray[0][0]+courtCIArray[0][1]+courtCIArray[0][2]+courtCIArray[0][3]+courtCIArray[0][4]+courtCIArray[0][5]+courtCIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[1][0]+courtCIArray[1][1]+courtCIArray[1][2]+courtCIArray[1][3]+courtCIArray[1][4]+courtCIArray[1][5]+courtCIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[2][0]+courtCIArray[2][1]+courtCIArray[2][2]+courtCIArray[2][3]+courtCIArray[2][4]+courtCIArray[2][5]+courtCIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[3][0]+courtCIArray[3][1]+courtCIArray[3][2]+courtCIArray[3][3]+courtCIArray[3][4]+courtCIArray[3][5]+courtCIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[4][0]+courtCIArray[4][1]+courtCIArray[4][2]+courtCIArray[4][3]+courtCIArray[4][4]+courtCIArray[4][5]+courtCIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[5][0]+courtCIArray[5][1]+courtCIArray[5][2]+courtCIArray[5][3]+courtCIArray[5][4]+courtCIArray[5][5]+courtCIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[6][0]+courtCIArray[6][1]+courtCIArray[6][2]+courtCIArray[6][3]+courtCIArray[6][4]+courtCIArray[6][5]+courtCIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[7][0]+courtCIArray[7][1]+courtCIArray[7][2]+courtCIArray[7][3]+courtCIArray[7][4]+courtCIArray[7][5]+courtCIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[8][0]+courtCIArray[8][1]+courtCIArray[8][2]+courtCIArray[8][3]+courtCIArray[8][4]+courtCIArray[8][5]+courtCIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[9][0]+courtCIArray[9][1]+courtCIArray[9][2]+courtCIArray[9][3]+courtCIArray[9][4]+courtCIArray[9][5]+courtCIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[10][0]+courtCIArray[10][1]+courtCIArray[10][2]+courtCIArray[10][3]+courtCIArray[10][4]+courtCIArray[10][5]+courtCIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[11][0]+courtCIArray[11][1]+courtCIArray[11][2]+courtCIArray[11][3]+courtCIArray[11][4]+courtCIArray[11][5]+courtCIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[12][0]+courtCIArray[12][1]+courtCIArray[12][2]+courtCIArray[12][3]+courtCIArray[12][4]+courtCIArray[12][5]+courtCIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[13][0]+courtCIArray[13][1]+courtCIArray[13][2]+courtCIArray[13][3]+courtCIArray[13][4]+courtCIArray[13][5]+courtCIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[14][0]+courtCIArray[14][1]+courtCIArray[14][2]+courtCIArray[14][3]+courtCIArray[14][4]+courtCIArray[14][5]+courtCIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[15][0]+courtCIArray[15][1]+courtCIArray[15][2]+courtCIArray[15][3]+courtCIArray[15][4]+courtCIArray[15][5]+courtCIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[16][0]+courtCIArray[16][1]+courtCIArray[16][2]+courtCIArray[16][3]+courtCIArray[16][4]+courtCIArray[16][5]+courtCIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[17][0]+courtCIArray[17][1]+courtCIArray[17][2]+courtCIArray[17][3]+courtCIArray[17][4]+courtCIArray[17][5]+courtCIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[18][0]+courtCIArray[18][1]+courtCIArray[18][2]+courtCIArray[18][3]+courtCIArray[18][4]+courtCIArray[18][5]+courtCIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[19][0]+courtCIArray[19][1]+courtCIArray[19][2]+courtCIArray[19][3]+courtCIArray[19][4]+courtCIArray[19][5]+courtCIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[20][0]+courtCIArray[20][1]+courtCIArray[20][2]+courtCIArray[20][3]+courtCIArray[20][4]+courtCIArray[20][5]+courtCIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[21][0]+courtCIArray[21][1]+courtCIArray[21][2]+courtCIArray[21][3]+courtCIArray[21][4]+courtCIArray[21][5]+courtCIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[22][0]+courtCIArray[22][1]+courtCIArray[22][2]+courtCIArray[22][3]+courtCIArray[22][4]+courtCIArray[22][5]+courtCIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[23][0]+courtCIArray[23][1]+courtCIArray[23][2]+courtCIArray[23][3]+courtCIArray[23][4]+courtCIArray[23][5]+courtCIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[24][0]+courtCIArray[24][1]+courtCIArray[24][2]+courtCIArray[24][3]+courtCIArray[24][4]+courtCIArray[24][5]+courtCIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[25][0]+courtCIArray[25][1]+courtCIArray[25][2]+courtCIArray[25][3]+courtCIArray[25][4]+courtCIArray[25][5]+courtCIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[26][0]+courtCIArray[26][1]+courtCIArray[26][2]+courtCIArray[26][3]+courtCIArray[26][4]+courtCIArray[26][5]+courtCIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[27][0]+courtCIArray[27][1]+courtCIArray[27][2]+courtCIArray[27][3]+courtCIArray[27][4]+courtCIArray[27][5]+courtCIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[28][0]+courtCIArray[28][1]+courtCIArray[28][2]+courtCIArray[28][3]+courtCIArray[28][4]+courtCIArray[28][5]+courtCIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[29][0]+courtCIArray[29][1]+courtCIArray[29][2]+courtCIArray[29][3]+courtCIArray[29][4]+courtCIArray[29][5]+courtCIArray[29][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[30][0]+courtCIArray[30][1]+courtCIArray[30][2]+courtCIArray[30][3]+courtCIArray[30][4]+courtCIArray[30][5]+courtCIArray[30][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[31][0]+courtCIArray[31][1]+courtCIArray[31][2]+courtCIArray[31][3]+courtCIArray[31][4]+courtCIArray[31][5]+courtCIArray[31][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[32][0]+courtCIArray[32][1]+courtCIArray[32][2]+courtCIArray[32][3]+courtCIArray[32][4]+courtCIArray[32][5]+courtCIArray[32][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[33][0]+courtCIArray[33][1]+courtCIArray[33][2]+courtCIArray[33][3]+courtCIArray[33][4]+courtCIArray[33][5]+courtCIArray[33][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[34][0]+courtCIArray[34][1]+courtCIArray[34][2]+courtCIArray[34][3]+courtCIArray[34][4]+courtCIArray[34][5]+courtCIArray[34][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[35][0]+courtCIArray[35][1]+courtCIArray[35][2]+courtCIArray[35][3]+courtCIArray[35][4]+courtCIArray[35][5]+courtCIArray[35][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[36][0]+courtCIArray[36][1]+courtCIArray[36][2]+courtCIArray[36][3]+courtCIArray[36][4]+courtCIArray[36][5]+courtCIArray[36][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", courtCIArray[37][0]+courtCIArray[37][1]+courtCIArray[37][2]+courtCIArray[37][3]+courtCIArray[37][4]+courtCIArray[37][5]+courtCIArray[37][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", nMovingPop[0]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[1]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[2]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[3]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[4]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[5]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[7]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[8]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[9]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[10]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[11]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[12]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[13]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[14]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[15]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[16]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[17]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[18]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", nMovingPop[19]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", truckCIArray[0][0]+truckCIArray[0][1]+truckCIArray[0][2]+truckCIArray[0][3]+truckCIArray[0][4]+truckCIArray[0][5]+truckCIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[1][0]+truckCIArray[1][1]+truckCIArray[1][2]+truckCIArray[1][3]+truckCIArray[1][4]+truckCIArray[1][5]+truckCIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[2][0]+truckCIArray[2][1]+truckCIArray[2][2]+truckCIArray[2][3]+truckCIArray[2][4]+truckCIArray[2][5]+truckCIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[3][0]+truckCIArray[3][1]+truckCIArray[3][2]+truckCIArray[3][3]+truckCIArray[3][4]+truckCIArray[3][5]+truckCIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[4][0]+truckCIArray[4][1]+truckCIArray[4][2]+truckCIArray[4][3]+truckCIArray[4][4]+truckCIArray[4][5]+truckCIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[5][0]+truckCIArray[5][1]+truckCIArray[5][2]+truckCIArray[5][3]+truckCIArray[5][4]+truckCIArray[5][5]+truckCIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[6][0]+truckCIArray[6][1]+truckCIArray[6][2]+truckCIArray[6][3]+truckCIArray[6][4]+truckCIArray[6][5]+truckCIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[7][0]+truckCIArray[7][1]+truckCIArray[7][2]+truckCIArray[7][3]+truckCIArray[7][4]+truckCIArray[7][5]+truckCIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[8][0]+truckCIArray[8][1]+truckCIArray[8][2]+truckCIArray[8][3]+truckCIArray[8][4]+truckCIArray[8][5]+truckCIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[9][0]+truckCIArray[9][1]+truckCIArray[9][2]+truckCIArray[9][3]+truckCIArray[9][4]+truckCIArray[9][5]+truckCIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[10][0]+truckCIArray[10][1]+truckCIArray[10][2]+truckCIArray[10][3]+truckCIArray[10][4]+truckCIArray[10][5]+truckCIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[11][0]+truckCIArray[11][1]+truckCIArray[11][2]+truckCIArray[11][3]+truckCIArray[11][4]+truckCIArray[11][5]+truckCIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[12][0]+truckCIArray[12][1]+truckCIArray[12][2]+truckCIArray[12][3]+truckCIArray[12][4]+truckCIArray[12][5]+truckCIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[13][0]+truckCIArray[13][1]+truckCIArray[13][2]+truckCIArray[13][3]+truckCIArray[13][4]+truckCIArray[13][5]+truckCIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[14][0]+truckCIArray[14][1]+truckCIArray[14][2]+truckCIArray[14][3]+truckCIArray[14][4]+truckCIArray[14][5]+truckCIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[15][0]+truckCIArray[15][1]+truckCIArray[15][2]+truckCIArray[15][3]+truckCIArray[15][4]+truckCIArray[15][5]+truckCIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[16][0]+truckCIArray[16][1]+truckCIArray[16][2]+truckCIArray[16][3]+truckCIArray[16][4]+truckCIArray[16][5]+truckCIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[17][0]+truckCIArray[17][1]+truckCIArray[17][2]+truckCIArray[17][3]+truckCIArray[17][4]+truckCIArray[17][5]+truckCIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[18][0]+truckCIArray[18][1]+truckCIArray[18][2]+truckCIArray[18][3]+truckCIArray[18][4]+truckCIArray[18][5]+truckCIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", truckCIArray[19][0]+truckCIArray[19][1]+truckCIArray[19][2]+truckCIArray[19][3]+truckCIArray[19][4]+truckCIArray[19][5]+truckCIArray[19][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", exp19);
        fprintf(fw, "%d\t", exp20to44);
        fprintf(fw, "%d\t", exp45to54);
        fprintf(fw, "%d\t", exp55to64);
        fprintf(fw, "%d\t", exp65to74);
        fprintf(fw, "%d\t", exp75to84);
        fprintf(fw, "%d\t", exp85);
        
        fprintf(fw, "%d\t", pre19);
        fprintf(fw, "%d\t", pre20to44);
        fprintf(fw, "%d\t", pre45to54);
        fprintf(fw, "%d\t", pre55to64);
        fprintf(fw, "%d\t", pre65to74);
        fprintf(fw, "%d\t", pre75to84);
        fprintf(fw, "%d\t", pre85);
        
        fprintf(fw, "%d\t", asy19);
        fprintf(fw, "%d\t", asy20to44);
        fprintf(fw, "%d\t", asy45to54);
        fprintf(fw, "%d\t", asy55to64);
        fprintf(fw, "%d\t", asy65to74);
        fprintf(fw, "%d\t", asy75to84);
        fprintf(fw, "%d\t", asy85);

        fprintf(fw, "%d\t", mil19);
        fprintf(fw, "%d\t", mil20to44);
        fprintf(fw, "%d\t", mil45to54);
        fprintf(fw, "%d\t", mil55to64);
        fprintf(fw, "%d\t", mil65to74);
        fprintf(fw, "%d\t", mil75to84);
        fprintf(fw, "%d\t", mil85);

        fprintf(fw, "%d\t", mod19);
        fprintf(fw, "%d\t", mod20to44);
        fprintf(fw, "%d\t", mod45to54);
        fprintf(fw, "%d\t", mod55to64);
        fprintf(fw, "%d\t", mod65to74);
        fprintf(fw, "%d\t", mod75to84);
        fprintf(fw, "%d\t", mod85);

        fprintf(fw, "%d\t", sev19);
        fprintf(fw, "%d\t", sev20to44);
        fprintf(fw, "%d\t", sev45to54);
        fprintf(fw, "%d\t", sev55to64);
        fprintf(fw, "%d\t", sev65to74);
        fprintf(fw, "%d\t", sev75to84);
        fprintf(fw, "%d\t", sev85);

        fprintf(fw, "%d\t", cle19);
        fprintf(fw, "%d\t", cle20to44);
        fprintf(fw, "%d\t", cle45to54);
        fprintf(fw, "%d\t", cle55to64);
        fprintf(fw, "%d\t", cle65to74);
        fprintf(fw, "%d\t", cle75to84);
        fprintf(fw, "%d\t", cle85);
        
        //cumulative infections
        fprintf(fw, "%d\t", zoneMinCIArray[0][0]+zoneMinCIArray[0][1]+zoneMinCIArray[0][2]+zoneMinCIArray[0][3]+zoneMinCIArray[0][4]+zoneMinCIArray[0][5]+zoneMinCIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[1][0]+zoneMinCIArray[1][1]+zoneMinCIArray[1][2]+zoneMinCIArray[1][3]+zoneMinCIArray[1][4]+zoneMinCIArray[1][5]+zoneMinCIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[2][0]+zoneMinCIArray[2][1]+zoneMinCIArray[2][2]+zoneMinCIArray[2][3]+zoneMinCIArray[2][4]+zoneMinCIArray[2][5]+zoneMinCIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[3][0]+zoneMinCIArray[3][1]+zoneMinCIArray[3][2]+zoneMinCIArray[3][3]+zoneMinCIArray[3][4]+zoneMinCIArray[3][5]+zoneMinCIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[4][0]+zoneMinCIArray[4][1]+zoneMinCIArray[4][2]+zoneMinCIArray[4][3]+zoneMinCIArray[4][4]+zoneMinCIArray[4][5]+zoneMinCIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[5][0]+zoneMinCIArray[5][1]+zoneMinCIArray[5][2]+zoneMinCIArray[5][3]+zoneMinCIArray[5][4]+zoneMinCIArray[5][5]+zoneMinCIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[6][0]+zoneMinCIArray[6][1]+zoneMinCIArray[6][2]+zoneMinCIArray[6][3]+zoneMinCIArray[6][4]+zoneMinCIArray[6][5]+zoneMinCIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[7][0]+zoneMinCIArray[7][1]+zoneMinCIArray[7][2]+zoneMinCIArray[7][3]+zoneMinCIArray[7][4]+zoneMinCIArray[7][5]+zoneMinCIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[8][0]+zoneMinCIArray[8][1]+zoneMinCIArray[8][2]+zoneMinCIArray[8][3]+zoneMinCIArray[8][4]+zoneMinCIArray[8][5]+zoneMinCIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[9][0]+zoneMinCIArray[9][1]+zoneMinCIArray[9][2]+zoneMinCIArray[9][3]+zoneMinCIArray[9][4]+zoneMinCIArray[9][5]+zoneMinCIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[10][0]+zoneMinCIArray[10][1]+zoneMinCIArray[10][2]+zoneMinCIArray[10][3]+zoneMinCIArray[10][4]+zoneMinCIArray[10][5]+zoneMinCIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[11][0]+zoneMinCIArray[11][1]+zoneMinCIArray[11][2]+zoneMinCIArray[11][3]+zoneMinCIArray[11][4]+zoneMinCIArray[11][5]+zoneMinCIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[12][0]+zoneMinCIArray[12][1]+zoneMinCIArray[12][2]+zoneMinCIArray[12][3]+zoneMinCIArray[12][4]+zoneMinCIArray[12][5]+zoneMinCIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[13][0]+zoneMinCIArray[13][1]+zoneMinCIArray[13][2]+zoneMinCIArray[13][3]+zoneMinCIArray[13][4]+zoneMinCIArray[13][5]+zoneMinCIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[14][0]+zoneMinCIArray[14][1]+zoneMinCIArray[14][2]+zoneMinCIArray[14][3]+zoneMinCIArray[14][4]+zoneMinCIArray[14][5]+zoneMinCIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[15][0]+zoneMinCIArray[15][1]+zoneMinCIArray[15][2]+zoneMinCIArray[15][3]+zoneMinCIArray[15][4]+zoneMinCIArray[15][5]+zoneMinCIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[16][0]+zoneMinCIArray[16][1]+zoneMinCIArray[16][2]+zoneMinCIArray[16][3]+zoneMinCIArray[16][4]+zoneMinCIArray[16][5]+zoneMinCIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[17][0]+zoneMinCIArray[17][1]+zoneMinCIArray[17][2]+zoneMinCIArray[17][3]+zoneMinCIArray[17][4]+zoneMinCIArray[17][5]+zoneMinCIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[18][0]+zoneMinCIArray[18][1]+zoneMinCIArray[18][2]+zoneMinCIArray[18][3]+zoneMinCIArray[18][4]+zoneMinCIArray[18][5]+zoneMinCIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[19][0]+zoneMinCIArray[19][1]+zoneMinCIArray[19][2]+zoneMinCIArray[19][3]+zoneMinCIArray[19][4]+zoneMinCIArray[19][5]+zoneMinCIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[20][0]+zoneMinCIArray[20][1]+zoneMinCIArray[20][2]+zoneMinCIArray[20][3]+zoneMinCIArray[20][4]+zoneMinCIArray[20][5]+zoneMinCIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[21][0]+zoneMinCIArray[21][1]+zoneMinCIArray[21][2]+zoneMinCIArray[21][3]+zoneMinCIArray[21][4]+zoneMinCIArray[21][5]+zoneMinCIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[22][0]+zoneMinCIArray[22][1]+zoneMinCIArray[22][2]+zoneMinCIArray[22][3]+zoneMinCIArray[22][4]+zoneMinCIArray[22][5]+zoneMinCIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[23][0]+zoneMinCIArray[23][1]+zoneMinCIArray[23][2]+zoneMinCIArray[23][3]+zoneMinCIArray[23][4]+zoneMinCIArray[23][5]+zoneMinCIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[24][0]+zoneMinCIArray[24][1]+zoneMinCIArray[24][2]+zoneMinCIArray[24][3]+zoneMinCIArray[24][4]+zoneMinCIArray[24][5]+zoneMinCIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[25][0]+zoneMinCIArray[25][1]+zoneMinCIArray[25][2]+zoneMinCIArray[25][3]+zoneMinCIArray[25][4]+zoneMinCIArray[25][5]+zoneMinCIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[26][0]+zoneMinCIArray[26][1]+zoneMinCIArray[26][2]+zoneMinCIArray[26][3]+zoneMinCIArray[26][4]+zoneMinCIArray[26][5]+zoneMinCIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[27][0]+zoneMinCIArray[27][1]+zoneMinCIArray[27][2]+zoneMinCIArray[27][3]+zoneMinCIArray[27][4]+zoneMinCIArray[27][5]+zoneMinCIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[28][0]+zoneMinCIArray[28][1]+zoneMinCIArray[28][2]+zoneMinCIArray[28][3]+zoneMinCIArray[28][4]+zoneMinCIArray[28][5]+zoneMinCIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinCIArray[29][0]+zoneMinCIArray[29][1]+zoneMinCIArray[29][2]+zoneMinCIArray[29][3]+zoneMinCIArray[29][4]+zoneMinCIArray[29][5]+zoneMinCIArray[29][6]); //number of those infected with HCV upon entry

        fprintf(fw, "%d\t", zoneMedCIArray[0][0]+zoneMedCIArray[0][1]+zoneMedCIArray[0][2]+zoneMedCIArray[0][3]+zoneMedCIArray[0][4]+zoneMedCIArray[0][5]+zoneMedCIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[1][0]+zoneMedCIArray[1][1]+zoneMedCIArray[1][2]+zoneMedCIArray[1][3]+zoneMedCIArray[1][4]+zoneMedCIArray[1][5]+zoneMedCIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[2][0]+zoneMedCIArray[2][1]+zoneMedCIArray[2][2]+zoneMedCIArray[2][3]+zoneMedCIArray[2][4]+zoneMedCIArray[2][5]+zoneMedCIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[3][0]+zoneMedCIArray[3][1]+zoneMedCIArray[3][2]+zoneMedCIArray[3][3]+zoneMedCIArray[3][4]+zoneMedCIArray[3][5]+zoneMedCIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[4][0]+zoneMedCIArray[4][1]+zoneMedCIArray[4][2]+zoneMedCIArray[4][3]+zoneMedCIArray[4][4]+zoneMedCIArray[4][5]+zoneMedCIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[5][0]+zoneMedCIArray[5][1]+zoneMedCIArray[5][2]+zoneMedCIArray[5][3]+zoneMedCIArray[5][4]+zoneMedCIArray[5][5]+zoneMedCIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[6][0]+zoneMedCIArray[6][1]+zoneMedCIArray[6][2]+zoneMedCIArray[6][3]+zoneMedCIArray[6][4]+zoneMedCIArray[6][5]+zoneMedCIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[7][0]+zoneMedCIArray[7][1]+zoneMedCIArray[7][2]+zoneMedCIArray[7][3]+zoneMedCIArray[7][4]+zoneMedCIArray[7][5]+zoneMedCIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[8][0]+zoneMedCIArray[8][1]+zoneMedCIArray[8][2]+zoneMedCIArray[8][3]+zoneMedCIArray[8][4]+zoneMedCIArray[8][5]+zoneMedCIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[9][0]+zoneMedCIArray[9][1]+zoneMedCIArray[9][2]+zoneMedCIArray[9][3]+zoneMedCIArray[9][4]+zoneMedCIArray[9][5]+zoneMedCIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[10][0]+zoneMedCIArray[10][1]+zoneMedCIArray[10][2]+zoneMedCIArray[10][3]+zoneMedCIArray[10][4]+zoneMedCIArray[10][5]+zoneMedCIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[11][0]+zoneMedCIArray[11][1]+zoneMedCIArray[11][2]+zoneMedCIArray[11][3]+zoneMedCIArray[11][4]+zoneMedCIArray[11][5]+zoneMedCIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[12][0]+zoneMedCIArray[12][1]+zoneMedCIArray[12][2]+zoneMedCIArray[12][3]+zoneMedCIArray[12][4]+zoneMedCIArray[12][5]+zoneMedCIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[13][0]+zoneMedCIArray[13][1]+zoneMedCIArray[13][2]+zoneMedCIArray[13][3]+zoneMedCIArray[13][4]+zoneMedCIArray[13][5]+zoneMedCIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[14][0]+zoneMedCIArray[14][1]+zoneMedCIArray[14][2]+zoneMedCIArray[14][3]+zoneMedCIArray[14][4]+zoneMedCIArray[14][5]+zoneMedCIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[15][0]+zoneMedCIArray[15][1]+zoneMedCIArray[15][2]+zoneMedCIArray[15][3]+zoneMedCIArray[15][4]+zoneMedCIArray[15][5]+zoneMedCIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[16][0]+zoneMedCIArray[16][1]+zoneMedCIArray[16][2]+zoneMedCIArray[16][3]+zoneMedCIArray[16][4]+zoneMedCIArray[16][5]+zoneMedCIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[17][0]+zoneMedCIArray[17][1]+zoneMedCIArray[17][2]+zoneMedCIArray[17][3]+zoneMedCIArray[17][4]+zoneMedCIArray[17][5]+zoneMedCIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[18][0]+zoneMedCIArray[18][1]+zoneMedCIArray[18][2]+zoneMedCIArray[18][3]+zoneMedCIArray[18][4]+zoneMedCIArray[18][5]+zoneMedCIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[19][0]+zoneMedCIArray[19][1]+zoneMedCIArray[19][2]+zoneMedCIArray[19][3]+zoneMedCIArray[19][4]+zoneMedCIArray[19][5]+zoneMedCIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[20][0]+zoneMedCIArray[20][1]+zoneMedCIArray[20][2]+zoneMedCIArray[20][3]+zoneMedCIArray[20][4]+zoneMedCIArray[20][5]+zoneMedCIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[21][0]+zoneMedCIArray[21][1]+zoneMedCIArray[21][2]+zoneMedCIArray[21][3]+zoneMedCIArray[21][4]+zoneMedCIArray[21][5]+zoneMedCIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[22][0]+zoneMedCIArray[22][1]+zoneMedCIArray[22][2]+zoneMedCIArray[22][3]+zoneMedCIArray[22][4]+zoneMedCIArray[22][5]+zoneMedCIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[23][0]+zoneMedCIArray[23][1]+zoneMedCIArray[23][2]+zoneMedCIArray[23][3]+zoneMedCIArray[23][4]+zoneMedCIArray[23][5]+zoneMedCIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[24][0]+zoneMedCIArray[24][1]+zoneMedCIArray[24][2]+zoneMedCIArray[24][3]+zoneMedCIArray[24][4]+zoneMedCIArray[24][5]+zoneMedCIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[25][0]+zoneMedCIArray[25][1]+zoneMedCIArray[25][2]+zoneMedCIArray[25][3]+zoneMedCIArray[25][4]+zoneMedCIArray[25][5]+zoneMedCIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[26][0]+zoneMedCIArray[26][1]+zoneMedCIArray[26][2]+zoneMedCIArray[26][3]+zoneMedCIArray[26][4]+zoneMedCIArray[26][5]+zoneMedCIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[27][0]+zoneMedCIArray[27][1]+zoneMedCIArray[27][2]+zoneMedCIArray[27][3]+zoneMedCIArray[27][4]+zoneMedCIArray[27][5]+zoneMedCIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[28][0]+zoneMedCIArray[28][1]+zoneMedCIArray[28][2]+zoneMedCIArray[28][3]+zoneMedCIArray[28][4]+zoneMedCIArray[28][5]+zoneMedCIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedCIArray[29][0]+zoneMedCIArray[29][1]+zoneMedCIArray[29][2]+zoneMedCIArray[29][3]+zoneMedCIArray[29][4]+zoneMedCIArray[29][5]+zoneMedCIArray[29][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxCIArray[0][0]+zoneMaxCIArray[0][1]+zoneMaxCIArray[0][2]+zoneMaxCIArray[0][3]+zoneMaxCIArray[0][4]+zoneMaxCIArray[0][5]+zoneMaxCIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[1][0]+zoneMaxCIArray[1][1]+zoneMaxCIArray[1][2]+zoneMaxCIArray[1][3]+zoneMaxCIArray[1][4]+zoneMaxCIArray[1][5]+zoneMaxCIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[2][0]+zoneMaxCIArray[2][1]+zoneMaxCIArray[2][2]+zoneMaxCIArray[2][3]+zoneMaxCIArray[2][4]+zoneMaxCIArray[2][5]+zoneMaxCIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[3][0]+zoneMaxCIArray[3][1]+zoneMaxCIArray[3][2]+zoneMaxCIArray[3][3]+zoneMaxCIArray[3][4]+zoneMaxCIArray[3][5]+zoneMaxCIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[4][0]+zoneMaxCIArray[4][1]+zoneMaxCIArray[4][2]+zoneMaxCIArray[4][3]+zoneMaxCIArray[4][4]+zoneMaxCIArray[4][5]+zoneMaxCIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[5][0]+zoneMaxCIArray[5][1]+zoneMaxCIArray[5][2]+zoneMaxCIArray[5][3]+zoneMaxCIArray[5][4]+zoneMaxCIArray[5][5]+zoneMaxCIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[6][0]+zoneMaxCIArray[6][1]+zoneMaxCIArray[6][2]+zoneMaxCIArray[6][3]+zoneMaxCIArray[6][4]+zoneMaxCIArray[6][5]+zoneMaxCIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[7][0]+zoneMaxCIArray[7][1]+zoneMaxCIArray[7][2]+zoneMaxCIArray[7][3]+zoneMaxCIArray[7][4]+zoneMaxCIArray[7][5]+zoneMaxCIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[8][0]+zoneMaxCIArray[8][1]+zoneMaxCIArray[8][2]+zoneMaxCIArray[8][3]+zoneMaxCIArray[8][4]+zoneMaxCIArray[8][5]+zoneMaxCIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[9][0]+zoneMaxCIArray[9][1]+zoneMaxCIArray[9][2]+zoneMaxCIArray[9][3]+zoneMaxCIArray[9][4]+zoneMaxCIArray[9][5]+zoneMaxCIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[10][0]+zoneMaxCIArray[10][1]+zoneMaxCIArray[10][2]+zoneMaxCIArray[10][3]+zoneMaxCIArray[10][4]+zoneMaxCIArray[10][5]+zoneMaxCIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[11][0]+zoneMaxCIArray[11][1]+zoneMaxCIArray[11][2]+zoneMaxCIArray[11][3]+zoneMaxCIArray[11][4]+zoneMaxCIArray[11][5]+zoneMaxCIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[12][0]+zoneMaxCIArray[12][1]+zoneMaxCIArray[12][2]+zoneMaxCIArray[12][3]+zoneMaxCIArray[12][4]+zoneMaxCIArray[12][5]+zoneMaxCIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[13][0]+zoneMaxCIArray[13][1]+zoneMaxCIArray[13][2]+zoneMaxCIArray[13][3]+zoneMaxCIArray[13][4]+zoneMaxCIArray[13][5]+zoneMaxCIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[14][0]+zoneMaxCIArray[14][1]+zoneMaxCIArray[14][2]+zoneMaxCIArray[14][3]+zoneMaxCIArray[14][4]+zoneMaxCIArray[14][5]+zoneMaxCIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[15][0]+zoneMaxCIArray[15][1]+zoneMaxCIArray[15][2]+zoneMaxCIArray[15][3]+zoneMaxCIArray[15][4]+zoneMaxCIArray[15][5]+zoneMaxCIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[16][0]+zoneMaxCIArray[16][1]+zoneMaxCIArray[16][2]+zoneMaxCIArray[16][3]+zoneMaxCIArray[16][4]+zoneMaxCIArray[16][5]+zoneMaxCIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[17][0]+zoneMaxCIArray[17][1]+zoneMaxCIArray[17][2]+zoneMaxCIArray[17][3]+zoneMaxCIArray[17][4]+zoneMaxCIArray[17][5]+zoneMaxCIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[18][0]+zoneMaxCIArray[18][1]+zoneMaxCIArray[18][2]+zoneMaxCIArray[18][3]+zoneMaxCIArray[18][4]+zoneMaxCIArray[18][5]+zoneMaxCIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[19][0]+zoneMaxCIArray[19][1]+zoneMaxCIArray[19][2]+zoneMaxCIArray[19][3]+zoneMaxCIArray[19][4]+zoneMaxCIArray[19][5]+zoneMaxCIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[20][0]+zoneMaxCIArray[20][1]+zoneMaxCIArray[20][2]+zoneMaxCIArray[20][3]+zoneMaxCIArray[20][4]+zoneMaxCIArray[20][5]+zoneMaxCIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[21][0]+zoneMaxCIArray[21][1]+zoneMaxCIArray[21][2]+zoneMaxCIArray[21][3]+zoneMaxCIArray[21][4]+zoneMaxCIArray[21][5]+zoneMaxCIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[22][0]+zoneMaxCIArray[22][1]+zoneMaxCIArray[22][2]+zoneMaxCIArray[22][3]+zoneMaxCIArray[22][4]+zoneMaxCIArray[22][5]+zoneMaxCIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[23][0]+zoneMaxCIArray[23][1]+zoneMaxCIArray[23][2]+zoneMaxCIArray[23][3]+zoneMaxCIArray[23][4]+zoneMaxCIArray[23][5]+zoneMaxCIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[24][0]+zoneMaxCIArray[24][1]+zoneMaxCIArray[24][2]+zoneMaxCIArray[24][3]+zoneMaxCIArray[24][4]+zoneMaxCIArray[24][5]+zoneMaxCIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[25][0]+zoneMaxCIArray[25][1]+zoneMaxCIArray[25][2]+zoneMaxCIArray[25][3]+zoneMaxCIArray[25][4]+zoneMaxCIArray[25][5]+zoneMaxCIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[26][0]+zoneMaxCIArray[26][1]+zoneMaxCIArray[26][2]+zoneMaxCIArray[26][3]+zoneMaxCIArray[26][4]+zoneMaxCIArray[26][5]+zoneMaxCIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[27][0]+zoneMaxCIArray[27][1]+zoneMaxCIArray[27][2]+zoneMaxCIArray[27][3]+zoneMaxCIArray[27][4]+zoneMaxCIArray[27][5]+zoneMaxCIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[28][0]+zoneMaxCIArray[28][1]+zoneMaxCIArray[28][2]+zoneMaxCIArray[28][3]+zoneMaxCIArray[28][4]+zoneMaxCIArray[28][5]+zoneMaxCIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxCIArray[29][0]+zoneMaxCIArray[29][1]+zoneMaxCIArray[29][2]+zoneMaxCIArray[29][3]+zoneMaxCIArray[29][4]+zoneMaxCIArray[29][5]+zoneMaxCIArray[29][6]); //number of those infected with HCV upon entry
        
        //int finalSevExp=0, finalSevPre=0, finalSevAsy=0, finalSevMil=0, finalSevMod=0, finalSevSev=0, finalSevNon=0;
        fprintf(fw, "%d\t", finalSevExp);
        fprintf(fw, "%d\t", finalSevPre);
        fprintf(fw, "%d\t", finalSevAsy);
        fprintf(fw, "%d\t", finalSevMil);
        fprintf(fw, "%d\t", finalSevMod);
        fprintf(fw, "%d\t", finalSevSev);
        fprintf(fw, "%d\t", finalSevNon);
        
        fprintf(fw, "%d\t", zoneMinDIArray[0][0]+zoneMinDIArray[0][1]+zoneMinDIArray[0][2]+zoneMinDIArray[0][3]+zoneMinDIArray[0][4]+zoneMinDIArray[0][5]+zoneMinDIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[1][0]+zoneMinDIArray[1][1]+zoneMinDIArray[1][2]+zoneMinDIArray[1][3]+zoneMinDIArray[1][4]+zoneMinDIArray[1][5]+zoneMinDIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[2][0]+zoneMinDIArray[2][1]+zoneMinDIArray[2][2]+zoneMinDIArray[2][3]+zoneMinDIArray[2][4]+zoneMinDIArray[2][5]+zoneMinDIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[3][0]+zoneMinDIArray[3][1]+zoneMinDIArray[3][2]+zoneMinDIArray[3][3]+zoneMinDIArray[3][4]+zoneMinDIArray[3][5]+zoneMinDIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[4][0]+zoneMinDIArray[4][1]+zoneMinDIArray[4][2]+zoneMinDIArray[4][3]+zoneMinDIArray[4][4]+zoneMinDIArray[4][5]+zoneMinDIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[5][0]+zoneMinDIArray[5][1]+zoneMinDIArray[5][2]+zoneMinDIArray[5][3]+zoneMinDIArray[5][4]+zoneMinDIArray[5][5]+zoneMinDIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[6][0]+zoneMinDIArray[6][1]+zoneMinDIArray[6][2]+zoneMinDIArray[6][3]+zoneMinDIArray[6][4]+zoneMinDIArray[6][5]+zoneMinDIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[7][0]+zoneMinDIArray[7][1]+zoneMinDIArray[7][2]+zoneMinDIArray[7][3]+zoneMinDIArray[7][4]+zoneMinDIArray[7][5]+zoneMinDIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[8][0]+zoneMinDIArray[8][1]+zoneMinDIArray[8][2]+zoneMinDIArray[8][3]+zoneMinDIArray[8][4]+zoneMinDIArray[8][5]+zoneMinDIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[9][0]+zoneMinDIArray[9][1]+zoneMinDIArray[9][2]+zoneMinDIArray[9][3]+zoneMinDIArray[9][4]+zoneMinDIArray[9][5]+zoneMinDIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[10][0]+zoneMinDIArray[10][1]+zoneMinDIArray[10][2]+zoneMinDIArray[10][3]+zoneMinDIArray[10][4]+zoneMinDIArray[10][5]+zoneMinDIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[11][0]+zoneMinDIArray[11][1]+zoneMinDIArray[11][2]+zoneMinDIArray[11][3]+zoneMinDIArray[11][4]+zoneMinDIArray[11][5]+zoneMinDIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[12][0]+zoneMinDIArray[12][1]+zoneMinDIArray[12][2]+zoneMinDIArray[12][3]+zoneMinDIArray[12][4]+zoneMinDIArray[12][5]+zoneMinDIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[13][0]+zoneMinDIArray[13][1]+zoneMinDIArray[13][2]+zoneMinDIArray[13][3]+zoneMinDIArray[13][4]+zoneMinDIArray[13][5]+zoneMinDIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[14][0]+zoneMinDIArray[14][1]+zoneMinDIArray[14][2]+zoneMinDIArray[14][3]+zoneMinDIArray[14][4]+zoneMinDIArray[14][5]+zoneMinDIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[15][0]+zoneMinDIArray[15][1]+zoneMinDIArray[15][2]+zoneMinDIArray[15][3]+zoneMinDIArray[15][4]+zoneMinDIArray[15][5]+zoneMinDIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[16][0]+zoneMinDIArray[16][1]+zoneMinDIArray[16][2]+zoneMinDIArray[16][3]+zoneMinDIArray[16][4]+zoneMinDIArray[16][5]+zoneMinDIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[17][0]+zoneMinDIArray[17][1]+zoneMinDIArray[17][2]+zoneMinDIArray[17][3]+zoneMinDIArray[17][4]+zoneMinDIArray[17][5]+zoneMinDIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[18][0]+zoneMinDIArray[18][1]+zoneMinDIArray[18][2]+zoneMinDIArray[18][3]+zoneMinDIArray[18][4]+zoneMinDIArray[18][5]+zoneMinDIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[19][0]+zoneMinDIArray[19][1]+zoneMinDIArray[19][2]+zoneMinDIArray[19][3]+zoneMinDIArray[19][4]+zoneMinDIArray[19][5]+zoneMinDIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[20][0]+zoneMinDIArray[20][1]+zoneMinDIArray[20][2]+zoneMinDIArray[20][3]+zoneMinDIArray[20][4]+zoneMinDIArray[20][5]+zoneMinDIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[21][0]+zoneMinDIArray[21][1]+zoneMinDIArray[21][2]+zoneMinDIArray[21][3]+zoneMinDIArray[21][4]+zoneMinDIArray[21][5]+zoneMinDIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[22][0]+zoneMinDIArray[22][1]+zoneMinDIArray[22][2]+zoneMinDIArray[22][3]+zoneMinDIArray[22][4]+zoneMinDIArray[22][5]+zoneMinDIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[23][0]+zoneMinDIArray[23][1]+zoneMinDIArray[23][2]+zoneMinDIArray[23][3]+zoneMinDIArray[23][4]+zoneMinDIArray[23][5]+zoneMinDIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[24][0]+zoneMinDIArray[24][1]+zoneMinDIArray[24][2]+zoneMinDIArray[24][3]+zoneMinDIArray[24][4]+zoneMinDIArray[24][5]+zoneMinDIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[25][0]+zoneMinDIArray[25][1]+zoneMinDIArray[25][2]+zoneMinDIArray[25][3]+zoneMinDIArray[25][4]+zoneMinDIArray[25][5]+zoneMinDIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[26][0]+zoneMinDIArray[26][1]+zoneMinDIArray[26][2]+zoneMinDIArray[26][3]+zoneMinDIArray[26][4]+zoneMinDIArray[26][5]+zoneMinDIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[27][0]+zoneMinDIArray[27][1]+zoneMinDIArray[27][2]+zoneMinDIArray[27][3]+zoneMinDIArray[27][4]+zoneMinDIArray[27][5]+zoneMinDIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[28][0]+zoneMinDIArray[28][1]+zoneMinDIArray[28][2]+zoneMinDIArray[28][3]+zoneMinDIArray[28][4]+zoneMinDIArray[28][5]+zoneMinDIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMinDIArray[29][0]+zoneMinDIArray[29][1]+zoneMinDIArray[29][2]+zoneMinDIArray[29][3]+zoneMinDIArray[29][4]+zoneMinDIArray[29][5]+zoneMinDIArray[29][6]); //number of those infected with HCV upon entry

        fprintf(fw, "%d\t", zoneMedDIArray[0][0]+zoneMedDIArray[0][1]+zoneMedDIArray[0][2]+zoneMedDIArray[0][3]+zoneMedDIArray[0][4]+zoneMedDIArray[0][5]+zoneMedDIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[1][0]+zoneMedDIArray[1][1]+zoneMedDIArray[1][2]+zoneMedDIArray[1][3]+zoneMedDIArray[1][4]+zoneMedDIArray[1][5]+zoneMedDIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[2][0]+zoneMedDIArray[2][1]+zoneMedDIArray[2][2]+zoneMedDIArray[2][3]+zoneMedDIArray[2][4]+zoneMedDIArray[2][5]+zoneMedDIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[3][0]+zoneMedDIArray[3][1]+zoneMedDIArray[3][2]+zoneMedDIArray[3][3]+zoneMedDIArray[3][4]+zoneMedDIArray[3][5]+zoneMedDIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[4][0]+zoneMedDIArray[4][1]+zoneMedDIArray[4][2]+zoneMedDIArray[4][3]+zoneMedDIArray[4][4]+zoneMedDIArray[4][5]+zoneMedDIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[5][0]+zoneMedDIArray[5][1]+zoneMedDIArray[5][2]+zoneMedDIArray[5][3]+zoneMedDIArray[5][4]+zoneMedDIArray[5][5]+zoneMedDIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[6][0]+zoneMedDIArray[6][1]+zoneMedDIArray[6][2]+zoneMedDIArray[6][3]+zoneMedDIArray[6][4]+zoneMedDIArray[6][5]+zoneMedDIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[7][0]+zoneMedDIArray[7][1]+zoneMedDIArray[7][2]+zoneMedDIArray[7][3]+zoneMedDIArray[7][4]+zoneMedDIArray[7][5]+zoneMedDIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[8][0]+zoneMedDIArray[8][1]+zoneMedDIArray[8][2]+zoneMedDIArray[8][3]+zoneMedDIArray[8][4]+zoneMedDIArray[8][5]+zoneMedDIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[9][0]+zoneMedDIArray[9][1]+zoneMedDIArray[9][2]+zoneMedDIArray[9][3]+zoneMedDIArray[9][4]+zoneMedDIArray[9][5]+zoneMedDIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[10][0]+zoneMedDIArray[10][1]+zoneMedDIArray[10][2]+zoneMedDIArray[10][3]+zoneMedDIArray[10][4]+zoneMedDIArray[10][5]+zoneMedDIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[11][0]+zoneMedDIArray[11][1]+zoneMedDIArray[11][2]+zoneMedDIArray[11][3]+zoneMedDIArray[11][4]+zoneMedDIArray[11][5]+zoneMedDIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[12][0]+zoneMedDIArray[12][1]+zoneMedDIArray[12][2]+zoneMedDIArray[12][3]+zoneMedDIArray[12][4]+zoneMedDIArray[12][5]+zoneMedDIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[13][0]+zoneMedDIArray[13][1]+zoneMedDIArray[13][2]+zoneMedDIArray[13][3]+zoneMedDIArray[13][4]+zoneMedDIArray[13][5]+zoneMedDIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[14][0]+zoneMedDIArray[14][1]+zoneMedDIArray[14][2]+zoneMedDIArray[14][3]+zoneMedDIArray[14][4]+zoneMedDIArray[14][5]+zoneMedDIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[15][0]+zoneMedDIArray[15][1]+zoneMedDIArray[15][2]+zoneMedDIArray[15][3]+zoneMedDIArray[15][4]+zoneMedDIArray[15][5]+zoneMedDIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[16][0]+zoneMedDIArray[16][1]+zoneMedDIArray[16][2]+zoneMedDIArray[16][3]+zoneMedDIArray[16][4]+zoneMedDIArray[16][5]+zoneMedDIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[17][0]+zoneMedDIArray[17][1]+zoneMedDIArray[17][2]+zoneMedDIArray[17][3]+zoneMedDIArray[17][4]+zoneMedDIArray[17][5]+zoneMedDIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[18][0]+zoneMedDIArray[18][1]+zoneMedDIArray[18][2]+zoneMedDIArray[18][3]+zoneMedDIArray[18][4]+zoneMedDIArray[18][5]+zoneMedDIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[19][0]+zoneMedDIArray[19][1]+zoneMedDIArray[19][2]+zoneMedDIArray[19][3]+zoneMedDIArray[19][4]+zoneMedDIArray[19][5]+zoneMedDIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[20][0]+zoneMedDIArray[20][1]+zoneMedDIArray[20][2]+zoneMedDIArray[20][3]+zoneMedDIArray[20][4]+zoneMedDIArray[20][5]+zoneMedDIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[21][0]+zoneMedDIArray[21][1]+zoneMedDIArray[21][2]+zoneMedDIArray[21][3]+zoneMedDIArray[21][4]+zoneMedDIArray[21][5]+zoneMedDIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[22][0]+zoneMedDIArray[22][1]+zoneMedDIArray[22][2]+zoneMedDIArray[22][3]+zoneMedDIArray[22][4]+zoneMedDIArray[22][5]+zoneMedDIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[23][0]+zoneMedDIArray[23][1]+zoneMedDIArray[23][2]+zoneMedDIArray[23][3]+zoneMedDIArray[23][4]+zoneMedDIArray[23][5]+zoneMedDIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[24][0]+zoneMedDIArray[24][1]+zoneMedDIArray[24][2]+zoneMedDIArray[24][3]+zoneMedDIArray[24][4]+zoneMedDIArray[24][5]+zoneMedDIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[25][0]+zoneMedDIArray[25][1]+zoneMedDIArray[25][2]+zoneMedDIArray[25][3]+zoneMedDIArray[25][4]+zoneMedDIArray[25][5]+zoneMedDIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[26][0]+zoneMedDIArray[26][1]+zoneMedDIArray[26][2]+zoneMedDIArray[26][3]+zoneMedDIArray[26][4]+zoneMedDIArray[26][5]+zoneMedDIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[27][0]+zoneMedDIArray[27][1]+zoneMedDIArray[27][2]+zoneMedDIArray[27][3]+zoneMedDIArray[27][4]+zoneMedDIArray[27][5]+zoneMedDIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[28][0]+zoneMedDIArray[28][1]+zoneMedDIArray[28][2]+zoneMedDIArray[28][3]+zoneMedDIArray[28][4]+zoneMedDIArray[28][5]+zoneMedDIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMedDIArray[29][0]+zoneMedDIArray[29][1]+zoneMedDIArray[29][2]+zoneMedDIArray[29][3]+zoneMedDIArray[29][4]+zoneMedDIArray[29][5]+zoneMedDIArray[29][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", zoneMaxDIArray[0][0]+zoneMaxDIArray[0][1]+zoneMaxDIArray[0][2]+zoneMaxDIArray[0][3]+zoneMaxDIArray[0][4]+zoneMaxDIArray[0][5]+zoneMaxDIArray[0][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[1][0]+zoneMaxDIArray[1][1]+zoneMaxDIArray[1][2]+zoneMaxDIArray[1][3]+zoneMaxDIArray[1][4]+zoneMaxDIArray[1][5]+zoneMaxDIArray[1][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[2][0]+zoneMaxDIArray[2][1]+zoneMaxDIArray[2][2]+zoneMaxDIArray[2][3]+zoneMaxDIArray[2][4]+zoneMaxDIArray[2][5]+zoneMaxDIArray[2][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[3][0]+zoneMaxDIArray[3][1]+zoneMaxDIArray[3][2]+zoneMaxDIArray[3][3]+zoneMaxDIArray[3][4]+zoneMaxDIArray[3][5]+zoneMaxDIArray[3][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[4][0]+zoneMaxDIArray[4][1]+zoneMaxDIArray[4][2]+zoneMaxDIArray[4][3]+zoneMaxDIArray[4][4]+zoneMaxDIArray[4][5]+zoneMaxDIArray[4][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[5][0]+zoneMaxDIArray[5][1]+zoneMaxDIArray[5][2]+zoneMaxDIArray[5][3]+zoneMaxDIArray[5][4]+zoneMaxDIArray[5][5]+zoneMaxDIArray[5][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[6][0]+zoneMaxDIArray[6][1]+zoneMaxDIArray[6][2]+zoneMaxDIArray[6][3]+zoneMaxDIArray[6][4]+zoneMaxDIArray[6][5]+zoneMaxDIArray[6][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[7][0]+zoneMaxDIArray[7][1]+zoneMaxDIArray[7][2]+zoneMaxDIArray[7][3]+zoneMaxDIArray[7][4]+zoneMaxDIArray[7][5]+zoneMaxDIArray[7][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[8][0]+zoneMaxDIArray[8][1]+zoneMaxDIArray[8][2]+zoneMaxDIArray[8][3]+zoneMaxDIArray[8][4]+zoneMaxDIArray[8][5]+zoneMaxDIArray[8][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[9][0]+zoneMaxDIArray[9][1]+zoneMaxDIArray[9][2]+zoneMaxDIArray[9][3]+zoneMaxDIArray[9][4]+zoneMaxDIArray[9][5]+zoneMaxDIArray[9][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[10][0]+zoneMaxDIArray[10][1]+zoneMaxDIArray[10][2]+zoneMaxDIArray[10][3]+zoneMaxDIArray[10][4]+zoneMaxDIArray[10][5]+zoneMaxDIArray[10][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[11][0]+zoneMaxDIArray[11][1]+zoneMaxDIArray[11][2]+zoneMaxDIArray[11][3]+zoneMaxDIArray[11][4]+zoneMaxDIArray[11][5]+zoneMaxDIArray[11][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[12][0]+zoneMaxDIArray[12][1]+zoneMaxDIArray[12][2]+zoneMaxDIArray[12][3]+zoneMaxDIArray[12][4]+zoneMaxDIArray[12][5]+zoneMaxDIArray[12][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[13][0]+zoneMaxDIArray[13][1]+zoneMaxDIArray[13][2]+zoneMaxDIArray[13][3]+zoneMaxDIArray[13][4]+zoneMaxDIArray[13][5]+zoneMaxDIArray[13][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[14][0]+zoneMaxDIArray[14][1]+zoneMaxDIArray[14][2]+zoneMaxDIArray[14][3]+zoneMaxDIArray[14][4]+zoneMaxDIArray[14][5]+zoneMaxDIArray[14][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[15][0]+zoneMaxDIArray[15][1]+zoneMaxDIArray[15][2]+zoneMaxDIArray[15][3]+zoneMaxDIArray[15][4]+zoneMaxDIArray[15][5]+zoneMaxDIArray[15][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[16][0]+zoneMaxDIArray[16][1]+zoneMaxDIArray[16][2]+zoneMaxDIArray[16][3]+zoneMaxDIArray[16][4]+zoneMaxDIArray[16][5]+zoneMaxDIArray[16][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[17][0]+zoneMaxDIArray[17][1]+zoneMaxDIArray[17][2]+zoneMaxDIArray[17][3]+zoneMaxDIArray[17][4]+zoneMaxDIArray[17][5]+zoneMaxDIArray[17][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[18][0]+zoneMaxDIArray[18][1]+zoneMaxDIArray[18][2]+zoneMaxDIArray[18][3]+zoneMaxDIArray[18][4]+zoneMaxDIArray[18][5]+zoneMaxDIArray[18][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[19][0]+zoneMaxDIArray[19][1]+zoneMaxDIArray[19][2]+zoneMaxDIArray[19][3]+zoneMaxDIArray[19][4]+zoneMaxDIArray[19][5]+zoneMaxDIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[20][0]+zoneMaxDIArray[20][1]+zoneMaxDIArray[20][2]+zoneMaxDIArray[20][3]+zoneMaxDIArray[20][4]+zoneMaxDIArray[20][5]+zoneMaxDIArray[20][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[21][0]+zoneMaxDIArray[21][1]+zoneMaxDIArray[21][2]+zoneMaxDIArray[21][3]+zoneMaxDIArray[21][4]+zoneMaxDIArray[21][5]+zoneMaxDIArray[21][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[22][0]+zoneMaxDIArray[22][1]+zoneMaxDIArray[22][2]+zoneMaxDIArray[22][3]+zoneMaxDIArray[22][4]+zoneMaxDIArray[22][5]+zoneMaxDIArray[22][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[23][0]+zoneMaxDIArray[23][1]+zoneMaxDIArray[23][2]+zoneMaxDIArray[23][3]+zoneMaxDIArray[23][4]+zoneMaxDIArray[23][5]+zoneMaxDIArray[23][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[24][0]+zoneMaxDIArray[24][1]+zoneMaxDIArray[24][2]+zoneMaxDIArray[24][3]+zoneMaxDIArray[24][4]+zoneMaxDIArray[24][5]+zoneMaxDIArray[24][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[25][0]+zoneMaxDIArray[25][1]+zoneMaxDIArray[25][2]+zoneMaxDIArray[25][3]+zoneMaxDIArray[25][4]+zoneMaxDIArray[25][5]+zoneMaxDIArray[25][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[26][0]+zoneMaxDIArray[26][1]+zoneMaxDIArray[26][2]+zoneMaxDIArray[26][3]+zoneMaxDIArray[26][4]+zoneMaxDIArray[26][5]+zoneMaxDIArray[26][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[27][0]+zoneMaxDIArray[27][1]+zoneMaxDIArray[27][2]+zoneMaxDIArray[27][3]+zoneMaxDIArray[27][4]+zoneMaxDIArray[27][5]+zoneMaxDIArray[27][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[28][0]+zoneMaxDIArray[28][1]+zoneMaxDIArray[28][2]+zoneMaxDIArray[28][3]+zoneMaxDIArray[28][4]+zoneMaxDIArray[28][5]+zoneMaxDIArray[28][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", zoneMaxDIArray[29][0]+zoneMaxDIArray[29][1]+zoneMaxDIArray[29][2]+zoneMaxDIArray[29][3]+zoneMaxDIArray[29][4]+zoneMaxDIArray[29][5]+zoneMaxDIArray[29][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", courtDIArray[0][0]+courtDIArray[0][1]+courtDIArray[0][2]+courtDIArray[0][3]+courtDIArray[0][4]+courtDIArray[0][5]+courtDIArray[0][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[1][0]+courtDIArray[1][1]+courtDIArray[1][2]+courtDIArray[1][3]+courtDIArray[1][4]+courtDIArray[1][5]+courtDIArray[1][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[2][0]+courtDIArray[2][1]+courtDIArray[2][2]+courtDIArray[2][3]+courtDIArray[2][4]+courtDIArray[2][5]+courtDIArray[2][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[3][0]+courtDIArray[3][1]+courtDIArray[3][2]+courtDIArray[3][3]+courtDIArray[3][4]+courtDIArray[3][5]+courtDIArray[3][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[4][0]+courtDIArray[4][1]+courtDIArray[4][2]+courtDIArray[4][3]+courtDIArray[4][4]+courtDIArray[4][5]+courtDIArray[4][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[5][0]+courtDIArray[5][1]+courtDIArray[5][2]+courtDIArray[5][3]+courtDIArray[5][4]+courtDIArray[5][5]+courtDIArray[5][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[6][0]+courtDIArray[6][1]+courtDIArray[6][2]+courtDIArray[6][3]+courtDIArray[6][4]+courtDIArray[6][5]+courtDIArray[6][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[7][0]+courtDIArray[7][1]+courtDIArray[7][2]+courtDIArray[7][3]+courtDIArray[7][4]+courtDIArray[7][5]+courtDIArray[7][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[8][0]+courtDIArray[8][1]+courtDIArray[8][2]+courtDIArray[8][3]+courtDIArray[8][4]+courtDIArray[8][5]+courtDIArray[8][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[9][0]+courtDIArray[9][1]+courtDIArray[9][2]+courtDIArray[9][3]+courtDIArray[9][4]+courtDIArray[9][5]+courtDIArray[9][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[10][0]+courtDIArray[10][1]+courtDIArray[10][2]+courtDIArray[10][3]+courtDIArray[10][4]+courtDIArray[10][5]+courtDIArray[10][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[11][0]+courtDIArray[11][1]+courtDIArray[11][2]+courtDIArray[11][3]+courtDIArray[11][4]+courtDIArray[11][5]+courtDIArray[11][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[12][0]+courtDIArray[12][1]+courtDIArray[12][2]+courtDIArray[12][3]+courtDIArray[12][4]+courtDIArray[12][5]+courtDIArray[12][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[13][0]+courtDIArray[13][1]+courtDIArray[13][2]+courtDIArray[13][3]+courtDIArray[13][4]+courtDIArray[13][5]+courtDIArray[13][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[14][0]+courtDIArray[14][1]+courtDIArray[14][2]+courtDIArray[14][3]+courtDIArray[14][4]+courtDIArray[14][5]+courtDIArray[14][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[15][0]+courtDIArray[15][1]+courtDIArray[15][2]+courtDIArray[15][3]+courtDIArray[15][4]+courtDIArray[15][5]+courtDIArray[15][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[16][0]+courtDIArray[16][1]+courtDIArray[16][2]+courtDIArray[16][3]+courtDIArray[16][4]+courtDIArray[16][5]+courtDIArray[16][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[17][0]+courtDIArray[17][1]+courtDIArray[17][2]+courtDIArray[17][3]+courtDIArray[17][4]+courtDIArray[17][5]+courtDIArray[17][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[18][0]+courtDIArray[18][1]+courtDIArray[18][2]+courtDIArray[18][3]+courtDIArray[18][4]+courtDIArray[18][5]+courtDIArray[18][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[19][0]+courtDIArray[19][1]+courtDIArray[19][2]+courtDIArray[19][3]+courtDIArray[19][4]+courtDIArray[19][5]+courtDIArray[19][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[20][0]+courtDIArray[20][1]+courtDIArray[20][2]+courtDIArray[20][3]+courtDIArray[20][4]+courtDIArray[20][5]+courtDIArray[20][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[21][0]+courtDIArray[21][1]+courtDIArray[21][2]+courtDIArray[21][3]+courtDIArray[21][4]+courtDIArray[21][5]+courtDIArray[21][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[22][0]+courtDIArray[22][1]+courtDIArray[22][2]+courtDIArray[22][3]+courtDIArray[22][4]+courtDIArray[22][5]+courtDIArray[22][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[23][0]+courtDIArray[23][1]+courtDIArray[23][2]+courtDIArray[23][3]+courtDIArray[23][4]+courtDIArray[23][5]+courtDIArray[23][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[24][0]+courtDIArray[24][1]+courtDIArray[24][2]+courtDIArray[24][3]+courtDIArray[24][4]+courtDIArray[24][5]+courtDIArray[24][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[25][0]+courtDIArray[25][1]+courtDIArray[25][2]+courtDIArray[25][3]+courtDIArray[25][4]+courtDIArray[25][5]+courtDIArray[25][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[26][0]+courtDIArray[26][1]+courtDIArray[26][2]+courtDIArray[26][3]+courtDIArray[26][4]+courtDIArray[26][5]+courtDIArray[26][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[27][0]+courtDIArray[27][1]+courtDIArray[27][2]+courtDIArray[27][3]+courtDIArray[27][4]+courtDIArray[27][5]+courtDIArray[27][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[28][0]+courtDIArray[28][1]+courtDIArray[28][2]+courtDIArray[28][3]+courtDIArray[28][4]+courtDIArray[28][5]+courtDIArray[28][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[29][0]+courtDIArray[29][1]+courtDIArray[29][2]+courtDIArray[29][3]+courtDIArray[29][4]+courtDIArray[29][5]+courtDIArray[29][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[30][0]+courtDIArray[30][1]+courtDIArray[30][2]+courtDIArray[30][3]+courtDIArray[30][4]+courtDIArray[30][5]+courtDIArray[30][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[31][0]+courtDIArray[31][1]+courtDIArray[31][2]+courtDIArray[31][3]+courtDIArray[31][4]+courtDIArray[31][5]+courtDIArray[31][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[32][0]+courtDIArray[32][1]+courtDIArray[32][2]+courtDIArray[32][3]+courtDIArray[32][4]+courtDIArray[32][5]+courtDIArray[32][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[33][0]+courtDIArray[33][1]+courtDIArray[33][2]+courtDIArray[33][3]+courtDIArray[33][4]+courtDIArray[33][5]+courtDIArray[33][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[34][0]+courtDIArray[34][1]+courtDIArray[34][2]+courtDIArray[34][3]+courtDIArray[34][4]+courtDIArray[34][5]+courtDIArray[34][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[35][0]+courtDIArray[35][1]+courtDIArray[35][2]+courtDIArray[35][3]+courtDIArray[35][4]+courtDIArray[35][5]+courtDIArray[35][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[36][0]+courtDIArray[36][1]+courtDIArray[36][2]+courtDIArray[36][3]+courtDIArray[36][4]+courtDIArray[36][5]+courtDIArray[36][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", courtDIArray[37][0]+courtDIArray[37][1]+courtDIArray[37][2]+courtDIArray[37][3]+courtDIArray[37][4]+courtDIArray[37][5]+courtDIArray[37][6]); //number of those infected with HCV upon entry
        
        fprintf(fw, "%d\t", truckDIArray[0][0]+truckDIArray[0][1]+truckDIArray[0][2]+truckDIArray[0][3]+truckDIArray[0][4]+truckDIArray[0][5]+truckDIArray[0][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[1][0]+truckDIArray[1][1]+truckDIArray[1][2]+truckDIArray[1][3]+truckDIArray[1][4]+truckDIArray[1][5]+truckDIArray[1][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[2][0]+truckDIArray[2][1]+truckDIArray[2][2]+truckDIArray[2][3]+truckDIArray[2][4]+truckDIArray[2][5]+truckDIArray[2][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[3][0]+truckDIArray[3][1]+truckDIArray[3][2]+truckDIArray[3][3]+truckDIArray[3][4]+truckDIArray[3][5]+truckDIArray[3][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[4][0]+truckDIArray[4][1]+truckDIArray[4][2]+truckDIArray[4][3]+truckDIArray[4][4]+truckDIArray[4][5]+truckDIArray[4][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[5][0]+truckDIArray[5][1]+truckDIArray[5][2]+truckDIArray[5][3]+truckDIArray[5][4]+truckDIArray[5][5]+truckDIArray[5][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[6][0]+truckDIArray[6][1]+truckDIArray[6][2]+truckDIArray[6][3]+truckDIArray[6][4]+truckDIArray[6][5]+truckDIArray[6][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[7][0]+truckDIArray[7][1]+truckDIArray[7][2]+truckDIArray[7][3]+truckDIArray[7][4]+truckDIArray[7][5]+truckDIArray[7][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[8][0]+truckDIArray[8][1]+truckDIArray[8][2]+truckDIArray[8][3]+truckDIArray[8][4]+truckDIArray[8][5]+truckDIArray[8][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[9][0]+truckDIArray[9][1]+truckDIArray[9][2]+truckDIArray[9][3]+truckDIArray[9][4]+truckDIArray[9][5]+truckDIArray[9][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[10][0]+truckDIArray[10][1]+truckDIArray[10][2]+truckDIArray[10][3]+truckDIArray[10][4]+truckDIArray[10][5]+truckDIArray[10][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[11][0]+truckDIArray[11][1]+truckDIArray[11][2]+truckDIArray[11][3]+truckDIArray[11][4]+truckDIArray[11][5]+truckDIArray[11][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[12][0]+truckDIArray[12][1]+truckDIArray[12][2]+truckDIArray[12][3]+truckDIArray[12][4]+truckDIArray[12][5]+truckDIArray[12][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[13][0]+truckDIArray[13][1]+truckDIArray[13][2]+truckDIArray[13][3]+truckDIArray[13][4]+truckDIArray[13][5]+truckDIArray[13][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[14][0]+truckDIArray[14][1]+truckDIArray[14][2]+truckDIArray[14][3]+truckDIArray[14][4]+truckDIArray[14][5]+truckDIArray[14][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[15][0]+truckDIArray[15][1]+truckDIArray[15][2]+truckDIArray[15][3]+truckDIArray[15][4]+truckDIArray[15][5]+truckDIArray[15][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[16][0]+truckDIArray[16][1]+truckDIArray[16][2]+truckDIArray[16][3]+truckDIArray[16][4]+truckDIArray[16][5]+truckDIArray[16][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[17][0]+truckDIArray[17][1]+truckDIArray[17][2]+truckDIArray[17][3]+truckDIArray[17][4]+truckDIArray[17][5]+truckDIArray[17][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[18][0]+truckDIArray[18][1]+truckDIArray[18][2]+truckDIArray[18][3]+truckDIArray[18][4]+truckDIArray[18][5]+truckDIArray[18][6]); //number of those infected with HCV upon entry
         fprintf(fw, "%d\t", truckDIArray[19][0]+truckDIArray[19][1]+truckDIArray[19][2]+truckDIArray[19][3]+truckDIArray[19][4]+truckDIArray[19][5]+truckDIArray[19][6]); //number of those infected with HCV upon entry
        fprintf(fw, "%d\t", totalTestPS);
        fprintf(fw, "%d\t", totalTestHS);
        fprintf(fw, "%d\t", totalPSTP);
        fprintf(fw, "%d\t", totalHSTP);
        fprintf(fw, "%d\t", totalPSTN);
        fprintf(fw, "%d\t", totalHSTN);
        fprintf(fw, "%d\t", totalPSFP);
        fprintf(fw, "%d\t", totalHSFP);
        fprintf(fw, "%d\t", totalPSFN);
        fprintf(fw, "%d\t", totalHSFN);
        
        for(mcCtr1=0; mcCtr1<27; mcCtr1++){ //correctional centres
            for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
                for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                    for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                        fprintf(fw, "%d\t", MinCellArray[mcCtr1][mcCtr2][mcCtr3][mcCtr4]);
                    }
                }
            }
        }
    
        for(mcCtr1=0; mcCtr1<11; mcCtr1++){ //correctional centres
            for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
                for(mcCtr3=0; mcCtr3<4; mcCtr3++){ //min med max
                    for(mcCtr4=0; mcCtr4<19; mcCtr4++){ //min med max
                        fprintf(fw, "%d\t", MedCellArray[mcCtr1][mcCtr2][mcCtr3][mcCtr4]);
                    }
                }
            }
        }

        for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
            for(mcCtr2=0; mcCtr2<4; mcCtr2++){ //min med max
                for(mcCtr3=0; mcCtr3<5; mcCtr3++){ //min med max
                    for(mcCtr4=0; mcCtr4<20; mcCtr4++){ //min med max
                        fprintf(fw, "%d\t", MaxCellArray[mcCtr1][mcCtr2][mcCtr3][mcCtr4]);
                    }
                }
            }
        }
        
        fprintf(fw, "%d\t", overcap);
        fprintf(fw, "%d\t", infectedRedZone);
        
        for(mcCtr1=0; mcCtr1<27; mcCtr1++){ //correctional centres
            //for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
              //  for(mcCtr3=0; mcCtr3<6; mcCtr3++){ //min med max
                //    for(mcCtr4=0; mcCtr4<13; mcCtr4++){ //min med max
                        fprintf(fw, "%d\t", lockdownMinZone[mcCtr1]);
                  //  }
                //}
            //}
        }

        for(mcCtr1=0; mcCtr1<11; mcCtr1++){ //correctional centres
            //for(mcCtr2=0; mcCtr2<2; mcCtr2++){ //min med max
               //for(mcCtr3=0; mcCtr3<4; mcCtr3++){ //min med max
                    //for(mcCtr4=0; mcCtr4<19; mcCtr4++){ //min med max
                        fprintf(fw, "%d\t", lockdownMedZone[mcCtr1]);
                    //}
                //}
            //}
        }

        for(mcCtr1=0; mcCtr1<18; mcCtr1++){ //correctional centres
            //for(mcCtr2=0; mcCtr2<4; mcCtr2++){ //min med max
                //for(mcCtr3=0; mcCtr3<5; mcCtr3++){ //min med max
                    //for(mcCtr4=0; mcCtr4<20; mcCtr4++){ //min med max
                        fprintf(fw, "%d\t", lockdownMaxZone[mcCtr1]);
                    //}
                //}
            //}
        }
        
        fprintf(fw, "\n");
        
        printf("DAYS:                %d\n", currDay);
        
        
        //refresh population age group breakdown
        for(agCtrR=0; agCtrR<ageGroups; agCtrR++){ //age groups
            for(agCtrC=0; agCtrC<3; agCtrC++){ //populations
                iAgeArray[agCtrR][agCtrC]=0;
                psAgeArray[agCtrR][agCtrC]=0;
                hsAgeArray[agCtrR][agCtrC]=0;
                evAgeArray[agCtrR][agCtrC]=0;
                fvAgeArray[agCtrR][agCtrC]=0;
            }
        }
        
        //refresh zone count
        for(zaCtrR=0; zaCtrR<30; zaCtrR++){ //age groups
            for(zaCtrC=0; zaCtrC<3; zaCtrC++){ //populations
                zoneArray[zaCtrR][zaCtrC]=0;
                zoneCovidArray[zaCtrR][zaCtrC]=0; //cumulative count
            }
        }
        
        //refresh age severity
        for(asCtrR=0; asCtrR<ageGroups; asCtrR++){ //age groups
            for(asCtrC=0; asCtrC<ctrSeverity; asCtrC++){ //populations
                iAgeSeverityArray[asCtrR][asCtrC]=0;
            }
        }
        
        for(tInfCtrR=0; tInfCtrR<20; tInfCtrR++){ //correctional centres
            nMovingPop[tInfCtrR]=0;
        }
        
        for(zInfCtrR=0; zInfCtrR<30; zInfCtrR++){ //correctional centres
            for(zInfCtrC=0; zInfCtrC<7; zInfCtrC++){ //min med max
                zoneMinDIArray[zInfCtrR][zInfCtrC]=0; //daily new
                zoneMedDIArray[zInfCtrR][zInfCtrC]=0;
                zoneMaxDIArray[zInfCtrR][zInfCtrC]=0;
                zoneMinDDArray[zInfCtrR][zInfCtrC]=0; //daily deaths
                zoneMedDDArray[zInfCtrR][zInfCtrC]=0;
                zoneMaxDDArray[zInfCtrR][zInfCtrC]=0;
            }
        }

        for(cInfCtrR=0; cInfCtrR<38; cInfCtrR++){ //refresh court daily counts
            for(cInfCtrC=0; cInfCtrC<7; cInfCtrC++){ //min med max
                courtDIArray[cInfCtrR][cInfCtrC]=0;
                courtDDArray[cInfCtrR][cInfCtrC]=0;
            }
        }
        
        for(tInfCtrR=0; tInfCtrR<20; tInfCtrR++){ //refresh trucks daily infections
            for(tInfCtrC=0; tInfCtrC<7; tInfCtrC++){
                truckDIArray[tInfCtrR][tInfCtrC]=0;
            }
        }
        
        //totalAge0=0, totalAge1=0, totalAge2=0, totalAge3=0, totalAge4=0, totalAge5=0, totalAge6=0;
        //totalAD0=0, totalAD1=0, totalAD2=0, totalAD3=0, totalAD4=0, totalAD5=0, totalAD6=0;
        newCases=0; infectedOutside=0;
        infectedRedZone=0,
        releasedI=0; releasedICov=0; releasedPS=0; releasedPSCov=0; releasedHS=0; releasedHSCov=0; releasedEV=0; releasedEVCov=0; releasedFV=0; releasedFVCov=0;
        nMinHospitalF=0; nMedHospitalF=0; nMaxHospitalF=0; nMinHospitalC=0; nMedHospitalC=0; nMaxHospitalC=0;
        nInmatesTested=0, nPStaffTested=0, nHStaffTested=0, nEVisitorsTested=0, nFVisitorsTested=0;
        nInmatesDetected=0;  nPStaffDetected=0; nHStaffDetected=0; nEVisitorsDetected=0; nFVisitorsDetected=0;
        nInmatesThermal=0; nPStaffThermal=0; nHStaffThermal=0; nEVisitorsThermal=0; nFVisitorsThermal=0;  nInmatesRapidTested=0; nPStaffRapidTested=0;  nHStaffRapidTested=0; nEVisitorsRapidTested=0; nFVisitorsRapidTested=0;
        newInfectedInmates=0; newInfectedPS=0; newInfectedHS=0; newInfectedEV=0; newInfectedFV=0;
        newInfectedInmatesMin=0; newInfectedPSMin=0; newInfectedHSMin=0; newInfectedEVMin=0; newInfectedFVMin=0;
        newInfectedInmatesMed=0; newInfectedPSMed=0; newInfectedHSMed=0; newInfectedEVMed=0; newInfectedFVMed=0;
        newInfectedInmatesMax=0; newInfectedPSMax=0; newInfectedHSMax=0; newInfectedEVMax=0; newInfectedFVMax=0;
        prisonHCV=0; communityHCV=0; prisonOpd=0; prisonOpdNotOST=0;
        released=0; releasedHCVAb=0; releasedHCVRNA=0; releasedHCVCom=0; releasedHCVPris=0; //set counter back to zero
        totalClrN=0;
        deadHCV=0; deadHCVCom=0; deadHCVPris=0;
        dead=0; totalPrisonPop=0; nReinfected=0;
        E0=0; E1=0; E2=0; E3=0; E4=0; E5=0; E6=0; E7=0; E8=0; E9=0; E10=0; E11=0; E12=0; E13=0; E14=0; E15=0; E16=0; E17=0; E18=0; E19=0; E20=0; E21=0; E22=0; E23=0; E24=0; E25=0; E26=0;
        R0=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; R11=0; R12=0;
        //printf("Total : %d \n", locArray[1][0]+locArray[2][0]+locArray[3][0]);
        
        overcap=0;
        
        totalTestPS=0; totalTestHS=0; totalPSTP=0; totalHSTP=0; totalPSTN=0; totalHSTN=0; totalPSFP=0; totalHSFP=0; totalPSFN=0; totalHSFN=0;
        
        if(nTail!=NULL){
            traverse(&nHead, &nTail);
        }
        
        currDay++;
        //printf("\n");
    }

	printf("DONE!\n");
    fclose(fw);
    gsl_rng_free(r);
    return 0;
}

