//
//  individual.h
//  HCV Model
//
//  Created by Neil Bretana on 1/10/2014.
//  Copyright (c) 2014 Neil Bretana. All rights reserved.
//
#ifndef HCV_Model_individual_h
#define HCV_Model_individual_h

struct individual
{
    int ID;
    int indivType; //0: inmate, 1: prison staff, 2: healthcare staff, 3: essential visitor, 4: family visitor
    int gender;
    int age; //18-100
    int atsi; //0-no; 1-yes
    int comorbidity; //0-no; 1-yes
    int vulnerable; //0-no; 1-yes depending on ATSi, age, and co-mornidity
    
    int location; //0 community, 1 min, 2 med, 3 max
    int locType; //0 normal prison, 1 reception, 2 hospital, 3etc.
    int locArea; //0 metropolitan, 1 regional (more for prison staff)
    int locPrison; //set up data dictionary for each specific prison
    int locUnit;
    int locPod;
    int locCell;
    
    int court; // 0 or 1?
    int courtNumber; //which specific court?
    int courtCell;
    int moving; //0 or 1;
    ////LAYER: LABEL prison wings, pods, cells
    int truckNumber;
    
    //testing
    //int tested;
    //int testResultTrue; //0 if false positive, 1 if true positive
    //int testDay;
    
    //vaccinated
    double vaccineProtection;
    double vaccineTransmissionReduction;
    int vaxDose1;
    int vaxDose1Day;
    int vaxDose2;
    int vaxDose2Day;
    
    //strategy: cohorting
    int cohorting; //0-no 1-yes
    int cohortNumber;
    int cohortMaxDays;
    int timeStartCohort; //first day of cohorting

    int timeOfImprisonment; //Record time of imprisonment
    
    int COVIDAb; //0-no, 1-yes //clearer
    int COVID; //0-no, 1-yes
    //int symptoms; //0-no, 1-yes
    int infectionNumber; //Tells the number of times an individual has been infected
    int severity;
    ////add 0:exposed, 1:asymptomatic -> 2:mild -> 3:mod -> 4:severe, 5:clearer, 6: no infection
    //add 0:exposed -> 1:pre-clinical -> 2:asymptomatic -> 3:mild -> 4:mod -> 5:severe, 6:clearer, 7: no infection
    int timeOfProgression;
    int timeOfInfection; //Record timepoint when it got infected
    int placeInfected; //0 (N/A), 1 (min), 2 (med), 3 (max), 4 (community)
    
    int tested; // 0 or 1
    int testDaysBeforeResult; // days befor result
    int testResultTrue; //0 if false positive, 1 if true positive
    int testDate; //date
    int testAccuracy;
    int TP;
    int FP;
    int FN;
    int TN;
    
    int PCRtested; //0 not yet, 1 done
    int PCRdate;
    int PCRresult;
    
    int detected; //0 or 1 //was covid detected during recent test?
    int ppe; //0 or 1 standard mask, 2 N95
    double ppeprotection;
    
    int hospitalised; //0-no, 1-yes field hospital, 2 community hospital
    
    int zone; // prison 0-30 depending on security setting (location)
    
    int group; //0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 This is for the population groups
    //COVID-19
    // Co-Morbidity+; COVID+; ATSI 0
    // Co-Morbidity+; COVID+; Non-ATSI 1
    // Co-Morbidity+; COVID-; ATSI; Antibody+ 2
    // Co-Morbidity+; COVID-; ATSI; Antibody- 3
    // Co-Morbidity+; COVID-; Non-ATSI; Antibody+ 4
    // Co-Morbidity+; COVID-; Non-ATSI; Antibody- 5
    // Co-Morbidity-; COVID+; ATSI; 6
    // Co-Morbidity-; COVID+; Non-ATSI; 7
    // Co-Morbidity-; COVID-; ATSI; Antibody+
    // Co-Morbidity-; COVID-; ATSI; Antibody- 9
    // Co-Morbidity-; COVID-; Non-ATSI; Antibody+
    // Co-Morbidity-; COVID-; Non-ATSI; Antibody- 11
    
    //HCV
    //[row][0] IDU+; HCV+; ATSI
    //[row][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE
    int risk; //0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 This is for the injecting risk group
    //Risk groups: 0: non-injecting;
    //1: injecting less than daily; opioid; not sharing
    //2: injecting less than daily; opioid; sharing less daily
    //3: injecting less than daily; opioid; sharing daily more
    //4: injecting less than daily; non-opioid; no sharing
    //5: injecting less than daily; non-opioid; sharing less daily
    //6: injecting less than daily; non-opioid; sharing daily or more
    //7: injecting daily or more; opioid; not sharing
    //8: injecting daily or more; opioid; sharing less daily
    //9: injecting daily or more; opioid; sharing daily more
    //10: injecting daily or more; non-opioid; no sharing
    //11: injecting daily or more; non-opioid; sharing less daily
    //12: injecting daily or more; non-opioid; sharing daily or more
    
    int isolated; //0 not isolated, 1 for 1 out isolation, 2 for 2 out isolation
    int timeIsolated;
    
    struct individual *prevIndiv;
    struct individual *nextIndiv; //nextIndiv
    
    //HCV dump characteristics NOT IN USE
    int injecting; //0-no; 1-yes;
    int sharing; //0-no; 1-yes;
    int injFreq; //0 - not applicable; 1 less than daily; 2 daily or more;
    int shaFreq; //0 - not applicable; 1 less than daily; 2 daily or more;
    int metavir; //0-4; 5-clearer; 6-non-infected//never-infected
    int cirrhosis; //0 and 1 link to metavir 4
    int viralLoad;
    int injOpd; //Injecting Opioids 0 1
    int her; //0 1
    int herDose;
    int herDaysMissed;
    int met; //0 1
    int metDose;
    int metDaysMissed;
    int bup;
    int bupDose;
    int bupDaysMissed;
    int OST; //0 no OST, 1 = old OST 2006 to 2018, 2 = optimised OST 2018 onwards
    int timeStartOST; //record time of OST start
    int bleach; //0 1
    int syringeEx; //0 1
    int DAA; //0 1
    int timeStartDAA; //record time of DAA start
    int everIDU; //0 1
};

typedef struct individual sIndiv;

//Methods

#endif
