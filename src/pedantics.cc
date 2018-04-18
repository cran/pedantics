// File: RPEDERR.cpp

// Version: 1.0

// Code last modified: 6 September 2007

// Created by: Michael Morrissey

// Part of software package PEDANTICS



#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <time.h>
#include <math.h>      
#include "R.h" 
using namespace std;



///////////////////////////////////////////
//////  SECTION 1: VARIABLES  /////////////
///////////////////////////////////////////


// definition of characters to identify missing pedigree data in input and
// output files:
string unknownParentageIndicator = "?";
string unknownParentageOutput = "*";

int performLHCheck = 1;

// if this indicator is set to 1, sex will be ignored when parents are selected:
int monoecy = 0;

// names of input and output files:
char inputFile[80] = "none";
char pedigreeFile[80] = "none";
char mainFile[80] = "none";

int pedigreeFileType = 0;

// controls whether or not the pedigree file is reordered by cohort.  Some
// software for quantitative genetic analyses may require this in order to run:
int reorderPedFile = 0;


class fish {
// an object to contain all individual data (only superficially a fish):
  public:
    string ID;
    int cohort;  // or birth year
    int sampled;  // if 0, then individual is a dummy individual that only
                  // exists so that some other individuals can have unsampled
                  // parents
    int sex;  // female=1, male=0
    int founder;  // founder=1, nonfounder=0

    // alphanumeric identities through which the investigator recognizes
    // individuals:
    string assignedMotherID; // these two are used by rpederr
    string assignedFatherID;

    string realMotherID; // these two are used by fpederr
    string realFatherID;

    // numerical identities for the parents in the pedigree with which
    // quantitative genetic analyses are performed, and with which variance
    // components are simulated:
    int assignedMother;
    int assignedFather;
    int realMother;
    int realFather;

    // the probabilities of erroneous parentage assignment:
    double probMotherError;
    double probFatherError;

    // when the reverse approach is applied, these are the probabilities that
    // the true parent is an unsampled individual, given the occurrence of an
    // erroneous parentage assignment.  When the forward approach is applied,
    // these are the probabilities parentage assignments are made.
    double probMotherSampled;  // these two are used by rpederr
    double probFatherSampled;
    double probMotherAssigned;  // these two are used by fpederr
    double probFatherAssigned;

    // the first and last cohorts to which an indivdual could potentially have
    // contributed offspring:
    int yearMat;
    int yearDeath;

    // indicators to control whether or not the assigned parents are to be
    // printed to the output files.  This simplifies some code, because
    // assigned parents can simply be assigned to all individuals, but not
    // necessarily all reported.
    bool asMother;
    bool asFather;

    // when the fish is actually a bird:
    int brood;
    int northing; int easting;
    int EPPcandidate;       //will be 0 for no and 1 for yes
    int broodDadAssigned;
    int broodMumAssigned;
};

// array to contain all individuals in pedigree, with a counter to track
// the actual number of individuals:
const int maxPedigreeSize = 150000;
fish pedigree[maxPedigreeSize];
int pedigreeSize = 0;
int males = 0;
int females = 0;
int unknownSexes = 0;
int founders = 0;
int minCohort = 99999;
int maxCohort = -99999;
int assignedMothers = 0;
int assignedFathers = 0;
int unsampled = 0;
int maxLifeSpan = 0;
int maxReproLifeSpan = 0;

///////////// some bird-specific variables:
double propEPPbroods = 0;
double propEPPchicksGivenEPPbrood = 0;
double propEPPbroodsTwoFathers = 0;
double EPPlambda = 0;
const int maxBroods = 10000;
const int maxEPPcandidates = 500;        // to sire a given brood
const int maxEPPsires = 10000;           // study-wide
const int maxBroodSize = 20;
int broodAffinities[maxBroods][maxBroodSize];
int potentialEPPbroods[maxBroods];
double broodNorthing[maxBroods];
double broodEasting[maxBroods];
int broodEPPfatherNums[maxBroods][maxEPPcandidates];   // numeric locations if vector 'pedigree', not alphanumeric ID
double broodEPPfatherProbs[maxBroods][maxEPPcandidates];
int broods;
int broodSize[maxBroods];
int broodYears[maxBroods];
double EPPsireYearLocs[5][maxEPPsires];
int numberCandidateEPPfathersForBrood[maxBroods];
int EPPcandidates;
double EPPbeta;
double EPPgamma;




/////////////////////////////////////////////////////////////////////
/////////  SECTION 2: FUNCTIONS COMMON TO ALL FUNCTIONS  ////////////
/////////////////////////////////////////////////////////////////////

double rndDouble() {
// Returns a random number between 0 and 1.  Not used directly, but
// rather supplies random numbers for all of the other rnd functions:

//    float a = (float) rand()/RAND_MAX;  return a;
   float a = (float) unif_rand();  return a;
}

bool rndError(double probError) {
// returns TRUE with probability probError:

    bool error = false;  double d = rndDouble();
    if(d <= probError) error = true;  return error;
}


int rndParent(int sex, int sampled, int reproductiveCohort) {
// returns the identity of an individual of a specified sex, that either was
// or was not sampled (as indicated) and that was able to contribute to a
// specified cohort.  If an appropriate parent cannot be found in 100,000
// attempts, all parents are systematically scanned to make sure that at
// least one exists.  If none exists, this method returns a value of -99
// for the parental identity.  This value acts as an error message for the
// calling method:

    bool parentFound = false;  bool parentAvailable = false;
    int p;  int counter = 0;

    if(monoecy == 0) {
        while(parentFound==false) {
            counter++;
            p = (int) (rndDouble()*pedigreeSize);
            if(pedigree[p].sex==sex
               && pedigree[p].sampled==sampled
               && pedigree[p].yearMat <= reproductiveCohort
               && pedigree[p].yearDeath >= reproductiveCohort)
                    parentFound = true;
            if(counter==100000 && parentAvailable==false){
                for(int i = 0; i < pedigreeSize; i++){
                    if(pedigree[p].sex==sex
                        && pedigree[p].sampled==sampled
                        && pedigree[p].yearMat <= reproductiveCohort
                        && pedigree[p].yearDeath >= reproductiveCohort)
                            parentAvailable = true;
                 }
            }
            if(parentAvailable==false && counter > 100000){
                parentFound = true;
                p = -99;
            }
        }
    }else{
        while(parentFound==false) {
            counter++;
            p = (int) (rndDouble()*pedigreeSize);
            if(pedigree[p].sampled==sampled
               && pedigree[p].yearMat <= reproductiveCohort
               && pedigree[p].yearDeath >= reproductiveCohort)
                    parentFound = true;
            if(counter==100000 && parentAvailable==false){
                for(int i = 0; i < pedigreeSize; i++){
                    if(pedigree[p].sampled==sampled
                        && pedigree[p].yearMat <= reproductiveCohort
                        && pedigree[p].yearDeath >= reproductiveCohort)
                            parentAvailable = true;
                 }
            }
            if(parentAvailable==false && counter > 100000){
                parentFound = true;
                p = -99;
            }
        }
    }
    return p;
}


string IDconvert(int ID){
  string s;
  stringstream sconv;
  sconv << ID;
  s = sconv.str();
  if(ID<0) s = unknownParentageIndicator;
  return s;
}



//////////////////////////////////////////////////////////////////
///////  SECTION 3: FUNCTIONS NOT USED BY ALL MODULES  ///////////
//////////////////////////////////////////////////////////////////

bool findAssignedParents(fish& fishy) {
// MODULE: rpederr...rpederr_bird can probably use this too

// finds the assigned parents of an individual in the assumed pedigree.
// If no parent is assigned in the assumed pedigree, one is chosen from
// those available for the individual's cohort.  Returns TRUE if parents
// were successfully located or simulated:

    bool foundMother = false;  bool foundFather = false;

    for(int p = 0; p < pedigreeSize; p++) {
        if(pedigree[p].ID == fishy.assignedMotherID) {
                fishy.assignedMother = p;
                foundMother = true;
        }
        if(pedigree[p].ID == fishy.assignedFatherID) {
                fishy.assignedFather = p;
                foundFather = true;
        }
    }

    if(fishy.assignedMotherID == unknownParentageIndicator) {
        bool motherSampled = rndError(fishy.probMotherSampled);
        int sam = 0; if(motherSampled==true) sam = 1;
        fishy.assignedMother = rndParent(1, sam, fishy.cohort);
        foundMother = true;
    }

    if(fishy.assignedFatherID == unknownParentageIndicator) {
        bool fatherSampled = rndError(fishy.probFatherSampled);
        int sam = 0; if(fatherSampled==true) sam = 1;
        fishy.assignedFather = rndParent(0, sam, fishy.cohort);
        foundFather = true;
    }

    bool foundParents = false;
    if(foundMother==true && foundFather==true
            && fishy.assignedMother!=-99 && fishy.assignedFather!=-99)
                 foundParents = true;
    return foundParents;
}



bool assignTrueParents(fish& fishy) {
// MODULE: rpederr

// This method first assigns an individual's parents to have the same identity
// their assumed parents.  Then, for mothers and fathers separately, errors
// are generated probabilistically.  If, for the mother and/or the father,
// the method rndParent was unable to find an appropriate parent, it will have
// returned a value of -99.  This will trigger the method to return FALSE.

    fishy.realMother = fishy.assignedMother;
    fishy.realFather = fishy.assignedFather;

    if( (rndError(fishy.probMotherError) == true)
               & (fishy.assignedMotherID != unknownParentageIndicator) ) {
        bool motherSampled = rndError(fishy.probMotherSampled);
        int sam = 0; if(motherSampled==true) sam = 1;
        fishy.realMother = rndParent(1, sam, fishy.cohort);
    }
    if( (rndError(fishy.probFatherError) == true)
               & (fishy.assignedFatherID != unknownParentageIndicator) ) {
        bool fatherSampled = rndError(fishy.probFatherSampled);
        int sam = 0; if(fatherSampled==true) sam = 1;
        fishy.realFather = rndParent(0, sam, fishy.cohort);
    }
    if(fishy.assignedMotherID!=unknownParentageIndicator){
        fishy.asMother = true;}else{fishy.asMother = false;}
    if(fishy.assignedFatherID!=unknownParentageIndicator){
        fishy.asFather = true;}else{fishy.asFather = false;}

    bool foundParents = true;
    if( (fishy.realMother==-99) | (fishy.realFather==-99) ) foundParents = false;
    return foundParents;
}


void copyDataFromRrpederr(int *ID, int *assumedMothers, int *assumedFathers, int *founder, int *sex, 
                 int *sampled, double *fatherError, double *motherError, double *fatherSamp, 
                 double *motherSamp, int *cohort, int *yearMat, int *yearDeath){
// MODULE: rpederr

  string newstring;
  stringstream stringconvert;

  for(int f=0; f<pedigreeSize; f++){
    pedigree[f].ID = IDconvert(ID[f]);
    pedigree[f].assignedFatherID = IDconvert(assumedFathers[f]);
    pedigree[f].assignedMotherID = IDconvert(assumedMothers[f]);
    pedigree[f].founder = founder[f];
    pedigree[f].sex = sex[f];
    pedigree[f].sampled = sampled[f];
    pedigree[f].probFatherError = fatherError[f];
    pedigree[f].probMotherError = motherError[f];
    pedigree[f].probFatherSampled = fatherSamp[f];
    pedigree[f].probMotherSampled = motherSamp[f];
    pedigree[f].cohort = cohort[f];
    pedigree[f].yearMat = yearMat[f];
    pedigree[f].yearDeath = yearDeath[f];
  }
}

void fillOutputArraysrpederr(int *realMothers, int *realFathers, int *completeAssumedMothers, int *completeAssumedFathers){
// MODULE: rpederr
  for(int f=0; f<pedigreeSize; f++){
    realMothers[f] = atoi(pedigree[pedigree[f].realMother].ID.c_str());
    if(pedigree[f].founder==1) realMothers[f] = -1;
    realFathers[f] = atoi(pedigree[pedigree[f].realFather].ID.c_str());
    if(pedigree[f].founder==1) realFathers[f] = -1;
    completeAssumedMothers[f] = atoi(pedigree[pedigree[f].assignedMother].ID.c_str());
    if(pedigree[f].founder==1) completeAssumedMothers[f] = -1;
    completeAssumedFathers[f] = atoi(pedigree[pedigree[f].assignedFather].ID.c_str());
    if(pedigree[f].founder==1) completeAssumedFathers[f] = -1;
  }
}

bool findRealParents(fish& fishy) {
// MODULE: fpederr
// finds the assigned parents of an individual in the assumed pedigree.
// If no parent is assigned in the assumed pedigree, one is chosen from
// those available for the individual's cohort.  Returns TRUE if parents
// were successfully located or simulated:

    bool foundMother = false;  bool foundFather = false;

    for(int p = 0; p < pedigreeSize; p++) {
        if(pedigree[p].ID == fishy.realMotherID) {
                fishy.realMother = p;
                foundMother = true;
        }
        if(pedigree[p].ID == fishy.realFatherID) {
                fishy.realFather = p;
                foundFather = true;
        }
    }

    bool foundParents = false;
    if(foundMother==true && foundFather==true)
                 foundParents = true;
    return foundParents;
}

bool assignAssignedParents(fish& fishy) {
// MODULE: fpederr
// This method first assigns an individual's parents to have the same identity
// their assumed parents.  Then, for mothers and fathers separately, errors
// are generated probabilistically.  If, for the mother and/or the father,
// the method rndParent was unable to find an appropriate parent, it will have
// returned a value of -99.  This will trigger the method to return FALSE.

    fishy.assignedMother = fishy.realMother;
    fishy.assignedFather = fishy.realFather;

    fishy.asMother = rndError(fishy.probMotherAssigned);
    fishy.asFather = rndError(fishy.probFatherAssigned);

    if(rndError(fishy.probMotherError) == true)
        fishy.assignedMother = rndParent(1, 1, fishy.cohort);

    if(rndError(fishy.probFatherError) == true)
        fishy.assignedFather = rndParent(0, 1, fishy.cohort);

    bool foundParents = true;
    if( (fishy.assignedMother==-99) | (fishy.assignedFather==-99) ) foundParents=false;
    return foundParents;
}


void copyDataFromRfpederr(int *ID, int *realMothers, int *realFathers, int *founder, int *sex, 
                 int *sampled, double *fatherError, double *motherError, double *fatherAssi, double *motherAssi, int *cohort, 
                 int *yearMat, int *yearDeath){
// MODULE: fpederr

  string newstring;
  stringstream stringconvert;

  for(int f=0; f<pedigreeSize; f++){
    pedigree[f].ID = IDconvert(ID[f]);
    pedigree[f].realFatherID = IDconvert(realFathers[f]);
    pedigree[f].realMotherID = IDconvert(realMothers[f]);
    pedigree[f].founder = founder[f];
    pedigree[f].sex = sex[f];
    pedigree[f].sampled = sampled[f];
    pedigree[f].probFatherError = fatherError[f];
    pedigree[f].probMotherError = motherError[f];
    pedigree[f].probFatherAssigned = fatherAssi[f];
    pedigree[f].probMotherAssigned = motherAssi[f];
    pedigree[f].cohort = cohort[f];
    pedigree[f].yearMat = yearMat[f];
    pedigree[f].yearDeath = yearDeath[f];
  }
}

void fillOutputArraysfpederr(int *assumedMothers, int *assumedFathers, int *completeAssumedMothers, int *completeAssumedFathers){
// MODULE: fpederr
  for(int f=0; f<pedigreeSize; f++){
    assumedMothers[f] = atoi(pedigree[pedigree[f].assignedMother].ID.c_str());
    if(pedigree[f].founder==1) assumedMothers[f] = -1;
    if(pedigree[f].asMother==false) assumedMothers[f] = -1;
    assumedFathers[f] = atoi(pedigree[pedigree[f].assignedFather].ID.c_str());
    if(pedigree[f].founder==1) assumedFathers[f] = -1;
    if(pedigree[f].asFather==false) assumedFathers[f] = -1;
    completeAssumedMothers[f] = atoi(pedigree[pedigree[f].assignedMother].ID.c_str());
    if(pedigree[f].founder==1) completeAssumedMothers[f] = -1;
    completeAssumedFathers[f] = atoi(pedigree[pedigree[f].assignedFather].ID.c_str());
    if(pedigree[f].founder==1) completeAssumedFathers[f] = -1;
  }
}



extern "C" {


  void RPEDERR_R(int *ID, int *assumedMothers, int *assumedFathers, int *founder, int *sex, 
                 int *sampled, double *fatherError, double *motherError, double *fatherSamp, 
                 double *motherSamp, int *cohort, int *yearMat, int *yearDeath, int *pedSize, 
                 int *realMothers, int *realFathers, int *completeAssumedMothers, 
                 int *completeAssumedFathers, int *monoecyTRUE){

    pedigreeSize = *pedSize;
    monoecy = *monoecyTRUE;

    copyDataFromRrpederr(ID, assumedMothers, assumedFathers, founder, sex, 
                 sampled, fatherError, motherError, fatherSamp, motherSamp, cohort, 
                 yearMat, yearDeath);

    Rprintf("\nSorting assigned parents...");
    for(int f = 0; f < pedigreeSize; f++) {
         if(pedigree[f].founder == 0) {
              bool found = findAssignedParents(pedigree[f]);
              if(found==false)
                 Rprintf("\nParent(s) not found for individual at pedigree position %i\n", f);
         }
    }

    Rprintf("Done.\n");

    Rprintf("\nAssigning real parents...");

    for(int f = 0; f < pedigreeSize; f++) {
        if(pedigree[f].founder == 0) {
                bool found = assignTrueParents(pedigree[f]);
                if(found==false)
                 Rprintf("\nTrue parent could not be simulated for individual at pedigree position %i\n", f);
        }
    }

    Rprintf("Done.\n");

    fillOutputArraysrpederr(realMothers, realFathers, completeAssumedMothers, completeAssumedFathers);
  }

  void FPEDERR_R(int *ID, int *realMothers, int *realFathers, int *founder, int *sex, 
                 int *sampled, double *fatherError, double *motherError, double *fatherAssi, double *motherAssi, int *cohort, 
                 int *yearMat, int *yearDeath, int *pedSize, int *assumedMothers, 
                 int *assumedFathers, int *completeAssumedMothers, int *completeAssumedFathers, int *monoecyTRUE){

        pedigreeSize = *pedSize;
        monoecy = *monoecyTRUE;

        copyDataFromRfpederr(ID, realMothers, realFathers, founder, sex, 
                 sampled, fatherError, motherError, fatherAssi, motherAssi, cohort, 
                 yearMat, yearDeath);

    Rprintf("\nSorting real parents...");

    for(int f = 0; f < pedigreeSize; f++) {
         if(pedigree[f].founder == 0) {
              bool found = findRealParents(pedigree[f]);
              if(found==false)
                 Rprintf("\nParent(s) not found for individual at pedigree position %i\n", f);
         }
    }

    Rprintf("Done.\n");

    Rprintf("\nAssigning parents for assumed pedigree...");

    for(int f = 0; f < pedigreeSize; f++) {
        if(pedigree[f].founder == 0) {
                bool found = assignAssignedParents(pedigree[f]);
                if(found==false)
                 Rprintf("\nTrue parent could not be simulated for individual at pedigree position %i\n", f);
        }
    }

    Rprintf("Done.\n");


    fillOutputArraysfpederr(assumedMothers, assumedFathers, completeAssumedMothers, completeAssumedFathers);

  }


  void CALC_SIB_NUMBERS(int *n, int *id, int *mum, int *dad, int *full, int *maternal, int *paternal){
    int f=0; int m=0; int p=0;
    for(int x = 0; x<(*n-1); x++){
      for(int y = (x+1); y<*n; y++){
        if((mum[x]!=-1)&(dad[x]!=-1)&(mum[y]!=-1)&(dad[y]!=-1)&(mum[x]==mum[y])&(dad[x]==dad[y])) f++;
        if((mum[x]!=-1)&(mum[y]!=-1)&(mum[x]==mum[y])) m++;
        if((dad[x]!=-1)&(dad[y]!=-1)&(dad[x]==dad[y])) p++;
      }
    }
    *full=f;  *maternal=m;  *paternal=p;
  }
}
