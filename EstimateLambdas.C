//**************************************************
//
// Code to estimate lambda parameters using THERMINATOR.  
// This code runs over THERMINATOR events.  
// It counts all primary and secondary lambdas in each event, 
// taking into account some reconstruction cuts.
// It stores each event multiplicty count in a TH1D.
// From the TH1D, the code then computes lambda = <N_ij>/<Sum_ij N_ij>, 
// where N_ij is the number of pairs particles i and j
// The code also estimates the average multiplicity of each species of lambda.
//
// Run the code via .x EstimateLambdas()
//
//**************************************************


#include <cstdio>
#include <cassert>
#include <vector>

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"

bool globalDebug = false;
class ParticleCoor;

//Was having a problem including a therminator file.
//This fixes it.
#define _CXX_VER_ "g++(4.9.1)" 

enum Particle {kProt = 2212, kAntiProt = -2212, 
	       kPiPlus = 211, kPiMinus = -211, 
	       kLam = 3122, kAntiLam = -3122};
enum PairType {kLamLam, kALamALam, kLamALam};

void EstimateLambdas(const PairType pairType, int nFiles=2)
{
  // Main function. Generates list of file names
  // and pass them to RunSimpleLambdaEstimate
  if( (pairType < kLamLam) || (pairType > kLamALam) ) {
    cout<<"Bad input pair type\n";
    return;
  }
  TString filePath = "/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/";
  vector<TString> inputFileNames;

  cout<<"\nAnalyzing the following files:\n";
  for(int iFile = 0; iFile < nFiles; iFile++)
  {
    TString fileName = filePath;
    fileName += "event";
    if     (iFile <  10) fileName+= "00";
    else if(iFile < 100) fileName+= "0";
    fileName += iFile;
    fileName += ".root";
    cout<<fileName<<endl;
    inputFileNames.push_back(fileName);
  }

  // Load in the necessary therminator particle class
  gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
  gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  RunSimpleLambdaEstimate(inputFileNames, pairType);
}


void RunSimpleLambdaEstimate(/*TString inputFileName = "event000.root"*/ vector<TString> inputFileNames, const PairType pairType)
{
  // Outputs and saves histograms of 
  // lambda parameters and average particle multiplicities

  //Make vectors to hold particle yield histograms for each event
  vector<TH1D*> eventParticles1;
  vector<TH1D*> secondParticleCollection; //Only use this if we want to estimate particle-antiparticle lambda parameters
  vector<TH1D*> &eventParticles2 = eventParticles1;
  if(pairType == kLamALam) {
    &eventParticles2 = secondParticleCollection; //We do have separate particles and antiparticles, so use a second collection
  }
  int nParticleTypes = 5;
  GenerateYieldHistograms(nParticleTypes, inputFileNames, eventParticles1, eventParticles2, pairType);

  //Now that we have found all the strange particles, calculate lambda parameters for each pair type
  TH2D* hLambdaPars = GenerateLambdaParHisto(nParticleTypes, eventParticles1, eventParticles2, pairType);
  if(pairType == kLamLam) {
    SetLambdaHistAxisLabels(hLambdaPars->GetXaxis(),kLam);
    SetLambdaHistAxisLabels(hLambdaPars->GetYaxis(),kLam);
  }
  else if(pairType == kALamALam) {
    SetLambdaHistAxisLabels(hLambdaPars->GetXaxis(),kALam);
    SetLambdaHistAxisLabels(hLambdaPars->GetYaxis(),kALam);
  }
  else {
    SetLambdaHistAxisLabels(hLambdaPars->GetXaxis(),kLam);
    SetLambdaHistAxisLabels(hLambdaPars->GetYaxis(),AkLam);
  }

  //Draw and save the lambda par histo
  // TCanvas *c1 = new TCanvas("LambdaPars","Lambda Parameters");
  // hLambdaPars->DrawCopy("colzTEXTe");

  // cout<<"Debug:: All Finished"<<endl;
  TString outfileName = "LambdaPars.root";

  TFile outFile(outfileName,"update");
  outFile.cd();
  cout<<"Writing lambda histogram. Name:\t"<<hLambdaPars->GetName()<<endl;
  hLambdaPars->Write(0,TObject::kOverwrite);
 
  //Make, draw, and save a histo of average particle yields
  TH1D *hAvgYields = NULL;
  if( (pairType == kLamLam) || (pairType == kLamALam) ) {
    hAvgYields = ComputeAverageYields(eventParticles1,kLam);
    SetLambdaHistAxisLabels(hAvgYields->GetXaxis(),kLam);
  }
  else {
    hAvgYields = ComputeAverageYields(eventParticles1,kALam);
    SetLambdaHistAxisLabels(hAvgYields->GetXaxis(),kALam);
  }
  cout<<"Writing AvgYield1"<<endl;
  hAvgYields->Write(0,TObject::kOverwrite);

  if(pairType == kLamALam) {
    TH1D *hAvgYieldsAnti = ComputeAverageYields(eventParticles2,kALam);
    SetLambdaHistAxisLabels(hAvgYieldsAnti->GetXaxis(),kALam);
    cout<<"Writing AvgYield2"<<endl;
    hAvgYieldsAnti->Write(0,TObject::kOverwrite);
  }
  cout<<"All finished!"<<endl;
}

void UseDaughtersToFindParents(const TTree *thermTree, int &nextEntry, const PairType thisPairType, const int finalEntry, vector<int> &lambdaIDs, vector<int> &antiLambdaIDs)
{
  // Returns the particle IDs of all lambdas whose daughters
  // have been found.
  // Code loops over all the particles in an event.
  // Find all the protons and pions that are daughters of lambdas.
  // Check if each particle passes particle cuts.
  // If both proton and pion daughter pass cuts, then add
  // lambda particle ID to vector.

  vector<int> protParentIDs;
  vector<int> antiprotParentIDs;
  vector<int> piplusParentIDs;
  vector<int> piminusParentIDs;

  bool wantLams = false;
  bool wantALams = false;
  if( (thisPairType == 1) || (thisPairType == 3) ) wantLams == true;
  if( (thisPairType == 2) || (thisPairType == 3) ) wantALams == true;

  ParticleCoor *particleEntry = new ParticleCoor();
  TBranch *thermBranch = thermTree->GetBranch("particle");
  thermBranch->SetAddress(particleEntry);

  //Get first particle in this event
  int nBytesInEntry = thermTree->GetEntry(nextEntry);
  assert(nBytesInEntry > 0);
  int currentEventID = particleEntry->eventid;
  cout<<"Event ID\t"<<currentEventID<<endl;
  int iPart = nextEntry;

  //////////////////
  // Loop over particles.  Stop when we reach a new event
  while(particleEntry->eventID == currentEventID)
  {
    int pid = particleEntry->pid;
    int parentPid = particleEntry->fatherpid;
    int parentID = particleEntry->fathereid;
    if( (pid == kProt) && wantLams && (parentPid == kLam) ) {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	protParentIDs.push_back(parentID);
      }
    }
    else if( (pid == kAntiProt) && wantALams && (parentPid == kAntiLam) ) 
    {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	antiprotParentIDs.push_back(parentID);
      }
    }
    else if( (pid == kPiPlus) && wantALams && (parentPid == kAntiLam) )
    {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	antiprotParentIDs.push_back(parentID);
      }
    }
    else if( (pid == kPiMinus) && wantLams && (parentPid == kLam) ) {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	protParentIDs.push_back(parentID);
      }
    }

    iPart++;
    if(iPart == finalEntry) break;
    nBytesInEntry = thermTree->GetEntry(iPart);
    assert(nBytesInEntry > 0);
  } //end while loop

  ///////////////
  // Check for (anti)lambdas in both daughter lists

  if(wantLams) {
    for(int iProt; iProt < protParentIDs.size(); iProt++){
      const int v0ID = protParentIDs[iProt];
      for(int iPi; iPi < piminusParentIDs.size(); iPi++){
	if(v0ID == piminusParentIDs[iPi]) {
	  lambdaIDs.push_back(v0ID);
	  break;
	}
      }
    }
  }
  if(wantALams) {
    for(int iProt; iProt < antiprotParentIDs.size(); iProt++){
      const int v0ID = antiprotParentIDs[iProt];
      for(int iPi; iPi < piplusParentIDs.size(); iPi++){
	if(v0ID == piplusParentIDs[iPi]) {
	  antiLambdaIDs.push_back(v0ID);
	  break;
	}
      }
    }
  }

  //Make sure we didnt collect any of the wrong particle type
  if(!wantLams) assert(lambdaIDs.size() == 0);
  if(!wantALams) assert(antiLambdaIDs.size() == 0);


  //Set nextEntry to the particle ID of the first particle in the 
  //next event
  nextEntry = iPart;
  
  cout<<"Deleting particleEntry"<<endl;
  delete particleEntry;
}


void GenerateYieldHistograms(const int nParticleTypes, vector<TString> &inputFileNames, vector<TH1D*> &eventParticles1, vector<TH1D*> &eventParticles2, const PairType pairType)
{
  
  cout<<"Starting particle collection"<<endl;

  // //Load the ParticleCoor therminator class
  // gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
  // gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");

  int iEntry = 0;
  int eventCounter = 0;

  // Loop over each file in the collection
  for(int iFile = 0; iFile < inputFileNames.size(); iFile++)
  {
    // Read in the therminator event file and get the particle branch
    TFile thermTFile(inputFileNames[iFile], "READ");
    assert(NULL!=&thermTFile);
    TTree *thermTree = (TTree*)thermTFile.Get("particles");
    assert(NULL!=thermTree);


    int nextEventFirstEntry = 0;

    int nThermEntries = thermTree->GetEntries();


    while(nextEventFirstEntry < nThermEntries) {
      // For each event, find the daughters and use them to find the parents
      eventCounter++;
      cout<<"Processing event \t"<<eventCounter<<endl;
      vector<int> lambdaIDs;
      vector<int> antiLambdaIDs;
      UseDaughtersToFindParents(thermTree, nextEventFirstEntry, pairType, lambdaIDs, antiLambdaIDs);
    
      //Only push back histogram vectors if we find V0s.
      if(!(lamddaIDs.size() == 0)) {
	TH1D *lambdaYields = FillLambdaYieldHist(thermTree,lambdaIDs,eventCounter,kLam, nParticleTypes);
	if(lambdaYields && lambdaYields->GetEntries > 0) {
	  eventParticles1.push_back(lambdaYields);
	}
      }
      if(!(antiLamddaIDs.size() == 0)) {
	TH1D *antiLambdaYields = FillLambdaYieldHist(thermTree,antiLambdaIDs,eventCounter,kALam,nParticleTypes);
	if(antiLambdaYields && antiLambdaYields->GetEntries > 0) {
	  if(pairType == kLamALam) {
	    eventParticles2.push_back(antiLambdaYields);
	  }
	  else eventParticles1.push_back(antiLambdaYields);
	}
      }
      
    } // end event (while) loop
  } // end file loop

  cout<<"Finished particle collection"<<endl;
  cout<<"Total events used:\t"<<eventCounter<<endl;
}

TH1D *FillLambdaYieldHist(const TTree *thermTree, const vector<int> &v0IDs, const int eventCounter, const Particle part, const int nParticleTypes)
{
  cout<<"Filling yield histogram for event "<<eventCounter<<endl;
  int strangePDGs[5] = {3122, //Lambda
			3212, //Sigma0
			// 3224, //Sigma*+
			// 3214, //Sigma*0
			// 3114, //Sigma*-
			3322, //Xi0
			3312, //Xi-
			// 3324, //Xi*0
			// 3314, //Xi*-
			3334}; //Omega-

  // Build hist name
  TString histName = "h";
  if(part < 0) histName += "Anti";
  histName += "Particles";
  histName += eventCounter;

  // Create histogram
  TH1D *hParticles = new TH1D(histName,"Particles Per Type", nParticleTypes, 0, nParticleTypes);
  hParticles->SetDirectory(0);
  // hParticles->Sumw2();

  // Loop over each v0ID, check if particle passes cuts, and fill histo appropriately
  ParticleCoor *particleEntry = new ParticleCoor;
  TBranch *thermBranch = thermTree->GetBranch("particle");
  thermBranch->SetAddress(particleEntry);

  for(int iID = 0; iID < v0IDs.size(); iID++) 
  {

    int currentID = v0IDs[iID];
    int nBytesInEntry = thermTree->GetEntry(currentID);
    assert(nBytesInEntry > 0);
    int pid = particleEntry->pid;
    assert(pid == part);

    //Make sure the particle passes reconstruction cuts
    if(!CheckIfPassLambdaCuts(particleEntry)) continue;
    
    //Bin the particle according to its parent info
    for(int iPar = 0; iPar < nParticleTypes; iPar++)
    {
      const int fpid = particleEntry->fatherpid;
      if(fpid == strangePDGs[iPar]) {
	hParticles->Fill(iPar);
	break;
      }
      // If we reach the end of the last loop iteration,
      // that means that the lambda parent isn't in
      // our short list of PDG codes.  It must be a 
      // short-lived resonance.  We'll just bin that
      // as a primary (anti)lambda.
      if(nParticleTypes-1 == iPar) {
	hParticles->Fill(0);
      }
    } // End loop over parent type
  } // End loop over V0s
  return hParticles;
}



// void GenerateEventParticleHistograms(int nParticleTypes, vector<TString> &inputFileNames, vector<TH1D*> &eventParticles1, vector<TH1D*> &eventParticles2, bool isAntipart1, bool isAntipart2)
// {

//   cout<<"Starting particle collection"<<endl;

//   //Load this therminator class
//   gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
//   gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  
//   bool isPartAntipartPairs = false;
//   if(isAntipart1 != isAntipart2) isPartAntipartPairs = true;


//   //Make a list of the PID codes for each particle we care about
//   //We will treat all resonance -> lambda decays as primary lambdas,
//   //so we don't need to include those PDG codes here.

//   int strangePDGs[5] = {3122, //Lambda
// 			3212, //Sigma0
// 			// 3224, //Sigma*+
// 			// 3214, //Sigma*0
// 			// 3114, //Sigma*-
// 			3322, //Xi0
// 			3312, //Xi-
// 			// 3324, //Xi*0
// 			// 3314, //Xi*-
// 			3334}; //Omega-

  
//   ParticleCoor *particleEntry = new ParticleCoor();
//   int nFiles = inputFileNames.size();
//   int eventCounter = -1;
//   for(int iFile = 0; iFile < nFiles; iFile++)
//   {

//     //If debugging, do calculations using a single event
//     if(globalDebug) if(1 == eventCounter) break;

//     //Read in the therminator event file and get the particle branch
//     TFile thermTFile(inputFileNames[iFile], "READ");
//     assert(NULL!=&thermTFile);
//     TTree *thermTree = (TTree*)thermTFile.Get("particles");
//     assert(NULL!=thermTree);
//     TBranch *thermBranch = thermTree->GetBranch("particle");
//     thermBranch->SetAddress(particleEntry);

//     UInt_t eventID = 0; //This will be set to each event ID (they are not number 1, 2, ...)

//     TH1D *currentHist1 = NULL;
//     TH1D *currentHist2 = NULL;

//     //Loop over all the particles and find the lambdas
//     int nThermEntries = thermTree->GetEntries();
//     for(int i=0; i<nThermEntries; i++)
//     {
//       //Get the next particle
//       int nBytesInEntry = thermTree->GetEntry(i);
//       assert(nBytesInEntry>0);

//       //Check if this entry is a (anti)lambda (primary or secondary)
//       int pid = particleEntry->pid;
//       if( (!isAntipart1 && 3122==pid) || 
// 	  (isAntipart1 && -3122==pid) || 
// 	  (isPartAntipartPairs && -3122==pid) )
//       {
// 	//Make sure the particle passes the reconstruction cuts
// 	if(!CheckIfPassLambdaCuts(particleEntry)) continue;
// 	//Check if this is a new event.  If so, we'll need to make a new multiplicity histogram
// 	if( (particleEntry->eventid != eventID) || !currentHist1){ 

// 	  //We have reached a new event (or the first event)
// 	  eventID = particleEntry->eventid;
// 	  eventCounter++;
// 	  cout<<"Processing event \t"<<eventCounter<<"\tEvent ID: \t"<<particleEntry->eventid<<endl;

// 	  //If debugging, do calculations using a single event
// 	  if(globalDebug) if(1 == eventCounter) break;

// 	  //Make a new particle histogram and add it to the vector
// 	  TString histName = "h";
// 	  if(isAntipart1) histName += "Anti";
// 	  histName += "Particles";
// 	  histName += eventCounter;
// 	  TH1D *hParticles1 = new TH1D(histName,"Particles Per Type", nParticleTypes, 0, nParticleTypes);
// 	  hParticles1->SetDirectory(0); //Disassociate it from the TFile
// 	  hParticles1->Sumw2();
// 	  eventParticles1.push_back(hParticles1);
// 	  //Now use this histogram in future particle binning
// 	  currentHist1 = hParticles1;

// 	  //If we are doing particle antiparticle pairs, make a second
// 	  //set of hists for holding antiparticle yields.
// 	  if(isPartAntipartPairs) {
// 	    TString histName2 = "hAntiParticles";
// 	    histName2 += eventCounter;
// 	    TH1D *hParticles2 = new TH1D(histName2,"Particles Per Type", nParticleTypes, 0, nParticleTypes);
// 	    hParticles2->SetDirectory(0); //Disassociate it from the TFile
// 	    hParticles2->Sumw2();
// 	    eventParticles2.push_back(hParticles2);
// 	    currentHist2 = hParticles2;
// 	  }
// 	}

// 	//Check parent info
// 	for(int iPar = 0; iPar < nParticleTypes; iPar++)
// 	{
// 	  int fpid = particleEntry->fatherpid;
// 	  if( !isAntipart1 && fpid == strangePDGs[iPar] ) {
// 	    currentHist1->Fill(iPar);
// 	    break;
// 	  }
// 	  else if( isAntipart1 && fpid == (-1 * strangePDGs[iPar]) ) {
// 	    currentHist1->Fill(iPar);
// 	    break;
// 	  }
// 	  if(isPartAntipartPairs && fpid == (-1 * strangePDGs[iPar]) ) {
// 	    currentHist2->Fill(iPar);
// 	    break;
// 	  }

// 	  // If we reach the end of the last loop iteration,
// 	  // that means that the lambda parent isn't in
// 	  // our short list of PDG codes.  It must be a 
// 	  // short-lived resonance.  We'll just bin that
// 	  // as a primary (anti)lambda.
// 	  if(nParticleTypes-1 == iPar) 
// 	  {
// 	    if( !isAntipart1 && pid > 0 ) currentHist1->Fill(0);
// 	    else if( isAntipart1 && pid < 0 ) currentHist1->Fill(0);
// 	    else if( isPartAntipartPairs && pid < 0) currentHist2->Fill(0);
// 	  }
// 	} // End loop over iPar
//       } // end is(lambda)?
//     } // end loop over particle entries
//   } // end loop over files
//   cout<<"Finished particle collection"<<endl;
//   cout<<"Total events used:\t"<<eventCounter+1<<endl;
// }


TH2D *GenerateLambdaParHisto(int nParticleTypes, const vector<TH1D*> &eventParticles1, const vector<TH1D*> &eventParticles2, const PairType pairType){
  cout<<"Generating lambda parameter histograms."<<endl;
  // Make a histogram containing lambda parameters for each
  // pair type
  TString histName = "LambdaParams";
  if(pairType == kLamLam) histName += "Particle";
  else if(pairType == kALamALam) histName += "Antiparticle";
  else histName += "ParticleAntiparticle";
  TH2D *hLambdaPars = new TH2D(histName, "Lambda Parameters by Pair Type", nParticleTypes, 0, nParticleTypes, nParticleTypes, 0, nParticleTypes);
  hLambdaPars->SetDirectory(0);

  // Loop over each type of pairs, and set the lambda
  cout<<"Calculating lambda parameters"<<endl;
  double sumLambdaPars = 0.;
  for(int iPart1 = 0; iPart1 < nParticleTypes; iPart1++){//First type of particle in pair
    int iter2 = iPart1;
    if(pairType == kLamALam) iter2 = 0;
    for(int iPart2 = iter2; iPart2 < nParticleTypes; iPart2++){
      double lambdaPar = -1.;
      double lambdaParError = 1.;
      ComputeLambda(iPart1, iPart2, eventParticles1, eventParticles2, lambdaPar, lambdaParError);
      hLambdaPars->SetBinContent(iPart1+1, iPart2+1, lambdaPar);
      hLambdaPars->SetBinError(iPart1+1, iPart2+1, lambdaParError);
      sumLambdaPars += lambdaPar;
      // cout<<iPart1<<"\t"<<iPart2<<"\t"<<lambdaPar<<"\t"<<lambdaParError<<"\n"<<endl;      
    }
  }
  //The sum of lambda pars should add up to unity.
  cout<<"Sum of lambda pars:\t"<<sumLambdaPars<<endl;
  cout<<"Finished generating lambda parameter histograms."<<endl;
  return hLambdaPars;
}

void ComputeLambda(const int part1, const int part2, const vector<TH1D*> &eventParticles1, const vector<TH1D*> &eventParticles2, double &lambda, double &lambdaError)
{
  //Determine whether we are computing lambda for particle-particle (anti-anti) or particle-antiparticle, then calculate the corresponding lambda parameter and lambda error bars.
  bool identical = false;
  if(&eventParticles1 == &eventParticles2) ComputeLambdaIdentical(part1, part2, eventParticles1, lambda, lambdaError);
  else ComputeLambdaNonIdentical(part1, part2, eventParticles1, eventParticles2, lambda, lambdaError);
}

 
void ComputeLambdaIdentical(const int part1, const int part2, const vector<TH1D*> &eventParticles1, double &lambda, double &lambdaError)
{
  
  // For identical particles (all particle or all antiparticle).
  // Calculate lambda parameter for a given pair type via 
  // <N pairs>/<N total pairs>

  int nSpecPairs = 0;
  int nTotalPairs = 0;
  int nEvents = eventParticles1.size();
  int nEffectiveEvents = 0;

  for(int iEv = 0; iEv < nEvents; iEv++){ 
    double nTotalYield = eventParticles1[iEv]->Integral();
    // cout<<"Event "<<iEv<<":\tTotal yield\t"<<nTotalYield<<endl;
    if(1 >= nTotalYield) continue;
    nEffectiveEvents++;
    // cout<<"Getting yield of particle 1"<<endl;
    double nPart1 = eventParticles1[iEv]->GetBinContent(part1+1);
    // cout<<"Part1\t"<<nPart1<<endl;
    // Sum pairs for each event
    if(part1 == part2) { //pairs of identical (parent) particles
      nSpecPairs += nPart1*(nPart1 - 1)/2;

    }
    else { //pairs of non-identical (parent) particles
      // cout<<"Getting yield of particle 2"<<endl;
      double nPart2 = eventParticles1[iEv]->GetBinContent(part2+1);
      nSpecPairs += nPart1*nPart2;
      // cout<<"Part2\t"<<nPart2<<endl;
    }
    nTotalPairs += nTotalYield * (nTotalYield - 1)/2;
  } //end looping over events
    
  //Compute average pairs and average total pairs
  // cout<<"Calculating num"<<endl;
  double avgSpecPairs = (1.*nSpecPairs) / (1.*nEffectiveEvents);
  // cout<<"Calculating den"<<endl;
  double avgTotalPairs = (1.*nTotalPairs) / (1.*nEffectiveEvents);

  // Find lambda (ratio of avgSpec/avgTotal) and error
  cout<<"Calculating lambda now"<<endl;
  lambda = avgSpecPairs/avgTotalPairs;
  cout<<"Lambda is\t"<<lambda<<endl
      <<"Now calculate identical pair error"<<endl;
  lambdaError = ComputeLambdaIdenticalError(part1, part2, eventParticles1, avgTotalPairs, lambda);
  return;
}

void ComputeLambdaNonIdentical(const int part1, const int part2, const vector<TH1D*> &eventParticles1, const vector<TH1D*> &eventParticles2, double &lambda, double &lambdaError)
{
  // For mix of particles-antiparticles.
  // Calculate lambda parameter for a given pair type via 
  // <N pairs>/<N total pairs>

  cout<<"Computing lambda parameter for non-identical particles:\t"<<part1<<"\t"<<part2<<endl;

  int nSpecPairs = 0;
  int nTotalPairs = 0;
  int nEvents = eventParticles1.size();
  int nEffectiveEvents = 0;

  for(int iEv = 0; iEv < nEvents; iEv++){
    // Sum pairs for each event

    // Ignore events where you have no particles or no antiparticles
    const double nTotalParticles = eventParticles1[iEv]->Integral();
    const double nTotalAntiParticles = eventParticles2[iEv]->Integral();
    // cout<<"nTotalParticles:\t"<<nTotalParticles<<endl;
    // cout<<"nTotalAntiParticles:\t"<<nTotalAntiParticles<<endl;
    if(0 == nTotalParticles) continue;
    if(0 == nTotalAntiParticles) continue;
    nEffectiveEvents++;
    const double nPart1 = eventParticles1[iEv]->GetBinContent(part1+1);
    const double nPart2 = eventParticles2[iEv]->GetBinContent(part2+1);
    // cout<<"const Part1\t"<<nPart1<<endl;
    // cout<<"const Part2\t"<<nPart2<<endl;
    // cout<<"non-const Part1\t"<<eventParticles1[iEv]->GetBinContent(part1+1)<<endl;
    // cout<<"non-const Part2\t"<<eventParticles2[iEv]->GetBinContent(part2+1)<<endl;
    nSpecPairs += nPart1*nPart2;
    nTotalPairs += nTotalParticles * nTotalAntiParticles;
  } //end looping over events
    
  //Compute average pairs and average total pairs
  double avgSpecPairs = (1.*nSpecPairs) / (1.*nEffectiveEvents);
  double avgTotalPairs = (1.*nTotalPairs) / (1.*nEffectiveEvents);
  // cout<<"AvgSpecPairs:\t"<<avgSpecPairs<<",\tAvgTotalPairs"<<avgTotalPairs<<endl;

  // Find lambda (ratio of avgSpec/avgTotal) and error
  cout<<"Calculating lambda now"<<endl;
  lambda = avgSpecPairs/avgTotalPairs;
  cout<<"Lambda is\t"<<lambda<<endl
      <<"Now calculate non-identical pair error"<<endl;
  lambdaError = ComputeLambdaNonIdenticalError(part1, part2, eventParticles1, eventParticles2, avgTotalPairs, lambda);
  return;
}

double ComputeLambdaIdenticalError(const int part1, const int part2, const vector<TH1D*> &eventParticles, const double avgTotalPairs, const double lambda)
{
  // Compute and return the mean error for a given lambda parameter
  cout<<"ComputeLambdaIdenticalError"<<endl;
  // Different formulae for identical particles vs non-identical particles
  // (See ComputeLambdaNonIdenticalError for particle-antiparticle error)
  double errSqr = 0.;
  double nEvents = eventParticles.size();
  int nEffectiveEvents = 0;
  cout<<"Number of events:\t"<<nEvents<<endl;

  // Sum in quadrature the error coming from each particle yield in each event
  // See lab notebook for details of error calculation
  if(part1 == part2){
    for(int iEv = 0; iEv < nEvents; iEv++){
      // cout<<"Event "<<iEv<<"\tP1==P2\n";
      double nTotalYield = eventParticles[iEv]->Integral();
      // cout<<"Total event yield\t"<<nTotalYield;
      if(1 >= nTotalYield) continue;
      nEffectiveEvents++;

      double p1Yield = eventParticles[iEv]->GetBinContent(part1+1);
      // cout<<"p1Yield\t"<<p1Yield<<endl;
      errSqr += pow(lambda,2) * nTotalYield * pow( (0.5 - nTotalYield) ,2);
      // cout<<"After mathing1"<<endl;
      errSqr += p1Yield * (p1Yield - 0.5) * (p1Yield - 0.5 + lambda * (0.5 - nTotalYield) );
      // cout<<"After mathing2"<<endl;
    }
  }
  else{
    for(int iEv = 0; iEv < nEvents; iEv++){
      // cout<<"Event "<<iEv<<"\tP1!=P2\n";
      double nTotalYield = eventParticles[iEv]->Integral();
      // cout<<"Total event yield\t"<<nTotalYield;
      if(1 >= nTotalYield) continue;
      nEffectiveEvents++;
      double p1Yield = eventParticles[iEv]->GetBinContent(part1+1);
      double p2Yield = eventParticles[iEv]->GetBinContent(part2+1);
      // cout<<"p1Yield\t"<<p1Yield<<"\n"
	  // <<"p2Yield\t"<<p2Yield<<"\n";
      errSqr += p1Yield * p2Yield * (p1Yield + p2Yield - 2*lambda * (nTotalYield - 0.5));
      // cout<<"After mathing1"<<endl;
      errSqr += pow(lambda,2) * nTotalYield * pow( (nTotalYield - 0.5), 2);
      // cout<<"After mathing2"<<endl;
    }
  }
  errSqr /= pow(nEffectiveEvents,2);
  errSqr /= pow(avgTotalPairs,2);
  cout<<"Err = "<<sqrt(errSqr)<<endl;
  return sqrt(errSqr);
}

double ComputeLambdaNonIdenticalError(const int part1, const int part2, const vector<TH1D*> &eventParticles1, const vector<TH1D*> &eventParticles2, const double avgTotalPairs, const double lambda)
{
  //Compute and return the mean error for a given lambda parameter
  cout<<"ComputeLambdaNonIdenticalError"<<endl;
  //Different formulae for identical particles vs non-identical particles
  double errSqr = 0.;
  double nEvents = eventParticles1.size();
  int nEffectiveEvents = 0;
  
  // Sum in quadrature the error coming from each particle yield in each event
  // See lab notebook for details of error calculation

  for(int iEv = 0; iEv < nEvents; iEv++){

    // Ignore events where you have no particles or no antiparticles
    const double nTotalYieldParticles = eventParticles1[iEv]->Integral();
    const double nTotalYieldAntiParticles = eventParticles2[iEv]->Integral();

    if(0 == nTotalYieldParticles) continue;
    if(0 == nTotalYieldAntiParticles) continue;
    nEffectiveEvents++;

    const double p1Yield = eventParticles1[iEv]->GetBinContent(part1+1);
    const double p2Yield = eventParticles2[iEv]->GetBinContent(part2+1);
    
    errSqr += p1Yield * p2Yield * (p1Yield + p2Yield - lambda * (nTotalYieldParticles + nTotalYieldAntiParticles));
    errSqr += pow(lambda,2) * nTotalYieldAntiParticles * nTotalYieldParticles * (nTotalYieldParticles + nTotalYieldAntiParticles);
  }

  errSqr /= pow(nEffectiveEvents,2);
  errSqr /= pow(avgTotalPairs,2);
  cout<<"Err = "<<sqrt(errSqr)<<endl;
  return sqrt(errSqr);
}


TH1D *ComputeAverageYields(const vector<TH1D*> &eventParticles, const Particle part)
{
  //Find the average yields and std deviation for each particle type
  //Current error bars show the error of the avg calculation, not the
  //std deviation of the yields
  cout<<"Generating average yield histogram"<<endl;
  int nEvents = eventParticles.size();
  TString histName = "AvgYields";
  if(part == kALam) histName += "Antiparticle";
  else histName += "Particle";
  TH1D *hAvgYields = eventParticles[0]->Clone(histName);
  TH1D *hAvgSquaredYields = eventParticles[0]->Clone("AvgSqrYields");
  hAvgSquaredYields->Multiply(eventParticles[0]);
  for(int iEv = 1; iEv < nEvents; iEv++){
    TH1D *hYield = eventParticles[iEv]->Clone();
    // cout<<"Omegas in this event:\t"<<hYield->GetBinContent(hYield->GetNbinsX())<<endl;
    hAvgYields->Add(hYield);
    hYield->Multiply(hYield);
    hAvgSquaredYields->Add(hYield);
    
    delete hYield; //Do I need to remove it from a directory first?
  }
  // Now that we've summed the yields or yields^2, let's divide
  // by N_events to get the averages
  hAvgYields->Scale(1./(1.*nEvents));
  hAvgSquaredYields->Scale(1./(1.*nEvents));

  // TH1D *hStdDeviationYields = CalculateStdDeviationOfYields(hAvgYields, hStdDeviationYields);
  
  TH1D *hAvgYieldsSquared = hAvgYields->Clone("AvgYieldsSqr");
  hAvgYieldsSquared->Multiply(hAvgYieldsSquared);
  TH1D *hVariance = hAvgSquaredYields; //Just renaming for clarity
  hVariance->Add(hAvgYieldsSquared,-1.);

  //Now take the square root bin by bin and use that as std deviation
  for(int iBin = 1; iBin <= hVariance->GetNbinsX(); iBin++){
    double rootVal = hVariance->GetBinContent(iBin);
    rootVal = sqrt(rootVal);
    hAvgYields->SetBinError(iBin,rootVal);
  }
  
  cout<<"Finished generating average yield histogram"<<endl;
  delete hAvgSquaredYields;
  delete hAvgYieldsSquared;
  return hAvgYields;
}

bool CheckIfPassLambdaCuts(const ParticleCoor *particle)
{
  //Implement all sorts of reconstruction/detection cuts here
  double etaCut = 0.8;
  if( abs( particle->GetEtaP() ) >= etaCut) return false;

  //If we make it here, the particle passed all the cuts
  return true;
}

bool CheckIfPassDaughterCuts(const ParticleCoor *particle, const int pid)
{
  // Check if a daughter track passes various cuts.
  // Cuts can vary based on the type (i.e. pid) of particle

  // Check generic cuts



  // Check pid specific cuts
  if(pid == kProt) {
    // Do proton cuts
  }
  else if(pid == kAntiProt) {

  }
  else if(pid == kPiPlus) {

  }
  else if(pid == kPiMinus) {

  }

  // If we made it this far, the particle has passed all cuts
  return true;
}

void SetLambdaHistAxisLabels(TAxis *axis, const Particle part)
{
  // Use LaTeX to set the bin labels to their
  // corresponding particles
  if(part == kLam){
    axis->SetBinLabel(1,"#Lambda");
    axis->SetBinLabel(2,"#Sigma^{0}");
    axis->SetBinLabel(3,"#Xi^{0}");
    axis->SetBinLabel(4,"#Xi^{-}");
    axis->SetBinLabel(5,"#Omega");
  }
  else {
    axis->SetBinLabel(1,"#bar{#Lambda}");
    axis->SetBinLabel(2,"#bar{#Sigma}^{0}");
    axis->SetBinLabel(3,"#bar{#Xi}^{0}");
    axis->SetBinLabel(4,"#bar{#Xi}^{+}");
    axis->SetBinLabel(5,"#bar{#Omega}");
  }
  axis->SetLabelSize(0.06);
}
