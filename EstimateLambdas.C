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
// Run the code via .x EstimateLambdas(pairType, nFiles)
// pairType is 1, 2, or 3 (LamLam, ALamALam, LamALam)
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
	       kLam = 3122, kALam = -3122};
enum PairType {kLamLam=0, kALamALam=1, kLamALam=2};

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
  
  // // Finally, before we run, set the seed for the random
  // // number generator that handles efficiency generation
  // gRandom->SetSeed(1);
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
  // cout<<"Address of eventParticles1\t"<<&eventParticles1<<endl;
  // cout<<"Address of eventParticles2\t"<<&eventParticles2<<endl;
  // cout<<"Address of secondParticleCollection\t"<<&secondParticleCollection<<endl;
  if(kLamALam==pairType) {
    cout<<"I am changing the particle collection.  Pair Type:\t"<<pairType<<endl;
    &eventParticles2 = secondParticleCollection; //We do have separate particles and antiparticles, so use a second collection
  }
  int nParticleTypes = 5;

  // cout<<"Address of eventParticles1\t"<<&eventParticles1<<endl;
  // cout<<"Address of eventParticles2\t"<<&eventParticles2<<endl;
  // cout<<"Address of secondParticleCollection\t"<<&secondParticleCollection<<endl;
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
    SetLambdaHistAxisLabels(hLambdaPars->GetYaxis(),kALam);
  }
  // This will make the text larger when we use Draw("colztexte")
  hLambdaPars->SetMarkerSize(1.5);

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

void UseDaughtersToFindParents(const TTree *thermTree, Int_t &nextEntry, const Int_t finalEntry, const PairType thisPairType, vector<int> &lambdaIDs, vector<int> &antiLambdaIDs)
{
  // Returns the particle IDs of all lambdas whose daughters
  // have been found.
  // Code loops over all the particles in an event.
  // Find all the protons and pions that are daughters of lambdas.
  // Check if each particle passes particle cuts.
  // If both proton and pion daughter pass cuts, then add
  // lambda particle ID to vector.

  vector<Int_t> protParentIDs;
  vector<Int_t> antiprotParentIDs;
  vector<Int_t> piplusParentIDs;
  vector<Int_t> piminusParentIDs;
  // cout<<"Using daughters to find parents.  Pair type is\t"<<thisPairType<<endl;
  bool wantLams = false;
  bool wantALams = false;
  if( (thisPairType == kLamLam) || (thisPairType == kLamALam) ) {wantLams = true;}
  if( (thisPairType == kALamALam) || (thisPairType == kLamALam) ) {wantALams = true;}

  ParticleCoor *particleEntry = new ParticleCoor();
  TBranch *thermBranch = thermTree->GetBranch("particle");
  thermBranch->SetAddress(particleEntry);

  //Get first particle in this event
  int nBytesInEntry = thermTree->GetEntry(nextEntry);
  assert(nBytesInEntry > 0);
  UInt_t currentEventID = particleEntry->eventid;
  if(globalDebug) cout<<"Event ID\t"<<currentEventID<<endl;
  Int_t iPart = nextEntry;
  //////////////////
  // Loop over particles.  Stop when we reach a new event
  if(globalDebug) cout<<"About to loop over particles"<<endl;


  while(particleEntry->eventid == currentEventID)
  {
    const Int_t pid = particleEntry->pid;
    const Int_t parentPid = particleEntry->fatherpid;
    const Int_t parentID = particleEntry->fathereid;
    if( (pid == kProt) && wantLams && (parentPid == kLam) ) {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	protParentIDs.push_back(parentID+nextEntry);
      }
    }
    else if( (pid == kAntiProt) && wantALams && (parentPid == kALam) ) 
    {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	antiprotParentIDs.push_back(parentID+nextEntry);
      }
    }
    else if( (pid == kPiPlus) && wantALams && (parentPid == kALam) )
    {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	piplusParentIDs.push_back(parentID+nextEntry);
      }
    }
    else if( (pid == kPiMinus) && wantLams && (parentPid == kLam) ) {
      if(CheckIfPassDaughterCuts(particleEntry,pid)) {
	piminusParentIDs.push_back(parentID+nextEntry);
      }
    }

    iPart++;
    if(iPart == finalEntry) break;
    nBytesInEntry = thermTree->GetEntry(iPart);
    assert(nBytesInEntry > 0);
  } //end while loop

  ///////////////
  // Check for (anti)lambdas in both daughter lists
  if(globalDebug) cout<<"Comparing daughter ID lists"<<endl;
  // cout<<"wantLams\t"<<wantLams<<"\twantALams\t"<<wantALams<<endl;
  // int lambdaCount = 0;
  if(wantLams) {
    if(globalDebug) cout<<"Number of proton daughters\t"<<protParentIDs.size()<<endl;
    if(globalDebug) cout<<"Number of pion daughters\t"<<piminusParentIDs.size()<<endl;
    for(int iProt = 0; iProt < protParentIDs.size(); iProt++){
      const Int_t v0ID = protParentIDs[iProt];
      for(int iPi = 0; iPi < piminusParentIDs.size(); iPi++){
	if(v0ID == piminusParentIDs[iPi]) {
	  lambdaIDs.push_back(v0ID);
	  // cout<<"LambdaCount:\t"<<++lambdaCount<<endl;

	  break;
	}
      }
    }
  }
  // int antilambdaCount = 0;
  if(wantALams) {
    if(globalDebug) cout<<"Number of antiproton daughters\t"<<antiprotParentIDs.size()<<endl;
    if(globalDebug) cout<<"Number of piplus daughters\t"<<piplusParentIDs.size()<<endl;
    for(int iProt = 0; iProt < antiprotParentIDs.size(); iProt++){
      const Int_t v0ID = antiprotParentIDs[iProt];
      // cout<<"Looking at antiproton number\t"<<iProt+1<<endl;
      for(int iPi = 0; iPi < piplusParentIDs.size(); iPi++){
	if(v0ID == piplusParentIDs[iPi]) {
	  antiLambdaIDs.push_back(v0ID);
	  // cout<<"AntilambdaCount:\t"<<++antilambdaCount<<endl;
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


    Int_t nextEventFirstEntry = 0;

    Int_t nThermEntries = thermTree->GetEntries();


    while(nextEventFirstEntry < nThermEntries) {
      // For each event, find the daughters and use them to find the parents
      eventCounter++;
      cout<<"Processing event \t"<<eventCounter<<endl;
      vector<Int_t> lambdaIDs;
      vector<Int_t> antiLambdaIDs;
      UseDaughtersToFindParents(thermTree, nextEventFirstEntry, nThermEntries, pairType, lambdaIDs, antiLambdaIDs);
      if((kLamLam==pairType) || (kLamALam == pairType)) cout<<"Lambda candidates\t"<<lambdaIDs.size()<<endl;
      if((kALamALam==pairType) || (kLamALam == pairType)) cout<<"AntiLambda candidates\t"<<antiLambdaIDs.size()<<endl;
      
      //Only push back histogram vectors if we find V0s.
      if(!(lambdaIDs.size() == 0)) {
	TH1D *lambdaYields = FillLambdaYieldHist(thermTree,lambdaIDs,eventCounter,kLam, nParticleTypes);
	if(lambdaYields && lambdaYields->GetEntries() > 0) {
	  eventParticles1.push_back(lambdaYields);
	}
      }
      if(!(antiLambdaIDs.size() == 0)) {
	TH1D *antiLambdaYields = FillLambdaYieldHist(thermTree,antiLambdaIDs,eventCounter,kALam,nParticleTypes);
	if(antiLambdaYields && antiLambdaYields->GetEntries() > 0) {
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

TH1D *FillLambdaYieldHist(const TTree *thermTree, const vector<Int_t> &v0IDs, const int eventCounter, const Particle part, const int nParticleTypes)
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
    Int_t currentID = v0IDs[iID];
    int nBytesInEntry = thermTree->GetEntry(currentID);
    assert(nBytesInEntry > 0);
    Int_t pid = particleEntry->pid;
    // cout<<"I should be a\t"<<part<<"\tMy pid is\t"<<pid<<"\teid:\t"<<currentID<<endl;
    // if(kALam == pid) cout<<"AntiLambda!!! My eid is \t"<<currentID<<endl;
    assert(pid == part);

    //Make sure the particle passes reconstruction cuts
    if(!CheckIfPassLambdaCuts(particleEntry)) continue;
    
    //Bin the particle according to its parent info
    for(int iPar = 0; iPar < nParticleTypes; iPar++)
    {
      const Int_t fpid = particleEntry->fatherpid;
      if((kLam == part) && (fpid == strangePDGs[iPar])) {
	hParticles->Fill(iPar);
	break;
      }
      else if((kALam == part) && (fpid == -1*strangePDGs[iPar])) {
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
  cout<<"Computing identical lambda params"<<endl
      // <<"eventParticle object:\t"<<eventParticles1<<"\n"
      <<"eventParticle address:\t"<<&eventParticles1<<"\n";
  int nSpecPairs = 0;
  int nTotalPairs = 0;
  int nEvents = eventParticles1.size();
  cout<<"Number of events in collection is\t"<<nEvents<<endl;
  int nEffectiveEvents = 0;
  // cout<<"Looping over events"<<endl;
  for(int iEv = 0; iEv < nEvents; iEv++){ 
    // cout<<"Event\t"<<iEv<<endl;
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
  assert(nEffectiveEvents > 0);
  cout<<"effective events: \t"<<nEffectiveEvents<<endl;
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
  assert(nEffectiveEvents > 0);
  cout<<"effective events: \t"<<nEffectiveEvents<<endl;
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

bool CheckIfPassDaughterCuts(const ParticleCoor *particle, const int pid)
{
  // Check if a daughter track passes various cuts.
  // Cuts can vary based on the type (i.e. pid) of particle

  // Check generic cuts
  if(particle->GetPt() < 0.16) return false;
  if(abs(particle->GetEtaP()) > 0.8) return false;

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

bool CheckIfPassLambdaCuts(const ParticleCoor *particle)
{
  //Implement all sorts of reconstruction/detection cuts here
  double etaCut = 0.8;
  if( abs( particle->GetEtaP() ) >= etaCut) return false;
  if( particle->GetPt() < 0.4 ) return false;

  // //If we make it here, the particle passed all the cuts
  // return true;

  //Finally, let's simulate efficiency
  return SimulateV0Efficiency();
}

bool SimulateV0Efficiency(/*const double pT*/)
{
  // Simulate the efficiency of reconstructing V0 efficiency.
  //This could be pT dependent. But for now, just roll the
  //dice and keep the particle if it falls under the limit.

  //Overall efficiency is ~15%.  But that includes the 64%
  //branching ratio, which is already accounted for by only
  //taking proton+pion daughters. So, factoring 64% out of 
  //15% gives us 23%.
  double efficiency = 0.23;
  if(gRandom->Rndm() < efficiency) return true;
  else return false;
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
