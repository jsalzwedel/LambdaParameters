//**************************************************
//
// Code to estimate lambda parameters using THERMINATOR.  
// This code runs over THERMINATOR events.  
// It counts all primary and secondary lambdas in each event, 
// taking into account some reconstruction cuts.
// It stores each event multiplicty count in a TH1I.
// From the TH1I, the code then computes lambda = <N_ij>/<Sum_ij N_ij>, 
// where N_ij is the number of pairs particles i and j
// The code also estimates the average multiplicity of each species of lambda.
//
// Run the code via .x simpleLambdaEstimate()
//
//**************************************************


#include <cstdio>
#include <cassert>
#include <vector>

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"

bool debug = false;
class ParticleCoor;

//Was having a problem including a therminator file.
//This fixes it.
#define _CXX_VER_ "g++(4.9.1)" 
void simpleLambdaEstimate(TString inputFileName = "event000.root")
{
  // Main function, which outputs and saves histograms of 
  // lambda parameters and average particle multiplicities


  //Load this therminator class
  gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
  gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  
  
  //Read in the therminator event file and get the particle branch
  TFile *thermTFile = new TFile(inputFileName, "READ");
  assert(NULL!=thermTFile);
  TTree *thermTree = (TTree*)thermTFile->Get("particles");
  assert(NULL!=thermTree);
  TBranch *thermBranch = thermTree->GetBranch("particle");
  ParticleCoor *particleEntry = new ParticleCoor();
  thermBranch->SetAddress(particleEntry);


  //Make a vector to hold a pair count histogram for each event
  std::vector<TH1I*> eventParticles;
  UInt_t eventID = 0; //This will be set to each event ID (they are not number 1, 2, ...)
  int eventCounter = -1;
  TH1I *currentHist = NULL;

  //Make a list of the PID codes for each particle we care about
  //We will treat all resonance -> lambda decays as primary lambdas,
  //so we don't need to include those PDG codes here.
  int nParticleTypes = 5;
  int strangePDGs[50] = {3122, //Lambda
			 3212, //Sigma0
			 // 3224, //Sigma*+
			 // 3214, //Sigma*0
			 // 3114, //Sigma*-
			 3322, //Xi0
			 3312, //Xi-
			 // 3324, //Xi*0
			 // 3314, //Xi*-
			 3334}; //Omega-

  //Loop over all the particles and find the lambdas
  int nThermEntries = thermTree->GetEntries();
  for(int i=0; i<nThermEntries; i++)
  {
    //Get the next particle
    int nBytesInEntry = thermTree->GetEntry(i);
    assert(nBytesInEntry>0);
    if(3122==particleEntry->pid) //Is it a lambda? (primary or secondary are allowed here)
    {
      //Make sure the particle passes the reconstruction cuts
      if(!CheckIfPassParticleCuts(particleEntry)) continue;

      //Check if this is a new event.  If so, we'll need to make a new multiplicity histogram
      if((particleEntry->eventid != eventID) || !currentHist){ 
	//We have reached a new event (or the first event)
	cout<<"Previous event: \t"<<eventID<<"\t current event: \t"<<particleEntry->eventid<<endl;
	eventID = particleEntry->eventid;
	eventCounter++;
	//If debugging, do calculations using a single event
	if(debug) if(1 == eventCounter) break;

	//Make a new pair histogram and add it to the vector
	TString histName = "hParticles";
	histName += eventCounter;
	TH1I *hParticles = new TH1I(histName,"Particles Per Type", nParticleTypes, 0, nParticleTypes);
	hParticles->Sumw2();
	eventParticles.push_back(hParticles);

	//Now use this histogram in future particle binning
	currentHist = hParticles;
      }

      //Check parent info
      for(int iPar = 0; iPar < nParticleTypes; iPar++)
      {
	if(particleEntry->fatherpid == strangePDGs[iPar]) {
	  currentHist->Fill(iPar);
	  break;
	}

	// If we reach the end of the last loop iteration,
	// that means that the lambda parent isn't in
	// our short list of PDG codes.  It must be a 
	// short-lived resonance.  We'll just bin that
	// as a primary lambda.
	if(nParticleTypes-1 == iPar) currentHist->Fill(0); // Count all the resonances decays as primary lambdas
      }
    } // end is(lambda)?
  } // end thermEntries loop
  cout<<"Finished looping over particles"<<endl;
  cout<<"Total events used:\t"<<eventCounter+1<<endl;

  //Now that we have found all the strange particles, calculate lambda parameters for each pair type
  TH2D* hLambdaPars = GenerateLambdaParHisto(nParticleTypes, eventParticles);

  //Draw and save the lambda par histo
  hLambdaPars->DrawCopy("colzTEXT");
  TString outfileName = "LambdaPars.root";
  TFile outFile(outfileName,"recreate");
  outFile.cd();
  hLambdaPars->Write();
 
  //Make, draw, and save a histo of average particle yields
  TH1D *hAvgYields = ComputeAverageYields(eventParticles);
  hAvgYields->Write();
  TCanvas *c2 = new TCanvas("yields","Avg Yields");
  hAvgYields->DrawCopy("ptext");
  
}


TH2D *GenerateLambdaParHisto(int nParticleTypes, const vector<TH1I*> &eventParticles){
  // Make a histogram containing lambda paramters for each
  // pair type
  double sumLambdaPars = 0.;
  TH2D *hLambdaPars = new TH2D("LambdaPars", "Lambda Parameters by Pair Type", nParticleTypes, 0, nParticleTypes, nParticleTypes, 0, nParticleTypes);
  // Loop over each type of pairs, and set the lambda
  for(int iPart1 = 0; iPart1 < nParticleTypes; iPart1++){//First type of particle in pair
    for(int iPart2 = iPart1; iPart2 < nParticleTypes; iPart2++){
      double lambdaPar = ComputeLambda(iPart1, iPart2, eventParticles);
      hLambdaPars->SetBinContent(iPart1+1, iPart2+1, lambdaPar);
      sumLambdaPars += lambdaPar;
    }
  }
  //The sum of lambda pars should add up to unity.
  cout<<"Sum of lambda pars:\t"<<sumLambdaPars<<endl;
  return hLambdaPars;
}

double ComputeLambda(const int part1, const int part2, const vector<TH1I*> &eventParticles)
{
  
  // Calculate lambda parameter for a given pair type via <N pairs>/<N total pairs>
  // Should I account somehow for events that have a zero of a certain type of particles?
  int nSpecPairs = 0;
  int nTotalPairs = 0;
  int nEvents = eventParticles.size();

  for(int iEv = 0; iEv < nEvents; iEv++){ //Sum pairs for each event
    if(part1 == part2) { //pairs of identical (parent) particles
      int nPart1 = eventParticles[iEv]->GetBinContent(part1+1);
      nSpecPairs += nPart1*(nPart1 - 1)/2;
    }
    else { //pairs of non-identical (parent) particles
      int nPart1 = eventParticles[iEv]->GetBinContent(part1+1);
      int nPart2 = eventParticles[iEv]->GetBinContent(part2+1);
      nSpecPairs += nPart1*nPart2;
    }
    int totalPart = eventParticles[iEv]->GetEntries();
    nTotalPairs += totalPart * (totalPart - 1)/2;
  } //end looping over events
    
  //Compute average pairs and average total pairs
  double avgSpecPairs = (1.*nSpecPairs) / (1.*nEvents);
  double avgTotalPairs = (1.*nTotalPairs) / (1.*nEvents);

  //Return lambda parameter (ratio of avgSpec/avgTotal)
  return avgSpecPairs/avgTotalPairs;
}



TH1D *ComputeAverageYields(const vector<TH1I*> &eventParticles)
{
  //Find the average yields and std deviation for each particle type
  int nEvents = eventParticles.size();
  TH1D *hAvgYields = eventParticles[0]->Clone("AvgYields");
  for(int iEv = 1; iEv < nEvents; iEv++){
    hAvgYields->Add(eventParticles[iEv]);
  }
  hAvgYields->Scale(1./(1.*nEvents));
  // Compute std deviation...
  return hAvgYields;
}



bool CheckIfPassParticleCuts(ParticleCoor *particle)
{
  //Implement all sorts of reconstruction/detection cuts here
  double etaCut = 0.8;
  if( abs( particle->GetEtaP() ) >= etaCut) return false;

  //If we make it here, the particle passed all the cuts
  return true;
}
