#include <cstdio>
#include <cassert>
#include <vector>

#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
// #include "TH11"
// #include "TH2D"

#include "/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include/ParticleCoor.h"
// #include "ParticleCoor.h"

bool debug = false;

void simpleLambdaEstimate(char inTFileName[] = "event000.root")
{
  TFile *thermTFile = new TFile(inTFileName, "READ");
  assert(NULL!=thermTFile);

  TTree *thermTree = (TTree*)thermTFile->Get("particles");
  assert(NULL!=thermTree);

  TBranch *thermBranch = thermTree->GetBranch("particle");
  ParticleCoor *particleEntry = new ParticleCoor();
  thermBranch->SetAddress(particleEntry);

  //      thermTree->SetBranchAddress("particle", &thermBranch);

  //Make a vector to hold a pair count histogram for each event
  std::vector<TH1I*> eventParticles;
  UInt_t eventID = 0;
  int eventCounter = -1;
  TH1I *currentHist = NULL;

  int strangePDGs[10] = {3122, //Lambda
			 3212, //Sigma0
			 3224, //Sigma*+
			 3214, //Sigma*0
			 3114, //Sigma*-
			 3322, //Xi0
			 3312, //Xi-
			 3324, //Xi*0
			 3314, //Xi*-
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
      //We have reached a new event (or first event)
      if((particleEntry->eventid != eventID) || !currentHist){ 
	cout<<"Previous event: \t"<<eventID<<"\t current event: \t"<<particleEntry->eventid<<endl;
	eventID = particleEntry->eventid;
	eventCounter++;
	if(debug) if(1 == eventCounter) break; //Calc lam with one evt
	//Make a new pair histogram and add it to the vector
	TString histName = "hParticles";
	histName += eventCounter;
	TH1I *hParticles = new TH1I(histName,"Particles Per Type", 10, 0, 10);
	eventParticles.push_back(hParticles);
	//Now use this histogram
	currentHist = hParticles;
      }
      //Check parent info
      for(int iPar = 0; iPar < 10; iPar++)
      {
	if(particleEntry->fatherpid == strangePDGs[iPar]) {
	  currentHist->Fill(iPar);
	  break;
	}
      }
      // if(particleEntry->fatherpid == 3334) cout<<"Found an Omega\n";
    }
  } // end thermEntreis loop

  //Now that we have found all the strange particles, calculate lambda parameters for each pair type
  double totalLambda = 0.;
  TH2D *hLambdaPars = new TH2D("LambdaPars", "Lambda Parameters by Pair Type", 10, 0, 10, 10, 0, 10);
  for(int iPart1 = 0; iPart1 < 10; iPart1++){//First type of particle in pair
    for(int iPart2 = iPart1; iPart2 < 10; iPart2++){
      double lambdaPar = ComputeLambda(iPart1, iPart2, eventParticles);
      hLambdaPars->SetBinContent(iPart1+1, iPart2+1, lambdaPar);
      totalLambda += lambdaPar;
    }
  }
  cout<<"Total events used:\t"<<eventCounter+1<<endl;
  cout<<"Sum of lambda pars:\t"<<totalLambda<<endl;

  //Draw and save the lambda par histo
  hLambdaPars->DrawCopy("colzTEXT");
  TString outfileName = "LambdaPars.root";
  TFile outFile(outfileName,"recreate");
  outFile.cd();
  hLambdaPars->Write();
  TH1D *hAvgYields = ComputeAverageYields(eventParticles);
  hAvgYields->Write();
  TCanvas *c2 = new TCanvas("yields","Avg Yields");
  hAvgYields->DrawCopy("text");
  
}

double ComputeLambda(const int part1, const int part2, const vector<TH1I*> &eventParticles)
{
  // Very naive way of calculating <pairs>/<total pairs>
  // Should probably account somehow for events that have particles==0

  int nSpecPairs = 0;
  int nTotalPairs = 0;
  int nEvents = eventParticles.size();
  // cout<<"Number Events:\n";

  for(int iEv = 0; iEv < nEvents; iEv++){ //Sum pairs for each event
    if(part1 == part2) { //particles with indentical pairs
      int nPart1 = eventParticles[iEv]->GetBinContent(part1+1);
      nSpecPairs += nPart1*(nPart1 - 1)/2;
      // cout<<"Part1\t"<<part1<<"\tPart2\t"<<part2
      // 	  <<"\tPairs:\t"<<nPart1*(nPart1 - 1)/2<<endl;
    }
    else { //pairs of particles with non-identical parents
      int nPart1 = eventParticles[iEv]->GetBinContent(part1+1);
      int nPart2 = eventParticles[iEv]->GetBinContent(part2+1);
      nSpecPairs += nPart1*nPart2;
      // cout<<"Part1\t"<<part1<<"\tPart2\t"<<part2
      // 	  <<"\tPairs:\t"<<nPart1*nPart2<<endl;
    }
    int totalPart = eventParticles[iEv]->GetEntries();
    nTotalPairs += totalPart * (totalPart - 1)/2;
    // cout<<"Events:\t"<<nEvents<<"\tTotalPairs\t"
    // 	<<totalPart * (totalPart - 1)/2<<endl;
  } //end looping over events
    
  //Compute average pairs and average total pairs
  double avgSpecPairs = (1.*nSpecPairs) / (1.*nEvents);
  double avgTotalPairs = (1.*nTotalPairs) / (1.*nEvents);
  // cout<<"avgSpecPairs\t"<<avgSpecPairs
  //     <<"\tavgTotalPairs\t"<<avgTotalPairs<<endl;
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
