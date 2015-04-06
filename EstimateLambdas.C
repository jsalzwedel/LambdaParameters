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

bool debug = false;
class ParticleCoor;

//Was having a problem including a therminator file.
//This fixes it.
#define _CXX_VER_ "g++(4.9.1)" 

void EstimateLambdas(int nFiles=2)
{
  // Main function. Generates list of file names
  // and pass them to RunSimpleLambdaEstimate
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
  RunSimpleLambdaEstimate(inputFileNames);
}


void RunSimpleLambdaEstimate(/*TString inputFileName = "event000.root"*/ vector<TString> inputFileNames)
{
  // Outputs and saves histograms of 
  // lambda parameters and average particle multiplicities

  int nFiles = inputFileNames.size();

  //Load this therminator class
  gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
  gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  
  //Make a vector to hold a pair count histogram for each event
  std::vector<TH1D*> eventParticles;
  int eventCounter = -1;


  ParticleCoor *particleEntry = new ParticleCoor();


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

  
  
  for(int iFile = 0; iFile < nFiles; iFile++)
  {

    //If debugging, do calculations using a single event
    if(debug) if(1 == eventCounter) break;

    //Read in the therminator event file and get the particle branch
    TFile thermTFile(inputFileNames[iFile], "READ");
    assert(NULL!=&thermTFile);
    TTree *thermTree = (TTree*)thermTFile.Get("particles");
    assert(NULL!=thermTree);
    TBranch *thermBranch = thermTree->GetBranch("particle");
    thermBranch->SetAddress(particleEntry);

    UInt_t eventID = 0; //This will be set to each event ID (they are not number 1, 2, ...)

    TH1D *currentHist = NULL;

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
	if( (particleEntry->eventid != eventID) || !currentHist){ 
	  //We have reached a new event (or the first event)
	  eventID = particleEntry->eventid;
	  eventCounter++;
	  cout<<"Processing event \t"<<eventCounter<<"\tEvent ID: \t"<<particleEntry->eventid<<endl;
	  //If debugging, do calculations using a single event
	  if(debug) if(1 == eventCounter) break;

	  //Make a new pair histogram and add it to the vector
	  TString histName = "hParticles";
	  histName += eventCounter;
	  TH1D *hParticles = new TH1D(histName,"Particles Per Type", nParticleTypes, 0, nParticleTypes);
	  hParticles->SetDirectory(0); //Disassociate it from the TFile
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
    } // end loop over particle entries
  } // end loop over files
  cout<<"Finished looping over particles"<<endl;
  cout<<"Total events used:\t"<<eventCounter+1<<endl;

  //Now that we have found all the strange particles, calculate lambda parameters for each pair type
  TH2D* hLambdaPars = GenerateLambdaParHisto(nParticleTypes, eventParticles);
  SetLambdaHistAxisLabels(hLambdaPars->GetXaxis());
  SetLambdaHistAxisLabels(hLambdaPars->GetYaxis());

  //Draw and save the lambda par histo
  hLambdaPars->DrawCopy("colzTEXTe");
  TString outfileName = "LambdaPars.root";
  TFile outFile(outfileName,"recreate");
  outFile.cd();
  hLambdaPars->Write();
 
  //Make, draw, and save a histo of average particle yields
  TH1D *hAvgYields = ComputeAverageYields(eventParticles);
  SetLambdaHistAxisLabels(hAvgYields->GetXaxis());
  hAvgYields->Write();
  TCanvas *c2 = new TCanvas("yields","Avg Yields");
  hAvgYields->DrawCopy("ptexte");
  
}


TH2D *GenerateLambdaParHisto(int nParticleTypes, const vector<TH1D*> &eventParticles){
  // Make a histogram containing lambda paramters for each
  // pair type
  double sumLambdaPars = 0.;
  TH2D *hLambdaPars = new TH2D("LambdaPars", "Lambda Parameters by Pair Type", nParticleTypes, 0, nParticleTypes, nParticleTypes, 0, nParticleTypes);
  // Loop over each type of pairs, and set the lambda
  for(int iPart1 = 0; iPart1 < nParticleTypes; iPart1++){//First type of particle in pair
    for(int iPart2 = iPart1; iPart2 < nParticleTypes; iPart2++){
      double lambdaPar = -1.;
      double lambdaParError = 1.;
      ComputeLambda(iPart1, iPart2, eventParticles, lambdaPar, lambdaParError);
      hLambdaPars->SetBinContent(iPart1+1, iPart2+1, lambdaPar);
      hLambdaPars->SetBinError(iPart1+1, iPart2+1, lambdaParError);
      sumLambdaPars += lambdaPar;
    }
  }
  //The sum of lambda pars should add up to unity.
  cout<<"Sum of lambda pars:\t"<<sumLambdaPars<<endl;
  return hLambdaPars;
}

void ComputeLambda(const int part1, const int part2, const vector<TH1D*> &eventParticles, double &lambda, double &lambdaError)
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

  // Find lambda (ratio of avgSpec/avgTotal) and error
  lambda = avgSpecPairs/avgTotalPairs;
  lambdaError = ComputeLambdaError(part1, part2, eventParticles, avgTotalPairs, lambda);
  return;
}

double ComputeLambdaError(const int part1, const int part2, const vector<TH1D*> &eventParticles, const double avgTotalPairs, const double lambda)
{
  //Compute and return the mean error for a given lambda parameter
  
  //Different formulae for identical particles vs non-identical particles
  double errSqr = 0.;
  double nEvents = eventParticles.size();
  
  // Sum in quadrature the error coming from each particle yield in each event
  // See lab notebook for details of error calculation
  if(part1 == part2){
    for(int iEv = 0; iEv < nEvents; iEv++){
      double nTotalYield = eventParticles[iEv]->Integral();
      double p1Yield = eventParticles[iEv]->GetBinContent(part1+1);
      errSqr += pow(lambda,2) * nTotalYield * pow( (0.5 - nTotalYield) ,2);
      errSqr += p1Yield * (p1Yield - 0.5) * (p1Yield - 0.5 + lambda * (0.5 - nTotalYield) );
    }
  }
  else{
    for(int iEv = 0; iEv < nEvents; iEv++){
      double nTotalYield = eventParticles[iEv]->Integral();
      double p1Yield = eventParticles[iEv]->GetBinContent(part1+1);
      double p2Yield = eventParticles[iEv]->GetBinContent(part2+1);
      errSqr += pow(lambda,2) * nTotalYield * pow( (0.5 - nTotalYield) ,2);
      errSqr += pow(p1Yield,2) * ( p1Yield + lambda * (0.5 - nTotalYield) );
      errSqr += pow(p2Yield,2) * ( p2Yield + lambda * (0.5 - nTotalYield) );
    }
  }
  errSqr /= pow(nEvents,2);
  errSqr /= pow(avgTotalPairs,2);

  return sqrt(errSqr);
}


TH1D *ComputeAverageYields(const vector<TH1D*> &eventParticles)
{
  //Find the average yields and std deviation for each particle type
  //Current error bars show the error of the avg calculation, not the
  //std deviation of the yields
  int nEvents = eventParticles.size();
  TH1D *hAvgYields = eventParticles[0]->Clone("AvgYields");
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
  
  delete hAvgSquaredYields;
  delete hAvgYieldsSquared;
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

void SetLambdaHistAxisLabels(TAxis *axis)
{
  // Use LaTeX to set the bin labels to their
  // corresponding particles
  axis->SetBinLabel(1,"#Lambda");
  axis->SetBinLabel(2,"#Sigma^{0}");
  axis->SetBinLabel(3,"#Xi^{0}");
  axis->SetBinLabel(4,"#Xi^{-}");
  axis->SetBinLabel(5,"#Omega");
  axis->SetLabelSize(0.06);
}