//********************************************************************
// Generate transform matrices for use in residual correlation studies
//
//
//
//
//
//********************************************************************

#include <vector>
#include <map>

enum ParticlePDG {kProt   = 2212,  kAntiProt = -2212, 
		  kPiPlus = 211,   kPiMinus  = -211, 
		  kLam    = 3122,  kALam     = -3122,
		  kSigma  = 3212,  kASigma   = -3212,
		  kXiC    = 3312,  kAXiC     = -3312,
		  kXi0    = 3322,  kAXi0     = -3322,
		  kOmega  = 3334,  kAOmega   = -3334};
class ParticleCoor;

//Was having a problem including a therminator file.
//This fixes it.
#define _CXX_VER_ "g++(4.9.1)" 

Bool_t debug = kFALSE;


void GenerateTransformMatrix(const Int_t nFiles, Int_t parent1PDG, Int_t parent2PDG, Bool_t isLocal)
{
  // Main function.  Specify how many input files to use and this
  // will generate a list of file names to use and pass them on.
  TStopwatch fullTimer;
  // Load in the necessary therminator particle class
  if(isLocal) {
    gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
    gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  }
  else {
      gROOT->LoadMacro("./ParticleCoor.cxx");
  }
  Bool_t isIdentical = kFALSE;
  if(parent1PDG == parent2PDG) isIdentical = kTRUE;
  vector<TString> fileNames = GetTFileNames(nFiles, isLocal);
  gInterpreter->GenerateDictionary("vector<vector<Int_t> >","vector");
  vector< vector<Int_t> > parent1IDs;
  vector< vector<Int_t> > parent2IDs;
  vector< vector<Int_t> > lambdaIDs;


  TStopwatch inputTimer;
  // Loop over each file and gather the track IDs
  for (Int_t iFile = 0; iFile < nFiles; iFile++)
  {
    
    // Read in the file and get the particle branch
    TFile inFile(fileNames[iFile], "read");
    assert(NULL != &inFile);
    TTree *thermTree = (TTree*) inFile.Get("particles");
    assert(NULL != thermTree);

    cout<<"Now using file "<<iFile<<endl;

    Int_t nThermEntries = thermTree->GetEntries();
    Int_t nextEntry = 0;
    
    Int_t nEvents = 0;
    // These will hold all the IDs in a file
    vector<Int_t> parent1IDsFile;
    vector<Int_t> parent2IDsFile;
    vector<Int_t> lambdaIDsFile;
    // Find lambdas and parents for each event, and bin in kstar
    while(nextEntry < nThermEntries) {
      nEvents++;
      // GetIDsOfLambdas will run through each particle starting from 
      // nextEntry and ending when it hits a new event.  That entry 
      // will then be set to nextEntry.
      Int_t firstEntryInEvent = nextEntry;

      // Get all the IDs in an event
      vector<Int_t> lambdaIDsEvent = GetIDsOfLambdas(parent1PDG, parent2PDG, thermTree, nextEntry);
      vector<Int_t> parent1IDsEvent = GetParentIDs(lambdaIDsEvent, parent1PDG, thermTree, firstEntryInEvent);
      vector<Int_t> parent2IDsEvent = GetParentIDs(lambdaIDsEvent, parent2PDG, thermTree, firstEntryInEvent);
      
      // Append the event IDs into the file IDs list
      lambdaIDsFile.insert(lambdaIDsFile.end(), lambdaIDsEvent.begin(), lambdaIDsEvent.end());
      parent1IDsFile.insert(parent1IDsFile.end(), parent1IDsEvent.begin(), parent1IDsEvent.end());
      parent2IDsFile.insert(parent2IDsFile.end(), parent2IDsEvent.begin(), parent2IDsEvent.end());
      if(debug) {
  	cout<<"# of 1st particle\t"<<parent1IDsEvent.size()<<endl;
  	cout<<"# of 2nd particle\t"<<parent2IDsEvent.size()<<endl;
      }
    }
    if(debug) cout<<"Events: \t"<<nEvents<<endl;
    
    // Add the file IDs vector to the list
    lambdaIDs.push_back(lambdaIDsFile);
    parent1IDs.push_back(parent1IDsFile);
    parent2IDs.push_back(parent2IDsFile);

    assert(lambdaIDs.size() == iFile+1);
  }

  cout<<"Real time elapsed during file reading:\t"<<inputTimer.RealTime()<<endl;
  // Loop over the files and use the track IDs to calculate kstar and
  // fill the transform matrix.  We use event mixing for more
  // statistics.
  cout<<"Preparing to fill transform matrix"<<endl;
  TH2D *hTransform = SetupMatrixHist(parent1PDG, parent2PDG);
  for (Int_t iFile1 = 0; iFile1 < nFiles; iFile1++)
  {
    // Read in the file and get the particle branch
    TStopwatch fileTimer;
    cout<<"File1: \t"<<iFile1<<endl;
    TFile inFile1(fileNames[iFile1], "read");
    assert(NULL != &inFile1);
    TTree *thermTree1 = (TTree*) inFile1.Get("particles");
    assert(NULL != thermTree1);

    Int_t file2Start = 0;
    if(isIdentical) file2Start = iFile1;
    for (Int_t iFile2 = file2Start; iFile2 < nFiles; iFile2++)
    {
      TStopwatch smallLoopTimer;
      cout<<"\t\tFile2: \t"<<iFile2; // Stay on the same line for the stopwatch

      TFile inFile2(fileNames[iFile2], "read");
      assert(NULL != &inFile2);
      TTree *thermTree2;
      if(iFile2 == iFile1) thermTree2 = thermTree1;
      else thermTree2 = (TTree*) inFile2.Get("particles");
      assert(NULL != thermTree2);
      
      FillTransformMatrix(hTransform, parent1IDs, parent2IDs, lambdaIDs, iFile1, iFile2, thermTree1, thermTree2, isIdentical);
      cout<<"\tSmall loop real time\t"<<smallLoopTimer.RealTime()<<endl;
    }
    
    // Save a work-in-progress copy of the histogram
    TFile tempOutFile("WIPTransform.root","Update");
    hTransform->Write(hTransform->GetName(),TObject::kOverwrite);
    tempOutFile.Close();
    cout<<"File loop real Time:\t"<<fileTimer.RealTime()<<endl
	<<"File loop CPU Time: \t"<<fileTimer.CpuTime()<<endl;
  }


  
  TString outFileName = "TransformMatrices.root";
  TFile outFile(outFileName, "update");
  outFile.cd();
  hTransform->Write(0,TObject::kOverwrite);
  cout<<"Transform matrix "<<hTransform->GetName()
      <<" written to "<<outFileName<<endl; 
  cout<<"Total real Time:\t"<<fullTimer.RealTime()<<endl
      <<"Total CPU Time: \t"<<fullTimer.CpuTime()<<endl;
  
}



TH2D *SetupMatrixHist(Int_t parent1PDG, Int_t parent2PDG)
{
  // Make the transform matrix histogram
  // Setup the axis labels
  // Make the plot pretty
  // Take as inputs the pdg codes or axis labels 

  cout<<"Setting up the transform matrix histogram.\n";
  
  TString part1Name = GetNameFromPDG(parent1PDG, true);
  TString part2Name = GetNameFromPDG(parent2PDG, true);
  TString histName = "TransformMatrix";
  histName += part1Name;
  histName += part2Name;
  Int_t kstarBins = 800;
  Double_t maxKstar = 2.;
  TH2D *hTransform = new TH2D(histName, histName, 
			      kstarBins, 0., maxKstar,
			      kstarBins, 0., maxKstar);
  
  // Prepare the axis labels
  TString xLabel = "k^{*}_{";
  if(parent1PDG < 0) xLabel += "#bar{#Lambda}";
  else xLabel += "#Lambda";
  if(parent2PDG < 0) xLabel += "#bar{#Lambda}";
  else xLabel += "#Lambda";
  xLabel += "}";
  TString yLabel = "k^{*}_{#";
  yLabel += GetNameFromPDG(parent1PDG,false);
  yLabel += "#";
  yLabel += GetNameFromPDG(parent2PDG,false);
  yLabel += "}";

  hTransform->GetXaxis()->SetTitle(xLabel);
  hTransform->GetYaxis()->SetTitle(yLabel);
  hTransform->SetNdivisions(505,"XY");
  return hTransform;
}

TString GetNameFromPDG(Int_t pdgCode, Bool_t isHistTitle)
{
  //Take a pdg code and return a TString of the particle's name
  TString pdgName = "";
  if((pdgCode < 0) && isHistTitle) pdgName += "Anti";
  else if(pdgCode < 0) pdgName += "bar{#";
  if(fabs(pdgCode) == kLam) pdgName += "Lambda";
  if(fabs(pdgCode) == kSigma) pdgName += "Sigma";
  if((fabs(pdgCode) == kXiC) && isHistTitle) pdgName += "XiC";
  else if(fabs(pdgCode) == kXiC)  pdgName += "Xi^{-}";
  if((fabs(pdgCode) == kXi0)  && isHistTitle) pdgName += "Xi0";
  else if(fabs(pdgCode) == kXi0)  pdgName += "Xi^{0}";
  if(fabs(pdgCode) == kOmega) pdgName += "Omega";

  if((pdgCode < 0) && !isHistTitle) pdgName += "}";
  return pdgName;
}

vector<TString> GetTFileNames(const Int_t nFiles, Bool_t isLocal)
{
  // Use the number of files to generate a list of event file names
  cout<<"Getting the TFile names of "<<nFiles<<" files\n";
  vector<TString> fileNames;
  
  for(Int_t i = 0; i < nFiles; i++){
    TString name;
    if(isLocal) name += "~/Analysis/lambda/AliAnalysisLambda/therminator2/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
    else name += "/home/jsalzwedel/Model/lhyqid3v_LHCPbPb_2760_b2/event";
    if(i < 10) name += "00";
    else if(i < 100) name += "0";
    name += i;
    name += ".root";
    fileNames.push_back(name);
  }

  return fileNames;
}

vector<Int_t> GetIDsOfLambdas(const Int_t parent1PDG, const Int_t parent2PDG, const TTree *thermTree, Int_t &nextEntry)
{
   // Find daughters of lambdas/antilambdas.  Check if they pass cuts.
   // If they do, check if their parent lambda passes cuts.  Return
   // list of track IDs of lambdas that pass cuts.
   if(debug) cout<<"Using charged tracks to find lambda parent IDs\n";

   // Determine which of proton, antiproton, pi+, pi- we need
   vector<Int_t> daughterPDGs = DetermineRelevantDaughterPDGs(parent1PDG, parent2PDG);

   // Get the particle branch
   ParticleCoor *particleEntry = new ParticleCoor();
   TBranch *thermBranch = thermTree->GetBranch("particle");
   thermBranch->SetAddress(particleEntry);
   // cout<<"Test1"<<endl;
   // thermTree->GetEntry(1);
   // cout<<particleEntry<<"\t"<<particleEntry->eid<<endl;
   // thermTree->GetEntry(2);
   // cout<<particleEntry<<"\t"<<particleEntry->eid<<endl;
   // cout<<"Test2"<<endl;

   //Get first particle in this event
   const Int_t startingEntry = nextEntry;
   const Int_t finalEntry = thermTree->GetEntries();
   int nBytesInEntry = thermTree->GetEntry(startingEntry);
   assert(nBytesInEntry > 0);
   UInt_t currentEventID = particleEntry->eventid;

   // Vector to store lambda/antilambda event IDs for daughters that
   // pass cuts
   vector<Int_t> lambdaCandidateIDs;
   // cout<<"This event is "<<currentEventID<<endl;
   // cout<<"The final entry is "<<finalEntry<<endl;
   //loop over all the particles in the event

   while(particleEntry->eventid == currentEventID)
   {
     // cout<<"In event "<<particleEntry->eventid<<". Current entry is "<<nextEntry<<endl;
     const Int_t pid = particleEntry->pid;
     const Int_t parentPid = particleEntry->fatherpid;
     const Int_t parentID = particleEntry->fathereid;
     // cout<<"My pid is "<<pid<<endl;
     // If particle doesn't have a lambda/antilambda parent, ignore it
     if((parentPid == kLam) || (parentPid == kALam)) 
     {
       //check if particle matches one of the relevant daughterPDGs
       for(Int_t iPDG = 0; iPDG < daughterPDGs.size(); iPDG++)
       {
	 if(pid == daughterPDGs[iPDG]) {
	   // if the daughter passes all cuts, add its parent ID to list
	   if(CheckIfPassDaughterCuts(particleEntry)) {
	     lambdaCandidateIDs.push_back(parentID+startingEntry);
	     // cout<<"Pid:\t"<<pid<<"\t\tparentPid:\t"<<parentPid<<"\t\tParentID:\t"<<parentID<<endl;
	     // cout<<"Found a candidate"<<endl;
	   }
	   break;
	 }
       }
     }
     nextEntry++;
     if(finalEntry == nextEntry) {
       if(debug) cout<<"Found end of file"<<endl;
       break;
     } //Reached the end of the file
    
    //Get the next entry
    nBytesInEntry = thermTree->GetEntry(nextEntry);
    assert(nBytesInEntry > 0);
  } // end loop over each track in event


 // Vector to store IDs of (anti)lambdas that pass all cuts
  vector<Int_t> finalV0IDs;
  Int_t nCandidates = lambdaCandidateIDs.size();
  if(debug) cout<<"Number of daughters found:\t"<<nCandidates<<endl;
  for(Int_t iV01=0; iV01 < nCandidates; iV01++)
  {
    // Check which V0 entries appear twice in the list (once for
    // the positive daughter, once for the negative daughter).
    // If it appears twice AND it passes cuts, add it to final V0
    // list.
    const Int_t v0ID = lambdaCandidateIDs[iV01];

    //debug
    Bool_t foundMate = kFALSE;
    // end debug

    for(Int_t iV02 = iV01+1; iV02 < nCandidates; iV02++)
    {
      const Int_t v0ID2 = lambdaCandidateIDs[iV02];
      // cout<<"\t\t"<<iV02<<"\t"<<v0ID2<<endl;
      // cout<<iV01<<"\t"<<iV02<<endl;
      if(v0ID == v0ID2){
	assert(!foundMate);
	foundMate = kTRUE;
    	//Now get the V0 and check if it passes cuts
    	nBytesInEntry = thermTree->GetEntry(v0ID);

	// cout<<"V0ID:\t"<<v0ID<<endl;
	// cout<<"Particle's pid:\t"<<particleEntry->pid<<"\tID:\t"<<particleEntry->eid<<endl;
	// cout<<"Bytes In Entry:\t"<<nBytesInEntry<<endl;
    	assert(nBytesInEntry > 0);
	// cout<<"fabs(pid) = "<<
	assert(fabs(particleEntry->pid) == kLam);
	
    	if(CheckIfPassLambdaCuts(particleEntry, parent1PDG, parent2PDG)) {
    	  // cout<<"Candidate passed cut"<<endl;
    	  finalV0IDs.push_back(v0ID);
	}
    	// else cout<<"Candidate failed cut"<<endl;
    	break;
      }
    }
  } // End loop over V0 candidates

  if(debug){
    for(int i = 0; i < finalV0IDs.size(); i++) {
      cout<<"ID:\t"<<finalV0IDs[i]<<endl;
    } 
  }

  delete particleEntry;
  return finalV0IDs;
}

  // Determine which of proton, antiproton, pi+, pi- we need
vector<Int_t> DetermineRelevantDaughterPDGs(const Int_t parent1PDG, const Int_t parent2PDG)
{
  // Use particle/antiparticle info of lambda/antilambdas to determine
  // which daughters we need to look for
  vector<Int_t> daughterPDGs;
  if((parent1PDG > 0) || (parent2PDG > 0)) { //want lambda daughters
    daughterPDGs.push_back(kProt);
    daughterPDGs.push_back(kPiMinus);
  }
  if((parent1PDG < 0) || (parent2PDG < 0)) { // want antilambda daughters
    daughterPDGs.push_back(kAntiProt); 
    daughterPDGs.push_back(kPiPlus); 
  }
  return daughterPDGs;
}

bool CheckIfPassDaughterCuts(const ParticleCoor *particle)
{
  // Check if a daughter track passes various cuts.
  // Cuts can vary based on the type (i.e. pid) of particle

  // Check generic cuts
  if(particle->GetPt() < 0.16) return false;
  if(fabs(particle->GetEtaP()) > 0.8) return false;
  
  // Check pid specific cuts
  const Int_t pid = particle->pid;
  if(pid == kProt) {
    if(particle->GetPt() < 0.5) return false;
  }
  else if(pid == kAntiProt) {
    if(particle->GetPt() < 0.3) return false;
  }
  else if(pid == kPiPlus) {
  }
  else if(pid == kPiMinus) {
  }
  // If we made it this far, the particle has passed all cuts
  return true;
}

bool CheckIfPassLambdaCuts(const ParticleCoor *particle, const Int_t parent1PDG, const Int_t parent2PDG)
{
  //Implement all sorts of reconstruction/detection cuts here
  double etaCut = 0.8;
  if( fabs( particle->GetEtaP() ) >= etaCut) return false;
  if( particle->GetPt() < 0.4 ) return false;

  // // Make sure the V0 has one of the parents that we care about
  // const Int_t parentPID = particle->fatherpid;
  // if(parentPID == parent1PDG) return true;
  // else if(parentPID == parent2PDG) return true;

  // // If we want "primary" lambdas, check that we also accept
  // // lambdas coming from resonance decays
  // if((fabs(parent1PDG) == 3122) || (fabs(parent2PDG) == 3122)){
  //   //We want primary (anti)lambdas, and (anti)lambdas that come
  //   //from resonanaces.  Let's make sure this particle
  //   //doesn't come from some other source.
    
  //   // loop over other PDG types, return false if this matches that
  //   // PDG value
  // }


  // //If we make it here, the particle passed all the cuts
  return true;
}

vector<Int_t> GetParentIDs(const vector<Int_t> &v0IDs, Int_t parentPDG, TTree *thermTree, Int_t firstEventEntry)
{
  // Use the list of lambdas to find lambdas with the correct parent.
  // Save the parent IDs in the vector.

  // Make a vector to hold parentIDs.  This will be the same size as the
  // LambdaID vector.  If the lambda doesn't have an appropriate parent,
  // just leave the entry as -1.  Otherwise, put in the parent's ID number
  vector<Int_t> parentIDs (v0IDs.size(),-555);

  // Prepare to get particles from the TBranch
  ParticleCoor *particleEntry = new ParticleCoor;
  TBranch *thermBranch = thermTree->GetBranch("particle");
  thermBranch->SetAddress(particleEntry);
  
  // Loop over all the lambdas and find the lambdas with a parent
  // that matches parentPDG
  Int_t nV0s = v0IDs.size();
  // cout<<"Found "<<nV0s<<" in the last step"<<endl;
  for(Int_t iID = 0; iID < nV0s; iID++){
    // Get a V0 and check that it exists and is a (anti)lambda
    Int_t currentID = v0IDs[iID];
    int nBytesInEntry = thermTree->GetEntry(currentID);
    assert(nBytesInEntry > 0);
    assert(fabs(particleEntry->pid) == kLam);
    // Check parentage
    const Int_t fatherPID = particleEntry->fatherpid;
    if(debug) cout<<"FatherPID is \t"<<fatherPID<<". FatherEID is\t"<<particleEntry->fathereid<<endl;
    if(fatherPID == parentPDG) {
      // We found one of the parents we want.  Save its absolute ID.
      if(-1 == particleEntry->fathereid) {
	// a fathereid of -1 means that the daughter was actually a primary particle
	parentIDs[iID] = currentID;
	if(debug) {
	  cout<<"For ID "<<currentID<<". I am a "<<particleEntry->pid
	      <<". My parent's ID is\t"<<parentIDs[iID]
	      <<". My parent's pid is\t"<<fatherPID<<"\tPrimary"<<endl;
	}
      }
      else {
	parentIDs[iID] = particleEntry->fathereid + firstEventEntry;
	if(debug) {
	  cout<<"For ID "<<currentID<<". I am a "<<particleEntry->pid
	      <<". My parent's ID is\t"<<parentIDs[iID]
	      <<". My parent's pid is\t"<<fatherPID<<"\tEMWeak"<<endl;
	}
      }
    }
    // Special case where we want "primary" lambdas.  Here, we also
    // accept lambdas where the parent is a short-lived resonance.
    else if(particleEntry->pid == parentPDG) {
      // Check that it isn't one of the standard weak/EM decays
      bool isElectroWeakDecay = false;
      for(int iEMW = 0; iEMW < 4; iEMW++)
      {
	if( fabs(fatherPID) == fabs(GetElectroWeakPDG(iEMW)) ) {
	  isElectroWeakDecay = true;
	  break;
	}
      }
      if(!isElectroWeakDecay) {
	// This came from a resonance.  Use the lambda's eid.
	parentIDs[iID] = particleEntry->eid + firstEventEntry;
	if(debug){
	  cout<<"For ID "<<currentID<<". I am a "<<particleEntry->pid
	      <<". My parent's ID is\t"<<parentIDs[iID]
	      <<". My parent's pid is\t"<<fatherPID<<"\tResonance"<<endl;
	}
      }
    }
  } // end loop over v0 candidates

  //debug
  if(debug){
    for(int i = 0; i < parentIDs.size(); i++) {
      cout<<"ID:\t"<<parentIDs[i]<<endl;
    } 
  }
  return parentIDs;
}

Int_t GetElectroWeakPDG(Int_t index)
{
  assert ((-1 < index) && (index < 4));
  if     (0 == index) return kSigma;
  else if(1 == index) return kXiC; 
  else if(2 == index) return kXi0; 
  else                return kOmega;
}


void FillTransformMatrix(TH2D *hTransform, vector<vector<Int_t> > &parent1IDs, vector<vector<Int_t> > &parent2IDs, vector<vector<Int_t> > &lambdaIDs, Int_t iFile1, Int_t iFile2, TTree *thermTree1, TTree *thermTree2, Bool_t isIdentical)
{
  if(debug) cout<<"Filling transform matrix\n";
  
  // Prepare to get particles from the TBranch
  ParticleCoor *parent1 = new ParticleCoor;
  ParticleCoor *parent2 = new ParticleCoor;
  ParticleCoor *daughter1 = new ParticleCoor;
  ParticleCoor *daughter2 = new ParticleCoor;

  TBranch *thermBranch1 = thermTree1->GetBranch("particle");
  TBranch *thermBranch2 = thermTree2->GetBranch("particle");
    // thermBranch->SetAddress(particleEntry);
  
  // Loop over both sets of parents and make pairs
  for(Int_t iPart1 = 0; iPart1 < parent1IDs[iFile1].size(); iPart1++)
  {
    // Get the parent particle
    Int_t parID1 = parent1IDs[iFile1][iPart1];
    if(parID1 == -555) continue; 
    thermBranch1->SetAddress(parent1);
    assert(thermTree1->GetEntry(parID1) > 0);

    // Get the daughter particle
    Int_t dauID1 = lambdaIDs[iFile1][iPart1];
    thermBranch1->SetAddress(daughter1);
    assert(thermTree1->GetEntry(dauID1) > 0);
    
    if(debug) cout<<"Particle1:\t"<<iPart1<<endl;
  
    Int_t part2StartValue = 0;
    if(isIdentical && (iFile1 == iFile2)) part2StartValue = iPart1 + 1;
    for(Int_t iPart2 = part2StartValue; iPart2 < parent2IDs[iFile2].size(); iPart2++)
    {
      //Check to make sure we arent using the same parent particle
      // if((iFile1 == iFile2) && (iPart1 == iPart2)) continue;

      // Get the parent particle
      Int_t parID2 = parent2IDs[iFile2][iPart2];
      if(parID2 == -555) continue; 
      thermBranch2->SetAddress(parent2);
      assert(thermTree2->GetEntry(parID2) > 0);

      // Get the daughter particle
      Int_t dauID2 = lambdaIDs[iFile2][iPart2];
      thermBranch2->SetAddress(daughter2);
      assert(thermTree2->GetEntry(dauID2) > 0);

      // We should not have the same daughter or parent IDs
      if(iFile1 == iFile2) {
	assert(dauID1 != dauID2);
	assert(parID1 != parID2);
      }
      
      // Calculate kstar and fill the histogram
      Double_t kstarParents = CalcKstar(parent1, parent2);
      Double_t kstarDaughters = CalcKstar(daughter1, daughter2);
      hTransform->Fill(kstarDaughters, kstarParents);
    } // end particle 2 loop
  } // end particle 1 loop
}  

Double_t CalcKstar(ParticleCoor *part1, ParticleCoor *part2)
{
  // Calculate the relative momentum difference (kstar) for the pair

  // Get the momenta and put them into a four-vector
  Double_t p1[4];
  Double_t p2[4];
  part1->GetMomentum(&p1[3], &p1[0], &p1[1], &p1[2]);
  part2->GetMomentum(&p2[3], &p2[0], &p2[1], &p2[2]);
  TLorentzVector lv1(p1);
  TLorentzVector lv2(p2);

  // Now calculate kstar
  TLorentzVector lvSum = lv1 + lv2; 
  TLorentzVector lvDiff = lv1 - lv2;
  TLorentzVector lvKstar = lvDiff - 
    ((lvDiff*lvSum) / (lvSum*lvSum)) * lvSum;
  lvKstar *= 0.5;
  Double_t kstar = fabs(lvKstar.M());
  // cout<<"kstar is:\t"<<kstar<<endl;
  return kstar;
}


TH2D* BasicCharacterLoadTransform(TString fileName = "TransformMatricesSigSig32.root",
				  TString histName = "TransformMatrixSigmaSigma")
{
  // function for plotting nice transform matrices

  //Read in the file
  // cout<<"Input root file name"<<endl;
  // TString fileName = "TransformMatricesSigSig32.root";
  // cin>>fileName;
  cout<<"Reading from "<<fileName<<endl;
  TFile inFile(fileName);

  //Read in the histogram
  // cout<<"Input 2D hist name"<<endl;
  // TString histName = "TransformMatrixSigmaSigma";
  // cin>>histName;
  cout<<"Loading histogram "<<histName<<endl;
  TH2D *h = (TH2D*) inFile.Get(histName);
  h->SetDirectory(0);
  return h;
}

TH2D* SpecialCharacterLoadTransform(TString fileName = "WIPTransform.root",
				  TString histName = "TransformMatrixXi^{0}Lambda")
{
  TFile file(fileName);
  TList *list = (TList *)file.GetListOfKeys();
  TKey *key = (TKey *)list->FindObject(histName);
  TObject *obj = ((TKey *)key)->ReadObj();
  TH2D *h = (TH2D*)obj;
  h->SetDirectory(0);
  return h;
}

void PlotTransformMatrix(TH2D *h, TString histName, Int_t rebinN)
{
  //Given a transform matrix, plot it

  //Make it pretty and draw it.
  gStyle->SetPalette(55);
  TCanvas *cTransform = new TCanvas("transform",histName,650,650);
  cTransform->SetLogz();
  
  cout<<"Found histogram at "<<h<<endl;
  cout<<"Changing axis range"<<endl;
  h->SetAxisRange(0.,0.5,"X");
  h->SetAxisRange(0.,0.5,"Y");
  h->DrawCopy("colz");

  TCanvas *cTransformRebin = new TCanvas("transformRebin",histName,650,650);
  cTransformRebin->SetLogz();
  cout<<"Rebinning by "<<rebinN<<endl;
  
  // cin>>rebinN;

  h->RebinX(rebinN);
  h->RebinY(rebinN);
  h->SetAxisRange(0.,0.5,"X");
  h->SetAxisRange(0.,0.5,"Y");
  h->DrawCopy("colz");
}

void RunTransformPlotting(TString loadType = "special",
			  TString histName = "TransformMatrixXi^{0}Lambda",
			  TString fileName = "WIPTransform.root",
			  Int_t rebinN = 4)
{
  TH2D *h;

  if(loadType == "special") h = SpecialCharacterLoadTransform(fileName,histName);
  else if(loadType == "basic") h = BasicCharacterLoadTransform(fileName,histName);
  else {cout<<"Not a valid loadType"<<endl; return;}
  
  if(!h) {cout<<"Could not find histogram"<<endl; return;};
  PlotTransformMatrix(h, histName, rebinN);

}




void CheckForDuplicateTracks(Int_t nFiles, Bool_t isLocal)
{
  // Worried that there might be duplicate tracks/events in the
  // Therminator events.  Use a hash table to see if any track
  // appears more than once.

  if(isLocal) {
    gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
    gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  }
  else {
    gROOT->LoadMacro("./ParticleCoor.cxx");
  }
  //vector<TString> fileNames = GetTFileNames(nFiles, isLocal);
  vector<TString> fileNames = GetDebugTFileNames(nFiles, isLocal);

  union Digest
  {
    UChar_t foo[16];
    ULong64_t key;
  };

  std::map<ULong64_t, int> hashes;
  Digest d;
  Int_t nTotalDupes = 0;
  Int_t nTotalParticles = 0;
  vector<Int_t> problemFiles;
  Int_t mostDupes = 0;
  
  for (Int_t iFile = 0; iFile < nFiles; iFile++)
  {
    // Read in the file and get the particle branch
    TFile inFile(fileNames[iFile], "read");
    assert(NULL != &inFile);
    TTree *thermTree = (TTree*) inFile.Get("particles");
    assert(NULL != thermTree);

    cout<<"Now using file "<<fileNames[iFile]<<endl;
    
    Int_t nThermEntries = thermTree->GetEntries();
    // Get the particle branch
    ParticleCoor *particleEntry = new ParticleCoor();
    TBranch *thermBranch = thermTree->GetBranch("particle");
    thermBranch->SetAddress(particleEntry);


    Int_t fileDupes = 0;
    Int_t fileOkays = 0;
    for(Int_t iPart = 0; iPart < nThermEntries; iPart++)
    {
      // Get a particle entry
      int nBytesInEntry = thermTree->GetEntry(iPart);
      assert(nBytesInEntry > 0);

      TMD5 hash;
      hash.Update((UChar_t *)particleEntry,sizeof(ParticleCoor));
      hash.Final(d.foo);
      hashes[d.key]++;
      if(hashes[d.key] > 1) {
	//cout<<"File:\t"<<iFile<<"\t"
	//    <<"Particle:\t"<<iPart<<"\t"
	//    <<"Duplicate Key:\t"<<d.key<<"\t"
	//    <<"\tNumber of dupes:\t"<<hashes[d.key]<<endl;
	nTotalDupes++;
	fileDupes++;

	// Find the max number of problem dupes
	if(hashes[d.key] > mostDupes) mostDupes = hashes[d.key];
	
	// Count the number of problem files problem file
	if( (problemFiles.size() > 0) && (problemFiles.Final() != iFile) ) {
	  problemFiles.push_back(iFile);
	  cout<<"New problem file:\t"<<iFile<<endl;
	}
      }
      else {
	fileOkays++;
      }
      nTotalParticles++;
    }
    cout<<"Dupes in this file:\t"<<fileDupes<<endl;
    cout<<"Okays in this file:\t"<<fileOkays<<endl;
  }
  cout<<"End of analysis.\n"
      <<"Total files:\t"<<nFiles<<"\n"
      <<"Total particles:\t"<<nTotalParticles<<"\n"
      <<"Total duplicates:\t"<<nTotalDupes<<"\n"
      <<"Most duplicates of one track:\t"<<mostDupes<<"\n";
  cout<<"\nProblem files:\n";
  for(Int_t iFile = 0; iFile < problemFiles.size(); iFile++){
    cout<<problemFiles[iFile]<<endl;
  }
}


vector<TString> GetDebugTFileNames(const Int_t nFiles, Bool_t isLocal)
{
  // Use the number of files to generate a list of event file names
  cout<<"Getting the TFile names of "<<nFiles<<" files\n";
  vector<TString> fileNames;

  TString nameBase;
  if(isLocal) nameBase += "~/Analysis/lambda/AliAnalysisLambda/therminator2/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
  else nameBase += "/home/jsalzwedel/Model/lhyqid3v_LHCPbPb_2760_b2/event";


  //Weird stuff starts happening with the 28th file.  Let's put that and 29 near the beginning of the lineup.
  fileNames.push_back(nameBase + "028.root");
  fileNames.push_back(nameBase + "001.root");
  fileNames.push_back(nameBase + "029.root");

  
  for(Int_t i = 2; i < nFiles+3; i++){
    if((i == 28) || (i == 29)) continue;
    TString name;
    if(isLocal) name += "~/Analysis/lambda/AliAnalysisLambda/therminator2/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
    else name += "/home/jsalzwedel/Model/lhyqid3v_LHCPbPb_2760_b2/event";
    if(i < 10) name += "00";
    else if(i < 100) name += "0";
    name += i;
    name += ".root";
    fileNames.push_back(name);
  }

  return fileNames;
}
