//********************************************************************
// Generate transform matrices for use in residual correlation studies
//
//
//
//
//
//********************************************************************

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

void GenerateTransformMatrix(const Int_t nFiles, Int_t parent1PDG, Int_t parent2PDG)
{
  // Main function.  Specify how many input files to use and this
  // will generate a list of file names to use and pass them on.

  // Load in the necessary therminator particle class
  gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
  gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");

  TH2D *hTransform = SetupMatrixHist(parent1PDG, parent2PDG);
  vector<TString> fileNames = GetTFileNames(nFiles);

  // Loop over each file
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
    

    // Find lambdas and parents for each event, and bin in kstar
    while(nextEntry < nThermEntries) {
      // GetIDsOfLambdas will run through each particle starting from 
      // nextEntry and ending when it hits a new event.  That entry 
      // will then be set to nextEntry.
      Int_t firstEntryInEvent = nextEntry;
      vector<Int_t> lambdaIDs  = GetIDsOfLambdas(parent1PDG, parent2PDG, thermTree, nextEntry);
      vector<Int_t> parent1IDs = GetParentIDs(lambdaIDs, parent1PDG, thermTree, firstEntryInEvent);
      vector<Int_t> parent2IDs = GetParentIDs(lambdaIDs, parent2PDG, thermTree, firstEntryInEvent);
      FillTransformMatrix(hTransform, parent1IDs, parent2IDs, lambdaIDs, thermTree);
    }
  }
  
  TString outFileName = "TransformMatrices.root";
  TFile outFile(outFileName, "update");
  outFile.cd();
  hTransform->Write(0,TObject::kOverwrite);
  cout<<"Transform matrix "<<hTransform->GetName()
      <<" written to "<<outFileName<<endl; 
}



TH2D *SetupMatrixHist(Int_t parent1PDG, Int_t parent2PDG)
{
  // Make the transform matrix histogram
  // Setup the axis labels
  // Make the plot pretty
  // Take as inputs the pdg codes or axis labels 

  cout<<"Setting up the transform matrix histogram.\n";
  
  TString part1Name = GetNameFromPDG(parent1PDG);
  TString part2Name = GetNameFromPDG(parent2PDG);
  TString histName = "TransformMatrix";
  histName += part1Name;
  histName += part2Name;
  Int_t kstarBins = 800;
  Double_t maxKstar = 2.;
  TH2D *hTransform = new TH2D(histName, histName, 
			      kstarBins, 0., maxKstar,
			      kstarBins, 0., maxKstar);

  // ***************  Antiparticle labels ****************


  TString xLabel = "k^{*}_{";
  if(parent1PDG < 0) xLabel += "#bar{#Lambda}";
  else xLabel += "#Lambda";
  if(parent2PDG < 0) xLabel += "#bar{#Lambda}";
  else xLabel += "#Lambda";
  xLabel += "}";
  TString yLabel = "k^{*}_{#";
  yLabel += part1Name;
  yLabel += "#";
  yLabel += part2Name;
  yLabel += "}";

  hTransform->GetXaxis()->SetTitle(xLabel);
  hTransform->GetYaxis()->SetTitle(yLabel);
  hTransform->SetNdivisions(505,"XY");
  return hTransform;
}

TString GetNameFromPDG(Int_t pdgCode)
{
  //Take a pdg code and return a TString of the particle's name
  TString pdgName = "test";
  //...
  return pdgName;
}

vector<TString> GetTFileNames(const Int_t nFiles)
{
  // Use the number of files to generate a list of event file names
  cout<<"Getting the TFile names of "<<nFiles<<" files\n";
  vector<TString> fileNames;
  
  for(Int_t i = 0; i < nFiles; i++){
    TString name = "~/Analysis/lambda/AliAnalysisLambda/therminator2/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
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
  cout<<"Using charged tracks to find lambda parent IDs\n";

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
  cout<<"This event is "<<currentEventID<<endl;
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
	    // cout<<"Found a candidate"<<endl;
	  }
	  break;
	}
      }
    }
    nextEntry++;
    if(finalEntry == nextEntry) {cout<<"Found end of file"<<endl; break;} //Reached the end of the file
    
    //Get the next entry
    nBytesInEntry = thermTree->GetEntry(nextEntry);
    assert(nBytesInEntry > 0);
  } // end loop over each track in event


 // Vector to store IDs of (anti)lambdas that pass all cuts
  vector<Int_t> finalV0IDs;
  Int_t nCandidates = lambdaCandidateIDs.size();
  cout<<"Number of candidates:\t"<<nCandidates<<endl;
  for(Int_t iV01=0; iV01 < nCandidates; iV01++)
  {
    // Check which V0 entries appear twice in the list (once for
    // the positive daughter, once for the negative daughter).
    // If it appears twice AND it passes cuts, add it to final V0
    // list.
    const Int_t v0ID = lambdaCandidateIDs[iV01];
    // cout<<iV01<<"\t"<<v0ID<<endl;
    // const Int_t v0ID2 = lambdaCandidateIDs[iV01+1];
    // cout<<iV01<<"\t"<<iV01+1<<endl;
    for(Int_t iV02 = iV01+1; iV02 < nCandidates; iV02++)
    {
      const Int_t v0ID2 = lambdaCandidateIDs[iV02];
      // cout<<"\t\t"<<iV02<<"\t"<<v0ID2<<endl;
      // cout<<iV01<<"\t"<<iV02<<endl;
      if(v0ID == v0ID2){
    	//Now get the V0 and check if it passes cuts
    	nBytesInEntry = thermTree->GetEntry(v0ID);
    	assert(nBytesInEntry > 0);
    	assert(particleEntry->pid == abs(kLam));
	
    	if(CheckIfPassLambdaCuts(particleEntry, parent1PDG, parent2PDG)) {
    	  // cout<<"Candidate passed cut"<<endl;
    	  finalV0IDs.push_back(v0ID);
	}
    	// else cout<<"Candidate failed cut"<<endl;
    	break;
      }
    }
  } // End loop over V0 candidates

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
  if(abs(particle->GetEtaP()) > 0.8) return false;
  
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
  if( abs( particle->GetEtaP() ) >= etaCut) return false;
  if( particle->GetPt() < 0.4 ) return false;

  // // Make sure the V0 has one of the parents that we care about
  // const Int_t parentPID = particle->fatherpid;
  // if(parentPID == parent1PDG) return true;
  // else if(parentPID == parent2PDG) return true;

  // // If we want "primary" lambdas, check that we also accept
  // // lambdas coming from resonance decays
  // if((abs(parent1PDG) == 3122) || (abs(parent2PDG) == 3122)){
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
  cout<<"Finding parents of lambdas\n";


  // Make a vector to hold parentIDs.  This will be the same size as the
  // LambdaID vector.  If the lambda doesn't have an appropriate parent,
  // just leave the entry as -1.  Otherwise, put in the parent's ID number
  vector<Int_t> parentIDs (v0IDs.size(),-1);

  // Prepare to get particles from the TBranch
  ParticleCoor *particleEntry = new ParticleCoor;
  TBranch *thermBranch = thermTree->GetBranch("particle");
  thermBranch->SetAddress(particleEntry);
  
  // Loop over all the lambdas and find the lambdas with a parent
  // that matches parentPDG
  Int_t nV0s = v0IDs.size();
  cout<<"Found "<<nV0s<<" in the last step"<<endl;
  for(Int_t iID = 0; iID < nV0s; iID++){
    // Get a V0 and check that it exists and is a (anti)lambda
    Int_t currentID = v0IDs[iID];
    int nBytesInEntry = thermTree->GetEntry(currentID);
    assert(nBytesInEntry > 0);
    assert(abs(particleEntry->pid) == kLam);
    // Check parentage
    const Int_t fatherPID = particleEntry->fatherpid;
    if(fatherPID == parentPDG) {
      // We found one of the parents we want.  Save its absolute ID.
      parentIDs[iID] = particleEntry->fathereid + firstEventEntry;
    }
    // Special case where we want "primary" lambdas.  Here, we also
    // accept lambdas where the parent is a short-lived resonance.
    else if(particleEntry->pid == parentPDG) {
      // Check that it isn't one of the standard weak/EM decays
      bool isElectroWeakDecay = false;
      for(int iEMW = 0; iEMW < 4; iEMW++)
      {
	if( abs(fatherPID) == abs(GetElectroWeakPDG(iEMW)) ) {
	  isElectroWeakDecay = true;
	  break;
	}
      }
      if(!isElectroWeakDecay) {
	// This came from a resonance.  Use the lambda's eid.
	parentIDs[iID] = particleEntry->eid + firstEventEntry;
      }
    }
  } // end loop over v0 candidates
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


void FillTransformMatrix(TH2D *hTransform, vector<Int_t> &parent1IDs, vector<Int_t> &parent2IDs, vector<Int_t> &lambdaIDs, TTree *thermTree)
{
  cout<<"Filling transform matrix\n";

  // Prepare to get particles from the TBranch
  ParticleCoor *parent1 = new ParticleCoor;
  ParticleCoor *parent2 = new ParticleCoor;
  ParticleCoor *daughter1 = new ParticleCoor;
  ParticleCoor *daughter2 = new ParticleCoor;

  TBranch *thermBranch = thermTree->GetBranch("particle");
    // thermBranch->SetAddress(particleEntry);
  
  // Loop over both sets of parents and make pairs
  for(Int_t iPart1 = 0; iPart1 < parent1IDs.size(); iPart1++)
  {
    // Get the parent particle
    Int_t parID1 = parent1IDs[iPart1];
    if(parID1 == -1) continue; 
    thermBranch->SetAddress(parent1);
    assert(thermTree->GetEntry(parID1) > 0);

    // Get the daughter particle
    Int_t dauID1 = lambdaIDs[iPart1];
    thermBranch->SetAddress(daughter1);
    assert(thermTree->GetEntry(dauID1) > 0);
    

    for(Int_t iPart2 = 0; iPart2 < parent2IDs.size(); iPart2++)
    {
      //Check to make sure we arent using the same parent particle
      if(iPart1 == iPart2) continue;

      // Get the parent particle
      Int_t parID2 = parent2IDs[iPart2];
      if(parID2 == -1) continue; 
      thermBranch->SetAddress(parent2);
      assert(thermTree->GetEntry(parID2) > 0);

      // Get the daughter particle
      Int_t dauID2 = lambdaIDs[iPart2];
      thermBranch->SetAddress(daughter2);
      assert(thermTree->GetEntry(dauID2) > 0);

      // We should not have the same daughter or parent IDs
      assert(dauID1 != dauID2);
      assert(parID1 != parID2);
      
      // Calculate kstar and fill the histogram
      Double_t kstarParents = CalcKstar(parent1, parent2);
      Double_t kstarDaughters = CalcKstar(daughter1, daughter2);
      hTransform->Fill(kstarDaughters, kstarParents);
    }

  }
  
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
  cout<<"kstar is:\t"<<kstar<<endl;
  return kstar;
}
