//********************************************************************
// Generate transform matrices for use in residual correlation studies
//
//
//
//
//
//********************************************************************

enum ParticlePID {kProt = 2212, kAntiProt = -2212, 
		  kPiPlus = 211, kPiMinus = -211, 
		  kLam = 3122, kALam = -3122};



void GenerateTransformMatrix(const Int_t nFiles)
{
  // Main function.  Specify how many input files to use and this
  // will generate a list of file names to use and pass them on.
  
  Int_t parent1PDG; //Pass this in somehow
  Int_t parent2PDG;

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

    // do stuff ...
    Int_t nThermEntries = thermTree->GetEntries();
    Int_t nextEntry = 0;

    // Find lambdas and parents for each event, and bin in kstar
    while(nextEntry < nThermEntries) {
      // GetIDsOfLambdas will run through each particle starting from 
      // nextEntry and ending when it hits a new event.  That entry 
      // will then be set to nextEntry.
      vector<Int_t> lambdaIDs  = GetIDsOfLambdas(parent1PDG, parent2PDG, thermTree, nextEntry);
      vector<Int_t> parent1IDs = GetParentIDs(&lambdaIDs, parent1PDG, thermTree);
      vector<Int_t> parent2IDs = GetParentIDs(&lambdaIDs, parent2PDG, thermTree);
      FillTransformMatrix(hTransform, &parent1IDs, &parent2IDs, &lambdaIDs, thermTree);
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
  
  TString part1Name = GetNameFromPDG(parent1PDG);
  TString part2Name = GetNameFromPDG(parent2PDG);
  TString histName = "TransformMatrix";
  histName += part1Name;
  histName += part2Name;
  Int_t kstarBins = 800;
  Double_t maxKstar = 2.;
  TH2D *hTransform = new TH2D(histName, histName, 
			      kstarBins, 0., maxKstar,
			      kstarBins, 0. maxKstar);
  TString xLabel = "k^{*}_{#Lambda#Lambda}";
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
  TString pdgName;
  //...
  return pdgName;
}

vector<TString> GetTFileNames(const Int_t nFiles)
{
  // Use the number of files to generate a list of event file names
  vector<TString> fileNames;
  // ...
  return fileNames;
}

vector<Int_t> GetIDsOfLambdas(const Int_t parent1PDG, const Int_T parent2PDG, const TTree *thermTree, Int_t &firstEntry)
{


  // Get the particle branch
  ParticleCoor *particleEntry = new ParticleCoor();
  TBranch *thermBranch = thermTree->GetBranch("particle");
  thermBranch->SetAddress(particleEntry);

  //Get first particle in this event
  int nBytesInEntry = thermTree->GetEntry(firstEntry);
  assert(nBytesInEntry > 0);
  UInt_t currentEventID = particleEntry->eventid;

  // Determine which of proton, antiproton, pi+, pi- we need
  vector<Int_t> daughterPIDs;
  if((parent1PDG > 0) || (parent2PDG > 0)) { //want lambdas
    daughterPIDs.push_back(kProt);
    daughterPIDs.push_back(kPiMinus);
  }
  if((parent1PDG < 0) || (parent2PDG < 0)) { // want antilambdas
    daughterPIDs.push_back(kAntiProt); 
    daughterPIDs.push_back(kPiPlus); 
  }
  
  //loop over all the particles in the event
  while(particleEntry->eventid == currentEventID)
  {


    //check if particle matcheS one of daughterPIDs
    //check if it has a lambda parent
    //check if the daughter passes cuts
    //add parent track ID to list

  }

  for(/* each parent trac ID*/)
  {
    //check if it appears twice in the list (once for each daughter)
    //check if it passes cuts
    //add track ID to final list
  }

  delete particleEntry;

  return lambdaIDs;
}

vector<Int_t> GetParentIDs(vector<Int_t> &lambdaIDs, Int_t parentPDG, TTree *thermTree)
{

}

void FillTransformMatrix(TH2D *hTransform, vector<Int_t> &parent1IDs, vector<Int_t> &parent2IDs, vector<Int_t> &lambdaIDs, TTree *thermTree);

