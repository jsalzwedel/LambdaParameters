//********************************************************************
// Generate transform matrices for use in residual correlation studies
//
//
//
//
//
//********************************************************************





void GenerateTransformMatrix(const Int_t nFiles)
{
  // Main function.  Specify how many input files to use and this
  // will generate a list of file names to use and pass them on.
  
  Int_t parent1PDG; //Pass this in somehow
  Int_t parent2PDG;

  TH2D *hTransform = SetupMatrixHist(parent1PDG, parent2PDG);
  
  vector<TString> fileNames = GetTFileNames(nFiles);

  for (Int_t iFile = 0; iFile < nFiles; iFile++)
  {
    TFile inFile(fileNames[iFile], "read");
    assert(NULL != &inFile);

    TTree *thermTree = (TTree*) inFile.Get("particles");
    assert(NULL != thermTree);

    // do stuff ...
    for(/* each particular event */) {
      vector<Int_t> lambdaIDs = GetPassingLambdaIDs(thermTree);
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


