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

  TH2D *hTransform = SetupMatrixHist();
  
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



TH2D *SetupMatrixHist()
{
  // Make the transform matrix histogram
  // Setup the axis labels
  // Make the plot pretty
  // Take as inputs the pdg codes or axis labels 

  TH2D *hTransform;
  // ...
  return hTransform;
}

vector<TString> GetTFileNames(const Int_t nFiles)
{
  // Use the number of files to generate a list of event file names
  vector<TString> fileNames;
  // ...
  return fileNames;
}


