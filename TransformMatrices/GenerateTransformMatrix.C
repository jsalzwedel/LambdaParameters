//********************************************************************
// Generate transform matrices for use in residual correlation studies
//
//
//
//
//
//********************************************************************





void GenerateTransformMatrix(const int nFiles)
{
  // Main function.  Specify how many input files to use and this
  // will generate a list of file names to use and pass them on.
  
  TH2D *hTransform = SetupMatrixHist();
  
  vector<TString> fileNames = GetTFileNames(nFiles);

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

vector<TString> GetTFileNames(const int nFiles)
{
  // Use the number of files to generate a list of event file names
  vector<TString> fileNames;
  // ...
  return fileNames;
}
