//********************************************************************
// Look at a bunch of THERMINATOR root files and check if any of the 
// particles are identical (i.e. have the same hash).
//********************************************************************

#include <vector>
#include <map>

class ParticleCoor;

//Was having a problem including a therminator file.
//This fixes it.
#define _CXX_VER_ "g++(4.9.1)" 


void CheckForDuplicateTracks(Int_t nFiles, bool useDebugFileList = kFALSE)
{
  // Worried that there might be duplicate tracks/events in the
  // Therminator events.  Use a hash table to see if any track
  // appears more than once.  A "duplicate" track is one that has
  // previously appeared in any event/file (including the current 
  // event/file).


  // Load in the ParticleCoor class (change this to the correct
  // directory)
  gROOT->LoadMacro("./ParticleCoor.cxx");


  // Set useDebugFileList to false to look at a generic list of files.
  // Set useDebugFileList to true if you want to compare a few specific files.
  // Bool_t useDebugFileList = kFALSE;
  vector<TString> fileNames = GetTFileNames(nFiles, useDebugFileList);

  // We'll use this Digest union when we make a hash of ParticleCoor objects
  union Digest
  {
    UChar_t foo[16];
    ULong64_t key;
  };

  // Make a hash table
  std::map<ULong64_t, int> hashes;
  Digest d;

  // Some simple counting variables
  Int_t nTotalDupes = 0;
  Int_t nTotalParticles = 0;
  vector<Int_t> problemFiles;
  Int_t mostDupes = 0;
  
  for (Int_t iFile = 0; iFile < nFiles; iFile++)
  {
    // Read in the file and get the particle branch
    cout<<"Now reading in file "<<iFile<<" located at "<<fileNames[iFile]<<endl;
    TFile inFile(fileNames[iFile], "read");
    assert(NULL != &inFile);
    TTree *thermTree = (TTree*) inFile.Get("particles");
    assert(NULL != thermTree);

    
    
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

      // Make a hash of this ParticleCoor and add it to the hash table
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
	
	// Count the number of problem files
	if( (problemFiles.size() > 0) && (problemFiles.back() != iFile) ) {
	  problemFiles.push_back(iFile);
	  cout<<"New problem file:\t"<<iFile<<endl;
	}
      }
      else {
	fileOkays++;
      }
      nTotalParticles++;
    }
    cout<<"Repeat particles found in this file:\t"<<fileDupes<<endl;
    cout<<"New particles found in this file:\t"<<fileOkays<<endl;
  }
  cout<<"End of analysis.\n"
      <<"Total files:\t"<<nFiles<<"\n"
      <<"Total particles:\t"<<nTotalParticles<<"\n"
      <<"Total duplicates:\t"<<nTotalDupes<<"\n"
      <<"Most duplicate copies of one track:\t"<<mostDupes<<"\n";
  cout<<"\nProblem files:\n";
  for(Int_t iFile = 0; iFile < problemFiles.size(); iFile++){
    cout<<problemFiles[iFile]<<endl;
  }
}

vector<TString> GetTFileNames(const Int_t nFiles, Bool_t useDebugFileList)
{
  // Use the number of files to generate a list of event file names
  // 
  cout<<"Getting the TFile names of "<<nFiles<<" files\n";
  vector<TString> fileNames;

  // TString nameBase = "/home/jsalzwedel/Model/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
  // TString nameBase = "/home/jsalzwedel/Model/lhyqid3v_LHCPbPb_2760_b2/event";
  // TString nameBase = "/home/jsalzwedel/therminator/therminator2Old/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
  TString nameBase = "/home/jsalzwedel/therminator/therminator2Fixed/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";

  if(useDebugFileList)
  {
    // It looks like files 28 and 12 contain some duplicate particles.
    // Let's look at them specifically.
    // Likewise with 29 and 13 seem to have a similar problem.
    fileNames.push_back(nameBase + "028.root");
    fileNames.push_back(nameBase + "012.root");
  }
  
  for(Int_t i = 0; i < nFiles; i++){
    if(useDebugFileList && ((i == 28) || (i == 12))) continue;
    TString name = nameBase;
    char nameTemp[1024];
    sprintf(nameTemp, "%s%03d.root", nameBase.Data(), i);
    TString name(nameTemp);
    fileNames.push_back(name);
  }

  return fileNames;
}
