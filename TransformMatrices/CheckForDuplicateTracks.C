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


void CheckForDuplicateTracks(Int_t nFiles, Bool_t isLocal)
{
  // Worried that there might be duplicate tracks/events in the
  // Therminator events.  Use a hash table to see if any track
  // appears more than once.  A "duplicate" track is one that has
  // previously appeared in any event/file (including the current 
  // event/file).

  if(isLocal) {
    gInterpreter->AddIncludePath("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/include");
    gROOT->LoadMacro("/home/jai/Analysis/lambda/AliAnalysisLambda/therminator2/build/src/ParticleCoor.cxx");
  }
  else {
    gROOT->LoadMacro("./ParticleCoor.cxx");
  }
  // vector<TString> fileNames = GetTFileNames(nFiles, isLocal);
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
      <<"Most duplicates of one track:\t"<<mostDupes<<"\n";
  cout<<"\nProblem files:\n";
  for(Int_t iFile = 0; iFile < problemFiles.size(); iFile++){
    cout<<problemFiles[iFile]<<endl;
  }
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



vector<TString> GetDebugTFileNames(const Int_t nFiles, Bool_t isLocal)
{
  // Use the number of files to generate a list of event file names
  // 
  cout<<"Getting the TFile names of "<<nFiles<<" files\n";
  vector<TString> fileNames;

  TString nameBase;
  if(isLocal) nameBase += "~/Analysis/lambda/AliAnalysisLambda/therminator2/events/lhyquid3v-LHCPbPb2760b2.3Ti512t0.60Tf140a0.08b0.08h0.24x2.3v2/event";
  else nameBase += "/home/jsalzwedel/Model/lhyqid3v_LHCPbPb_2760_b2/event";


  // It looks like file 28 and 12 contain some duplicate particles.
  // Likewise with 29 and 13.  Let's look at them specifically.
  fileNames.push_back(nameBase + "028.root");
  fileNames.push_back(nameBase + "012.root");
  fileNames.push_back(nameBase + "029.root");
  fileNames.push_back(nameBase + "013.root");

  
  for(Int_t i = 0; i < nFiles+3; i++){
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
