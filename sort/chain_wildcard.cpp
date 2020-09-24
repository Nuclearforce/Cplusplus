#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <cstdio>
#include "TSystem.h"

void chain_wildcard(){
    TChain DataChain("data");
    TTreeReader myReader(&DataChain);
    TTreeReaderArray<unsigned int> TRenergy(myReader, "energy");
    TTreeReaderArray<unsigned int> TRtime(myReader, "time");
    TTreeReaderArray<unsigned int> TRmodchan(myReader, "modchan");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0*.root");
    DataChain.ls();
}
