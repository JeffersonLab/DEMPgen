/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Convert .txt file into ROOT tree and store in .root file.
  Includes routines to handle data given by Goloskokov & Kroll,
  and data from http://rprmodel.ugent.be/calc/
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"

using namespace std;

void GKConvert(string in_filepath, string out_filepath){
  ifstream in_filestream;

  in_filestream.open(in_filepath.c_str());

  if (! in_filestream){
    cout << "File not Found" << endl;
    exit(1);
  }

  //Skip to Data
  string line;
  int i;
  for (i = 0; i < 18; i++) {
    getline(in_filestream, line);
  }

  //Create a tree
  TFile out_TFile(out_filepath.c_str(), "UPDATE");

  gDirectory->Delete("GK_Raw;*");
  //This ensures only the most recent tree in in the file

  TTree GK_Raw("GK_Raw", "Asymmetry Data");

  //Initialize variables for each data column
  Double_t w, q2, tp, t, asy, asysfi, asy2fi, asyfpfs, asy3f;

  //Create branches for each variable
  GK_Raw.Branch("w", &w, "w/D");
  GK_Raw.Branch("q2", &q2, "q2/D");
  GK_Raw.Branch("tp", &tp, "tp/D");
  GK_Raw.Branch("t",&t,"t/D");
  GK_Raw.Branch("asy",&asy,"asy/D");
  GK_Raw.Branch("asysfi",&asysfi,"asysfi/D");
  GK_Raw.Branch("asy2fi",&asy2fi,"asy2fi/D");
  GK_Raw.Branch("asyfpfs",&asyfpfs,"asyfpfs/D");
  GK_Raw.Branch("asy3f",&asy3f,"asy3f/D");

  //Fill the tree
  while (!in_filestream.eof()) {
    in_filestream >> w >> q2 >> tp >> t >> asy >> asysfi >> asy2fi
                  >> asyfpfs >> asy3f;
    //cout << w << "\t" << q2 << "\t" << tp <<"\t";
    //cout << t << "\t" << asy << "\t" << asysfi << "\t";
    //cout << asy2fi << "\t" << asyfpfs << "\t" << asy3f;
    //cout << endl;
    GK_Raw.Fill();
  }

  GK_Raw.Write();
}

void VR_SigL(string in_filepath, string out_filepath){
  ifstream in_filestream;

  in_filestream.open(in_filepath.c_str());

  if (! in_filestream){
    cout << "File not Found" << endl;
    exit(1);
  }

  //Skip to Data
  string line;
  int i;
  for (i = 0; i < 3; i++) {
    getline(in_filestream, line);
  }

  //Create a tree
  TFile out_TFile(out_filepath.c_str(), "UPDATE");

  gDirectory->Delete("VR_SigL_Raw;*");
  //This ensures only the most recent tree in in the file

  TTree VR_SigL_Raw("VR_SigL_Raw", "VR Model");

  //Initialize variables for each data column
  Double_t Q2, W, Negt, sigL;

  //Create branches for each variable
  VR_SigL_Raw.Branch("Q2", &Q2, "Q2/D");
  VR_SigL_Raw.Branch("W", &W, "W/D");
  VR_SigL_Raw.Branch("Negt", &Negt, "-t/D");
  VR_SigL_Raw.Branch("dsig_L",&sigL,"dsig_L/D");

  //Fill the tree
  while (!in_filestream.eof()) {
    in_filestream >> Q2 >> W >> Negt >> sigL;
    VR_SigL_Raw.Fill();
  }

  VR_SigL_Raw.Write();
}


void AsyPars(string in_filepath, string out_filepath,
             string asyname, int nPars)
{
  ifstream in_filestream;

  in_filestream.open(in_filepath.c_str());

  if (! in_filestream){
    cout << "File not Found" << endl;
    exit(1);
  }


  //Create a tree
  TFile out_TFile(out_filepath.c_str(), "UPDATE");

  TTree t1(asyname.c_str(), "Asymmetry");

  double pars[nPars];
  double Qsq;

  char llist[100] = "pars[%d]/D";
  char temp[100];

  sprintf(temp, llist, nPars);

  t1.Branch("Qsq",&Qsq,"Qsq/D");
  t1.Branch("pars",&pars,temp);

  while(!in_filestream.eof()){
    in_filestream >> Qsq;
    for (int i=0; i<nPars; i++){
      in_filestream >> pars[i];
    }
    t1.Fill();
  }

  t1.Write();
}

void AsyAll()
{
  AsyPars("../data/input/asyout.txt",
          "test.root",
          "asy",
          4);
  AsyPars("../data/input/asysfiout.txt",
          "test.root",
          "asysfi",
          3);
  AsyPars("../data/input/asy2fiout.txt",
          "test.root",
          "asy2fi",
          3);
  AsyPars("../data/input/asyfpfsout.txt",
          "test.root",
          "asyfpfs",
          3);
  AsyPars("../data/input/asy3fout.txt",
          "test.root",
          "asy3f",
          3);
}
