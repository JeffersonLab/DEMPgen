// Original file from Love Preet, University of Regina
// This file converts nEvents from an input HEPMC3 file (raw from DEMPGen or afterburned) into an output root tree.
// Modified on - 04/03/23 - by Stephen JD Kay University of Regina

#define HEPMC3_Converter_cxx

// Include relevant stuff
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <TSystem.h>
#include <TTree.h>

// SJDK - 04/03/23 - Read in a few things as arguments
void HEPMC3_Converter(Long_t nEvents = -1, string InputHEPMC3File = "", string OutputRootFile = "", string FileType = ""){
  
  // 04/04/23 - SJDK - Explicitly set batch mode so it doesn't spam stuff to screen
  gROOT->SetBatch(kTRUE);

  // SJDK - 04/03/23 - Explicitly require 3 input arguments
  if(nEvents == -1){
    cout << "Enter a the number of events to convert: ";
    cin >> nEvents;
    if( nEvents<=0 ) return;
  }

  if(InputHEPMC3File == ""){
    cout << "Enter a HEPMC3 file to convert: ";
    cin >> InputHEPMC3File;
  } 

  if(OutputRootFile == "") {
    cout << "Enter a root file to write to: ";
    cin >> OutputRootFile;
  }
  if(FileType == ""){
    cout << "File type not specified, enter raw or AB. Raw for DEMPGen output, AB for afterburned files." << endl << "Defaulting to raw." << endl;
    FileType = "raw";
  }
  else if (FileType != "raw" && FileType != "AB"){
    cout << "Invalid file type specified, enter only raw or AB. Defaulting to raw." << endl;
    FileType = "raw";  
  } 

  Long_t nSkip; // SJDK - 04/04/23 - Number of lines to skip, depends upon file type
  if (FileType == "raw"){
    nSkip = 2;
  }
  else if (FileType == "AB"){
    nSkip = 21;
  }
  
  TString TInputHEPMC3File = InputHEPMC3File;
  TString TOutputRootFile = OutputRootFile;
  
  if (gSystem->AccessPathName(TInputHEPMC3File) == kTRUE){
    cerr << "!!!!! ERROR !!!!! " << endl << TInputHEPMC3File <<  " not found" << endl <<  "!!!!! ERRROR !!!!!" << endl;
    exit;
  }
  
  fstream HEPMC3In;
  HEPMC3In.open(TInputHEPMC3File, ios::in);
  
  //.............................................................................................................................................  
  // Stored the whole file in a vector
  //.............................................................................................................................................     
  string s;     
  vector<string> v;
  Long_t nLines = 0;
  Long_t nEvents_read = 0;

  for (Long_t i = 0; i <nSkip; i++){ // Skip the first 21 lines - NOTE - Only valid for afterburner! 
    getline(HEPMC3In,s);
  }
 
  while(getline(HEPMC3In,s)){
    nLines++;
    if (HEPMC3In.eof()) break;
    v.push_back(s);  // store all the file to the vector
  }

  if (FileType == "raw"){
    nEvents_read = nLines/9; // 04/04/23 - Number of events in the file is the number of lines read divided by 9 (each event is 9 lines)
  }

  else if (FileType == "AB"){
    nEvents_read = (nLines-2)/9; // 04/04/23 - Number of events in the file is the number of lines read (-2, AB has two tailing lines) divided by 9 (each event is 9 lines)
  }
  
  HEPMC3In.close();

  // 04/04/23 - Check number of events doesn't exceed number in the file!
  if (nEvents > nEvents_read){
    cerr << "!!!!! ERROR !!!!! " << endl << "Requested to process more events - ("<< nEvents << ") than are actually in the file - (" << nEvents_read << ")" << endl << "Double check and try again!" << endl <<  "!!!!! ERRROR !!!!!" << endl;
    exit;
  }
  else if (nEvents != nEvents_read){
    cout << endl << "!!!!! NOTICE !!!!!" << endl << "Requested number of events to process is SMALLER than the number of events in the file, this is fine, but re-run if you want all events." << endl;
    cout << "There are - " << nEvents_read << " in your input HEPMC3 file" << endl << "!!!!! NOTICE !!!!!" << endl << endl;
  }

  //.............................................................................................................................................  
  // Getting the values from the vector 
  //.............................................................................................................................................       
  vector<double> weight; 
  string w1,w2; // columns in the weight line
  double w3,w4; // colums in the weight line
  
  vector<double> e_px; 
  vector<double> e_py; 
  vector<double> e_pz; 
  vector<double> e_E; 
  string l_e; // columns in the outgoing electron line
  double id_e,vid_e,pdg_e,px_e,py_e,pz_e,en_e,m_e,s_e; // columns in the outgoing electron line
  
  vector<double> k_px; 
  vector<double> k_py; 
  vector<double> k_pz; 
  vector<double> k_E; 
  string l_k; // columns in the outgoing kaon line
  double id_k,vid_k,pdg_k,px_k,py_k,pz_k,en_k,m_k,s_k; // columns in the outgoing kaon line
  
  vector<double> l_px; 
  vector<double> l_py; 
  vector<double> l_pz; 
  vector<double> l_E; 
  string l_l; // columns in the outgoing lambda line
  double id_l,vid_l,pdg_l,px_l,py_l,pz_l,en_l,m_l,s_l; // columns in the outgoing lambda line
  if(nEvents <= nEvents_read){ // Only process if less than or equal to the # events in the file
    for(int i=0; i<nEvents; i++){          // Loop over the number of events

      //.............................................................................................................................................     
      int w=(2+(i*9));  // Accessing the weights line from the whole file 
      //  cout<<v[2]<<" "<<endl;
  
      stringstream weight_l; // Extract the doubles from the vector using stringstream (accessing the weight line for each event)
   
      weight_l << v[w]; // store the vector to string stream
        
      while (! weight_l.eof()) {
        
	weight_l >>w1>>w3>>w2>>w4;
	weight.push_back(w4);
           
      }
     
      int o_e=(6+(i*9));  // Accessing the outgoing particle(i.e.electron) line from the whole file 
  
      stringstream o_e_l; // Extract the doubles from the vector using stringstream (accessing the outgoing particle line for each event)
   
        
      o_e_l << v[o_e]; // store the vector to string stream
        
      //cout << o_e_l.str() <<endl;
        
      while (! o_e_l.eof()) {
          
	o_e_l>>l_e>>id_e>>vid_e>>pdg_e>>px_e>>py_e>>pz_e>>en_e>>m_e>>s_e;
	//cout<<" "<<px_e<<endl;
           
	e_px.push_back(px_e);
	e_py.push_back(py_e);
	e_pz.push_back(pz_e);
	e_E.push_back(en_e);
      }
      //.............................................................................................................................................            
      int o_k=(7+(i*9));  // Accessing the outgoing particle(i.e.kaon) line from the whole file 
  
      stringstream o_k_l; // Extract the doubles from the vector using stringstream (accessing the outgoing particle line for each event)
   
      o_k_l << v[o_k]; // store the vector to string stream
        
      while (! o_k_l.eof()) {
        
	o_k_l >>l_k>>id_k>>vid_k>>pdg_k>>px_k>>py_k>>pz_k>>en_k>>m_k>>s_k;
	k_px.push_back(px_k);
	k_py.push_back(py_k);
	k_pz.push_back(pz_k);
	k_E.push_back(en_k);
      }
      //.............................................................................................................................................            
      int o_l=(8+(i*9));  // Accessing the outgoing particle(i.e.lambda) line from the whole file 
  
      stringstream o_l_l; // Extract the doubles from the vector using stringstream (accessing the outgoing particle line for each event)
   
      o_l_l << v[o_l]; // store the vector to string stream
        
      while (! o_l_l.eof()) {
        
	o_l_l >>l_l>>id_l>>vid_l>>pdg_l>>px_l>>py_l>>pz_l>>en_l>>m_l>>s_l;
	l_px.push_back(px_l);
	l_py.push_back(py_l);
	l_pz.push_back(pz_l);
	l_E.push_back(en_l);
      }
    }//-> End of for loop over the events

    //.............................................................................................................................................  
    // Storing the weights in Root TTree
    //.............................................................................................................................................       
  
    double w_gp,e_px_gp,e_py_gp,e_pz_gp,e_E_gp;
    double k_px_gp,k_py_gp,k_pz_gp,k_E_gp;
    double l_px_gp,l_py_gp,l_pz_gp,l_E_gp;
    
    TFile *OutputFile = new TFile(TOutputRootFile,"RECREATE");
   
    TTree *tree = new TTree("Truth_Events","Truth_Events");
   
    tree->Branch("w_gp", &w_gp); // weights of generated particles
   
    tree->Branch("e_px_gp", &e_px_gp); // momentum of generated particles (i.e. electrons)
    tree->Branch("e_py_gp", &e_py_gp);
    tree->Branch("e_pz_gp", &e_pz_gp);
    tree->Branch("e_E_gp", &e_E_gp);
   
    tree->Branch("k_px_gp", &k_px_gp); // momentum of generated particles (i.e. kaons)
    tree->Branch("k_py_gp", &k_py_gp);
    tree->Branch("k_pz_gp", &k_pz_gp);
    tree->Branch("k_E_gp", &k_E_gp);
   
    tree->Branch("l_px_gp", &l_px_gp); // momentum of generated particles (i.e. lambdas)
    tree->Branch("l_py_gp", &l_py_gp);
    tree->Branch("l_pz_gp", &l_pz_gp);
    tree->Branch("l_E_gp", &l_E_gp);
 
   
    for (int j=0; j<nEvents; j++) {
          
      w_gp = weight[j];
         
      e_px_gp = e_px[j];
      e_py_gp = e_py[j];
      e_pz_gp = e_pz[j];
      e_E_gp = e_E[j];
          
      k_px_gp = k_px[j];
      k_py_gp = k_py[j];
      k_pz_gp = k_pz[j];
      k_E_gp = k_E[j];
          
      l_px_gp = l_px[j];
      l_py_gp = l_py[j];
      l_pz_gp = l_pz[j];
      l_E_gp = l_E[j];
           
      tree->Fill(); 
    }
          
    //OutputFile->Print(); // 04/04/23 - Don't print to screen to prevent spam when running
    OutputFile->Write();
    OutputFile->Close();

    cout << "Created output root file - " << TOutputRootFile << " - sucessfully." <<endl;
    cout << "Bye!" << endl;
    exit;
  } 
} //-> Main void End
