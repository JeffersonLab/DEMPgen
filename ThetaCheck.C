// Stephen JD Kay, University of York, 21/02/25
// A short script to read in DEMPgen output and find the scattered electron theta ranges that (roughly) correspond to fixed Q2 ranges.
#include <string>

void ThetaCheck(TString infile=""){

  // If no input file provide as argument, promot for one
  if(infile == ""){
    cout << "Enter a filename to analyse: ";
    cin >> infile;
  }
  
  // Check input file exists, exit if not
  if(gSystem->AccessPathName(infile) == kTRUE){
    cerr << "!!!!! ERROR !!!!!" << endl << infile << " not found" << endl << "!!!!! ERROR !!!!!" << endl;
    exit(0);
  }

  TString ofile_name, ofile_tmp, BeamE;
  TObjArray *tmp_Name_arr;

  // Check the input file is a .root file as we would expect
  if(infile.Contains(".root") == false){
    cerr << "!!!!!!!!!!! - ERROR - !!!!!!!!!!!!" << endl;
    cerr << "Input files should be a root file!" << endl;
    cerr << "!!!!!!!!!!! - ERROR - !!!!!!!!!!!!" << endl;
    exit(1);
  }
  else{
    cout << "Opening and finding Q2 boundaries " << infile << endl;
  }

  // If the input file name is a file path containing /, extract only the actual file name for further use. Assign the temporary name of the output to be the input, minus its .root extension
  if(infile.Contains("/")){
    tmp_Name_arr = infile.Tokenize("/");
    ofile_tmp = (((TObjString *)(tmp_Name_arr->At(tmp_Name_arr->GetLast())))->String()).ReplaceAll(".root","");
  }
  else{
    ofile_tmp = infile;
    ofile_tmp.ReplaceAll(".root","");
  }

  TChain *mychain = new TChain("Events");
  mychain->Add(infile);

  // Initialise reader
  TTreeReader tree_reader(mychain);
  
  // Get particle info
  TTreeReaderArray<double> ScatElec_pX(tree_reader, "scat_e_px");
  TTreeReaderArray<double> ScatElec_pY(tree_reader, "scat_e_py");
  TTreeReaderArray<double> ScatElec_pZ(tree_reader, "scat_e_pz");
  //TTreeReaderArray<double> ScatElec_E(tree_reader, "scat_e_E");
  TTreeReaderArray<double> KinQ2(tree_reader, "Q2");

  double Theta, Q2;
  TVector3 ScatElec;
  
  // Set output file name
  ofile_name = Form("%s_ThetaQ2.root", ofile_tmp.Data());
  TH1D* H1_Theta = new TH1D("H1_Theta", "#theta_{e'}; #theta_{e'} (Deg)", 240, 120, 180);
  TH2D* H2_Q2_Theta = new TH2D("H2_Q2_Theta", "#theta_{e'}(Q^{2});Q^{2} (GeV^{2}); #theta_{e'} (Deg)", 400, 0, 40, 240, 120, 180); 
  TH2D* H2_Q2_Theta_R1 = new TH2D("H2_Q2_Theta_R1", "#theta_{e'}(Q^{2}), 3 < Q^{2} < 10;Q^{2} (GeV^{2}); #theta_{e'} (Deg)", 400, 0, 40, 240, 120, 180); 
  TH2D* H2_Q2_Theta_R2 = new TH2D("H2_Q2_Theta_R2", "#theta_{e'}(Q^{2}) 10 < Q^{2} < 20;Q^{2} (GeV^{2}); #theta_{e'} (Deg)", 400, 0, 40, 240, 120, 180); 
  TH2D* H2_Q2_Theta_R3 = new TH2D("H2_Q2_Theta_R3", "#theta_{e'}(Q^{2}) 20 < Q^{2} < 35;Q^{2} (GeV^{2}); #theta_{e'} (Deg)", 400, 0, 40, 240, 120, 180); 
  // Set and open output file for the histograms
  TFile *ofile = TFile::Open(ofile_name,"RECREATE");

  while (tree_reader.Next()){
    ScatElec.SetXYZ(ScatElec_pX[0], ScatElec_pY[0], ScatElec_pZ[0]);
    Theta=ScatElec.Theta()*TMath::RadToDeg();
    Q2 = KinQ2[0];
    H1_Theta->Fill(Theta);
    H2_Q2_Theta->Fill(Q2, Theta);
    if (Q2 > 3 && Q2 < 10){
      H2_Q2_Theta_R1->Fill(Q2, Theta);
    }
    else if (Q2 > 10 && Q2 < 20){
      H2_Q2_Theta_R2->Fill(Q2, Theta);
    }
    if (Q2 > 20 && Q2 < 35){
      H2_Q2_Theta_R3->Fill(Q2, Theta);
    }
    
  }

  H1_Theta->Write();
  H2_Q2_Theta->Write();
  H2_Q2_Theta_R1->Write();
  H2_Q2_Theta_R2->Write();
  H2_Q2_Theta_R3->Write();

  cout << "Min elec theta = " << H2_Q2_Theta_R1->GetYaxis()->GetBinLowEdge(H2_Q2_Theta_R1->FindFirstBinAbove(0., 2)) << " deg, max elec theta = " << H2_Q2_Theta_R1->GetYaxis()->GetBinLowEdge(H2_Q2_Theta_R1->FindLastBinAbove(0., 2)) << " deg for 3 < Q2 < 10" << endl;
  cout << "Min elec theta = " << H2_Q2_Theta_R2->GetYaxis()->GetBinLowEdge(H2_Q2_Theta_R2->FindFirstBinAbove(0., 2)) << " deg, max elec theta = " << H2_Q2_Theta_R2->GetYaxis()->GetBinLowEdge(H2_Q2_Theta_R2->FindLastBinAbove(0., 2)) << " deg for 10 < Q2 < 20" << endl;
  cout << "Min elec theta = " << H2_Q2_Theta_R3->GetYaxis()->GetBinLowEdge(H2_Q2_Theta_R3->FindFirstBinAbove(0., 2)) << " deg, max elec theta = " << H2_Q2_Theta_R3->GetYaxis()->GetBinLowEdge(H2_Q2_Theta_R3->FindLastBinAbove(0., 2)) << " deg for 20 < Q2 < 35" << endl;
  
  cout << H2_Q2_Theta_R1->GetEntries()/H2_Q2_Theta_R2->GetEntries() << " times more events in range 3 < Q2 < 10 than 10 < Q2 < 20" << endl;
  cout << H2_Q2_Theta_R1->GetEntries()/H2_Q2_Theta_R3->GetEntries() << " times more events in range 3 < Q2 < 10 than 20 < Q2 < 35" << endl;
  
  ofile->Write();
  ofile->Close();
  delete ofile;
   
}
