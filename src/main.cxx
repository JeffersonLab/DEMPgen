///*--------------------------------------------------*/
/// main.cxx:
///*--------------------------------------------------*/
/// Old Generator name: DEMP_EvGen
/// Original author: Rory Evens
/// Date: 2015-2018
///
///*--------------------------------------------------*/
/// New Generator name: DEMPgen
/// Modifier: Wenliang (Bill) Li
/// Date: Feb 24 2020
/// Email: wenliang.billlee@gmail.com
///
/// Comment: March 04, 2020: 

#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

//C++lib includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "json/json.h"
#include "json/json-forwards.h"

//Root includes
#include "TRandom.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TMath.h"
#include "TVector3.h"

//Project includes
#include "Particle.hxx"
#include "CustomRand.hxx"
#include "ScatteredParticleGen.hxx"
#include "ProductGen.hxx"
#include "TreeBuilder.hxx"
#include "Constants.hxx"
#include "DEMPEvent.hxx"
#include "SigmaCalc.hxx"
#include "TargetGen.hxx"
#include "MatterEffects.hxx"
#include "FSI.hxx"

#include "eic_evgen/eic.h"

using namespace std;
using namespace constants;

TFile * WorkFile;

Json::Value obj; //Declared here for global access

string get_date(void);

int Gen_seed;

int main(int argc, char** argv){
 
   // Parsing of config file.
   //ifstream ifs("../Config.json");
   ifstream ifs(argv[1]);
   Json::Reader reader;
   reader.parse(ifs, obj);
 
   unsigned long long int nEvents = obj["n_events"].asUInt64();
   cout << "Generating "<< nEvents << " events."<<endl;

   TString file_name = obj["file_name"].asString();

   if (file_name == "") {
	 file_name.Form("%i", nEvents);	
	 file_name = file_name + "_" + get_date();
   } 

   int gen_seed = obj["generator_seed"].asInt();
   TString particle = obj["particle"].asString();

   //TString HBeamPart = obj["hbeam_part"].asString(); // Work in progress

//	cout << obj["experiment"].asString() << endl;
//	cout << "File name: " << file_name << "  " << get_date() << endl;
//	cout << particle << endl;
//	exit(0);

   if (obj["experiment"].asString() == "eic") {
 
		cout << "EIC is used " << endl;
 
  	 	int target_direction = obj["Targ_dir"].asInt();
  		int kinematics_type = obj["Kinematics_type"].asInt();
		double EBeam = obj["ebeam"].asDouble();
		double HBeam = obj["hbeam"].asDouble();
		TString hadron = obj["hadron"].asString(); // SJDK 08/02/22 - Add the hadron type as an argument
//		bool = obj["pi0_particle"].asBool()
 		eic(obj);
 
   } else if (obj["experiment"].asString() == "solid") {
 
     Gen_seed = gen_seed;

     if (obj["ionisation"].asBool())
       cout << "Ionisation Enabled" << endl;
     if (obj["bremsstrahlung"].asBool())
       cout << "Bremsstrahlung Enabled" << endl;
     if (obj["fermi_momentum"].asBool())
       cout << "Fermi Momentum Enabled" << endl;
     if (obj["multiple_scattering"].asBool())
       cout << "Multiple Scattering Enabled" << endl;
     if (obj["final_state_interaction"].asBool())
       cout << "FSI Enabled" << endl;
 
   
     MatterEffects* ME = new MatterEffects();
      
     WorkFile = new TFile("../data/output/test.root");
   
     // Initilization of DEMPEvent objects for different reference frames
     DEMPEvent* VertEvent = new DEMPEvent("Vert");
     // VertEvent contains particles with kinematic properties at the vertex.
     // VertEvent is viewed in the laboratory reference frame.
     // This event is not to be modified after initial generation of the event.
     DEMPEvent* RestEvent = new DEMPEvent("RF");
     // Event in rest frame of the target neutron.
     DEMPEvent* CofMEvent = new DEMPEvent("CofM");
     // Event in the rest frame of the system's center of mass.
     DEMPEvent* TConEvent = new DEMPEvent("TCon"); 
     // Event rotated to follow the Trento Conventions
     DEMPEvent* LCorEvent = new DEMPEvent("Lab");
     // LCorEvent is the event after all corrections and effects have been applied,
     // in the laboratory reference frame.
   
     SigmaCalc* Sig = new SigmaCalc(VertEvent, CofMEvent, RestEvent, TConEvent);
   
     // Retrieval of pointers to particles in VertEvent.
     // For clarity: these are the same objects as are referenced by VertEvent,
     // not copies. Operations on these objects affect the original. 
     Particle* VertBeamElec = VertEvent->BeamElec;
     Particle* VertTargNeut = VertEvent->TargNeut;
     Particle* VertScatElec = VertEvent->ScatElec;
     Particle* VertProdPion = VertEvent->ProdPion;
     Particle* VertProdProt = VertEvent->ProdProt;
   
     Particle* Photon = VertEvent->VirtPhot;
     Photon->SetName("VirtPhot");

//     cout << "aaa  " << VertScatElec->Px() << "  " << VertScatElec->Py() << "  " 
//        << VertScatElec->Pz() << "  " << VertScatElec->E() << "  " << VertScatElec->GetMass() << endl;
   
     Particle* FSIProt = new Particle(proton_mass_mev, "FSIProt", pid_prot);
   
     Particle* InTotal = new Particle();
     Particle* OutTotal = new Particle();
   
     VertBeamElec->SetThetaPhiE(0, 0, obj["beam_energy"].asDouble());
   
     double elecERange[2] = {obj["scat_elec_Emin"].asDouble(),
                             obj["scat_elec_Emax"].asDouble()};
     double elecThetaRange[2] = {obj["scat_elec_thetamin"].asDouble()/DEG,
                                 obj["scat_elec_thetamax"].asDouble()/DEG};
     double elecPhiRange[2] = {0, 360/DEG};
   
     TargetGen * NeutGen = new TargetGen(neutron_mass_mev, obj["fermi_momentum"].asBool());
   
     ScatteredParticleGen * ElecGen =
     new ScatteredParticleGen(electron_mass_mev,
                                elecERange,
                                elecThetaRange,
                                elecPhiRange);
   
     FSI* FSIobj = new FSI();
   
     /*
       Particle * Photon = new Particle();
       VertEvent->VirtPhot = Photon;
     */

     ProductGen * ProtonPionGen = new ProductGen(Photon,
                                                 VertTargNeut);
   
     int nSuccess = 0;
     int nFail = 0;
     int nNeg = 0;
     int nCut = 0;
   
     int FSIfail = 0;
   
     int event_status = 0;
   

	 file_name = "RootFiles/Solid_DEMP_" + file_name + ".root";

     TreeBuilder * Output = new TreeBuilder(file_name.Data(), "t1");
   
     Output->AddEvent(VertEvent);
     //Output->AddEvent(CofMEvent);
     //Output->AddEvent(RestEvent);
     Output->AddEvent(LCorEvent);
   
     //if (obj["final_state_interaction"].asBool())
     Output->AddParticle(FSIProt);
   
     Output->AddParticle(Photon);
   
     // These parameters are calculated using multiple reference frames (DEMPEvent objects),
     // and need to be added to the output seperately.
     double sigma_l;
     double sigma_t;
     double sigma_tt;
     double sigma_lt;
   
     double sigma_uu;
     double sigma_ut;
   
     double sigma_k0;
     double sigma_k1;
     double sigma_k2;
     double sigma_k3;
     double sigma_k4;
   
     double sigma_k5 = 0;
   
     double sigma;
   
     double weight;
     double epsilon;
   
     double targetthickness, airthickness, targwindowthickness;
   
     targwindowthickness = 0.0120*Window_Density / Window_Thickness;
   
     Output -> AddDouble(&sigma_l,"sigma_l");
     Output -> AddDouble(&sigma_t,"sigma_t");
     Output -> AddDouble(&sigma_tt,"sigma_tt");
     Output -> AddDouble(&sigma_lt,"sigma_lt");
     Output -> AddDouble(&sigma_uu,"sigma_uu");
     Output -> AddDouble(&sigma_ut,"sigma_ut");
   
     Output -> AddDouble(&sigma_k0,"AutPhiMinusPhiS");
     Output -> AddDouble(&sigma_k1,"AutPhiS");
     Output -> AddDouble(&sigma_k2,"Aut2PhiMinusPhiS");
     Output -> AddDouble(&sigma_k3,"AutPhiPlusPhiS");
     Output -> AddDouble(&sigma_k4,"Aut3PhiMinusPhiS");
     Output -> AddDouble(&sigma_k5,"Aut2PhiPlusPhiS");
   
     Output -> AddDouble(&sigma, "sigma");
   
     Output -> AddDouble(&weight,"weight");
     Output -> AddDouble(&epsilon, "epsilon");
   
     //if(obj["final_state_interaction"].asBool()){
     Output -> AddDouble(FSIobj->PhaseShiftWeight, "PhaseShiftWeight");
     Output -> AddDouble(FSIobj->WilliamsWeight, "WilliamsWeight");
     Output -> AddDouble(FSIobj->DedrickWeight, "DedrickWeight");
     Output -> AddDouble(FSIobj->CatchenWeight, "CatchenWeight");
       //}
   
     // Event global variables to add to output
     Output -> AddDouble(VertEvent->qsq_GeV, "Vert_Qsq_GeV");
     Output -> AddDouble(VertEvent->w_GeV, "Vert_w_GeV");
     Output -> AddDouble(VertEvent->t_GeV, "Vert_t_GeV");
     Output -> AddDouble(VertEvent->t_para_GeV, "Vert_t_para_GeV");
     Output -> AddDouble(VertEvent->t_prime_GeV, "Vert_t_prime_GeV");
     Output -> AddDouble(VertEvent->negt, "Vert_negt_GeV");
     Output -> AddDouble(VertEvent->x_B, "Vert_x_B");
     Output -> AddDouble(VertEvent->Phi_deg, "Vert_Phi");
     Output -> AddDouble(VertEvent->Phi_s_deg, "Vert_Phi_s");
     Output -> AddDouble(VertEvent->Theta_deg, "Vert_Theta");
   
     Output -> AddDouble(LCorEvent->qsq_GeV, "Lab_Qsq_GeV");
     Output -> AddDouble(LCorEvent->w_GeV, "Lab_w_GeV");
     Output -> AddDouble(LCorEvent->t_GeV, "Lab_t_GeV");
     Output -> AddDouble(LCorEvent->t_para_GeV, "Lab_t_para_GeV");
     Output -> AddDouble(LCorEvent->t_prime_GeV, "Lab_t_prime_GeV");
     Output -> AddDouble(LCorEvent->negt, "Lab_negt_GeV");
     Output -> AddDouble(LCorEvent->x_B, "Lab_x_B");
     Output -> AddDouble(LCorEvent->Phi_deg, "Lab_Phi");
     Output -> AddDouble(LCorEvent->Phi_s_deg, "Lab_Phi_s");
     Output -> AddDouble(LCorEvent->Theta_deg, "Lab_Theta");
   
     Output -> AddDouble(VertEvent->Vertex_x, "Vertex_x");
     Output -> AddDouble(VertEvent->Vertex_y, "Vertex_y");
     Output -> AddDouble(VertEvent->Vertex_z, "Vertex_z");
   
     cout << "Starting Main Loop." << endl;
     // Main loop of the generator
     for (int i=0; i<nEvents; i++){
   
       if (i%100 == 0)
         cout <<"\r"<< (double)i/nEvents * 100 << "%";
   
       // Generate vertex location
       *VertEvent->Vertex_x = gRandom->Uniform(-0.25, 0.25);
       *VertEvent->Vertex_y = gRandom->Uniform(-0.25,0.25);
       *VertEvent->Vertex_z = gRandom->Uniform(-370,-330);
   
       // Reset inc. electron to beam energy
   
       VertBeamElec->SetThetaPhiE(0, 0, obj["beam_energy"].asDouble());
   
       // Perform matter effects on inc electron
       // These effects occur before the reaction, so affect the vertex values
       targetthickness = ((*VertEvent->Vertex_z+370.0) * Helium_Density)/
         (ME->X0(Helium_Z, Helium_A));
       if (obj["ionisation"].asBool()){
         ME->IonLoss(VertBeamElec, Window_A, Window_Z, Window_Density, targwindowthickness);
         ME->IonLoss(VertBeamElec, Helium_A, Helium_Z, Helium_Density, targetthickness);
       }
       if (obj["bremsstrahlung"].asBool()){
         ME->BremLoss(VertBeamElec,
                      targetthickness*ME->b(Helium_Z)
                      /ME->X0(Helium_Z, Helium_A));
         ME->BremLoss(VertBeamElec,
                      targwindowthickness*ME->b(Window_Z)
                      /ME->X0(Window_Z, Window_A));
       }
       if (obj["multiple_scattering"].asBool()){
         ME->MultiScatter(VertBeamElec,
                          targetthickness / ME->X0(Helium_Z, Helium_A));
         ME->MultiScatter(VertBeamElec,
                          targwindowthickness / ME->X0(Window_Z, Window_A));
       }
   
       // Generate target and scattered electron
       *VertTargNeut = *NeutGen->GetParticle();
       *VertScatElec = *ElecGen->GetParticle();

		/*--------------------------------------------------*/ 
        /// Test only 
//        VertScatElec->Px() = 15.934; 
//        VertScatElec->Py() = 1106.06; 
//        VertScatElec->Pz() = 2281.09; 
//        VertScatElec->E()  = 2535.16; 

        VertScatElec->SetPxPyPzE(15.934, 1106.06, 2281.09, 2535.16);

       *Photon = *VertBeamElec - *VertScatElec;

       cout << "              " << VertBeamElec->Px() << "  " << VertBeamElec->Py() << "  " << VertBeamElec->Pz() << "  " << VertBeamElec->E() << "  " << VertBeamElec->GetMass() << endl;
 
       cout << "              " << VertScatElec->Px() << "  " << VertScatElec->Py() << "  " << VertScatElec->Pz() << "  " << VertScatElec->E() << "  " << VertScatElec->GetMass() << endl;

       cout << "asdasdabbbbb  " << Photon->Px() << "  " << Photon->Py() << "  " 
       << Photon->Pz() << "  " << Photon->E() << "  " << Photon->GetMass() << endl; 



  
       // Solve for remaining particles
       event_status = ProtonPionGen->Solve();
       if (event_status == 0)
         nSuccess ++;
       if (event_status == 1){
         nFail ++;
         continue;
       }
       *VertProdPion = *ProtonPionGen->ProdPion();
       *VertProdProt = *ProtonPionGen->ProdProton();
       //    cout<<VertProdPion->GetPid() << endl;
      
       VertEvent->Update();
       // VertEvent and its components are not to be modified beyond this point.
   
       // Copy VertEvent and boost to CofM frame.
       *CofMEvent = *VertEvent;
       CofMEvent->Boost(-VertEvent->CoM());
       CofMEvent->Update();
   
       // Copy VertEvent and boost to rest frame
       *RestEvent = *VertEvent;
       RestEvent->Boost(-(VertEvent->TargNeut->Vect()*(1/VertEvent->TargNeut->E())));
       RestEvent->Update();
   
       // Copy RestEvent and rotate to Trento coords.
       *TConEvent = *RestEvent;
       TConEvent->Rotate(RestEvent->VirtPhot->Theta(),-RestEvent->ScatElec->Phi());
       TConEvent->Update();
   
       // Copy VertEvent onto LCorEvent, revert beamelec back to beam value
       *LCorEvent = *VertEvent;
       LCorEvent->BeamElec->SetThetaPhiE(0, 0, obj["beam_energy"].asDouble());
       LCorEvent->Update();
   
       // Get cross-sections
       sigma_l = Sig->sigma_l();
       sigma_t = Sig->sigma_t();
       sigma_lt = Sig->sigma_lt();
       sigma_tt = Sig->sigma_tt();
       sigma_uu = Sig->sigma_uu();
       sigma_ut = Sig->sigma_ut();
   
       sigma_k0 = Sig->Sigma_k(0);
       sigma_k1 = Sig->Sigma_k(1);
       sigma_k2 = Sig->Sigma_k(2);
       sigma_k3 = Sig->Sigma_k(3);
       sigma_k4 = Sig->Sigma_k(4);
   
       sigma = Sig->sigma();
   
       epsilon = Sig->epsilon();

       // Cuts
   
       if (obj["Qsq_cut"].asBool()){
           if (*VertEvent->qsq_GeV<obj["Qsq_min"].asDouble()){
             nSuccess --;
             nCut++;
             continue;
           }
         }
   
       if (obj["w_cut"].asBool()){
           if (*VertEvent->w_GeV<obj["w_min"].asDouble()){
             nSuccess --;
             nCut++;
             continue;
           }
         }
   
       if (obj["t_cut"].asBool()){
           if (*VertEvent->t_GeV<obj["t_min"].asDouble()){
             nSuccess --;
             nCut++;
             continue;
           }
         }
   
       if (sigma<=0){
         //cout << "NEGATIVE EVENT" << endl;
         nNeg++;
       }
   
       if (obj["weight_cut"].asBool()){
         if (sigma<=0) {
           nSuccess --;
           // cout << "NEGATIVE EVENT" << endl;
           // cout << "Sigma_l: \t" << sigma_l << endl;
           // cout << "Sigma_t: \t" << sigma_t << endl;
           // cout << "Sigma_lt: \t" << sigma_lt << endl;
           // cout << "Sigma_tt: \t" << sigma_tt << endl;
           // cout << "Sigma_uu: \t" << sigma_uu << endl;
           // cout << "Sigma_ut: \t" << sigma_ut << endl;
           // cout << "Sigma: \t" << sigma << endl;
           continue;
         }
       }
      
       weight = Sig->weight(nEvents);
   
       // Final State Interaction.
       
       if (obj["final_state_interaction"].asBool()){
         *FSIobj->VertInPion = *VertProdPion;
         FSIfail = FSIobj->Generate();
         if (FSIfail == 1){
           cout << "FSI Generation Failure!" << endl;
           // No instances of this so far.
           continue;
         }
         *FSIProt = *FSIobj->VertOutProt;
         *LCorEvent->ProdPion = *FSIobj->VertOutPion;
         //if (FSIfail == 0) cout << "FSI Generation Successful" << endl;
         FSIobj -> CalculateWeights();
       }
   
       // Final Cons Law Check
   
       *InTotal = (*(VertEvent->BeamElec)+
                   *(VertEvent->TargNeut)
                   );
       *OutTotal = (*(VertEvent->ScatElec)+
                    *(LCorEvent->ProdPion)+
                    *(VertEvent->ProdProt)
                    );
       if (obj["final_state_interaction"].asBool()){
         *InTotal += *(FSIobj->VertTargProt);
         *OutTotal += *(FSIobj->VertOutProt);
       }
   
       if ((*InTotal-*OutTotal).Px() > 1.0)
         cout << "Px Violation" << endl;
       if ((*InTotal-*OutTotal).Py() > 1.0)
         cout << "Px Violation" << endl;
       if ((*InTotal-*OutTotal).Pz() > 1.0)
         cout << "Px Violation" << endl;
       if ((*InTotal-*OutTotal).E() > 1.0)
         cout << "Px Violation" << endl;
   
       //Matter Effects~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
       targetthickness = ((-330.0 - *VertEvent->Vertex_z) * Helium_Density)/
         (ME->X0(Helium_Z, Helium_A));
       // Assuming path along z only inside the target.
       // Shape of target is not really being considered here.
       // A similar approximation is made for through-air distances.
       // Exact geometry of the detector is not considered, only large angle
       // or small angle.
       // If more detailed analysis is desired, geant can do it better.
       // Turn off effects and use GEMC.
   
   
       //cout << targetthickness << endl;
   
       if (obj["ionisation"].asBool()){
   
         //Scattered Electron
         ME->IonLoss(LCorEvent->ScatElec, Window_A, Window_Z,
                     Window_Density, targwindowthickness);
         ME->IonLoss(LCorEvent->ScatElec, Helium_A, Helium_Z,
                     Helium_Density, targetthickness);
         airthickness = 450 * Air_Density;
         if (LCorEvent->ScatElec->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density;
         ME->IonLoss(LCorEvent->ScatElec, Air_A, Air_Z,
                     Air_Density, airthickness);
   
         //Proton
         airthickness = 450 * Air_Density;
         if (LCorEvent->ProdProt->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density;
   
         ME->IonLoss(LCorEvent->ProdProt, Helium_A, Helium_Z,
                     Helium_Density, targetthickness);
   
         ME->IonLoss(LCorEvent->ProdProt, Window_A, Window_Z,
                     Window_Density, targwindowthickness);
   
         ME->IonLoss(LCorEvent->ProdProt, Air_A, Air_Z,
                     Air_Density, airthickness);
   
         //Pion
         airthickness = 450 * Air_Density;
         if (LCorEvent->ProdPion->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density;
   
         ME->IonLoss(LCorEvent->ProdPion, Helium_A, Helium_Z,
                     Helium_Density, targetthickness);
   
         ME->IonLoss(LCorEvent->ProdPion, Window_A, Window_Z,
                     Window_Density, targwindowthickness);
   
         ME->IonLoss(LCorEvent->ProdPion, Air_A, Air_Z,
                     Air_Density, airthickness);
   
       }
   
       if (obj["bremsstrahlung"].asBool()){
   
         //Scattered Electron
         ME->BremLoss(LCorEvent->ScatElec,
                      targetthickness*ME->b(Helium_Z)
                      /ME->X0(Helium_Z, Helium_A));
         ME->BremLoss(LCorEvent->ScatElec,
                      targwindowthickness*ME->b(Window_Z)
                      /ME->X0(Window_Z, Window_A));
         airthickness = 450 * Air_Density / ME->X0(Air_Z, Air_A);
         if (LCorEvent->ScatElec->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density / ME->X0(Air_Z, Air_A);
         ME->BremLoss(LCorEvent->ScatElec, airthickness*ME->b(Air_Z));
   
       }
   
       if (obj["multiple_scattering"].asBool()){
         //Scattered Electron
         ME->MultiScatter(LCorEvent->ScatElec,
                          targetthickness / ME->X0(Helium_Z, Helium_A));
         ME->MultiScatter(LCorEvent->ScatElec,
                          targwindowthickness / ME->X0(Window_Z, Window_A));
         airthickness = 450 * Air_Density / ME->X0(Air_Z,Air_A);
         if (LCorEvent->ScatElec->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density / ME->X0(Air_Z,Air_A);
         ME->MultiScatter(LCorEvent->ScatElec, airthickness);
   
         //Proton
         ME->MultiScatter(LCorEvent->ProdProt,
                          targetthickness / ME->X0(Helium_Z, Helium_A));
         ME->MultiScatter(LCorEvent->ProdProt,
                          targwindowthickness / ME->X0(Window_Z, Window_A));
         airthickness = 450 * Air_Density / ME->X0(Air_Z,Air_A);
         if (LCorEvent->ProdProt->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density / ME->X0(Air_Z,Air_A);
         ME->MultiScatter(LCorEvent->ProdProt, airthickness);
   
         //Pion
         ME->MultiScatter(LCorEvent->ProdPion,
                          targetthickness / ME->X0(Helium_Z, Helium_A));
         ME->MultiScatter(LCorEvent->ProdPion,
                          targwindowthickness / ME->X0(Window_Z, Window_A));
         airthickness = 450 * Air_Density / ME->X0(Air_Z,Air_A);
         if (LCorEvent->ProdPion->Theta() < 16) // Large Angle
           airthickness = 950 * Air_Density / ME->X0(Air_Z,Air_A);
         ME->MultiScatter(LCorEvent->ProdPion, airthickness);
   
       }
      
       Output->Fill();
   
     }
   
     cout << endl;
   
     if (nSuccess>0)
       Output->Save();
   
     cout << "Successful Events: \t" << nSuccess << endl;
     cout << "Failed Events: \t\t" << nFail << endl;
     cout << "Negative Events: \t\t" << nNeg << endl;
     cout << "Cut Events: \t\t" << nCut << endl;

     // Checks against old event generator:
   
     if(nEvents <0){
   
       int tests = -nEvents;
       int N = 10000;
   
       TFile * Check = new TFile("RootFiles/DEMP_Ee_11_Events_10000_File_0.root");
       TTree * t1 = (TTree*)Check->Get("t1");
   
       cout << "Running Debug/Check Values" << endl;
   
       Double_t Epsilon, Qsq_GeV, T_GeV, W_GeV, x, y, z, WSq_GeV, EventWeight, PhaseShiftWeight, PhaseSpaceWeight;
       Double_t Qsq_Corrected_GeV, T_Corrected_GeV, W_Corrected_GeV;
   
       double WilliamsWeight, DedrickWeight, CatchenWeight;
   
       Double_t ScatElec_Energy_Col_GeV,ScatElec_MomX_Col_GeV,ScatElec_MomY_Col_GeV,ScatElec_MomZ_Col_GeV,ScatElec_Mom_Col_GeV;
       Double_t ScatElec_Phi_Col,ScatElec_Theta_Col;
   
       Double_t ScatElec_Corrected_Energy_Col_GeV,ScatElec_Corrected_MomX_Col_GeV,ScatElec_Corrected_MomY_Col_GeV,ScatElec_Corrected_MomZ_Col_GeV,ScatElec_Corrected_Mom_Col_GeV;
       Double_t ScatElec_Corrected_Phi_Col,ScatElec_Corrected_Theta_Col;
   
       Double_t Pion_FSI_Energy_Col_GeV,Pion_FSI_MomX_Col_GeV,Pion_FSI_MomY_Col_GeV,Pion_FSI_MomZ_Col_GeV, Pion_FSI_Mom_Col_GeV;
       Double_t Pion_FSI_Phi_Col, Pion_FSI_Theta_Col;
   
       Double_t Pion_Energy_Col_GeV,Pion_MomX_Col_GeV,Pion_MomY_Col_GeV,Pion_MomZ_Col_GeV, Pion_Mom_Col_GeV;
       Double_t Pion_Phi_Col, Pion_Theta_Col;
   
       Double_t Pion_Corrected_Energy_Col_GeV,Pion_Corrected_MomX_Col_GeV,Pion_Corrected_MomY_Col_GeV,Pion_Corrected_MomZ_Col_GeV, Pion_Corrected_Mom_Col_GeV;
       Double_t Pion_Corrected_Phi_Col, Pion_Corrected_Theta_Col;
   
       Double_t RecoilProton_Energy_Col_GeV, RecoilProton_MomX_Col_GeV, RecoilProton_MomY_Col_GeV, RecoilProton_MomZ_Col_GeV, RecoilProton_Mom_Col_GeV;
       Double_t RecoilProton_Phi_Col, RecoilProton_Theta_Col;
   
       Double_t RecoilProton_Corrected_Energy_Col_GeV, RecoilProton_Corrected_MomX_Col_GeV, RecoilProton_Corrected_MomY_Col_GeV, RecoilProton_Corrected_MomZ_Col_GeV;
       Double_t RecoilProton_Corrected_Phi_Col, RecoilProton_Corrected_Theta_Col, RecoilProton_Corrected_Mom_Col_GeV;
   
       Double_t AsymPhiMinusPhi_S, AsymPhi_S, Asym2PhiMinusPhi_S, AsymPhiPlusPhi_S, Asym3PhiMinusPhi_S, Asym2PhiPlusPhi_S;
   
       double T_Para_GeV, T_Prime_GeV;
   
       t1->SetBranchAddress("Epsilon", &Epsilon );
       t1->SetBranchAddress("Qsq_GeV", &Qsq_GeV );
       t1->SetBranchAddress("T_GeV", &T_GeV );
       t1->SetBranchAddress("W_GeV", &W_GeV );
       t1->SetBranchAddress("T_Para_GeV", &T_Para_GeV);
       //t1->SetBranchAddress("T_Prime_GeV", &T_Prime_GeV);
   
       t1->SetBranchAddress("Qsq_Corrected_GeV", &Qsq_Corrected_GeV );
       t1->SetBranchAddress("T_Corrected_GeV", &T_Corrected_GeV );
       t1->SetBranchAddress("W_Corrected_GeV", &W_Corrected_GeV );
   
       t1->SetBranchAddress("ScatElec_Energy_Col_GeV", &ScatElec_Energy_Col_GeV );
       t1->SetBranchAddress("ScatElec_MomX_Col_GeV", &ScatElec_MomX_Col_GeV );
       t1->SetBranchAddress("ScatElec_MomY_Col_GeV", &ScatElec_MomY_Col_GeV );
       t1->SetBranchAddress("ScatElec_MomZ_Col_GeV", &ScatElec_MomZ_Col_GeV );
       t1->SetBranchAddress("ScatElec_Mom_Col_GeV", &ScatElec_Mom_Col_GeV );
       t1->SetBranchAddress("ScatElec_Theta_Col", &ScatElec_Theta_Col );
       t1->SetBranchAddress("ScatElec_Phi_Col", &ScatElec_Phi_Col );
   
       t1->SetBranchAddress("ScatElec_Corrected_Energy_Col_GeV", &ScatElec_Corrected_Energy_Col_GeV );
       t1->SetBranchAddress("ScatElec_Corrected_MomX_Col_GeV", &ScatElec_Corrected_MomX_Col_GeV );
       t1->SetBranchAddress("ScatElec_Corrected_MomY_Col_GeV", &ScatElec_Corrected_MomY_Col_GeV );
       t1->SetBranchAddress("ScatElec_Corrected_MomZ_Col_GeV", &ScatElec_Corrected_MomZ_Col_GeV );
       t1->SetBranchAddress("ScatElec_Corrected_Mom_Col_GeV", &ScatElec_Corrected_Mom_Col_GeV );
       t1->SetBranchAddress("ScatElec_Corrected_Theta_Col", &ScatElec_Corrected_Theta_Col );
       t1->SetBranchAddress("ScatElec_Corrected_Phi_Col", &ScatElec_Corrected_Phi_Col );
   
       // t1->SetBranchAddress("Pion_FSI_Energy_Col_GeV", &Pion_FSI_Energy_Col_GeV );
       // t1->SetBranchAddress("Pion_FSI_MomX_Col_GeV", &Pion_FSI_MomX_Col_GeV );
       // t1->SetBranchAddress("Pion_FSI_MomY_Col_GeV", &Pion_FSI_MomY_Col_GeV );
       // t1->SetBranchAddress("Pion_FSI_MomZ_Col_GeV", &Pion_FSI_MomZ_Col_GeV );
       // t1->SetBranchAddress("Pion_FSI_Mom_Col_GeV", &Pion_FSI_Mom_Col_GeV );
       // t1->SetBranchAddress("Pion_FSI_Theta_Col", &Pion_FSI_Theta_Col );
       // t1->SetBranchAddress("Pion_FSI_Phi_Col", &Pion_FSI_Phi_Col );
   
       t1->SetBranchAddress("Pion_Energy_Col_GeV", &Pion_Energy_Col_GeV );
       t1->SetBranchAddress("Pion_MomX_Col_GeV", &Pion_MomX_Col_GeV );
       t1->SetBranchAddress("Pion_MomY_Col_GeV", &Pion_MomY_Col_GeV );
       t1->SetBranchAddress("Pion_MomZ_Col_GeV", &Pion_MomZ_Col_GeV );
       t1->SetBranchAddress("Pion_Mom_Col_GeV", &Pion_Mom_Col_GeV );
       t1->SetBranchAddress("Pion_Theta_Col", &Pion_Theta_Col );
       t1->SetBranchAddress("Pion_Phi_Col", &Pion_Phi_Col );
   
       t1->SetBranchAddress("Pion_Corrected_Energy_Col_GeV", &Pion_Corrected_Energy_Col_GeV );
       t1->SetBranchAddress("Pion_Corrected_MomX_Col_GeV", &Pion_Corrected_MomX_Col_GeV );
       t1->SetBranchAddress("Pion_Corrected_MomY_Col_GeV", &Pion_Corrected_MomY_Col_GeV );
       t1->SetBranchAddress("Pion_Corrected_MomZ_Col_GeV", &Pion_Corrected_MomZ_Col_GeV );
       t1->SetBranchAddress("Pion_Corrected_Mom_Col_GeV", &Pion_Corrected_Mom_Col_GeV );
       t1->SetBranchAddress("Pion_Corrected_Theta_Col", &Pion_Corrected_Theta_Col );
       t1->SetBranchAddress("Pion_Corrected_Phi_Col", &Pion_Corrected_Phi_Col );
   
       t1->SetBranchAddress("RecoilProton_Energy_Col_GeV", &RecoilProton_Energy_Col_GeV );
       t1->SetBranchAddress("RecoilProton_MomX_Col_GeV", &RecoilProton_MomX_Col_GeV );
       t1->SetBranchAddress("RecoilProton_MomY_Col_GeV", &RecoilProton_MomY_Col_GeV );
       t1->SetBranchAddress("RecoilProton_MomZ_Col_GeV", &RecoilProton_MomZ_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Mom_Col_GeV", &RecoilProton_Mom_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Theta_Col", &RecoilProton_Theta_Col );
       t1->SetBranchAddress("RecoilProton_Phi_Col", &RecoilProton_Phi_Col );
   
       t1->SetBranchAddress("RecoilProton_Corrected_Energy_Col_GeV", &RecoilProton_Corrected_Energy_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Corrected_MomX_Col_GeV", &RecoilProton_Corrected_MomX_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Corrected_MomY_Col_GeV", &RecoilProton_Corrected_MomY_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Corrected_MomZ_Col_GeV", &RecoilProton_Corrected_MomZ_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Corrected_Mom_Col_GeV", &RecoilProton_Corrected_Mom_Col_GeV );
       t1->SetBranchAddress("RecoilProton_Corrected_Theta_Col", &RecoilProton_Corrected_Theta_Col );
       t1->SetBranchAddress("RecoilProton_Corrected_Phi_Col", &RecoilProton_Corrected_Phi_Col );
   
       t1->SetBranchAddress("EventWeight", &EventWeight );
       t1->SetBranchAddress("WilliamsWeight", &WilliamsWeight );
       t1->SetBranchAddress("DedrickWeight", &DedrickWeight );
       t1->SetBranchAddress("CatchenWeight", &CatchenWeight );
   
       double ZASig_T, ZASig_L, ZASig_LT, ZASig_TT, ZASigma_Lab, ZASigma_UU, RorySigma_UT;
       double Theta, Phi, Phi_S;
   
       t1->SetBranchAddress("ZASig_T", &ZASig_T);
       t1->SetBranchAddress("ZASig_L", &ZASig_L);
       t1->SetBranchAddress("ZASig_LT", &ZASig_LT);
       t1->SetBranchAddress("ZASig_TT", &ZASig_TT);
       t1->SetBranchAddress("ZASigma_Lab", &ZASigma_Lab);
       t1->SetBranchAddress("ZASigma_UU", &ZASigma_UU);
       t1->SetBranchAddress("RorySigma_UT", &RorySigma_UT);
   
       t1->SetBranchAddress("AsymPhiMinusPhi_S", &AsymPhiMinusPhi_S);
       t1->SetBranchAddress("AsymPhi_S", &AsymPhi_S);
       t1->SetBranchAddress("Asym2PhiMinusPhi_S", &Asym2PhiMinusPhi_S);
       t1->SetBranchAddress("AsymPhiPlusPhi_S", &AsymPhiPlusPhi_S);
       t1->SetBranchAddress("Asym3PhiMinusPhi_S", &Asym3PhiMinusPhi_S);
       t1->SetBranchAddress("Asym2PhiPlusPhi_S", &Asym2PhiPlusPhi_S);
   
       t1->SetBranchAddress("Phi", &Phi);
       t1->SetBranchAddress("PhiS", &Phi_S);
       t1->SetBranchAddress("Photon_Theta_Col", &Theta);
   
       double Jacobian_CM, Jacobian_CM_RF, Jacobian_CM_Col;
       double Flux_Factor_RF, Flux_Factor_Col, A;
   
       t1->SetBranchAddress("Jacobian_CM",&Jacobian_CM);
       t1->SetBranchAddress("Jacobian_CM_RF",&Jacobian_CM_RF);
       t1->SetBranchAddress("Jacobian_CM_Col",&Jacobian_CM_Col);
       t1->SetBranchAddress("A", &A);
   
       t1->SetBranchAddress("Flux_Factor_RF",&Flux_Factor_RF);
       t1->SetBranchAddress("Flux_Factor_Col",&Flux_Factor_Col);
   
       bool ALERT = false;
       for (int i=0; i<tests; i++){
   
         t1->GetEntry(i);
   
   
         VertScatElec->SetThetaPhiE(ScatElec_Theta_Col/DEG, ScatElec_Phi_Col/DEG,
                                    ScatElec_Energy_Col_GeV * 1000);
         *VertTargNeut = *NeutGen->GetParticle();
   
         /// Setting Photon
         *Photon = *VertBeamElec - *VertScatElec;
   
         ProtonPionGen->Solve(Pion_Theta_Col/DEG,Pion_Phi_Col/DEG);
   

         /// Setting Pion
         *VertProdPion = *ProtonPionGen->ProdPion();
         /// Setting Proton
         *VertProdProt = *ProtonPionGen->ProdProton();
   
         VertEvent->Update();
   
         *CofMEvent = *VertEvent;
         CofMEvent->Boost(-VertEvent->CoM());
         CofMEvent->Update();
   
         *RestEvent = *VertEvent;
         RestEvent->Boost(-(VertEvent->TargNeut->Vect()*(1/VertEvent->TargNeut->E())));
         RestEvent->Update();
   
         *TConEvent = *RestEvent;
         TConEvent->Rotate(RestEvent->VirtPhot->Theta(),-RestEvent->ScatElec->Phi());
         TConEvent->Update();
   
         double qsq_GeV = *VertEvent->qsq_GeV;
         double w_GeV = *VertEvent->w_GeV;
         double t_GeV = *VertEvent->t_GeV;
         double t_prime_GeV = *VertEvent->t_prime_GeV;
         double t_para_GeV = *VertEvent->t_para_GeV;
         double phi = *TConEvent->Phi;
         if (phi<0) phi+=2*TMath::Pi();
         double phi_s = *TConEvent->Phi_s;
         if (phi_s<0) phi_s+=2*TMath::Pi();
         double theta = *RestEvent->Theta;
         if (theta <0) theta+=2*TMath::Pi();
   
   
         sigma_l = Sig->sigma_l();
         sigma_t = Sig->sigma_t();
         sigma_lt = Sig->sigma_lt();
         sigma_tt = Sig->sigma_tt();
         sigma_uu = Sig->sigma_uu();
         sigma_ut = Sig->sigma_ut();
   
         sigma_k0 = Sig->Sigma_k(0);
         sigma_k1 = Sig->Sigma_k(1);
         sigma_k2 = Sig->Sigma_k(2);
         sigma_k3 = Sig->Sigma_k(3);
         sigma_k4 = Sig->Sigma_k(4);
   
         sigma = Sig->sigma();
   
         if (obj["final_state_interaction"].asBool()){
           *FSIobj->VertInPion = *VertProdPion;
           FSIobj->VertOutPion -> SetPxPyPzE(Pion_Corrected_MomX_Col_GeV*1000,
                                             Pion_Corrected_MomY_Col_GeV*1000,
                                             Pion_Corrected_MomZ_Col_GeV*1000,
                                             Pion_Corrected_Energy_Col_GeV*1000);
           //FSIfail = FSIobj->Generate();
           //if (FSIfail == 1){
           //  cout << "FSI Generation Failure!" << endl;
           //  continue;
           //}
           FSIobj->GenerateNoRandom();
           *FSIProt = *FSIobj->VertOutProt;
           *LCorEvent->ProdPion = *FSIobj->VertOutPion;
           FSIobj -> CalculateWeights();
         }
   
         int printw = 20;
   
         ALERT = false;
         if(TMath::Abs((sigma-ZASigma_Lab)/sigma)>0.01) {ALERT = true; cout << "SIGMA:\t"<<TMath::Abs((sigma-ZASigma_Lab)/sigma)<<endl;}
         if(TMath::Abs((sigma_l-ZASig_L)/sigma_l)>0.01) {ALERT = true; cout << "SIGMA_L:\t"<<TMath::Abs((sigma_l-ZASig_L)/sigma_l)<<endl;}
         if(TMath::Abs((sigma_t-ZASig_T)/sigma_t)>0.01) {ALERT = true; cout << "SIGMA_T:\t"<<TMath::Abs((sigma_t-ZASig_T)/sigma_t)<<endl;}
         if(TMath::Abs((sigma_lt-ZASig_LT)/sigma_lt)>0.01) {ALERT = true; cout << "SIGMA_LT:\t"<<TMath::Abs((sigma_lt-ZASig_LT)/sigma_lt)<<endl;}
         if(TMath::Abs((sigma_tt-ZASig_TT)/sigma_tt)>0.01) {ALERT = true; cout << "SIGMA_TT:\t"<<TMath::Abs((sigma_tt-ZASig_TT)/sigma_tt)<<endl;}
         if(TMath::Abs((sigma_uu-ZASigma_UU)/sigma_uu)>0.01) {ALERT = true; cout << "SIGMA_UU:\t"<<TMath::Abs((sigma_uu-ZASigma_UU)/sigma_uu)<<endl;}
         if(TMath::Abs((sigma_ut-RorySigma_UT)/sigma_ut)>0.01) {ALERT = true; cout << "SIGMA_UT:\t"<<TMath::Abs((sigma_ut-RorySigma_UT)/sigma_ut)<<endl;}
         if(TMath::Abs((Sig->weight(N)-EventWeight)/EventWeight)>0.01) {ALERT = true; cout << "WEIGHT:\t"<<TMath::Abs((Sig->weight(N)-EventWeight)/EventWeight)<<endl;}
   
         if(ALERT){
   
         cout<<left<<setw(printw)<<"Event:"<<left<<setw(printw)<<i<<endl;
         cout<<left<<setw(printw)<<"Kinematics:"<<endl;
   
         cout<<left<<setw(printw)<< "ElecE:"<<left<<setw(printw) << VertScatElec->E()/1000<<left<<setw(printw)<<ScatElec_Energy_Col_GeV << endl;
         cout<<left<<setw(printw)<< "ElecTh:"<<left<<setw(printw)<<VertScatElec->Theta()*DEG<<left<<setw(printw)<<ScatElec_Theta_Col << endl;
         cout<<left<<setw(printw)<< "ElecPhi:"<<left<<setw(printw)<<VertScatElec->Phi()*DEG <<left<<setw(printw)<<ScatElec_Phi_Col << endl;
   
         cout<<left<<setw(printw)<< "PionE:"<<left<<setw(printw)<<VertProdPion->E()/1000<<left<<setw(printw)<<Pion_Energy_Col_GeV << endl;
         cout<<left<<setw(printw)<< "PionTh:"<<left<<setw(printw)<<VertProdPion->Theta()*DEG<<left<<setw(printw)<<Pion_Theta_Col<< endl;
         cout<<left<<setw(printw)<< "PionPhi:"<<left<<setw(printw)<<VertProdPion->Phi()*DEG<<left<<setw(printw)<<Pion_Phi_Col<< endl;
   
         cout<<left<<setw(printw)<< "ProtE:"<<left<<setw(printw)<<VertProdProt->E()/1000<<left<<setw(printw)<<RecoilProton_Energy_Col_GeV << endl;
         cout<<left<<setw(printw)<< "ProtTh:"<<left<<setw(printw)<<VertProdProt->Theta()*DEG<<left<<setw(printw)<<RecoilProton_Theta_Col<< endl;
         cout<<left<<setw(printw)<< "ProtPhi:"<<left<<setw(printw)<<VertProdProt->Phi()*DEG<<left<<setw(printw)<<RecoilProton_Phi_Col << endl;
   
         cout<<left<<setw(printw)<< "Qsq:"<<left<<setw(printw) << qsq_GeV<<left<<setw(printw)<<Qsq_GeV << endl;
         cout<<left<<setw(printw)<< "W:"<<left<<setw(printw) << w_GeV<<left<<setw(printw)<<W_GeV << endl;
         cout<<left<<setw(printw)<< "t:"<<left<<setw(printw) << t_GeV<<left<<setw(printw)<<T_GeV << endl;
         cout<<left<<setw(printw)<< "t_para:"<<left<<setw(printw) << t_para_GeV<<left<<setw(printw)<<T_Para_GeV << endl;
         cout<<left<<setw(printw)<< "t':"<<left<<setw(printw) << t_prime_GeV<<left<<setw(printw)<<T_Prime_GeV << endl;
         cout<<left<<setw(printw)<<"Theta:"<<left<<setw(printw)<<theta*DEG<<left<<setw(printw)<<Theta<<endl;
         cout<<left<<setw(printw)<<"Phi"<<left<<setw(printw)<<phi*DEG<<left<<setw(printw)<<Phi<<endl;
         cout<<left<<setw(printw)<<"PhiS"<<left<<setw(printw)<<phi_s*DEG<<left<<setw(printw)<<Phi_S<<endl;
         cout<<left<<setw(printw)<< endl;
   
         cout<<left<<setw(printw)<<"Cross Sections: "<<endl;
         cout<<left<<setw(printw)<<"sigma:"<<left<<setw(printw)<< sigma<<left<<setw(printw)<<ZASigma_Lab << endl;
         cout<<left<<setw(printw)<<"sigma_l:"<<left<<setw(printw)<< sigma_l<<left<<setw(printw)<<ZASig_L << endl;
         cout<<left<<setw(printw)<<"sigma_t:"<<left<<setw(printw)<< sigma_t<<left<<setw(printw)<<ZASig_T << endl;
         cout<<left<<setw(printw)<<"sigma_lt:"<<left<<setw(printw)<< sigma_lt<<left<<setw(printw)<<ZASig_LT << endl;
         cout<<left<<setw(printw)<<"sigma_tt:"<<left<<setw(printw)<< sigma_tt<<left<<setw(printw)<<ZASig_TT << endl;
         cout<<left<<setw(printw)<<"sigma_uu:"<<left<<setw(printw)<< sigma_uu<<left<<setw(printw)<<ZASigma_UU << endl;
         cout<<left<<setw(printw)<<"sigma_ut:"<<left<<setw(printw)<< sigma_ut<<left<<setw(printw)<<RorySigma_UT << endl;
         cout<<left<<setw(printw)<<"AsyP-Ps:"<<left<<setw(printw)<< sigma_k0<<left<<setw(printw)<<AsymPhiMinusPhi_S << endl;
         cout<<left<<setw(printw)<<"AsyPs:"<<left<<setw(printw)<< sigma_k1<<left<<setw(printw)<<AsymPhi_S << endl;
         cout<<left<<setw(printw)<<"Asy2P-Ps:"<<left<<setw(printw)<< sigma_k2<<left<<setw(printw)<<Asym2PhiMinusPhi_S << endl;
         cout<<left<<setw(printw)<<"AsyP+Ps:"<<left<<setw(printw)<< sigma_k3<<left<<setw(printw)<<AsymPhiPlusPhi_S << endl;
         cout<<left<<setw(printw)<<"Asy3P-Ps:"<<left<<setw(printw)<< sigma_k4<<left<<setw(printw)<<Asym3PhiMinusPhi_S << endl;
         cout<<left<<setw(printw)<<"Asy2P_Ps:"<<left<<setw(printw)<< sigma_k5<<left<<setw(printw)<<Asym2PhiPlusPhi_S << endl;
         cout<<left<<setw(printw)<<"Epsilon:"<<left<<setw(printw)<< Sig->epsilon()<<left<<setw(printw)<<Epsilon << endl;
         cout<<left<<setw(printw)<<"Jacobian CM:"<<left<<setw(printw)<<Sig->jacobian_cm()<<left<<setw(printw)<<Jacobian_CM<<endl;
         cout<<left<<setw(printw)<<"Jacobian CM-Col:"<<left<<setw(printw)<<Sig->jacobian_cm_col()<<left<<setw(printw)<<Jacobian_CM_Col<<endl;
         cout<<left<<setw(printw)<<"Flux Factor Col:"<<left<<setw(printw)<<Sig->fluxfactor_col()<<left<<setw(printw)<<Flux_Factor_Col<<endl;
         cout<<left<<setw(printw)<<"A:"<<left<<setw(printw)<<Sig->jacobian_A()<<left<<setw(printw)<<A<<endl;
         cout<<left<<setw(printw)<<"Event Weight:"<<left<<setw(printw)<<Sig->weight(N)<<left<<setw(printw)<<EventWeight<<left<<setw(printw)<< endl;
         //cout<<left<<setw(printw)<<EventWeight/Sig->weight(100000) << endl;
         cout<<left<<setw(printw)<< endl;
   
         cout << (ZASigma_UU+RorySigma_UT)*Flux_Factor_Col*A*Jacobian_CM_Col << endl;
   
         cout<<left<<setw(printw)<< endl;
   
         // Sig->Asyms->at(0)->PrintPars();
   
         cout<<left<<setw(printw)<< endl;
         cout<<left<<setw(printw)<<"WilliamWeight:"<<left<<setw(printw)<<*FSIobj->WilliamsWeight<<left<<setw(printw)<<WilliamsWeight<<left<<setw(printw)<<endl;
         cout<<left<<setw(printw)<<"DedrickWeight:"<<left<<setw(printw)<<*FSIobj->DedrickWeight<<left<<setw(printw)<<DedrickWeight<<left<<setw(printw)<<endl;
         cout<<left<<setw(printw)<<"CatchenWeight:"<<left<<setw(printw)<<*FSIobj->CatchenWeight<<left<<setw(printw)<<CatchenWeight<<left<<setw(printw)<<endl;
   
   
         cout << "--------------------------------------------------------------------------------------------"<<endl;
         }
       }
     }
 
  	} else {
 
 	  cerr << endl;
       cerr << "/*--------------------------------------------------*/" << endl;
       cerr << "/*--------------------------------------------------*/" << endl;
       cerr << "/*                      Error!                      */" << endl;
       cerr << "/*            Not a valid experiment!               */" << endl;
       cerr << "/*            Option use EIC or solid!              */" << endl;
       cerr << "/*--------------------------------------------------*/" << endl;
       cerr << "/*--------------------------------------------------*/" << endl;
 	  cerr << endl;
 
 	  exit(0);
 
 	}
 
  return 0;
}



string get_date(void) {
 
  int MAX_DATE = 22;

  time_t now;
   char the_date[MAX_DATE];

   the_date[0] = '\0';

   now = time(NULL);

   if (now != -1)
   {
      strftime(the_date, MAX_DATE, "%Y_%h_%d_%H_%M_%S", localtime(&now));
   }

   return string(the_date);
}

	

