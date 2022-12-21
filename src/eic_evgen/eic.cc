///*--------------------------------------------------*/
/// eic.cc:
/// Original author: Dr. Ahmed Zafar
/// Date: 2015-2018
///
///*--------------------------------------------------*/
/// Modifier: Wenliang (Bill) Li
/// Date: Feb 24 2020
/// Email: wenliang.billlee@gmail.com
///
/// Comment: Feb 24, 2020: the main function is excuted in main.cc


#include "eic.h"

using std::setw;
using std::setprecision;
using std::cout;
using std::cin;
using std::endl;
using namespace std;

//---------------------------------------------------------
// g++ -o pim pimFermi.C `root-config --cflags --glibs`
//---------------------------------------------------------

//int main() {
// 
//  eic();
//  
//  return 0;
//}

void eic() {

    Int_t target_direction, kinematics_type;
    Double_t EBeam, HBeam;
 
   	cout << "Target Direction (1->Up, 2->Down): "; cin >> target_direction; cout << endl;
   	cout << "Kinematics type (1->FF, 2->TSSA): ";  cin >> kinematics_type;  cout << endl;
   	cout << "Enter the number of events: ";        cin >> fNEvents;         cout << endl;
   	cout << "Enter the file number: ";             cin >> fNFile;           cout << endl;
	cout << "Enter the electron beam energy: ";    cin >> EBeam;            cout << endl;
	cout << "Enter the hadron beam energy: ";      cin >> HBeam;            cout << endl;
 
//	eic(target_direction, kinematics_type, fNEvents);

}

/*--------------------------------------------------*/

void eic(int event_number, int target_direction, int kinematics_type, TString file_name, int fEIC_seed, TString particle, TString hadron, TString det_location, TString OutputType, double EBeam, double HBeam) {

   	TString targetname;
	TString charge;

   	if( target_direction == 1 ) targetname = "up";
  	if( target_direction == 2 ) targetname = "down";
	
	gKinematics_type = kinematics_type;
	gfile_name = file_name;

	fNFile = 1;

	fNEvents = event_number;

	fSeed = fEIC_seed;
	cout << EBeam << " elec " << HBeam << " hadrons" << endl; 
	fEBeam = EBeam;
	fPBeam = HBeam;

	pim* myPim = new pim(fSeed);
  	myPim->Initilize();
	// 09/02/22 - SJDK - Special case for the kaon, if hadron not specified, default to Lambda
	if (particle == "K+"){
	  if (hadron != "Lambda" && hadron != "Sigma0"){
	    hadron = "Lambda";
	  }
	  else{
	    hadron = ExtractParticle(hadron);
	  }
	  Reaction* r1 = new Reaction(particle, hadron);
	  r1->process_reaction();
	  delete r1;
	}
	else if (particle == "pi+" || particle == "Pion+" ||  particle == "Pi+"){
	  hadron = "Neutron";
	  particle = ExtractParticle(particle);
	  charge = ExtractCharge(particle);
	  Reaction* r1 = new Reaction(particle, hadron);
	  r1->process_reaction();
	  delete r1;
	}
	else if (particle == "pi0" || particle == "Pion0" || particle == "Pi0"){
	  hadron = "Proton";
	  particle = ExtractParticle(particle);
	  charge = ExtractCharge(particle);
	  //Reaction* r1 = new Reaction(particle);
	  Reaction* r1 = new Reaction(particle, hadron);
	  r1->process_reaction();
	  delete r1;
	}
	else{
	  particle = ExtractParticle(particle);
	  charge = ExtractCharge(particle);
	  Reaction* r1 = new Reaction(particle);
	  r1->process_reaction();
	  delete r1;
	}
}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
// SJDK 21/12/22 - Note that this is the one that actually gets used, reads in the .json file
void eic(Json::Value obj) {

   	TString targetname;  
 	TString charge;

	int target_direction = obj["Targ_dir"].asInt();
 	gKinematics_type     = obj["Kinematics_type"].asInt();

   	if( target_direction == 1 ) targetname = "up";
   	if( target_direction == 2 ) targetname = "down";
 
 	gfile_name = obj["file_name"].asString();
 
 	gPi0_decay = obj["pi0_decay"].asBool();

 	fNFile = 1;
 	fNEvents = obj["n_events"].asUInt64();

 	fSeed = obj["generator_seed"].asInt();

 	pim* myPim = new pim(fSeed);
   	myPim->Initilize();
 
//  	TDatime dsTime;
//  	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;

	TString particle = obj["particle"].asString();
	TString hadron = obj["hadron"].asString(); // 09/02/22 - SJDK - Added in hadron type argument for K+
	// SJDK - 08/02/22 - This is terrible, need to change this, particle should just be K+
	// Add a new flag which, hadron - where this is specified too, then add conditionals elsewhere based on this
	// New conditional, special case for Kaon
	
	particle = ExtractParticle(particle);
	charge = ExtractCharge(particle);
	if (particle == "K+"){
	  if (hadron != "Lambda" && hadron != "Sigma0"){
	    hadron = "Lambda";
	  }
	  else{
	    hadron = ExtractParticle(hadron);
	  }
	}
	// SJDK - 19/12/22 - Specify hadron to neutron/proton for pi+/pi0 production, for pi0 production, may want to adjust? 
	else if (particle == "pi+" || particle == "Pion+" || particle == "Pi+"){
	  hadron = "Neutron";
	}
	else if (particle == "pi0" || particle == "Pion0" || particle == "Pi0"){
	  hadron = "Proton";
	}
	else { // SJDK -09/02/22 - Note that in future this could be changed to get different hadrons in other reactions if desired
	  hadron = "";
	}

	// SJDK - 01/06/21
	// Set beam energies from .json read in
	fEBeam = obj["ebeam"].asDouble();
	fPBeam = obj["hbeam"].asDouble();

	// SJDK - 12/01/22
	// Set output type as a .json read in
	// Should be Pythia6, LUND or HEPMC3
	gOutputType = obj["OutputType"].asString();
	if (gOutputType == "Pythia6"){
	  cout << "Using Pythia6 output format for Fun4All" << endl;
	}
	else if (gOutputType == "LUND"){
	  cout << "Using LUND output format" << endl;
	}
	else if (gOutputType == "HEPMC3"){
	  cout << "Using HEPMC3 output format for EPIC" << endl;
	}
	else{
	  cout << "Output type not recognised!" << endl;
	  cout << "Setting output type to Pythia6 by default!" << endl;
	  gOutputType = "Pythia6";
	}

	///*--------------------------------------------------*/
	/// The detector selection is determined here
	/// The incidence proton phi angle is 

	gDet_location = obj["det_location"].asString();

	if (gDet_location == "ip8") {

		fProton_incidence_phi = 0.0;

	} else if (gDet_location == "ip6") {

		fProton_incidence_phi = fPi;

	} else {
		fProton_incidence_phi = 0.0;
		cout << "The interaction point not recognized!" << endl;
		cout << "Therefore default opition ip6 is used." << endl;
	}

	fScatElec_E_Lo = obj["Ee_Low"].asDouble();
	fScatElec_E_Hi = obj["Ee_High"].asDouble();
	fScatElec_Theta_I = obj["e_Theta_Low"].asDouble() *  fDEG2RAD;
	fScatElec_Theta_F = obj["e_Theta_High"].asDouble() *  fDEG2RAD;
	fEjectileX_Theta_I = obj["EjectileX_Theta_Low"].asDouble() * fDEG2RAD;
	fEjectileX_Theta_F = obj["EjectileX_Theta_High"].asDouble() * fDEG2RAD;

	if(particle != "pi0"){ // Default case now
	  Reaction* r1 = new Reaction(particle, hadron);
	  r1->process_reaction();
	  delete r1;
	}
	else{  // Treat pi0 slightly differently for now
	  Reaction* r1 = new Reaction(particle);
	  r1->process_reaction();
	  delete r1;
	}

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void SetEICSeed (int seed) {
	fSeed = seed;
}

///*--------------------------------------------------*/
///*--------------------------------------------------*/
///  Some utility functions


///*--------------------------------------------------*/
/// Extracting the particle type

TString ExtractParticle(TString particle) {

	/// Make the input particle case insansitive
	particle.ToLower();
	if (particle.Contains("on")) {
		particle.ReplaceAll("on", "");
	};
 	
	if (particle.Contains("plus")) {
		particle.ReplaceAll("plus", "+");
	}

	if (particle.Contains("minus")) {
		particle.ReplaceAll("minus", "-");
	}

	if (particle.Contains("zero")) {
		particle.ReplaceAll("zero", "0");
	}

	particle[0] = toupper(particle[0]);
	cout << "Particle: " << particle << endl;
	return particle;

}

///*--------------------------------------------------*/
/// Extracting the particle charge

TString ExtractCharge(TString particle) {

	TString charge;

	if (particle.Contains("+") || particle.Contains("plus")) {
		charge = "+";
	} else if (particle.Contains("-") || particle.Contains("minus")) {
		charge = "-";
	} else {
		charge = "0";
	}
	return charge;
}
