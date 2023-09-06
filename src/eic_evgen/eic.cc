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
using std::vector;
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
// 18/01/23 - SJDK- This function is never used since eic() is only called with a json object as the argument. Commented out for now, delete later?
/* void eic(int event_number, int target_direction, int kinematics_type, TString file_name, int fEIC_seed, TString particle, TString hadron, TString det_location, TString OutputType, double EBeam, double HBeam) {

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
*/

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
	// 21/12/22 - SJDK - Should do a check if these are defined or not, should crash if not defined or set defaults, see other quantities below
	TString particle = obj["particle"].asString();
	TString hadron = obj["hadron"].asString(); // 09/02/22 - SJDK - Added in hadron type argument for K+
	// SJDK - 08/02/22 - This is terrible, need to change this, particle should just be K+
	// Add a new flag which, hadron - where this is specified too, then add conditionals elsewhere based on this
	// New conditional, special case for Kaon	
	particle = ExtractParticle(particle);
	charge = ExtractCharge(particle);
	if(hadron == "Sigma" || hadron == "sigma"){ // SJDK - 31/01/23 - If hadron specified as Sigma, interpret this as Sigma0. Also correct for lower case
	  hadron = "Sigma0";
	}
	if (hadron == "lambda"){ // SJDK - 31/01/23 - Make Lambda selection case insensitive
	  hadron = "Lambda"; 
	}
	if (particle == "K+"){
	  if (hadron != "Lambda" && hadron != "Sigma0"){
	    hadron = "Lambda";
	    cout << "! WARNING !" << endl;
	    cout << "! WARNING !- K+ production specified but hadron not recognised, deaulting to Lambda - ! WARNING!" << endl;
	    cout << "! WARNING !" << endl;
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

	// SJDK 03/04/23 - Change to how Qsq range is set/chosen, could add as an override variable later too
	// Set min/max Qsq values depending upon particle type
	if (particle == "pi+" || particle == "Pion+" || particle == "Pi+"){
	  fQsq_Min = 3.0; fQsq_Max = 35.0;
	  fW_Min = 2.0; fW_Max = 10.2;
	}
	else if (particle == "pi0" || particle == "Pion0" || particle == "Pi0"){
	  fQsq_Min = 5.0; fQsq_Max = 1000.0;
	  fW_Min = 2.0; fW_Max = 10.0;
	}
	else if (particle == "K+"){
	  fQsq_Min = 1.0; fQsq_Max = 35.0;
	  fW_Min = 2.0; fW_Max = 10.0;
	}
	else{
	  fQsq_Min = 5.0; fQsq_Max = 35.0;
	  fW_Min = 2.0; fW_Max = 10.0;
	}

	// SJDK - 01/06/21
	// Set beam energies from .json read in
	if (obj.isMember("ebeam")){
	  fEBeam = obj["ebeam"].asDouble();
	}
	else{
	  fEBeam = 5;
	  cout << "Electron beam energy not specified in .json file, defaulting to 5 GeV." << endl;
	}
	if (obj.isMember("hbeam")){
	  fPBeam = obj["hbeam"].asDouble();
	}
	else{
	  fPBeam = 100;
	  cout << "Ion beam energy not specified in .json file, defaulting to 100 GeV." << endl;
	}

	if (obj.isMember("hbeam_part")){
	  gBeamPart = obj["hbeam_part"].asString();
	} 
    else {
      gBeamPart = "Proton";
    }

	// SJDK - 12/01/22
	// Set output type as a .json read in
	// Should be Pythia6, LUND or HEPMC3
	if (obj.isMember("OutputType")){
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
	    cout << "Setting output type to HEPMC3 by default!" << endl;
	    gOutputType = "HEPMC3";
	  }
	}
	else{
	  cout << "Output type not specified in .json file!" << endl;
	  cout << "Setting output type to HEPMC3 by default!" << endl;
	  gOutputType = "HEPMC3";
	}
	///*--------------------------------------------------*/
	/// The detector selection is determined here
	/// The incidence proton phi angle is 
	if (obj.isMember("det_location")){
	  gDet_location = obj["det_location"].asString();
	  if (gDet_location == "ip8") {
	    fProton_incidence_phi = 0.0;
	  } 
	  else if (gDet_location == "ip6") {
	    fProton_incidence_phi = fPi;
	  }
	  else {
	    fProton_incidence_phi = 0.0;
	    cout << "The interaction point requested is not recognized!" << endl;
	    cout << "Therefore ip6 is used by default." << endl;
	  }
	}
	else{ // 21/12/22 - This could probably be combined with the else statement above in some way
	    fProton_incidence_phi = 0.0;
	    cout << "The interaction points was not specified in the .json file!" << endl;
	    cout << "Therefore ip6 is used by default" << endl;
	}

	if (obj.isMember("Ee_Low")){
	  fScatElec_E_Lo = obj["Ee_Low"].asDouble();
	}
	else{
	  fScatElec_E_Lo = 0.5;
	  cout << "Minumum scattered electron energy not specified in .json file, defaulting to 0.5*EBeam." << endl;
	}

	if (obj.isMember("Ee_High")){
	  fScatElec_E_Hi = obj["Ee_High"].asDouble();
	}
	else{
	  fScatElec_E_Hi = 2.5;
	  cout << "Max scattered electron energy not specified in .json file, defaulting to 2.5*EBeam." << endl;
	}

	if (obj.isMember("e_Theta_Low")){
	  fScatElec_Theta_I = obj["e_Theta_Low"].asDouble() * fDEG2RAD;
	}
	else{
	  fScatElec_Theta_I = 60.0 * fDEG2RAD;
	  cout << "Min scattered electron theta not specified in .json file, defaulting to 60 degrees." << endl;
	}

	if (obj.isMember("e_Theta_High")){
	  fScatElec_Theta_F = obj["e_Theta_High"].asDouble() * fDEG2RAD;
	}
	else{
	  fScatElec_Theta_F = 175.0 * fDEG2RAD;
	  cout << "Max scattered electron theta not specified in .json file, defaulting to 175 degrees." << endl;
	}

	if (obj.isMember("EjectileX_Theta_Low")){
	  fEjectileX_Theta_I = obj["EjectileX_Theta_Low"].asDouble() * fDEG2RAD;
	}
	else{
	  fEjectileX_Theta_I = 0.0 * fDEG2RAD;
	  cout << "Min ejectile X theta not specified in .json file, defaulting to 0 degrees." << endl;
	}

	if (obj.isMember("EjectileX_Theta_High")){
	  fEjectileX_Theta_F = obj["EjectileX_Theta_High"].asDouble() * fDEG2RAD;
	}
	else{
	  fEjectileX_Theta_F = 60.0 * fDEG2RAD;
	  cout << "Max ejectile X theta not specified in .json file, defaulting to 60 degrees." << endl;
	}
        
	SigPar = ReadCrossSectionPar(particle, hadron);

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

vector<vector<vector<vector<double>>>> ReadCrossSectionPar(TString particle, TString hadron){
  
  string sigL_ParamFile, sigT_ParamFile;
 
  if (particle == "Pi+" && hadron == "Neutron"){
    //cout << "Add Pi+/Neutron case here" << endl;
  }
  else if (particle == "Pi-" && hadron == "Proton"){
    //cout << "Add Pi-/Proton case here" << endl;
  }
  else if (particle == "K+" && hadron == "Lambda"){
    //cout << "Add K+/Lambda case here" << endl;
    sigL_ParamFile = "../src/eic_evgen/CrossSection_Params/KPlusLambda_Param_sigL";
    sigT_ParamFile = "../src/eic_evgen/CrossSection_Params/KPlusLambda_Param_sigT"; // Shouldn't really have a relative path, should look at setting a DEMPGen variable and doing this in a better way later
  }
  else if (particle == "K+" && hadron == "Sigma0"){
    //cout << "Add K+/Sigma case here" << endl;
    sigL_ParamFile = "../src/eic_evgen/CrossSection_Params/KPlusSigma_Param_sigL";
    sigT_ParamFile = "../src/eic_evgen/CrossSection_Params/KPlusSigma_Param_sigT";
  }
  else if (particle == "Pi0"){
    //cout << "Add Pi0 case here" << endl;
  }
  else{
    cerr << " !!!!! " << endl << "Warning! Combination of specified ejectile and recoil particles not found!" << " !!!!! " << endl;
  }
 
  //....................................................................................................
  // Love's model parameters (Gojko, Stephen and Nishchey helped me to understand this part) 
  //....................................................................................................
  double ptmp;
  std::vector<std::vector<std::vector<std::vector<double>>>> p_vec;
  fstream file_vgl; // The parameterization file we will open and loop over

  for (int i = 0; i < 2; i++){
    if(i == 0){
      file_vgl.open(sigL_ParamFile, ios::in); 
    }
    if(i == 1){
      file_vgl.open(sigT_ParamFile, ios::in); 
    }
    p_vec.push_back(std::vector<std::vector<std::vector<double>>>());
    for(int j=0; j <9; j++){// Loop over all values of W - 2 to 10
      p_vec[i].push_back(std::vector<std::vector<double>>());
      for(int k=0; k<35; k++){ // Loop over all values of Q2 - 1 to 35 for each w
	p_vec[i][j].push_back(std::vector<double>());
      
	for(int l=0; l<13; l++){ //Loop over all columns at once
	  file_vgl>>ptmp;
	  p_vec[i][j][k].push_back(ptmp);
	}
      }
    }
    file_vgl.close();// Need to close the file at end of each loop over i 
  }       
  return p_vec;
}
