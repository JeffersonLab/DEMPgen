///*--------------------------------------------------*/
/// Reaction_routine.cc
///
/// Creator: Wenliang (Bill) Li
/// Date: Mar 12 2020
/// Email: wenliang.billlee@gmail.com
///
/// Comment: Mar 12, 2020: Subroutine Exclusive_Omega_Prodoction is created 
///	         modeled off the Exclusive_Omega_Prodoction routine written by 
///          A. Zafar
///

#include "reaction_routine.h"
#include "eic.h"
#include "particleType.h"

using namespace std;


/*--------------------------------------------------*/
/// Reaction 

Reaction::Reaction(TString particle_str) { 

 	rEjectile = particle_str;
	cout << "Produced particle is: " << GetParticle() << endl; 
	cout << "Generated process: e + p -> e' + p' + " << GetParticle() << endl; 
    	tTime.Start(); 
 
 	cout << "/*--------------------------------------------------*/" << endl;
 	cout << "Starting setting up process" << endl;
 	cout << endl;
 
    	TDatime dsTime;
   	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;

}

// SJDK 09/02/22 - New reaction where the particle and hadron are specified
Reaction::Reaction(TString particle_str, TString hadron_str) { 

 	rEjectile = particle_str;
	rRecoil = hadron_str;
	cout << "Produced particle is: " << GetParticle() << endl; 
	cout << "Produced hadron is: " << GetHadron() << endl;
	cout << "Generated process: e + p -> e'+ " << GetHadron() << " + " << GetParticle() << endl; 
    	tTime.Start(); 
 
 	cout << "/*--------------------------------------------------*/" << endl;
 	cout << "Starting setting up process" << endl;
 	cout << endl;
 
    	TDatime dsTime;
   	cout << "Start Time:   " << dsTime.GetHour() << ":" << dsTime.GetMinute() << endl;

}

Reaction::Reaction(){};

Reaction::~Reaction() {
 
//	ppiOut.close();
// 	ppiDetails.close();
 
 	cout << endl;
 	cout << "Ending the process" << endl;
 	cout << "/*--------------------------------------------------*/" << endl;
 
    	tTime.Stop();
    	tTime.Print();
 
    	TDatime deTime;
    	cout << "End Time:   " << deTime.GetHour() << ":" << deTime.GetMinute() << endl;
 
}

/*--------------------------------------------------*/
///

void Reaction::process_reaction() {
  if (rEjectile == "Pi0") {
    //		Pi0_Production* r1 = new Pi0_Production("Eta");
    Pi0_Production* rr1 = new Pi0_Production(rEjectile);
    rr1->process_reaction();
    delete rr1;
  }
  // SJDK - 19/12/22 - New generic DEMP reaction class, the intention is that this should be able to handle any case
  else{
    DEMP_Reaction* rr1 = new DEMP_Reaction(rEjectile, rRecoil);
    rr1->process_reaction();
    delete rr1;
  }

}
