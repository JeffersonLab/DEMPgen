#include "reaction_routine.h"
#include "eic.h"
#include "particleType.h"
#include <cstdlib>

using namespace std;

Pi0_Production::Pi0_Production() { 

  cout << "Program Start" << endl;

}

/*--------------------------------------------------*/
/// Pi0_Production 

Pi0_Production::Pi0_Production(TString particle_str) { 

  rParticle = particle_str;
  
  Init();

  if (gOutputType == "Pythia6"){
    Pi0_Decay_Pythia6_Out_Init();
  }
  else if (gOutputType == "HEPMC3"){
    Pi0_HEPMC3_Out_Init();
  }
}

/*--------------------------------------------------*/

Pi0_Production::~Pi0_Production() {

  DEMPOut.close();
  DEMPDetails.close();

}

/*--------------------------------------------------*/

void Pi0_Production::process_reaction() {

  if_pi0_decay = gPi0_decay;
  
  polar_out.open("test.txt", ofstream::out);
  
  for( long long int i = 0; i < rNEvents; i++ ) {
  
    rNEvent_itt = i;
    fNGenerated ++;
    
    Progress_Report();  // This is happens at each 10% of the total event is processed
    Processing_Event();
  	
  }
  
  polar_out.close();
  
  Detail_Output();

  system("./../../practice/outputs/bin/rootIOTree_example_write.exe LundFiles/eic_input_DEMPGen_EIC_test.dat LundFiles/Converted_Output_ROOT_File.root");

}

/*--------------------------------------------------*/

void Pi0_Production::Processing_Event() {

  // ----------------------------------------------------
  // Considering Fermi momentum for the proton
  // ----------------------------------------------------

  if( kCalcFermi ) {
    Consider_Proton_Fermi_Momentum();
  }

  // ----------------------------------------------------
  // Boost vector from collider (lab) frame to protons rest frame (Fix target)
  // ----------------------------------------------------
 
  beta_col_rf = r_lproton.BoostVector();        
  fGamma_Col_RF = 1.0/sqrt( 1 - pow( beta_col_rf.Mag() , 2 ) );

  // ---------------------------------------------------------------------
  // Specify the energy and solid angle of scatterd electron in Collider (lab) frame
  // ---------------------------------------------------------------------

  fScatElec_Theta_Col  = acos( fRandom->Uniform( cos( fScatElec_Theta_I ) , cos( fScatElec_Theta_F ) ) );
  fScatElec_Phi_Col    = fRandom->Uniform( 0 , 2.0 * fPi );
  fScatElec_Energy_Col = fRandom->Uniform( fScatElec_E_Lo * fElectron_Energy_Col , fScatElec_E_Hi * fElectron_Energy_Col );

  // fScatElec_Theta_Col  = fRandom->Uniform( 1.7 , 1.9 );
  // fScatElec_Phi_Col = fRandom->Uniform( 0-0.01 , 0+0.01 );
  // fScatElec_Energy_Col = fRandom->Uniform( 1.33 * fElectron_Energy_Col , 2.50 * fElectron_Energy_Col );

  // ----------------------------------------------------
  // Produced Particle X in Collider frame
  // ----------------------------------------------------
  
  /// The generic produced particle in the exclusive reaction is labelled as X 
  fX_Theta_Col      = acos( fRandom->Uniform( cos(fX_Theta_I), cos(fX_Theta_F ) ) );
  fX_Phi_Col        = fRandom->Uniform( 0 , 2.0 * fPi );

  // fX_Theta_Col      = fRandom->Uniform( 0.0 , 0.9 );
  // fX_Phi_Col        = fRandom->Uniform( -0.0000001 + fPi , 0.0000001 + fPi );

	// ---------------------------------------------------------------------
    // Specify the energy and solid angle of scatterd electron in Collider (lab) frame
    // ---------------------------------------------------------------------

  fScatElec_Mom_Col  = sqrt( pow( fScatElec_Energy_Col,2) - pow( fElectron_Mass , 2) );
  fScatElec_MomZ_Col = ( fScatElec_Mom_Col * cos(fScatElec_Theta_Col) );  
  fScatElec_MomX_Col = ( fScatElec_Mom_Col * sin(fScatElec_Theta_Col) * cos(fScatElec_Phi_Col) );
  fScatElec_MomY_Col = ( fScatElec_Mom_Col * sin(fScatElec_Theta_Col) * sin(fScatElec_Phi_Col) );

  r_lscatelec.SetPxPyPzE( fScatElec_MomX_Col, fScatElec_MomY_Col, fScatElec_MomZ_Col, fScatElec_Energy_Col);
 
  r_lscatelecg = r_lscatelec * fm;

  // ----------------------------------------------------
  // Photon in collider (lab) frame and Qsq
  // ----------------------------------------------------

  r_lphoton  = r_lelectron - r_lscatelec;
  r_lphotong = r_lelectrong - r_lscatelecg;

  fQsq_GeV = -1.* r_lphotong.Mag2();

  if ( fQsq_GeV < 5.0 ) {
     qsq_ev++;
     return;
  }

    // ----------------------------------------------------
    // W square, Invariant Mass (P_g + P_p)^2
    // ----------------------------------------------------
 
    TLorentzVector lwg;
    lwg = r_lprotong + r_lphotong;
    fW_GeV    = lwg.Mag();
    fWSq_GeV  = lwg.Mag2();

    if ( fWSq_GeV < 0 ) {
		w_neg_ev++;
		return;
    }    

  // 13/12/22 - SJDK - This is the start of the block that will need to be replaced by the ROOT function Rory used to determine the pion momentum
  // ---------------------------------------------------------
  // Pion momentum in collider frame, analytic solution starts
  // ---------------------------------------------------------

  if(!Solve()){
    return;
  }

  // ---------------------------------------------------------
  // Pion momentum in collider frame, analytic solution starts
  // ---------------------------------------------------------

  double fupx = sin( fX_Theta_Col ) * cos( fX_Phi_Col );
  double fupy = sin( fX_Theta_Col ) * sin( fX_Phi_Col );
  double fupz = cos( fX_Theta_Col );
 
  double fuqx = sin( r_lphoton.Theta() ) * cos( r_lphoton.Phi() );
  double fuqy = sin( r_lphoton.Theta() ) * sin( r_lphoton.Phi() );
  double fuqz = cos( r_lphoton.Theta() );
 
  double fa = -(r_lphoton.Vect()).Mag() * ( fupx * fuqx +  fupy * fuqy +  fupz * fuqz );
  double fb = pow ( (r_lphoton.Vect()).Mag() , 2 );
  double fc = r_lphoton.E() + fProton_Mass;

  fa = ( fa - std::abs( (r_lproton.Vect()).Mag() ) * ( ( ( r_lproton.X() / (r_lproton.Vect()).Mag() ) * fupx ) + 
  					       ( ( r_lproton.Y() / (r_lproton.Vect()).Mag() ) * fupy ) + 
  					       ( ( r_lproton.Z() / (r_lproton.Vect()).Mag() ) * fupz ) ) );
   
  double factor = ( pow( (r_lproton.Vect()).Mag() , 2 ) + 2.0 * (r_lphoton.Vect()).Mag() * (r_lproton.Vect()).Mag() *  
  	      ( ( ( r_lproton.X() / (r_lproton.Vect()).Mag() ) * fuqx ) + 
  		( ( r_lproton.Y() / (r_lproton.Vect()).Mag() ) * fuqy ) + 
  		( ( r_lproton.Z() / (r_lproton.Vect()).Mag() ) * fuqz ) ) );
   
  fb =  fb + factor;  
  fc = r_lphoton.E() + r_lproton.E();

  double e_X_rf      = lX_rf.E();
  double e_photon_rf = lphoton_rf.E();
  
  
  double ft = fc * fc - fb + fX_Mass * fX_Mass - fProton_Mass * fProton_Mass;
  
  double fu = lproton_rf.Dot(lproton_rf) + lX_rf.Dot(lX_rf) - lproton_rf.Dot(lX_rf);

  ft_min = -fQsq_GeV + pow(fX_Mass, 2) -2.*(e_X_rf*e_photon_rf -sqrt( (pow(e_X_rf, 2) - pow(fX_Mass, 2)) * (pow(e_photon_rf,2)+fQsq) ));
  
  fu_min = -fQsq_GeV + pow(f_Scat_hadron_Mass, 2) -2.*(e_X_rf*e_photon_rf -sqrt( (pow(e_X_rf, 2) - pow(f_Scat_hadron_Mass , 2))*(pow(e_photon_rf,2)+fQsq) ));

  double fQA = 4.0 * ( fa * fa - fc * fc );
  double fQB = 4.0 * fc * ft;

  double fQC = -4.0 * fa * fa * fX_Mass * fX_Mass - ft * ft;    
 
  fradical = fQB * fQB - 4.0 * fQA * fQC;
 
  fepi1 = ( -fQB - sqrt( fradical ) ) / ( 2.0 * fQA );
  fepi2 = ( -fQB + sqrt( fradical ) ) / ( 2.0 * fQA );
 
  ///---------------------------------------------------------
  /// Particle X momentum in collider frame, analytic solution
  /// And obtain recoiled proton in collider (lab) frame
  ///---------------------------------------------------------
       
  //r_lX.SetPxPyPzE( (sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * sin(fX_Theta_Col) * cos(fX_Phi_Col),
  //		  ( sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * sin(fX_Theta_Col) * sin(fX_Phi_Col),
  //		  ( sqrt( pow( fepi1 , 2) - pow(fX_Mass , 2) ) ) * cos(fX_Theta_Col),
  //		  fepi1 );
  //
  // r_l_scat_hadron.SetPxPyPzE( ( r_lproton + r_lelectron - r_lscatelec - r_lX).X(),
  // 		     ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Y(),
  // 		     ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Z(),
  // 		     sqrt( pow( ( ( ( r_lproton + r_lelectron - r_lscatelec - r_lX ).Vect() ).Mag()),2) +
  // 			   pow( f_Scat_Hadron_Mass ,2 ) ) );

  ///--------------------------------------------------
  /// Output with the Solve Function
  /// Setting the solution values to X and recoiled nucleon

  r_lX.SetPxPyPzE(r_lX_solved->Px(), r_lX_solved->Py(), r_lX_solved->Pz(), r_lX_solved->E());
  r_l_scat_hadron.SetPxPyPzE(r_l_scat_hadron_solved->Px(), r_l_scat_hadron_solved->Py(), r_l_scat_hadron_solved->Pz(),  r_l_scat_hadron_solved->E());

  ///--------------------------------------------------

  r_lX_g = r_lX * fm;
  r_l_scat_hadron_g = r_l_scat_hadron * fm;

  // ----------------------------------------------------------------------------------------------
  // Calculate w = (proton + photon)^2
  // ----------------------------------------------------------------------------------------------
   
  // cout << fW_GeV << endl;
  if ( fW_GeV < 3.0 || fW_GeV > 10.6 ) {
    w_ev++;
    return;
  }

  r_lw = r_lproton + r_lphoton;
  fW = r_lw.Mag();

  // ----------------------------------------------------------------------------------------------
  // Calculate w prime w' = (proton + photon - pion)^2                                             
  // ----------------------------------------------------------------------------------------------

  lwp = r_lprotong + r_lphotong - r_lX_g;
  fW_Prime_GeV = lwp.Mag();    

  fsini = r_lelectron + r_lproton;
  fsfin = r_lscatelec + r_lX + r_l_scat_hadron;
  
  fsinig = fsini * fm;
  fsfing = fsfin * fm; 
 
  fMandSConserve = std::abs( fsinig.Mag() - fsfing.Mag() );

  //*--------------------------------------------------*/ 
  // SJDK 15/06/21 - Added integer counters for conservation law check and for NaN check
  if (r_lX.E() != r_lX.E()){ // SJDK 15/06/21 - If the energy of the produced meson is not a number, return and add to counter
    fNaN++;
    return;
  }

  kSConserve = false;
  if( std::abs( fsinig.Mag() - fsfing.Mag() ) < fDiff ) {
    kSConserve = true;
  }
        
  if ( pd->CheckLaws( r_lelectron, r_lproton, r_lscatelec, r_lX, r_l_scat_hadron) != 1 ){
    fConserve++;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  //                                          Start                                         //
  // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
  ////////////////////////////////////////////////////////////////////////////////////////////
 
  lproton_rf = r_lproton;
  lproton_rf.Boost(-beta_col_rf);
  lproton_rfg = lproton_rf * fm;
  
  lelectron_rf = r_lelectron;
  lelectron_rf.Boost(-beta_col_rf);
  lelectron_rfg = lelectron_rf * fm;
  
  lscatelec_rf = r_lscatelec;
  lscatelec_rf.Boost(-beta_col_rf);
  lscatelec_rfg = lscatelec_rf * fm;
  
  lphoton_rf = r_lphoton;
  lphoton_rf.Boost(-beta_col_rf);
  lphoton_rfg = lphoton_rf * fm;
  
  lX_rf = r_lX;
  lX_rf.Boost(-beta_col_rf);
  lX_rfg = lX_rf * fm;
      
  l_scat_hadron_rf = r_l_scat_hadron;
  l_scat_hadron_rf.Boost(-beta_col_rf);
  l_scat_hadron_rf_g = l_scat_hadron_rf * fm;
  
  ///*--------------------------------------------------*/
  /// Doing pi0 decay 
  
  if (if_pi0_decay) {
      Pi0_decay(r_lX);
  } 


  ////////////////////////////////////////////////////////////////////////////////////////////
  //                                          End                                           //
  // Transformation of e', pi- and recoil proton to target's rest frmae without energy loss //
  ////////////////////////////////////////////////////////////////////////////////////////////

  /*--------------------------------------------------*/
  /*--------------------------------------------------*/
  /*--------------------------------------------------*/
  /*--------------------------------------------------*/
  /// Polar output

  if (fQsq_GeV > 5 && fQsq_GeV < 12 && fWSq_GeV > 8 && fWSq_GeV < 12) {

// 		cout << fQsq_GeV << "  " <<  fWSq_GeV << "  ";             
// 
// 		cout << r_lscatelecg.Vect().Theta()       << "  " <<  r_lscatelecg.Vect().Mag()    << "  " 
// 			 << r_l_scat_hadron_g.Vect().Theta() << "  " << r_l_scat_hadron_g.Vect().Mag() << "  " 
//              << r_lX_g.Vect().Theta()             << "  " << r_lX_g.Vect().Mag()                  
// 			 << endl; 
		
    polar_out << fQsq_GeV << "  " <<  fWSq_GeV << "  ";             
  
    polar_out << r_lscatelecg.Vect().Theta()  << "  " <<  r_lscatelecg.Vect().Mag()      << "  " 
  	      << r_l_scat_hadron_g.Vect().Theta() << "  " << r_l_scat_hadron_g.Vect().Mag() << "  " 
          << r_lX_g.Vect().Theta()             << "  " << r_lX_g.Vect().Mag()             << "  "                
          << l_photon_1.Vect().Theta()         << "  " << l_photon_1.Vect().Mag()         << "  "                
          << l_photon_2.Vect().Theta()         << "  " << l_photon_2.Vect().Mag()                  
  	      << endl; 

  }

//	exit(0);
//	return;


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/

  e_X_rf      = lX_rf.E();
  e_photon_rf = lphoton_rf.E();
  double e_p_rf      = lproton_rf.E();

  double e_photCM = (fWSq_GeV - fQsq_GeV - pow(fProton_Mass/1000, 2))/fW_GeV/2.;
  double e_pCM = (fWSq_GeV  + pow(fProton_Mass/1000, 2) - pow(f_Scat_hadron_Mass/1000, 2))/fW_GeV/2.;

  ft_min = -fQsq_GeV + pow(fX_Mass/1000, 2) -2.*(e_X_rf/1000*e_photon_rf/1000 -sqrt( (pow(e_X_rf/1000, 2) - pow(fX_Mass/1000, 2)) * (pow(e_photon_rf/1000,2)+fQsq_GeV) ));

  fu_min = -fQsq_GeV + pow(fProton_Mass/1000, 2) -2.*(e_pCM*e_photCM -sqrt( (pow(e_pCM/1000, 2) - pow(fProton_Mass/1000 , 2))*(pow(e_photCM/1000,2)+fQsq_GeV) ));

  // ----------------------------------------------------------------------------------------------
  // Cut on centre of mass angle of meson X
  // ----------------------------------------------------------------------------------------------

  // Theta_cm = r_lX.Vect().Angle(r_lphoton.Vect());

  // if (Theta_cm < 90*TMath::Pi()/180) {
  //   return;
  // }

  // -----------------------------------------------------------------------------------------
  // Calculate -t and -u
  // -----------------------------------------------------------------------------------------

  fBeta_CM_RF        = (lphoton_rf.Vect()).Mag() / ( lphoton_rf.E() + fProton_Mass );
  fGamma_CM_RF       = ( lphoton_rf.E() + fProton_Mass ) / fW;
  fX_Energy_CM       = ( pow( fW , 2) + pow(fX_Mass , 2) - pow(f_Scat_hadron_Mass , 2) ) / ( 2.0 * fW);    
  fX_Mom_CM          = sqrt( pow(fX_Energy_CM , 2) - pow(fX_Mass , 2));    
  fX_Energy_CM_GeV   = fX_Energy_CM / 1000.0;
  fX_Mom_CM_GeV      = fX_Mom_CM / 1000.0;

  // this equation is valid for parallel kinematics only!
  fT_Para = ( pow(((r_lphoton.Vect()).Mag() - (r_lX.Vect()).Mag()),2) - pow((r_lphoton.E() - r_lX.E()),2));
  fT_Para_GeV = fT_Para/1000000.0;

  lt = r_lphoton - r_lX;
  ltg = lt * fm;

  fT = -1.*lt.Mag2();
  fT_GeV = -1.*ltg.Mag2();

  tc = -1.*lt.Mag2();
  tc_GeV = -1.*ltg.Mag2();

  lu = r_lproton - r_lX;
  lug = lu * fm;

  uc = -1.*lu.Mag2();
  uc_GeV = -1.*lug.Mag2();

  // ----------------------------------------------------------------------------------------------
  // Cut on u-channel
  // ----------------------------------------------------------------------------------------------

  if ( abs(uc_GeV) > abs(tc_GeV) ) {
    return;
  }

  // if ( gKinematics_type == 1 && tc_GeV > 0.5 ) {
  //   t_ev++;
  //   return;
  // }

  // if ( gKinematics_type == 2 && tc_GeV > 1.3 ) {
  //   t_ev++;
  //   return;
  // }

  // if ( gKinematics_type == 1 && uc_GeV > 0.5 ) {
  //   t_ev++;
  //   return;
  // }

  // if ( gKinematics_type == 2 && uc_GeV > 1.3 ) {
  //   t_ev++;
  //   return;
  // }

  fx = fQsq_GeV / ( 2.0 * r_lprotong.Dot( r_lphotong ) );
  fy = r_lprotong.Dot( r_lphotong ) / r_lprotong.Dot( r_lelectrong );
  fz = r_lX.E()/r_lphoton.E();    

  // -------------------------------------------------------------------------------------------------------
  // Calculation of Phi  ( azimuthal angle of pion momentum w.r.t lepton plane in target's rest frame)
  // Calculation of PhiS ( azimuthal angle of target polarization w.r.t lepton plane in target's rest frame)
  // -------------------------------------------------------------------------------------------------------  

  v3Photon.SetX( lphoton_rfg.X() );     
  v3Photon.SetY( lphoton_rfg.Y() );     
  v3Photon.SetZ( lphoton_rfg.Z() );    

  v3Electron.SetX( lelectron_rfg.X() ); 
  v3Electron.SetY( lelectron_rfg.Y() ); 
  v3Electron.SetZ( lelectron_rfg.Z() );

  v3X.SetX( lX_rfg.X() ) ;        
  v3X.SetY( lX_rfg.Y() ) ;        
  v3X.SetZ( lX_rfg.Z() );

  v3S.SetX( -1 );                       
  v3S.SetY( 0 );                        
  v3S.SetZ( 0 );        

  v3PhotonUnit = v3Photon.Unit();    
  v3QxL        = v3Photon.Cross(v3Electron);
  v3QxP        = v3Photon.Cross(v3X);
  v3QxS        = v3Photon.Cross(v3S);
  v3LxP        = v3Electron.Cross(v3X);
  v3LxS        = v3Electron.Cross(v3S);
  v3PxL        = v3X.Cross(v3Electron);
  v3QUnitxL    = v3PhotonUnit.Cross(v3Electron);
  v3QUnitxP    = v3PhotonUnit.Cross(v3X);
  v3QUnitxS    = v3PhotonUnit.Cross(v3S);

  /*--------------------------------------------------*/
  // Get the Phi scattering angle with respect to the electron scattering plane
  fPhi   = Get_Phi_X_LeptonPlane_RF ();

  /*--------------------------------------------------*/
  // Get the Phi scattering angle between the target polarization plane and the electron scattering plane
  fPhiS  = Get_Phi_TargPol_LeptonPlane_RF();

  fTheta_X_Photon_RF       = fRAD2DEG * acos( ( v3Photon.Dot( v3X     ) ) / ( v3Photon.Mag()  * v3X.Mag()    ) );
  if ( fTheta_X_Photon_RF < 0 ) { fTheta_X_Photon_RF = 180.0 + fTheta_X_Photon_RF; }

  // -----------------------------------------------------------------------------------
  // If we have fermi momentum then epsilon should be in rest frame 
  // The theta angle of scattered angle used in expression of epsilon is the angle 
  // with respect to direction of incoming electron in the rest frame of target nucleon
  // epsilon=1./(1.+ 2.*(pgam_restg**2)/q2g * *(tand(thscat_rest/2.))**2)
  // -----------------------------------------------------------------------------------

  double fTheta_EEp = (lelectron_rf.Vect()).Angle(lscatelec_rf.Vect());

  fEpsilon = 1.0 / ( 1.0 + 2.0 * ( pow( (lphoton_rfg.Vect()).Mag(),2)/fQsq_GeV ) * pow( tan( fTheta_EEp / 2 ) , 2 ) );

  /// meson produced angle with respect to the Q-vector

  theta_X_rf = (lX_rf.Vect()).Angle(lphoton_rf.Vect());

  // ----------------------------------------------------
  // Virtual Photon flux factor in units of 1/(GeV*Sr)
  // ----------------------------------------------------

  fFlux_Factor_Col = (fAlpha/(2.0*pow(fPi,2))) * (r_lscatelecg.E() / r_lelectrong.E()) * 
    ( pow(fW_GeV,2) - pow(fProton_Mass_GeV,2) ) / (2.0*fProton_Mass_GeV*fQsq_GeV*(1.0 - fEpsilon));

  fFlux_Factor_RF = ( fAlpha / ( 2.0 * pow( fPi , 2 ) ) ) * ( lscatelec_rfg.E() / lelectron_rfg.E() ) *
    ( pow( fW_GeV , 2 ) - pow( fProton_Mass_GeV , 2 ) ) /
    ( 2.0 * fProton_Mass_GeV * fQsq_GeV * ( 1.0 - fEpsilon ) );

  // ----------------------------------------------------
  //  Jacobian  dt/dcos(theta*)dphi in units of GeV2/sr
  // ----------------------------------------------------

  fJacobian_CM = ( (lphoton_rfg.Vect()).Mag() - fBeta_CM_RF * lphoton_rfg.E() ) / ( fGamma_CM_RF * ( 1.0 - pow(fBeta_CM_RF,2) ) );

  fA = fJacobian_CM * fX_Mom_CM_GeV / fPi;

  double fttt = r_lprotong * r_lprotong - 2 * r_lprotong * r_l_scat_hadron_g  + r_l_scat_hadron_g * r_l_scat_hadron_g;
  double fuuu = r_lprotong * r_lprotong - 2 * r_lprotong * r_lX_g + r_lX_g * r_lX_g;

  // ----------------------------------------------------
  // Jacobian dOmega* / dOmega dimensionless
  // ----------------------------------------------------

  fJacobian_CM_RF  = ( pow((lX_rf.Vect()).Mag(),2)*fW) / 
    ( fX_Mom_CM * std::abs( ( fProton_Mass + lphoton_rf.E()) * (lX_rf.Vect()).Mag() - 
			    ( lX_rf.E() * (lphoton_rf.Vect()).Mag() * cos( lX_rf.Theta() ) ) ) );

  fJacobian_CM_Col = ( ( pow((r_lX.Vect()).Mag(),2) * fW ) /
		       ( fX_Mom_CM * std::abs( ( fProton_Mass + r_lphoton.E() ) * (r_lX.Vect()).Mag() -
					       ( r_lX.E() * (r_lphoton.Vect()).Mag() * cos( r_lX.Theta() ) ) ) ) );

  // ------------------------------------------------------------------------------------------
  // CKY sigma L and T starts
  // ------------------------------------------------------------------------------------------

  r_fSig = Get_CrossSection();

  // ------------------------------------------------------------------------------------------
  // CKY sigma L and T ends
  // ------------------------------------------------------------------------------------------

  fSigma_Col = r_fSig * fFlux_Factor_Col * fA * fJacobian_CM_Col;

  if ( ( fSigma_Col <= 0 ) || std::isnan( fSigma_Col ) ) { 
    fNSigmaNeg ++;
    return;
  }

  // ------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------
  //             Lab cross section     Phase Space   Conversion     Luminosity                Total events tried
  // Hz        = ub / ( sr^2 * GeV ) * GeV * sr^2 * ( cm^2 / ub ) * ( # / ( cm^2 * sec ) ) / ( # )

  fEventWeight = fSigma_Col * fPSF * fuBcm2 * fLumi / fNEvents;   // in Hz

  fNRecorded ++;
  fLundRecorded++;
  fRatio = fNRecorded / fNGenerated;

  if (gOutputType == "Pythia6"){
    Pi0_Decay_Pythia6_Output();
  }
  else if (gOutputType == "LUND"){
    if (if_pi0_decay){
      Pi0_Decay_Lund_Output();
    }
    else {
      Pi0_Lund_Output();
    }
  }
  else if (gOutputType == "HEPMC3"){	
    Pi0_HEPMC3_Output();
  }
	
}

/*--------------------------------------------------*/

void Pi0_Production::Detail_Output() {

  DEMPDetails << "Total events tried                           " << setw(50) << fNGenerated   << endl;
  DEMPDetails << "Total events recorded                        " << setw(50) << fNRecorded    << endl;

  DEMPDetails << "Seed used for the Random Number Generator    " << setw(50) << fSeed         << endl;
  DEMPDetails << "Number of lund events                        " << setw(50) << fLundRecorded << endl;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Pi0_Production::Pi0_Lund_Output() {

  DEMPOut << "5"
	 << " \t " << fPhi           // var 1
	 << " \t " << fPhiS          // var 2
	 << " \t " << fx             // var 3
	 << " \t " << "1"
	 << " \t " << fQsq_GeV       // var 4
	 << " \t " << fT_GeV         // var 5
	 << " \t " << fW_GeV 	     // var 6
	 << " \t " << fEpsilon       // var 7
	 << " \t " << fEventWeight   // var 8	   
	 << endl;
 
  ///*--------------------------------------------------*/
  // Initial State

  DEMPOut << setw(10) << "1"
	 << setw(10) << "-1"
	 << setw(10) << "0"
	 << setw(10) << "11"
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_lelectrong.X()
	 << setw(16) << r_lelectrong.Y()
	 << setw(16) << r_lelectrong.Z()
	 << setw(16) << r_lelectrong.E()
	 << setw(16) << fElectron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  DEMPOut << setw(10) << "2"
	 << setw(10) << "1"
	 << setw(10) << "0"
	 << setw(10) << "2212"
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_lprotong.X()
	 << setw(16) << r_lprotong.Y()
	 << setw(16) << r_lprotong.Z()
	 << setw(16) << r_lprotong.E()
	 << setw(16) << fProton_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  ///*--------------------------------------------------*/
  // Final State

  // Produced Particle X
  DEMPOut << setw(10) << "3"
	 << setw(10) << "1"
	 << setw(10) << "1"
	 << setw(10) << PDGtype(produced_X)
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_lX_g.X()
	 << setw(16) << r_lX_g.Y()
	 << setw(16) << r_lX_g.Z()
	 << setw(16) << r_lX_g.E()
	 << setw(16) << fX_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  // Scattered electron
  DEMPOut << setw(10) << "4"
	 << setw(10) << "-1"
	 << setw(10) << "1"
	 << setw(10) << "11"
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_lscatelecg.X()
	 << setw(16) << r_lscatelecg.Y()
	 << setw(16) << r_lscatelecg.Z()
	 << setw(16) << r_lscatelecg.E()
	 << setw(16) << fElectron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  // Recoiled neutron
  DEMPOut << setw(10) << "5"
	 << setw(10) << "1"
	 << setw(10) << "1"
	 << setw(10) << PDGtype(recoil_hadron)
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_l_scat_hadron_g.X()
	 << setw(16) << r_l_scat_hadron_g.Y()
	 << setw(16) << r_l_scat_hadron_g.Z()
	 << setw(16) << r_l_scat_hadron_g.E()
	 << setw(16) << f_Scat_hadron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Pi0_Production::Pi0_Decay_Lund_Output() {

  DEMPOut << "4"
	 << " \t " << fPhi           // var 1
	 << " \t " << fPhiS          // var 2
	 << " \t " << fx             // var 3
	 << " \t " << "1"
	 << " \t " << fQsq_GeV       // var 4
	 << " \t " << fT_GeV         // var 5
	 << " \t " << fW_GeV 	     // var 6
	 << " \t " << fEpsilon       // var 7
	 << " \t " << fEventWeight   // var 8	   
	 << endl;

  // Scattered electron
  DEMPOut << setw(10) << "2"
	 << setw(10) << "-1"
	 << setw(10) << "1"
	 << setw(10) << "11"
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_lscatelecg.X()
	 << setw(16) << r_lscatelecg.Y()
	 << setw(16) << r_lscatelecg.Z()
	 << setw(16) << r_lscatelecg.E()
	 << setw(16) << fElectron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  // Recoiled neutron
  DEMPOut << setw(10) << "3"
	 << setw(10) << "1"
	 << setw(10) << "1"
	 << setw(10) << PDGtype(recoil_hadron)
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << r_l_scat_hadron_g.X()
	 << setw(16) << r_l_scat_hadron_g.Y()
	 << setw(16) << r_l_scat_hadron_g.Z()
	 << setw(16) << r_l_scat_hadron_g.E()
	 << setw(16) << f_Scat_hadron_Mass_GeV
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  // Photon 1
  DEMPOut << setw(10) << "3"
	 << setw(10) << "1"
	 << setw(10) << "1"
	 << setw(10) << "22"
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << l_photon_1.X()
	 << setw(16) << l_photon_1.Y()
	 << setw(16) << l_photon_1.Z()
	 << setw(16) << l_photon_1.E()
	 << setw(16) << "0"
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

  // Photon 2
  DEMPOut << setw(10) << "4"
	 << setw(10) << "1"
	 << setw(10) << "1"
	 << setw(10) << "22"
	 << setw(10) << "0"
	 << setw(10) << "0"
	 << setw(16) << l_photon_2.X()
	 << setw(16) << l_photon_2.Y()
	 << setw(16) << l_photon_2.Z()
	 << setw(16) << l_photon_2.E()
	 << setw(16) << "0"
	 << setw(16) << fVertex_X
	 << setw(16) << fVertex_Y
	 << setw(16) << fVertex_Z
	 << endl;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Pi0_Production::Pi0_Decay_Pythia6_Out_Init() {

  print_itt = 0;

  //	DEMPOut << "PYTHIA EVENT FILE" << endl;
  DEMPOut << "SIMPLE Event FILE" << endl;
  DEMPOut << "============================================" << endl;
  DEMPOut << "I, ievent, nParticles" << endl;
  DEMPOut << "============================================" << endl;
  DEMPOut << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)" << endl;
  DEMPOut << "============================================" << endl;

}

/*--------------------------------------------------*/

void Pi0_Production::Pi0_Decay_Pythia6_Output() {

  DEMPOut << "0" << " \t\t\t "  << print_itt << " \t\t\t " << "1" << endl;           // var 1

  print_itt++;

  DEMPOut << "============================================" << endl;

  ///*--------------------------------------------------*/
  // Initial State
 
  DEMPOut  << "1"
	  << setw(6) << "21"
	  << setw(6) << "11"
	  << setw(6) << "0"
	  << setw(6) << "3"
	  << setw(6) << "4"

	  << setw(14) << r_lelectrong.X()
	  << setw(14) << r_lelectrong.Y()
	  << setw(14) << r_lelectrong.Z()
	  << setw(14) << r_lelectrong.E()
	  << setw(14) << fElectron_Mass_GeV
	  << setw(6) << fVertex_X
	  << setw(6) << fVertex_Y
	  << setw(6) << fVertex_Z
	  << endl;

  DEMPOut << "2"
	 << setw(6) << "21"
	 << setw(6) << "2212"
	 << setw(6) << "0"
	 << setw(6) << "5"
	 << setw(6) << "6"

	 << setw(14) << r_lprotong.X()
	 << setw(14) << r_lprotong.Y()
	 << setw(14) << r_lprotong.Z()
	 << setw(14) << r_lprotong.E()
	 << setw(14) << fProton_Mass_GeV
	 << setw(6) << fVertex_X
	 << setw(6) << fVertex_Y
	 << setw(6) << fVertex_Z
	 << endl;

  DEMPOut << "3"
	 << setw(6) << "21"
	 << setw(6) << "22"
	 << setw(6) << "1"
	 << setw(6) << "0"
	 << setw(6) << "0"

	 << setw(14) << r_lphotong.X()
	 << setw(14) << r_lphotong.Y()
	 << setw(14) << r_lphotong.Z()
	 << setw(14) << r_lphotong.E()
	 << setw(14) << r_lphotong.M()
	 << setw(6) << fVertex_X
	 << setw(6) << fVertex_Y
	 << setw(6) << fVertex_Z
	 << endl;

  ///*--------------------------------------------------*/
  // Final State

  // Scattered electron
  DEMPOut << "4"
	 << setw(6) << "1"
	 << setw(6) << "11"
	 << setw(6) << "1"
	 << setw(6) << "0"
	 << setw(6) << "0"

	 << setw(14) << r_lscatelecg.X()
	 << setw(14) << r_lscatelecg.Y()
	 << setw(14) << r_lscatelecg.Z()
	 << setw(14) << r_lscatelecg.E()
	 << setw(14) << fElectron_Mass_GeV
	 << setw(6) << fVertex_X
	 << setw(6) << fVertex_Y
	 << setw(6) << fVertex_Z
	 << endl;

  // Recoiled nucleon
  DEMPOut << "5"
	 << setw(6) << "1"
	 << setw(6) << PDGtype(recoil_hadron)
	 << setw(6) << "2"
	 << setw(6) << "0"
	 << setw(6) << "0"

	 << setw(14) << r_l_scat_hadron_g.X()
	 << setw(14) << r_l_scat_hadron_g.Y()
	 << setw(14) << r_l_scat_hadron_g.Z()
	 << setw(14) << r_l_scat_hadron_g.E()
	 << setw(14) << f_Scat_hadron_Mass_GeV
	 << setw(6) << fVertex_X
	 << setw(6) << fVertex_Y
	 << setw(6) << fVertex_Z
	 << endl;

  // Produced Particle X
  DEMPOut << "6"
	 << setw(6) << "1"
	 << setw(6) << PDGtype(produced_X)
	 << setw(6) << "2"
	 << setw(6) << "0"
	 << setw(6) << "0"

	 << setw(14) << r_lX_g.X()
	 << setw(14) << r_lX_g.Y()
	 << setw(14) << r_lX_g.Z()
	 << setw(14) << r_lX_g.E()
	 << setw(14) << fX_Mass_GeV
	 << setw(6) << fVertex_X
	 << setw(6) << fVertex_Y
	 << setw(6) << fVertex_Z
	 << endl;

  DEMPOut << "=============== Event finished ===============" << endl;

}

void Pi0_Production::Pi0_HEPMC3_Out_Init() {

  print_itt = 0;
  DEMPOut << "HepMC::Version 3.02.05" << endl;
  DEMPOut << "HepMC::Asciiv3-START_EVENT_LISTING" << endl;

}

/*--------------------------------------------------*/

void Pi0_Production::Pi0_HEPMC3_Output() {
  
  // HEPMC3 output for Athena simulations

  // First line - E - Event# - #Vertices - #Particles
  DEMPOut << "E" << " "  << print_itt <<  " " << "1" << " " << 5 << endl;
  print_itt++;
  // Second line, Units - U - ENERGY UNIT - DISTANCE UNIT
  DEMPOut << "U" << " " << "GEV" << " " << "MM" << endl;
  // Third line, optional attributes, the weight
  DEMPOut << "A" << " " << "0" << " " << "weight" << " " <<  fEventWeight << endl;
  // Beam particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (4, incoming beam particle)
  DEMPOut << "P" << " " << "1" << " " << "0" << " " << "11" << " " << r_lelectrong.X() << " " << r_lelectrong.Y() << " " << r_lelectrong.Z() << " " << r_lelectrong.E() << " " << fElectron_Mass_GeV << " " << "4" << endl;
  DEMPOut << "P" << " " << "2" << " " << "0" << " " << "2212" << " " << r_lprotong.X() << " " << r_lprotong.Y() << " " << r_lprotong.Z() << " " << r_lprotong.E() << " " << fProton_Mass_GeV << " " << "4" << endl;
  // Vertex line - V - 1 - 0 - [1,2]
  DEMPOut << "V" << " " << "-1" << " " << "0" << " " << "[1,2]" << endl;
  // Output particles, particle line - P - Particle ID - Parent Vertex ID - PDG id - px - py - pz - energy - particle mass - status (1, undecayed physical particle)
  // Scattered electron
  DEMPOut << "P" << " " << "3" << " " << "-1" << " " << "11" << " " << r_lscatelecg.X() << " "  << r_lscatelecg.Y() << " "  << r_lscatelecg.Z() << " " << r_lscatelecg.E() << " " << fElectron_Mass_GeV << " " << "1" << endl;
  // Produced meson
  DEMPOut << "P" << " " << "4" << " " << "-1" << " " << PDGtype(produced_X) << " " << r_lX_g.X() << " "  << r_lX_g.Y() << " "  << r_lX_g.Z() << " " << r_lX_g.E() << " " << fX_Mass_GeV << " " << "1" << endl;
  // Recoil nucleon
  DEMPOut << "P" << " " << "5" << " " << "-1" << " " << PDGtype(recoil_hadron) << " " << r_l_scat_hadron_g.X() << " "  << r_l_scat_hadron_g.Y() << " "  << r_l_scat_hadron_g.Z() << " " << r_l_scat_hadron_g.E() << " " << f_Scat_hadron_Mass_GeV << " " << "1" << endl;

}

/*--------------------------------------------------*/
/// Cross section model: Based on 12 GeV Hall C u-channel pi0 proposal
/// Author: Garth Huber
/// Date: April 1st, 2020

Double_t Pi0_Production::Get_CrossSection(){

  double_t sig_total;

  // use fit parameters for xB=0.36 data+model (nb/GeV2)
  Double_t Q2tab[4] = {1.75,  3.00,  4.00,  5.50};
  Double_t Wtab[4]  = {2.00,  2.46,  2.83,  3.26};

  // sigT is parameterized in the form p1+p2/(-t)
  Double_t p1T[4]   = {1577.,  168.,  67.4,  24.7};
  Double_t p2T[4]   = {-125., -11.1,  -4.6, -1.66};

  // sigL is simply constant vs t
  Double_t p1L[4]   = {0.,  2.77,  1.16,  0.43};

  // sigTT  p1+p2/(-t)
  Double_t p1TT[4]  = {-674., -141., -57.6, -21.3};
  Double_t p2TT[4]  = { 102.,  25.4,  10.4,  3.82};

  // sigLT  p1+p2/(-t)
  Double_t p1LT[4]  = { 156., 12.57,  5.17,  1.97};
  Double_t p2LT[4]  = {-60.0, -2.08, -0.85, -0.32};

  double sigThi, sigTlo, sigLhi, sigLlo;
  double sigTThi, sigTTlo, sigLThi, sigLTlo;
  double tmin, tprime;
  double umin, uprime;
  Int_t Q2count =0;

  double thetacm = theta_X_rf;

  double eps = fEpsilon;
  double phicm = fPhi;

  double Q2hi=0.;
  double Q2lo=0.;
  double delQ2=1.;

  double Whi, Wlo;

  double Q2tmp = fQsq_GeV;
  double m_p = fProton_Mass;
  double W_gev = fW_GeV;

  double pi = fPi;

  int ndat=4;

  tprime = ft-ft_min;
  uprime = fu-fu_min;

  if( Q2tmp < Q2tab[0] ) {

    Q2count = 0;
    Q2hi = Q2tab[1];
    Whi  = Wtab[1];
    Q2lo = Q2tab[0];
    Wlo  = Wtab[0];
    delQ2 = (Q2hi - Q2lo);

  } else {

    for(int Q2c=0; Q2c <= ndat-2; Q2c++){

      if( (Q2tmp >= Q2tab[Q2c] && Q2tmp < Q2tab[Q2c+1]) || Q2tmp >= Q2tab[Q2c+1] ) {

	Q2count = Q2c ;
	Q2hi = Q2tab[Q2count+1];
	Whi  = Wtab[Q2count+1];
	Q2lo = Q2tab[Q2count];
	Wlo  = Wtab[Q2count];
	delQ2 = (Q2hi - Q2lo);

      }
    }
  }

  ///*--------------------------------------------------*/
  // t-channel
  if (  thetacm > pi/2.) {

    sigThi = p1T[Q2count+1] + p2T[Q2count+1]/(tprime+fabs(tmin));
    sigTlo = p1T[Q2count]   + p2T[Q2count]/(tprime+fabs(tmin));
    sigLhi = p1L[Q2count+1];
    sigLlo = p1L[Q2count];

    sigTThi = p1TT[Q2count+1]+p2TT[Q2count+1]/(tprime+fabs(tmin));
    sigTTlo = p1TT[Q2count]  +p2TT[Q2count]/(tprime+fabs(tmin));
    sigLThi = p1LT[Q2count+1]+p2LT[Q2count+1]/(tprime+fabs(tmin));
    sigLTlo = p1LT[Q2count]  +p2LT[Q2count]/(tprime+fabs(tmin));

  } else {
    //
    ///*--------------------------------------------------*/
    // u-channel
    //   christian weiss recommends the following change for u-channel:
    //   switch u-slope for t-slope, then divide by 10, since back angle peak
    //   is ~10% of forward angle peak (at least for omega electroproduction)

    sigThi = (p1T[Q2count+1]+p2T[Q2count+1]/(uprime+fabs(tmin)))/10.;
    sigTlo = (p1T[Q2count]  +p2T[Q2count]/(uprime+fabs(tmin)))/10.;
    sigLhi = p1L[Q2count+1]/10.;
    sigLlo = p1L[Q2count]/10.;
    sigTThi = (p1TT[Q2count+1]+p2TT[Q2count+1]/(uprime+fabs(tmin)))/10.;
    sigTTlo = (p1TT[Q2count]  +p2TT[Q2count]/(uprime+fabs(tmin)))/10.;
    sigLThi = (p1LT[Q2count+1]+p2LT[Q2count+1]/(uprime+fabs(tmin)))/10.;
    sigLTlo = (p1LT[Q2count]  +p2LT[Q2count]/(uprime+fabs(tmin)))/10.;

  }

  double Wfac_hi = pow(pow(Whi, 2) - pow(m_p, 2), 2) / pow( pow(W_gev,2)- pow(m_p, 2), 2);
  double Wfac_lo = pow(pow(Wlo, 2) - pow(m_p, 2), 2) / pow( pow(W_gev,2)- pow(m_p, 2), 2);

  double sigThiW = sigThi*Wfac_hi;

  if (sigThiW  <0) {
    //       cout << pizero: sigThiW<0 ',sigThiW,uprime, abs(tmin)
    sigThiW=0.;
  }

  double sigTloW = sigTlo*Wfac_lo;

  if (sigTloW < 0) {
    //       write(6,*)' pizero: sigTloW<0 ',sigTloW,uprime,abs(tmin)
    sigTloW=0.;
  }

  double sigLhiW = sigLhi*Wfac_hi;

  if (sigLhiW < 0) {
    //       write(6,*)' pizero: sigLhiW<0 ',sigLhiW,uprime,abs(tmin)
    sigLhiW=0.;
  }

  double sigLloW = sigLlo*Wfac_lo;

  if (sigLloW < 0) {
    //       write(6,*)' pizero: sigLloW<0 ',sigLloW,uprime,abs(tmin)
    sigLloW=0.;
  }

  double sigTThiW = sigTThi*Wfac_hi;
  double sigTTloW;
  double sigLThiW;
  double sigLTloW;

  double sig, sigL, sigT, sigLT, sigTT;

  if (abs(sigTThiW) < sigThiW) {
    sigTThiW=Sign(sigThiW, sigTThiW);
  }
  sigTTloW = sigTTlo*Wfac_lo;

  if (abs(sigTTloW) < sigTloW) {
    sigTTloW = Sign(sigTloW,sigTTloW);
  }
  sigLThiW = sigLThi*Wfac_hi;

  if (abs(sigLThiW) < sigThiW) {
    sigLThiW = Sign(sigThiW,sigLThiW);
  }
  sigLTloW = sigLTlo*Wfac_lo;

  if (abs(sigLTloW) < sigTloW) {
    sigLTloW = Sign(sigTloW,sigLTloW);
  }

  /*--------------------------------------------------*/

  if( Q2count <= (ndat-2) && Q2tmp >= Q2tab[Q2count] && Q2tmp < Q2tab[Q2count+1] ) {

    sigT  = ( sigTloW*(Q2hi-Q2tmp) + sigThiW*(Q2tmp-Q2lo))/delQ2;
    sigL  = ( sigLloW*(Q2hi-Q2tmp) + sigLhiW*(Q2tmp-Q2lo))/delQ2;
    sigTT = (sigTTloW*(Q2hi-Q2tmp) + sigTThiW*(Q2tmp-Q2lo))/delQ2;
    sigLT = (sigLTloW*(Q2hi-Q2tmp) + sigLThiW*(Q2tmp-Q2lo))/delQ2;

  } else if (Q2tmp >= Q2tab[ndat-1]) {

    sigT  =  sigThiW + (sigThiW-sigTloW)/delQ2;
    sigL  =  sigLhiW + (sigLhiW-sigLloW)/delQ2;
    sigTT = sigTThiW + (sigTThiW-sigTTloW)/delQ2;
    sigLT = sigLThiW + (sigLThiW-sigLTloW)/delQ2;

  } else if (Q2tmp <= Q2tab[1] ) {

    sigT  =  sigTloW - (sigThiW-sigTloW)/delQ2;
    sigL  =  sigLloW - (sigLhiW-sigLloW)/delQ2;
    sigTT = sigTTloW - (sigTThiW-sigTTloW)/delQ2;
    sigLT = sigLTloW - (sigLThiW-sigLTloW)/delQ2;

  } else {

    //         write(6,*)' Q2tmp error ',Q2tmp,Q2count
  }
  //      endif

  sig = sigT + eps*sigL +eps*cos(2.*phicm)*sigTT +sqrt(2.*eps*(1.+eps))*cos(phicm)*sigLT;

  double sig_pi0gmh = sig/2./pi*1.e-09;  //dsig/dtdphicm in microbarns/MeV^2/rad

  return sig_total;
}

///*--------------------------------------------------*/
// Pi0 decay subroutine

void Pi0_Production::Pi0_decay(TLorentzVector pi0_vec) {

  TVector3 beta_col_rf_pi0;

  beta_col_rf_pi0 = pi0_vec.BoostVector();

  Double_t fGamma_Col_pi0_RF = 1.0/sqrt( 1 - pow( beta_col_rf_pi0.Mag() , 2 ) );

  TLorentzVector l_photon_1_rf;
  TLorentzVector l_photon_2_rf;

  //	TLorentzVector pi0_vec;
  TLorentzVector pi0_vec_rf;
  TLorentzVector pi0_vec_rf_g;

  pi0_vec_rf = pi0_vec;

  pi0_vec_rf.Boost(-beta_col_rf_pi0);
  pi0_vec_rf_g = pi0_vec_rf * fm;

  double photon_1_E = pi0_vec_rf_g.E()/2;
  double photon_2_E = pi0_vec_rf_g.E()/2;

  double photon_1_theta = fRandom->Uniform( 0, fPi );
  double photon_1_phi   = fRandom->Uniform( 0, 2.0 * fPi);

  double photon_1_Mom_Col  = photon_1_E;
  double photon_1_MomZ_Col = ( photon_1_Mom_Col * cos(photon_1_theta) );
  double photon_1_MomX_Col = ( photon_1_Mom_Col * sin(photon_1_theta) * cos(photon_1_phi) );
  double photon_1_MomY_Col = ( photon_1_Mom_Col * sin(photon_1_theta) * sin(photon_1_phi) );

  l_photon_1_rf.SetPxPyPzE( photon_1_MomX_Col, photon_1_MomY_Col, photon_1_MomZ_Col, photon_1_E );

  l_photon_2_rf.SetPxPyPzE( -photon_1_MomX_Col, -photon_1_MomY_Col, -photon_1_MomZ_Col, photon_1_E );

  l_photon_1 = l_photon_1_rf;
  l_photon_2 = l_photon_2_rf;

  TVector3 beta_col_rf_pi0_inverse = -pi0_vec.BoostVector();

  l_photon_1.Boost(beta_col_rf_pi0);
  l_photon_2.Boost(beta_col_rf_pi0);

  double photons_1_sep = (tan(l_photon_1.Vect().Angle(pi0_vec.Vect())) * 32)*100;
  double photons_2_sep = (tan(l_photon_2.Vect().Angle(pi0_vec.Vect())) * 32)*100;

  double photons_seperation = photons_1_sep + photons_2_sep;

  ///*--------------------------------------------------*/
  /// Max angle formula: sin(theta) = m_pi/(2*E_gamma)
  //
  // Decay at 7.8m
  // 2.6033 ± 0.0005 × 10−8 life time

  /// ZDC is at 30 down stream
  // Pi0 drift is 22.2 meters

  // 2.6*3 = 7.8m before decay

  TLorentzVector photon_1_diff = l_photon_1 - pi0_vec*fm;

  double angle = atan(l_photon_1.Px()/l_photon_1.Py());

}


