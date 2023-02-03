#include "KPlus_sig.h"
#include "reaction_routine.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"

double GetKPlus_CrossSection(double ft, double fw, double fqsq, double feps, TString fHadron) {

  double_t sig_total;
  
  double w,q2,t; // For the sake of testing, I redescribed the fw,fqsq,ft as w,q2,t as I alraedy had a file in my PC to calculate the crossection with these notations
  w=fw; q2=fqsq; t=ft;
     
  int w_1,w_2,q2_1,q2_2; // w_1, rounded value of the input W, w_2 is the next (or previous) value of W. q2_1 is the rounded input value of Q2, q2_2 is the next (or previous) value of Q2 in the array
  w_1 = round(w);  
  q2_1 = round(q2);

  //....................................................................................................................................................................
  // Below are set of if and else if conditions for boundary values of w and q2 
  //....................................................................................................................................................................
    
  if (w == w_1){
    w_1 = w;
    w_2 = w;
    
    if (q2 > q2_1){
      q2_2 = q2_1 + 1;
    }
    
    else if (q2 < q2_1){
      q2_2 = q2_1;
      q2_1 = q2_2 -1;
    }
     
    else if (q2 == q2_1 ){
      q2_1 = q2;
      q2_2 = q2;}     
  }
   
  //................................................................................. 
   
  else if (w > w_1){
    w_2 = w_1 + 1;
    
    if (q2 > q2_1){
      q2_2 = q2_1 + 1;
    }
    
    else if (q2 < q2_1){
      q2_2 = q2_1;
      q2_1 = q2_2 -1;
    }
     
    else if (q2 == q2_1 ){
      q2_1 = q2;
      q2_2 = q2;}
  }
    
  //................................................................................... 
  else if (w < w_1){
    w_2 =w_1;
    w_1 = w_2 - 1;
      
    if (q2 > q2_1){
      q2_2 = q2_1 + 1;
    }
    
    else if (q2 < q2_1){
      q2_2 = q2_1;
      q2_1 = q2_2 -1;
    }
     
    else if (q2 == q2_1 ){
      q2_1 = q2;
      q2_2 = q2;}
  }

  //.................................................................................................................................................................... 
  // Calcualtion of t-min values at each w and q2 (i.e. at all the four combinations of w's and q2's ) 
  //.................................................................................................................................................................... 
  // t1 is the t_min value UNLESS the first fit fails
  // If the first fit fails, t2 is the t_min value. If the first fit did NOT fail, t2 is the intersection point between fits 1 and 2
  // If the second fit fails, t3 is the t_min value. If the second fit did NOT fail, t3 is the intersection point between fits 2 and 3
  
  int x_1,x_2,y_1,y_2; // Redefination of w's and q2's as the array start from 0 but w starts from 2 and the q2 starts from 1 in the file.
  x_1 = w_1-2;
  x_2 = w_2-2;
  y_1 = q2_1-1;
  y_2 = q2_2-1;  // Till this point variables are same for sigT and sigL

  //....................................................................................................................................................................     
  // Variables that need to be defined globally somewhere ( These variables are only for sigT and I have used different variables for sigL)
  double t_1,t_2,t_3,t_4; // 4 corners (points) in t to interpolate over
  double sigT1,sigT2,sigT3,sigT4; // Value of SigT at each corner of the square to interpolate over
  double lsigT1,lsigT2,lsigT3,lsigT4, slb,slt,sll,slr,altb,allr,fsigTLa,fsigTa; // Logarithm of sigma values and slopes along the bottom, top, left, right and the average of the slope along the top/bottom and left/right. fsigTLa is the result in logarithmic form and fsigTa is the result
  double sigT11,sigT14; // Bottom left and top right sigma points, these are calculated if there are less than 3 (or 2) points to interpolate over
  //...................................................................................................................................................................
 
 
  if(SigPar[1][x_2][y_1][10] != -10001){ // t_3 is t_min for the third corner - this is always the top left corner
    t_3 = SigPar[1][x_2][y_1][10];
  }
  else if (SigPar[1][x_2][y_1][11] != -10001){
    t_3 = SigPar[1][x_2][y_1][11];
  }
  else if (SigPar[1][x_2][y_1][12] != -10001){
    t_3 = SigPar[1][x_2][y_1][12];
  }
   else {
    return -100;
  }
 
  //.................................................................................................................................................................... 
  // Calcualtion of sigT's at all the corners of the square 
  //....................................................................................................................................................................  
 
  if (t>= t_3 && t<2.0){ // t_3 corner is the one that we will loose at the end

    //Calculating the sigT1 at bottom left corner of the square..................................................................................
    if (t>=SigPar[1][x_1][y_1][10] && t<SigPar[1][x_1][y_1][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_1][2]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_1][3]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_1][4]);
      sigT1=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_1][y_1][11] && t<SigPar[1][x_1][y_1][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_1][5]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_1][6]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_1][7]);
      sigT1=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_1][y_1][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_1][y_1][8]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_1][9]);
      sigT1=parasigT->Eval(t); 
	 
      if (SigPar[1][x_1][y_1][8] == 0 && SigPar[1][x_1][y_1][9] == 0){
	sigT1=0;
      }
    }
	 
    else {
      sigT1=0;
    }

    //Calculating the sigT2 at bottom right corner of the square..................................................................................
    if (t>=SigPar[1][x_1][y_2][10] && t<SigPar[1][x_1][y_2][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_2][2]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_2][3]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_2][4]);
      sigT2=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_1][y_2][11] && t<SigPar[1][x_1][y_2][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_1][y_2][5]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_2][6]);
      parasigT->FixParameter(2, SigPar[1][x_1][y_2][7]);
      sigT2=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_1][y_2][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_1][y_2][8]); 
      parasigT->FixParameter(1, SigPar[1][x_1][y_2][9]);
      sigT2=parasigT->Eval(t); 
	 
      if (SigPar[1][x_1][y_2][8] == 0 && SigPar[1][x_1][y_2][9] == 0){
	sigT2=0;
      }
    }
	 
    else {
      sigT2=0;
    }

    //Calculating the sigT3 at top left corner of the square..................................................................................
    if (t>=SigPar[1][x_2][y_1][10] && t<SigPar[1][x_2][y_1][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_1][2]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_1][3]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_1][4]);
      sigT3=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_2][y_1][11] && t<SigPar[1][x_2][y_1][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_1][5]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_1][6]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_1][7]);
      sigT3=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_2][y_1][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_2][y_1][8]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_1][9]);
      sigT3=parasigT->Eval(t); 
	 
      if (SigPar[1][x_2][y_1][8] == 0 && SigPar[1][x_2][y_1][9] == 0){
	sigT3=0;
      }
    }
	 
    else {
      sigT3=0;
    }

    //Calculating the sigT4 at top right corner of the square..................................................................................

    if (t>=SigPar[1][x_2][y_2][10] && t<SigPar[1][x_2][y_2][11]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_2][2]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_2][3]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_2][4]);
      sigT4=parasigT->Eval(t);
    }
	 
    else if (t>=SigPar[1][x_2][y_2][11] && t<SigPar[1][x_2][y_2][12]){
      TF1* parasigT= new TF1("parasigT","pol2");
      parasigT->FixParameter(0, SigPar[1][x_2][y_2][5]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_2][6]);
      parasigT->FixParameter(2, SigPar[1][x_2][y_2][7]);
      sigT4=parasigT->Eval(t); 
    }
	 
    else if (t>=SigPar[1][x_2][y_2][12] && t<2.0){
      TF1* parasigT= new TF1("parasigT","expo");
      parasigT->FixParameter(0, SigPar[1][x_2][y_2][8]); 
      parasigT->FixParameter(1, SigPar[1][x_2][y_2][9]);
      sigT4=parasigT->Eval(t); 
	 
      if (SigPar[1][x_2][y_2][8] == 0 && SigPar[1][x_2][y_2][9] == 0){
	sigT4=0;
      }
    }
	 
    else {
      sigT4=0;
    }
    //....................................................................................................................................................................
    // Difterent if and else conditions to find the crossection values at the given points.
    //....................................................................................................................................................................
    
    if (sigT1 == sigT2 && sigT2 == sigT3 && sigT3 == sigT4 && sigT4 == sigT1){ // if the w and q2 will have whole number values
      fsigTa=sigT1;
    }   
    
    //...................................................................................

    else if (sigT1 != 0 && sigT2 != 0 && sigT3 != 0 && sigT4 != 0){ // if all the four corners are  present
   
      // Taking the log of claculated sigT values    
      lsigT1 = TMath::Log(sigT1); //log value of sigT1.
      lsigT2 = TMath::Log(sigT2); //log value of sigT2.
      lsigT3 = TMath::Log(sigT3); //log value of sigT3.
      lsigT4 = TMath::Log(sigT4); //log value of sigT4.
      // Calculated slopes of different lines
      slb = lsigT4-lsigT3; //interpolation from the third corner (i.e. top left corner)
      slt = lsigT2-lsigT1;
      sll = -(lsigT1-lsigT3);
      slr = -(lsigT2-lsigT4);
      // Taking averages of the slopes
      altb = ((slb +slt)/2 );
      allr = ((sll +slr)/2 );
        
      //  Applying taylor's series formula for the average slopes
      fsigTLa = lsigT3 + (q2-q2_1)*altb + (w-w_2)*allr;
      
      // Find the anti-log of the taylor's series formula value  
      fsigTa = exp(fsigTLa);
    }
    //...................................................................................

    else if (sigT1 == 0 && sigT2 == 0 && sigT4 == 0)  // if we loose the 1st, 2nd and the 4th corner simultaneously 
      {
	// In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) will always be there and for other two corners, find the value of the cross-section at the first and the fourth corner at the minimum value of t. After that we can interpolate them.

// First try to find t_1 and t_4 

 if(SigPar[1][x_1][y_1][10] != -10001){ // t_1 is t_min for the first corner - this is always the bottom left corner
    t_1 = SigPar[1][x_1][y_1][10];
  }
  else if (SigPar[1][x_1][y_1][11] != -10001){
    t_1 =  SigPar[1][x_1][y_1][11];
  } 
  else if ( SigPar[1][x_1][y_1][12] != -10001){
    t_1 =  SigPar[1][x_1][y_1][12];
  }
  else {
    return -100;
  }
 
//...................................................................................  
if(SigPar[1][x_2][y_2][10] != -10001){  // t_4 is t_min for the fourth corner - this is always the top right corner
    t_4 = SigPar[1][x_2][y_2][10];
  }
  else if (SigPar[1][x_2][y_2][11] != -10001){
    t_4 =SigPar[1][x_2][y_2][11];
  }
  else if (SigPar[1][x_2][y_2][12] != -10001){
    t_4 = SigPar[1][x_2][y_2][12];
  }
  else {
    return -100;
  }

//Calculating the sigT11 at bottom left corner of the square      
	if (t_1>=SigPar[1][x_1][y_1][10]&& t_1<SigPar[1][x_1][y_1][11]){ 
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0, SigPar[1][x_1][y_1][2]); 
	  parasigT->FixParameter(1, SigPar[1][x_1][y_1][3]);
	  parasigT->FixParameter(2, SigPar[1][x_1][y_1][4]);
	  sigT11=parasigT->Eval(t_1);
	}
	 
	else if (t_1>=SigPar[1][x_1][y_1][11] && t_1<SigPar[1][x_1][y_1][12]){
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0,SigPar[1][x_1][y_1][5]); 
	  parasigT->FixParameter(1,SigPar[1][x_1][y_1][6]);
	  parasigT->FixParameter(2,SigPar[1][x_1][y_1][7]);
	  sigT11=parasigT->Eval(t_1); 
	} 
	
	else if (t_1>=SigPar[1][x_1][y_1][12] && t_1<2.0){
	  TF1* parasigT= new TF1("parasigT","expo");
	  parasigT->FixParameter(0, SigPar[1][x_1][y_1][8]); 
	  parasigT->FixParameter(1, SigPar[1][x_1][y_1][9]);
	  sigT11=parasigT->Eval(t_1); 
	  
	  if (SigPar[1][x_1][y_1][8] == 0 && SigPar[1][x_1][y_1][9] == 0){
	    sigT11 =0;
	  }
	}
	 
	else {
	  sigT11=0;
	}

	//Calculating the sigT14 at bottom left corner of the square
   
	if (t_4>=SigPar[1][x_2][y_2][10]&& t_4<SigPar[1][x_2][y_2][11]){ 
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0, SigPar[1][x_2][y_2][2]); 
	  parasigT->FixParameter(1, SigPar[1][x_2][y_2][3]);
	  parasigT->FixParameter(2, SigPar[1][x_2][y_2][4]);
	  sigT14=parasigT->Eval(t_4);
	}
	 
	else if (t_4>=SigPar[1][x_2][y_2][11] && t_1<SigPar[1][x_2][y_2][12]){
	  TF1* parasigT= new TF1("parasigT","pol2");
	  parasigT->FixParameter(0,SigPar[1][x_2][y_2][5]); 
	  parasigT->FixParameter(1,SigPar[1][x_2][y_2][6]);
	  parasigT->FixParameter(2,SigPar[1][x_2][y_2][7]);
	  sigT14=parasigT->Eval(t_4); 
	} 
	
	else if (t_4>=SigPar[1][x_1][y_1][12] && t_1<2.0){
	  TF1* parasigT= new TF1("parasigT","expo");
	  parasigT->FixParameter(0, SigPar[1][x_2][y_2][8]); 
	  parasigT->FixParameter(1, SigPar[1][x_2][y_2][9]);
	  sigT14=parasigT->Eval(t_4); 
	  
	  if (SigPar[1][x_2][y_2][8] == 0 && SigPar[1][x_2][y_2][9] == 0){
	    sigT14 =0;
	  }
	}
	 
	else {
	  sigT14=0;
	}

	// Taking the log of claculated sigT values
	lsigT1 = TMath::Log(sigT11); //log value of sigT11.
	lsigT3 = TMath::Log(sigT3); //log value of sigT3.
	lsigT4 = TMath::Log(sigT14); //log value of sigT4.
	// Calculated slopes of different lines
	slb = lsigT4-lsigT3; //->interpolation from the third corner
	sll = -(lsigT1-lsigT3);
      
	// Applying taylor's series formula without averaging the slopes        
	fsigTLa = lsigT3 + (q2-q2_1)*slb + (w-w_2)*sll; //->interpolation from the third corner
  
	// Find the anti-log of the taylor's series formula value 
        fsigTa = exp(fsigTLa);
      }
    //...................................................................................

    else if (sigT1 == 0 && sigT2 == 0){  // if we loose the 1st and the 2nd corner simultaneously
      // In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) and the fourth corner (i.e. top right)will always be there and for other one corner, find the value of the cross-section at the first corner at the minimum value of t. After that we can interpolate them.

// First try to find t_1

 if(SigPar[1][x_1][y_1][10] != -10001){ // t_1 is t_min for the first corner - this is always the bottom left corner
    t_1 = SigPar[1][x_1][y_1][10];
  }
  else if (SigPar[1][x_1][y_1][11] != -10001){
    t_1 =  SigPar[1][x_1][y_1][11];
  } 
  else if ( SigPar[1][x_1][y_1][12] != -10001){
    t_1 =  SigPar[1][x_1][y_1][12];
  }
  else {
    return -100;
  } 
     
   //Calculating the sigT11 at bottom left corner of the square
      if (t_1>=SigPar[1][x_1][y_1][10]&& t_1<SigPar[1][x_1][y_1][11]){ 
	TF1* parasigT= new TF1("parasigT","pol2");
	parasigT->FixParameter(0, SigPar[1][x_1][y_1][2]); 
	parasigT->FixParameter(1, SigPar[1][x_1][y_1][3]);
	parasigT->FixParameter(2, SigPar[1][x_1][y_1][4]);
	sigT11=parasigT->Eval(t_1);
      }
	 
      else if (t_1>=SigPar[1][x_1][y_1][11] && t_1<SigPar[1][x_1][y_1][12]){
	TF1* parasigT= new TF1("parasigT","pol2");
	parasigT->FixParameter(0,SigPar[1][x_1][y_1][5]); 
	parasigT->FixParameter(1,SigPar[1][x_1][y_1][6]);
	parasigT->FixParameter(2,SigPar[1][x_1][y_1][7]);
	sigT11=parasigT->Eval(t_1); 
      } 
	
      else if (t_1>=SigPar[1][x_1][y_1][12] && t_1<2.0){
	TF1* parasigT= new TF1("parasigT","expo");
	parasigT->FixParameter(0, SigPar[1][x_1][y_1][8]); 
	parasigT->FixParameter(1, SigPar[1][x_1][y_1][9]);
	sigT11=parasigT->Eval(t_1); 	  
	if (SigPar[1][x_1][y_1][8] == 0 && SigPar[1][x_1][y_1][9] == 0){
	  sigT11 =0;
	}
      }
	 
      else {
	sigT11=0;
      }
      // Taking the log of claculated sigT values
       
      lsigT1 = TMath::Log(sigT11); //log value of sigT11.
      lsigT3 = TMath::Log(sigT3); //log value of sigT3.
      lsigT4 = TMath::Log(sigT4); //log value of sigT4.
      // Calculated slopes of different lines
      slb = lsigT3-lsigT4; //->interpolation from the third corner
      sll = -(lsigT1-lsigT3);
      // Applying taylor's series formula without averaging the slopes
      fsigTLa = lsigT3 + (q2-q2_1)*slb + (w-w_2)*sll; //->interpolation from the third corner
        
      // Find the anti-log of the taylor's series formula value         
      fsigTa = exp(fsigTLa);
    }          

    //..............................................................................

    else if (sigT2 == 0) { // if we loose 2nd corner, first we will always loose this corner as this correspond to highest -t value (interpolate from 3rd corner)
      // In this case, we will need atleat three corners to find the cross-section. And even after loosing second corner, we still have three corners to interpolate.

      // Taking the log of claculated sigT values         
      lsigT1 = TMath::Log(sigT1); //log value of sigT1.
      lsigT3 = TMath::Log(sigT3); //log value of sigT3.
      lsigT4 = TMath::Log(sigT4); //log value of sigT4.
    
      // Calculated slopes of different lines
      slb = lsigT3-lsigT4; //->interpolation from the third corner
      sll = -(lsigT1-lsigT3);
        
      // Applying taylor's series formula without averaging the slopes  
      fsigTLa = lsigT3 + (q2-q2_1)*slb + (w-w_2)*sll; //->interpolation from the third corner
      
      // Find the anti-log of the taylor's series formula value       
      fsigTa = exp(fsigTLa);
    }                 
  } // end of if statement over t
  //....................................................................................................................................................................
  else{
    //cerr<<" Invalid t-value "<<endl;
    return -100;
  }      

  //....................................................................................................................................................................
  //.................................................................................................................................................................... 
  // Calcualtion of t-min values at each w and q2 (i.e. at all the four combinations of w's and q2's ) [SigL Calculations started from here]
  //....................................................................................................................................................................
  //.................................................................................................................................................................... 
  // t1 is the t_min value UNLESS the first fit fails
  // If the first fit fails, t2 is the t_min value. If the first fit did NOT fail, t2 is the intersection point between fits 1 and 2
  // If the second fit fails, t3 is the t_min value. If the second fit did NOT fail, t3 is the intersection point between fits 2 and 3
 
  //....................................................................................................................................................................     
  // Variables that need to be defined globally somewhere ( These variables are only for sigL)
  double l_1,l_2,l_3,l_4; // 4 corners (points) in t to interpolate over
  double sigL1,sigL2,sigL3,sigL4; // Value of SigL at each corner of the square to interpolate over
  double lsigL1,lsigL2,lsigL3,lsigL4, stb,stt,stl,str,attb,atlr,fsigLLa,fsigLa; // Logarithm of sigma values and slopes along the bottom, top, left, right and the average of the slope along the top/bottom and left/right. fsigTLa is the result in logarithmic form and fsigTa is the result
  double sigL11,sigL14; // Bottom left and top right sigma points, these are calculated if there are less than 3 (or 2) points to interpolate over
  //...................................................................................................................................................................
       
  if(SigPar[0][x_2][y_1][9] != -10001){ // l_3 is t_min for the third corner - this is always the top left corner
    l_3 = SigPar[0][x_2][y_1][9];
  }
  else if (SigPar[0][x_2][y_1][10] != -10001){
    l_3 = SigPar[0][x_2][y_1][10];
  }
  else if (SigPar[0][x_2][y_1][11] != -10001){
    l_3 = SigPar[0][x_2][y_1][11];
  }
   else {
    return -100;
  } 
  
  //.................................................................................................................................................................... 
  // Calcualtion of sigL's at all the corners of the square 
  //....................................................................................................................................................................  
 
  if (t>= l_3 && t<2.0){ // l_3 corner is the one that we will loose at the end

    //Calculating the sigL1 at bottom left corner of the square..................................................................................
    if (t>=SigPar[0][x_1][y_1][9] && t<SigPar[0][x_1][y_1][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	sigL1=parasigL->Eval(t);
	if (SigPar[0][x_1][y_1][2] == 0 && SigPar[0][x_1][y_1][3] == 0){
	  sigL1=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	parasigL->FixParameter(2, SigPar[0][x_1][y_1][4]);
	sigL1=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_1][10] && t<SigPar[0][x_1][y_1][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_1][5]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_1][6]);
      sigL1=parasigL->Eval(t);
      if (SigPar[0][x_1][y_1][5] == 0 && SigPar[0][x_1][y_1][6] == 0){
	sigL1=0;
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_1][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_1][7]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_1][8]);
      sigL1=parasigL->Eval(t);
      if (SigPar[0][x_1][y_1][7] == 0 && SigPar[0][x_1][y_1][8] == 0){
	sigL1=0;
      }
    }
	 
    else {
      sigL1=0;
    }

    //Calculating the sigL2 at bottom right corner of the square..................................................................................

    if (t>=SigPar[0][x_1][y_2][9] && t<SigPar[0][x_1][y_2][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_2][3]);
	sigL2=parasigL->Eval(t);
	if (SigPar[0][x_1][y_2][2] == 0 && SigPar[0][x_1][y_2][3] == 0){
	  sigL2=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_1][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_2][3]);
	parasigL->FixParameter(2, SigPar[0][x_1][y_2][4]);
	sigL2=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_2][10] && t<SigPar[0][x_1][y_2][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_2][5]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_2][6]);
      sigL2=parasigL->Eval(t); 
      if (SigPar[0][x_1][y_2][5] == 0 && SigPar[0][x_1][y_2][6] == 0){
	sigL2=0;
      }
    }
	 
    else if (t>=SigPar[0][x_1][y_2][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_1][y_2][7]); 
      parasigL->FixParameter(1, SigPar[0][x_1][y_2][8]);
      sigL2=parasigL->Eval(t);
      if (SigPar[0][x_1][y_2][7] == 0 && SigPar[0][x_1][y_2][8] == 0){
	sigL2=0;
      }
    }
	 
    else {
      sigL2=0;
    }

    //Calculating the sigL3 at top left corner of the square..................................................................................
    if (t>=SigPar[0][x_2][y_1][9] && t<SigPar[0][x_2][y_1][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_2][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_1][3]);
	sigL3=parasigL->Eval(t); 
	 
	if (SigPar[0][x_2][y_1][2] == 0 && SigPar[0][x_2][y_1][3] == 0){
	  sigL3=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_2][y_1][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_1][3]);
	parasigL->FixParameter(2, SigPar[0][x_2][y_1][4]);
	sigL3=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_1][10] && t<SigPar[0][x_2][y_1][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_1][5]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_1][6]);
      sigL3=parasigL->Eval(t); 
      
      if (SigPar[0][x_2][y_1][5] == 0 && SigPar[0][x_2][y_1][6] == 0){
	sigL3=0;
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_1][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_1][7]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_1][8]);
      sigL3=parasigL->Eval(t); 
	 
      if (SigPar[0][x_2][y_1][7] == 0 && SigPar[0][x_2][y_1][8] == 0){
	sigL3=0;
      }
    }
	 
    else {
      sigL3=0;
    }

    //Calculating the sigL4 at top right corner of the square..................................................................................
    if (t>=SigPar[0][x_2][y_2][9] && t<SigPar[0][x_2][y_2][10]){
      if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	sigL4=parasigL->Eval(t); 
	 
	if (SigPar[0][x_2][y_2][2] == 0 && SigPar[0][x_2][y_2][3] == 0){
	  sigL4=0;
	}
      }
      
      else{
	TF1* parasigL= new TF1("parasigL","pol2");
	parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	parasigL->FixParameter(2, SigPar[0][x_2][y_2][4]);
	sigL4=parasigL->Eval(t);
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_2][10] && t<SigPar[0][x_2][y_2][11]){
      TF1* parasigL= new TF1("parasigL","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_2][5]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_2][6]);
      sigL4=parasigL->Eval(t); 
      
      if (SigPar[0][x_2][y_2][5] == 0 && SigPar[0][x_2][y_2][6] == 0){
	sigL4=0;
      }
    }
	 
    else if (t>=SigPar[0][x_2][y_2][11] && t<2.0){
      TF1* parasigL= new TF1("parasigT","expo");
      parasigL->FixParameter(0, SigPar[0][x_2][y_2][7]); 
      parasigL->FixParameter(1, SigPar[0][x_2][y_2][8]);
      sigL4=parasigL->Eval(t); 
	 
      if (SigPar[0][x_2][y_2][7] == 0 && SigPar[0][x_2][y_2][8] == 0){
	sigL4=0;
      }
    }
	 
    else {
      sigL4=0;
    }

    //....................................................................................................................................................................
    // Difterent if and else conditions to find the crossection values at the given points.
    //....................................................................................................................................................................
    
    if (sigL1 == sigL2 && sigL2 == sigL3 && sigL3 == sigL4 && sigL4 == sigL1){ // if the w and q2 will have whole number values
      fsigLa = sigL1;
    }   
    
    //...................................................................................

    else if (sigL1 != 0 && sigL2 != 0 && sigL3 != 0 && sigL4 != 0){ // if all the four corners are  present
    
      // Taking the log of claculated sigT values
      lsigL1 = TMath::Log(sigL1); //log value of sigT1.
      lsigL2 = TMath::Log(sigL2); //log value of sigT2.
      lsigL3 = TMath::Log(sigL3); //log value of sigT3.
      lsigL4 = TMath::Log(sigL4); //log value of sigT4.
     
      // Calculated slopes of different lines
       
      stb = lsigL4-lsigL3; //interpolation from the third corner (i.e. top left corner)
      stt = lsigL2-lsigL1;
      stl = -(lsigL1-lsigL3);
      str = -(lsigL2-lsigL4);
        
      // Taking averages of the slopes
        
      attb = ((stb +stt)/2 );
      atlr = ((stl +str)/2 );
        
      //  Applying taylor's series formula for the average slopes
        
      fsigLLa = lsigL3 + (q2-q2_1)*attb + (w-w_2)*atlr;
        
      // Find the anti-log of the taylor's series formula value
      fsigLa = exp(fsigLLa);
    }
    //...................................................................................

    else if (sigL1 == 0 && sigL2 == 0 && sigL4 == 0)  // if we loose the 1st, 2nd and the 4th corner simultaneously 
      {
	// In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) will always be there and for other two corners, find the value of the cross-section at the first and the fourth corner at the minimum value of t. After that we can interpolate them.
 
// First try to find t_1 and t_4

if(SigPar[0][x_1][y_1][9] != -10001){ // l_1 is t_min for the first corner - this is always the bottom left corner
    l_1 = SigPar[0][x_1][y_1][9];
  }
  else if (SigPar[0][x_1][y_1][10] != -10001){
    l_1 =  SigPar[0][x_1][y_1][10];
  }       
  else if ( SigPar[0][x_1][y_1][11] != -10001){
    l_1 =  SigPar[0][x_1][y_1][11];
  }
   else {
    return -100;
  } 
 //...................................................................................
 if(SigPar[0][x_2][y_2][9] != -10001){  // l_4 is t_min for the fourth corner - this is always the top right corner
    l_4 = SigPar[0][x_2][y_2][9];
  }
  else if (SigPar[0][x_2][y_2][10] != -10001){
    l_4 =SigPar[0][x_2][y_2][10];
  }
  else if (SigPar[0][x_2][y_2][11] != -10001){
    l_4 = SigPar[0][x_2][y_2][11];
  }
 else {
    return -100;
  } 
	//Calculating the sigL11 at bottom left corner of the square
	if (l_1>=SigPar[0][x_1][y_1][9] && l_1<SigPar[0][x_1][y_1][10]){
	  if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	    TF1* parasigL= new TF1("parasigL","expo");
	    parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	    sigL11=parasigL->Eval(l_1); 
	 
	    if (SigPar[0][x_1][y_1][2] == 0 && SigPar[0][x_1][y_1][3] == 0){
	      sigL11=0;
	    }
	  }
      
	  else{
	    TF1* parasigL= new TF1("parasigL","pol2");
	    parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	    parasigL->FixParameter(2, SigPar[0][x_1][y_1][4]);
	    sigL11=parasigL->Eval(l_1);
	  }
	}
	 
	else if (l_1>=SigPar[0][x_1][y_1][10] && l_1<SigPar[0][x_1][y_1][11]){
	  TF1* parasigL= new TF1("parasigL","expo");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][5]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][6]);
	  sigL11=parasigL->Eval(l_1); 
      
	  if (SigPar[0][x_1][y_1][5] == 0 && SigPar[0][x_1][y_1][6] == 0){
	    sigL11=0;
	  }
	}
	 
	else if (l_1>=SigPar[0][x_1][y_1][11] && t<2.0){
	  TF1* parasigL= new TF1("parasigT","expo");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][7]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][8]);
	  sigL11=parasigL->Eval(l_1); 
	 
	  if (SigPar[0][x_1][y_1][7] == 0 && SigPar[0][x_1][y_1][8] == 0){
	    sigL11=0;
	  }
	}
	 
	else {
	  sigL11=0;
	}   
    
	//Calculating the sigL14 at bottom left corner of the square 
	if (l_4>=SigPar[0][x_2][y_2][9] && l_4<SigPar[0][x_2][y_2][10]){
	  if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	    TF1* parasigL= new TF1("parasigL","expo");
	    parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	    sigL14=parasigL->Eval(l_4); 
	 
	    if (SigPar[0][x_2][y_2][2] == 0 && SigPar[0][x_2][y_2][3] == 0){
	      sigL14=0;
	    }
	  }
      
	  else{
	    TF1* parasigL= new TF1("parasigL","pol2");
	    parasigL->FixParameter(0, SigPar[0][x_2][y_2][2]); 
	    parasigL->FixParameter(1, SigPar[0][x_2][y_2][3]);
	    parasigL->FixParameter(2, SigPar[0][x_2][y_2][4]);
	    sigL14=parasigL->Eval(l_4);
	  }
	}
	 
	else if (l_4>=SigPar[0][x_2][y_2][10] && l_4<SigPar[0][x_2][y_2][11]){
	  TF1* parasigL= new TF1("parasigL","expo");
	  parasigL->FixParameter(0, SigPar[0][x_2][y_2][5]); 
	  parasigL->FixParameter(1, SigPar[0][x_2][y_2][6]);
	  sigL14=parasigL->Eval(l_4); 
      
	  if (SigPar[0][x_2][y_2][5] == 0 && SigPar[0][x_2][y_2][6] == 0){
	    sigL14=0;
	  }
	}
	 
	else if (l_4>=SigPar[0][x_2][y_2][11] && t<2.0){
	  TF1* parasigL= new TF1("parasigT","expo");
	  parasigL->FixParameter(0, SigPar[0][x_2][y_2][7]); 
	  parasigL->FixParameter(1, SigPar[0][x_2][y_2][8]);
	  sigL14=parasigL->Eval(l_4); 
	 
	  if (SigPar[0][x_2][y_2][7] == 0 && SigPar[0][x_2][y_2][8] == 0){
	    sigL14=0;
	  }
	}
	 
	else {
	  sigL14=0;
	}

	// Taking the log of claculated sigL values
	lsigL1 = TMath::Log(sigL11); //log value of sigL11.
	lsigL3 = TMath::Log(sigL3); //log value of sigL3.
	lsigL4 = TMath::Log(sigL14); //log value of sigL14.

	// Calculated slopes of different lines  
	stb = lsigL4-lsigL3; //->interpolation from the third corner
	stl = -(lsigL1-lsigL3);
	// Applying taylor's series formula without averaging the slopes
     	fsigLLa = lsigL3 + (q2-q2_1)*stb + (w-w_2)*stl; //->interpolation from the third corner
     
	// Find the anti-log of the taylor's series formula value 
        fsigLa = exp(fsigLLa);
      }
    //...................................................................................

    else if (sigL1 == 0 && sigL2 == 0){  // if we loose the 1st and the 2nd corner simultaneously
      // In this case, we will need atleat three corners to find the cross-section. The third corner (i.e. top left) and the fourth corner (i.e. top right)will always be there and for other one corner, find the value of the cross-section at the first corner at the minimum value of t. After that we can interpolate them.

// First try to find t_1

if(SigPar[0][x_1][y_1][9] != -10001){ // l_1 is t_min for the first corner - this is always the bottom left corner
    l_1 = SigPar[0][x_1][y_1][9];
  }
  else if (SigPar[0][x_1][y_1][10] != -10001){
    l_1 =  SigPar[0][x_1][y_1][10];
  }       
  else if ( SigPar[0][x_1][y_1][11] != -10001){
    l_1 =  SigPar[0][x_1][y_1][11];
  }
   else {
    return -100;
  } 
 //...................................................................................
      //Calculating the sigL11 at bottom left corner of the square
      if (l_1>=SigPar[0][x_1][y_1][9] && l_1<SigPar[0][x_1][y_1][10]){
	if(w_1 ==2 || w_2 ==2 || w_1 ==3 ||w_2 ==3){
	  TF1* parasigL= new TF1("parasigL","expo");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	  sigL11=parasigL->Eval(l_1); 
	 
	  if (SigPar[0][x_1][y_1][2] == 0 && SigPar[0][x_1][y_1][3] == 0){
	    sigL11=0;
	  }
	}
      
	else{
	  TF1* parasigL= new TF1("parasigL","pol2");
	  parasigL->FixParameter(0, SigPar[0][x_1][y_1][2]); 
	  parasigL->FixParameter(1, SigPar[0][x_1][y_1][3]);
	  parasigL->FixParameter(2, SigPar[0][x_1][y_1][4]);
	  sigL11=parasigL->Eval(l_1);
	}
      }
	 
      else if (l_1>=SigPar[0][x_1][y_1][10] && l_1<SigPar[0][x_1][y_1][11]){
	TF1* parasigL= new TF1("parasigL","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][5]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][6]);
	sigL11=parasigL->Eval(l_1);
	if (SigPar[0][x_1][y_1][5] == 0 && SigPar[0][x_1][y_1][6] == 0){
	  sigL11=0;
	}
      }
	 
      else if (l_1>=SigPar[0][x_1][y_1][11] && t<2.0){
	TF1* parasigL= new TF1("parasigT","expo");
	parasigL->FixParameter(0, SigPar[0][x_1][y_1][7]); 
	parasigL->FixParameter(1, SigPar[0][x_1][y_1][8]);
	sigL11=parasigL->Eval(l_1);	 
	if (SigPar[0][x_1][y_1][7] == 0 && SigPar[0][x_1][y_1][8] == 0){
	  sigL11=0;
	}
      }
	 
      else {
	sigL11=0;
      }
     
      // Taking the log of claculated sigL values
      lsigL1 = TMath::Log(sigL11); //log value of sigL11.
      lsigL3 = TMath::Log(sigL3); //log value of sigL3.
      lsigL4 = TMath::Log(sigL4); //log value of sigL4.
      // Calculated slopes of different lines 
      stb = lsigL3-lsigL4; //->interpolation from the third corner
      stl = -(lsigL1-lsigL3);
        
      // Applying taylor's series formula without averaging the slopes
      fsigLLa = lsigL3 + (q2-q2_1)*stb + (w-w_2)*stl; //->interpolation from the third corner
      
      // Find the anti-log of the taylor's series formula value        
      fsigLa = exp(fsigLLa);
    }          

    //..............................................................................
    else if (sigL2 == 0) { // if we loose 2nd corner, first we will always loose this corner as this correspond to highest -t value (interpolate from 3rd corner)
      // In this case, we will need atleat three corners to find the cross-section. And even after loosing second corner, we still have three corners to interpolate.      
      // Taking the log of claculated sigT values         
      lsigL1 = TMath::Log(sigL1); //log value of sigL1.
      lsigL3 = TMath::Log(sigL3); //log value of sigL3.
      lsigL4 = TMath::Log(sigL4); //log value of sigL4.

      // Calculated slopes of different lines
      stb = lsigL3-lsigL4; //->interpolation from the third corner
      stl = -(lsigL1-lsigL3);
        
      // Applying taylor's series formula without averaging the slopes
      fsigLLa = lsigL3 + (q2-q2_1)*stb + (w-w_2)*stl; //->interpolation from the third corner
   
      // Find the anti-log of the taylor's series formula value         
      fsigLa = exp(fsigLLa);
    }                 
  } // end of if statement over t
  //....................................................................................................................................................................
  else{
    //cerr<<" Invalid t-value "<<endl;
   return -100;
  } 
  //....................................................................................................................................................................
  //.................................................................................................................................................................... 

  //sig_total = 5222;
  sig_total = fsigTa +(feps*fsigLa); 

  return sig_total;

}
