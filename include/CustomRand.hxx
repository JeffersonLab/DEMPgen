/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Central class for different randomizers
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef CustomRand_H
#define CustomRand_H

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TRandom3.h"

class CustomRand
{
  TRandom3 * r;

  TF1 * PhiRand;
  TF1 * ThetaRand;
  TF1 * ERand;
  TF3 * VertexRand;

  bool PartSet = false;
  bool VertSet = false;
 public:
  double Phi();
  double Theta();
  double E();

  CustomRand(const char * fname,
             double ERange[2], double ThetaRange[2],
             double PhiRange[2]);
  CustomRand(const char * fname,
             double VRange[6]);

  CustomRand(const char * fname,
             double ERange[2], double ThetaRange[2],
             double PhiRange[2], int demp_seed);



  /* CustomRand(double ERange[2], double PhiRange[2], */
  /*            bouble PhiRange[2], */
  /*            double VRange[6]); */
};

#endif
