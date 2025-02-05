#include "TreeBuilder.hxx"
#include "Particle.hxx"
#include "Constants.hxx"
#include "DEMPEvent.hxx"

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include "json/json.h"

using namespace std;

TreeBuilder::TreeBuilder(const char * file_name, const char * name) {

  extern Json::Value obj;

  tree_name = name;

//  string str_copy = "RootFiles/" + obj["file_name"].asString() + ".root";
//  file_name = str_copy.c_str();
//
//  cout << obj["output_file"].asString() << endl;
//  cout << file_name << endl;
//  cout << file_name <<endl;


  File_Out = new TFile(file_name, "RECREATE");
  Tree_Out = new TTree(tree_name, tree_name);

  Tree_Out->SetAutoSave();

  nParticles = 0;

  nVars = 8;
}

void TreeBuilder::AddParticle(Particle * p)
{
  ParticleList.push_back(p);
  vector< double* >* TempVec = new vector<double*>(nVars);
  // Vector elements must be initialized to non-null
  // pointers
  for (int i = 0; i<nVars; i++)
    TempVec->at(i) = new double();

  PropList.push_back(TempVec);

  char * p_name = p->GetName();
  char b_name[100]; // Branch name
  char l_list[100]; // Leaf list

  strcpy(b_name, p_name);
  strcat(b_name, "_pid");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(0),l_list);

  strcpy(b_name, p_name);
  strcat(b_name, "_Px");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(1),l_list);

  strcpy(b_name, p_name);
  strcat(b_name, "_Py");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(2),l_list);

  strcpy(b_name, p_name);
  strcat(b_name, "_Pz");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(3),l_list);


  strcpy(b_name, p_name);
  strcat(b_name, "_E");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(4),l_list);

  strcpy(b_name, p_name);
  strcat(b_name, "_theta");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(5),l_list);

  strcpy(b_name, p_name);
  strcat(b_name, "_phi");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(6),l_list);

  strcpy(b_name, p_name);
  strcat(b_name, "_P");
  strcpy(l_list, b_name);
  strcat(l_list, "/D");
  Tree_Out->Branch(b_name,TempVec->at(7),l_list);

  nParticles++;

  //cout << ParticleList[0] << endl;
  //cout << TempVec->at(1) << endl;
  //cout << PropList.at(0)->at(1) << endl;
}

void TreeBuilder::Fill()
{
  Retrieve();
  Tree_Out->Fill();
}

// Particles have private data members that must be retrieved
// to be copied to an accessible address to be added to the
// TTree.
void TreeBuilder::Retrieve()
{
  for (int i = 0; i < nParticles; i++){
    //cout << PropList.at(i)->at(1) <<endl;
    *PropList.at(i)->at(0) = (double)ParticleList.at(i)->GetPid();
    *PropList.at(i)->at(1) = ParticleList.at(i)->Px()/1000;
    *PropList.at(i)->at(2) = ParticleList.at(i)->Py()/1000;
    *PropList.at(i)->at(3) = ParticleList.at(i)->Pz()/1000;
    *PropList.at(i)->at(4) = ParticleList.at(i)->E()/1000;
    *PropList.at(i)->at(5) = ParticleList.at(i)->Theta()*constants::DEG;
    *PropList.at(i)->at(6) = ParticleList.at(i)->Phi()*constants::DEG;
    *PropList.at(i)->at(7) = ParticleList.at(i)->P()/1000;
    //cout <<  &PropList.at(i)->at(1) << endl;
  }
}

void TreeBuilder::Save()
{
  Tree_Out->AutoSave();
  //File_Out->Write();
}

void TreeBuilder::AddDouble(double* x, const char* name)
{
  char b_name[100];
  char l_list[100];

  strcpy(b_name, name);
  strcpy(l_list, name);
  strcat(l_list,"/D");

  Tree_Out->Branch(b_name, x, l_list);

}

void TreeBuilder::AddEvent(DEMPEvent * event)
{
  this->AddParticle(event->BeamElec);
  this->AddParticle(event->ProdPion);
  this->AddParticle(event->ProdProt);
  this->AddParticle(event->ScatElec);
  this->AddParticle(event->TargNeut);
  //this->AddParticle(event->VirtPhot);

  /*
  this->AddDouble(event->qsq_GeV,"qsq_GeV");
  this->AddDouble(event->w_GeV,"w_GeV");
  this->AddDouble(event->t_GeV,"t_GeV");
  this->AddDouble(event->t_prime_GeV, "t_prime_GeV");
  this->AddDouble(event->x_B, "x_B");
  this->AddDouble(event->negt, "negt");
  
  this->AddDouble(event->Vertex_x,"Vertex_x");
  this->AddDouble(event->Vertex_y,"Vertex_y");
  this->AddDouble(event->Vertex_z,"Vertex_z");
  */
}
