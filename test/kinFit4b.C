#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

#include "TLorentzVector.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Modified version of the ExampleEtEtaPhi2CFit.C macro from CMS AN 2005/025.
// Jet error parametrization from CMS AN 2005/005.
//
// To run this macro in a root session do:
// root [0] gSystem->Load("libPhysicsToolsKinFitter.so");
// root [1] .x PhysicsTools/KinFitter/test/kinFit4b.C+
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


Double_t ErrEt(Float_t Et, Float_t Eta){
Double_t InvPerr2, a, b, c;
const int NYBINS = 5;
double YBND[NYBINS+1] = {0,0.5,1.0,1.5,2.0,2.5};
double PAR[NYBINS][3] = {
  {3.51, 0.826, 0.0364},{3.23, 0.851, 0.0367},{4.36, 0.871, 0.0415},
  {5.22, 0.713, 0.0229},{5.07, 0.610, 0.0207}};
int iy =0;
if (fabs(Eta) < 2.5){
if (fabs(Eta) < 0.5 && fabs(Eta) >0.) iy =0;
if (fabs(Eta) < 1. && fabs(Eta) >0.5) iy =1;
if (fabs(Eta) < 1.5 && fabs(Eta) >1.) iy =2;
if (fabs(Eta) < 2. && fabs(Eta) >1.5) iy =3;
if (fabs(Eta) < 2. && fabs(Eta) >2.5) iy =4;
///std::cout<< iy << " kin  iy "<<std::endl;

InvPerr2 = sqrt(pow(PAR[iy][0]/Et,2) + pow(PAR[iy][1],2)/Et + pow(PAR[iy][2],2));
//std::cout<< InvPerr2 << " kin  invPerr "<<std::endl;
}
else if(fabs(Eta) > 2.5){
  a = 4.8;
  b = 0.89;
  c = 0.043;
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
}

return InvPerr2;
}


Double_t ErrEta(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

Double_t ErrPhi(Float_t Et, Float_t Eta) {
  Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}

void print(TKinFitter *fitter)
{
  std::cout << "=============================================" << std ::endl;
  std::cout << "-> Number of measured Particles  : " << fitter->nbMeasParticles() << std::endl;
  std::cout << "-> Number of unmeasured particles: " << fitter->nbUnmeasParticles() << std::endl;
  std::cout << "-> Number of constraints         : " << fitter->nbConstraints() << std::endl;
  std::cout << "-> Number of degrees of freedom  : " << fitter->getNDF() << std::endl;
  std::cout << "-> Number of parameters A        : " << fitter->getNParA() << std::endl;
  std::cout << "-> Number of parameters B        : " << fitter->getNParB() << std::endl;
  std::cout << "-> Maximum number of iterations  : " << fitter->getMaxNumberIter() << std::endl;
  std::cout << "-> Maximum deltaS                : " << fitter->getMaxDeltaS() << std::endl;
  std::cout << "-> Maximum F                     : " << fitter->getMaxF() << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std ::endl;
  std::cout << "-> Status                        : " << fitter->getStatus() << std::endl;
  std::cout << "-> Number of iterations          : " << fitter->getNbIter() << std::endl;
  std::cout << "-> Chi2                             : " << fitter->getS() << std::endl;
  std::cout << "-> |constraints|                             : " << fitter->getF() << std::endl;
  std::cout << "=============================================" << std ::endl;
}

void kinFit4b()
{

  TLorentzVector v1(-77.92, 16.24, 117.64, 142.87);
  TLorentzVector v2( 15.41, 28.78,   6.06,  34.08);

  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  m1.Zero();
  m2.Zero();

  //In this example the covariant matrix depends on the transverse energy and eta of the jets
  m1(0,0) = ErrEt (v1.Et(), v1.Eta()); // et
  m1(1,1) = ErrEta(v1.Et(), v1.Eta()); // eta
  m1(2,2) = ErrPhi(v1.Et(), v1.Eta()); // phi
  m2(0,0) = ErrEt (v2.Et(), v2.Eta()); // et
  m2(1,1) = ErrEta(v2.Et(), v2.Eta()); // eta
  m2(2,2) = ErrPhi(v2.Et(), v2.Eta()); // phi
  TFitParticleEtEtaPhi *jet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &v1, &m1 );
  TFitParticleEtEtaPhi *jet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &v2, &m2 );

  //vec1 and vec2 must make a W boson
  TFitConstraintM *mCons1 = new TFitConstraintM( "WMassConstraint", "WMass-Constraint", 0, 0 , 125.6);
  mCons1->addParticles1( jet1, jet2 );

  //Definition of the fitter
  //Add three measured particles(jets)
  //Add two constraints
  TKinFitter* fitter = new TKinFitter("fitter", "fitter");
  fitter->addMeasParticle( jet1 );
  fitter->addMeasParticle( jet2 );
  fitter->addConstraint( mCons1 );

  //Set convergence criteria
  fitter->setMaxNbIter( 30 );
  fitter->setMaxDeltaS( 1e-2 );
  fitter->setMaxF( 1e-1 );
  fitter->setVerbosity(1);

  //Perform the fit
  std::cout << "Performing kinematic fit..." << std::endl;
  print(fitter);
  fitter->fit();
  std::cout << "Done." << std::endl;
  print(fitter);

  delete jet1;
  delete jet2;
  delete mCons1;
  delete fitter;

}

double chi2(TLorentzVector v1, TLorentzVector v2){

 double result = -1;

  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  m1.Zero();
  m2.Zero();

  //In this example the covariant matrix depends on the transverse energy and eta of the jets
  m1(0,0) = ErrEt (v1.Et(), v1.Eta()); // et
  m1(1,1) = ErrEta(v1.Et(), v1.Eta()); // eta
  m1(2,2) = ErrPhi(v1.Et(), v1.Eta()); // phi
  m2(0,0) = ErrEt (v2.Et(), v2.Eta()); // et
  m2(1,1) = ErrEta(v2.Et(), v2.Eta()); // eta
  m2(2,2) = ErrPhi(v2.Et(), v2.Eta()); // phi
  TFitParticleEtEtaPhi *jet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &v1, &m1 );
  TFitParticleEtEtaPhi *jet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &v2, &m2 );

  //vec1 and vec2 must make a W boson
  TFitConstraintM *mCons1 = new TFitConstraintM( "WMassConstraint", "WMass-Constraint", 0, 0 , 125.6);
  mCons1->addParticles1( jet1, jet2 );

  //Definition of the fitter
  //Add three measured particles(jets)
  //Add two constraints
  TKinFitter* fitter = new TKinFitter("fitter", "fitter");
  fitter->addMeasParticle( jet1 );
  fitter->addMeasParticle( jet2 );
  fitter->addConstraint( mCons1 );

  //Set convergence criteria
  fitter->setMaxNbIter( 30 );
  fitter->setMaxDeltaS( 1e-2 );
  fitter->setMaxF( 1e-1 );
  fitter->setVerbosity(1);

  //Perform the fit
  std::cout << "Performing kinematic fit..." << std::endl;
  print(fitter);
  fitter->fit();
  std::cout << "Done." << std::endl;
  print(fitter);

	result = fitter->getS(); 

  delete jet1;
  delete jet2;
  delete mCons1;
  delete fitter;

	return result;

}


