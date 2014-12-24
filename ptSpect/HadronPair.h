#ifndef HADRONPAIR_H
#define HADRONPAIR_H
#include "event/BelleEvent.h"
#include "belle.h"
#include "ptSpect/AnaDefs.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "ptSpect/ptSpect.h"
#include "ptSpect/ParticleInfo.h"
#include "ptSpect/AuxFunc.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  //  class ptSpect;
  // float ptSpect::getPhi(const Hep3Vector& axis, const Hep3Vector& input);
class HadronPair
{
 public:
  HadronPair()
    {
      hadCharge=AnaDef::NA;
      hadPType=AnaDef::UNKNOWN;

      hadCharge1=AnaDef::SH_ChargeUnknown;
      hadCharge2=AnaDef::SH_ChargeUnknown;

      hadPType1=AnaDef::SH_TypeUnknown;
      hadPType2=AnaDef::SH_TypeUnknown;

    };

  Particle* firstHadron;
  Particle* secondHadron;
  ///for ease of access
  float z1;
  float z2;

  float z;
  Hep3Vector P_h; //sum of the momenta
  Hep3Vector PDiff;

  double mass; //inv mass of two hadron system
  double kT;
  //only for mc
  double diffPhi;
  double diffTheta;
  float qT;
  AnaDef::TwoHadCharge hadCharge;
  AnaDef::TwoHadPType hadPType;

  AnaDef::SingleHadCharge hadCharge1;
  AnaDef::SingleHadCharge hadCharge2;

  AnaDef::SingleHadType hadPType1;
  AnaDef::SingleHadType hadPType2;


  //everything relative to thrust

  //set the various values based on the two hadrons
  void compute()
  {
      HepLorentzVector vPhoton=kinematics::firstElectronCM+kinematics::secondElectronCM;
      HepLorentzVector vR1=firstHadron->p();
      HepLorentzVector vR2=secondHadron->p();
      if(vR1.vect().mag()==0 || vR2.vect().mag() ==0)
	return;
      HepLorentzVector RSum=vR1+vR2;
      HepLorentzVector RSumBoosted=RSum;

      HepLorentzVector R1Boosted=vR1;
      HepLorentzVector R2Boosted=vR2;
      Hep3Vector rBoost=RSum.boostVector();

      vPhoton.boost(-rBoost);
      RSumBoosted.boost(-rBoost);
      R1Boosted.boost(-rBoost);
      R2Boosted.boost(-rBoost);
      qT=vPhoton.perp(R1Boosted.vect());



      //      cout <<"rsum boosted: "<< RSumBoosted.vect() << " r1 boosted: "<< R1Boosted.vect() <<" r2: " << R2Boosted.vect() <<endl;
      ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(firstHadron->userInfo());
      ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>(secondHadron->userInfo());
      z1=pinf.z;
      z2=pinf2.z;
      z=pinf.z+pinf2.z;
      double E1,E2;
      E1=firstHadron->p().e();
      E2=secondHadron->p().e();
      mass=sqrt((E1+E2)*(E1+E2)-(firstHadron->p().vect()+secondHadron->p().vect())*(firstHadron->p().vect()+secondHadron->p().vect()));
      P_h=firstHadron->p().vect()+secondHadron->p().vect();
      PDiff=firstHadron->p().vect()-secondHadron->p().vect();
      diffPhi=PDiff.phi();
      diffTheta=PDiff.theta();
      kT=firstHadron->p3().perp(secondHadron->p3());


      //set charges/particle types...
      hadCharge1=getHadCharge(firstHadron->lund());
      hadCharge2=getHadCharge(firstHadron->lund());

      hadPType1=getHadType(firstHadron->lund());
      hadPType2=getHadType(firstHadron->lund());

      hadCharge=getTwoHadCharge(hadCharge1,hadCharge2);
      hadPType=getTwoHadType(hadPType1,hadPType2);

  }


  AnaDef::SingleHadCharge getHadCharge(int lund)
    {
      if(lund==111)
	return AnaDef::Neut;
      if(lund==221)
	return AnaDef::Neut;
      if(lund==311)
	return AnaDef::Neut;
      if(lund==22)
	return AnaDef::Neut;
      
      if(lund<0)
	return AnaDef::Neg;
      else
	return AnaDef::Pos;


    }

  AnaDef::SingleHadType getHadType(int lund)
    {
      int llund=fabs(lund);
      if(llund==13 || llund==11)
	return AnaDef::MuonElectron;
      if(llund==211)
	return AnaDef::Pion;
      if(llund==321)

	return AnaDef::Kaon;
      if(llund==2212)
	return AnaDef::Proton;

      return AnaDef::SH_TypeUnknown;
    }

  AnaDef::TwoHadCharge getTwoHadCharge(AnaDef::SingleHadCharge c1, AnaDef::SingleHadCharge c2)
    {
      if(AnaDef::SH_ChargeUnknown==c1 || AnaDef::SH_ChargeUnknown==c2)
	return AnaDef::NA;
      if((AnaDef::Pos==c1 && AnaDef::Neg==c2)|| (AnaDef::Neg==c1 && AnaDef::Pos==c2))
	return AnaDef::PN;
      if(AnaDef::Pos==c1 && AnaDef::Pos==c2)
	return AnaDef::PP;
      if(AnaDef::Neg==c1 && AnaDef::Neg==c2)
	return AnaDef::NN;

      return AnaDef::NA;
    }

  AnaDef::TwoHadPType getTwoHadType(AnaDef::SingleHadType c1, AnaDef::SingleHadType c2)
    {
      if(AnaDef::SH_TypeUnknown==c1 || AnaDef::SH_TypeUnknown==c2)
	return AnaDef::UNKNOWN;

      if(AnaDef::Pion==c1&& AnaDef::Pion==c2)
	return AnaDef::PiPi;
      if(AnaDef::Kaon==c1&& AnaDef::Kaon==c2)
	return AnaDef::KK;

      if((AnaDef::Kaon==c1&& AnaDef::Pion==c2)&&(AnaDef::Pion==c1&&AnaDef::Kaon==c2))
	return AnaDef::PiK;

      return AnaDef::UNKNOWN;
    }


};
#if defined(BELLE_NAMESPACE)
}
#endif

#endif
