//#define DEBUG_EVENT 8035//please no output

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
 HadronPair():secondRun(false),thrustMethod(false)
    {
      hadCharge=AnaDef::NA;
      hadPType=AnaDef::UNKNOWN;

      hadCharge1=AnaDef::SH_ChargeUnknown;
      hadCharge2=AnaDef::SH_ChargeUnknown;

      hadPType1=AnaDef::SH_TypeUnknown;
      hadPType2=AnaDef::SH_TypeUnknown;

    };
  bool secondRun;
  bool thrustMethod;
  Particle* firstHadron;
  Particle* secondHadron;
  ///for ease of access
  float z1; 
  float z2;


  ///this is already the product, in principle we only need to save 6 values
  float p_PiPi;
  float p_PiK;
  float p_PiP;

  float p_KPi;
  float p_KK;
  float p_KP;

  float p_PPi;
  float p_PK;
  float p_PP;

  float p_PiPi2;
  float p_PiK2;
  float p_PiP2;

  float p_KPi2;
  float p_KK2;
  float p_KP2;

  float p_PPi2;
  float p_PK2;
  float p_PP2;

  float ep_PiPi;
  float ep_PiK;
  float ep_PiP;

  float ep_KPi;
  float ep_KK;
  float ep_KP;

  float ep_PPi;
  float ep_PK;
  float ep_PP;




  float z;
  Hep3Vector P_h; //sum of the momenta
  Hep3Vector PDiff;

  double mass; //inv mass of two hadron system
  float kT;


  double kT_PiPi;
  double kT_PiK;
  double kT_PiP;

  double kT_KPi;
  double kT_KK;
  double kT_KP;

  double kT_PPi;
  double kT_PK;
  double kT_PP;


  double qT_PiPi;
  double qT_PiK;
  double qT_PiP;

  double qT_KPi;
  double qT_KK;
  double qT_KP;

  double qT_PPi;
  double qT_PK;
  double qT_PP;

  float z1_PiPi;
  float z1_PiK;
  float z1_PiP;

  float z1_KPi;
  float z1_KK;
  float z1_KP;

  float z1_PPi;
  float z1_PK;
  float z1_PP;

  float z2_PiPi;
  float z2_PiK;
  float z2_PiP;

  float z2_KPi;
  float z2_KK;
  float z2_KP;


  float z2_PPi;
  float z2_PK;
  float z2_PP;


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


  //compute PID weights based on the original classification
  void computePIDWeights()
  {

  }

  float getProductUncert(float p1, float ep1,float p2, float ep2)
  {
    if(p2==0 || p1==0)
      {
	return 0.0;
      }
    return sqrt((ep1*ep1)/(p1*p1)+(ep2*ep2)/(p2*p2))*(p1*p2);

  }


  //everything relative to thrust

  float getQt(HepLorentzVector& v1, HepLorentzVector& v2)
  {
      HepLorentzVector vPhoton=kinematics::firstElectronCM+kinematics::secondElectronCM;
      HepLorentzVector vR1=v1;
      HepLorentzVector vR2=v2;
      if(vR1.vect().mag()==0 || vR2.vect().mag() ==0)
	return 0;
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
      return qT;
  }
  //get the QT based on the boosted vectors in the hadron pair. The boost will depend on the hadron type.
  float getQtOld(Hep3Vector& v1, Hep3Vector& v2, float z1, float z2)
  {
      HepLorentzVector vPhoton=kinematics::firstElectronCM+kinematics::secondElectronCM;
      HepLorentzVector vR1(v1,kinematics::Q*z1*0.5);
      HepLorentzVector vR2(v2,kinematics::Q*z2*0.5);
      if(vR1.vect().mag()==0 || vR2.vect().mag() ==0)
	return 0;
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
      return qT;
  }


  //set the various values based on the two hadrons
  void compute()
  {
    ////compute pid products, kT and z..
    ParticleInfo& pinf1=dynamic_cast<ParticleInfo&>(firstHadron->userInfo());
    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>(secondHadron->userInfo());

   if(DEBUG_EVENT==kinematics::evtNr)
      {
	cout <<"hadron pair compute p pion: " << pinf1.p_Pi <<" and " << pinf2.p_Pi<<endl;
	cout <<"dot product:  " << pinf1.boostedMoms[pionIdx].dot(pinf2.boostedMoms[pionIdx])<<endl;
      }
    //sub optimal, since the weights are not in arrays one can iterate over...
    p_PiPi=0.0;
    p_PiPi2=0.0;
    ep_PiPi=0.0;
    if(thrustMethod || (pinf1.boostedMoms[pionIdx].dot(pinf2.boostedMoms[pionIdx])<0))
      {
	p_PiPi=pinf1.p_Pi*pinf2.p_Pi;
	p_PiPi2=pinf1.p_Pi2*pinf2.p_Pi2;
	//	cout <<" looking at pion p: " << pinf1.p_Pi <<" and " << pinf2.p_Pi<<", alt: "<< pinf1.p_Pi2<<" and " <<pinf2.p_Pi2 <<" comb p: " << p_PiPi <<" and " << p_PiPi2<<endl;
	ep_PiPi=getProductUncert(pinf1.p_Pi,pinf1.ep_Pi,pinf2.p_Pi,pinf2.ep_Pi);
      }

    p_PiK=0.0;
    p_PiK2=0.0;
    ep_PiK=0.0;
    if(thrustMethod || (pinf1.boostedMoms[pionIdx].dot(pinf2.boostedMoms[kaonIdx])<0))
      {
	p_PiK=pinf1.p_Pi*pinf2.p_K;
	p_PiK2=pinf1.p_Pi2*pinf2.p_K2;
	ep_PiK=getProductUncert(pinf1.p_Pi,pinf1.ep_Pi,pinf2.p_K,pinf2.ep_K);
      }
    p_PiP=0.0;
    p_PiP2=0.0;
    ep_PiP=0.0;
    if(thrustMethod || (pinf1.boostedMoms[pionIdx].dot(pinf2.boostedMoms[protonIdx])<0))
      {
	p_PiP=pinf1.p_Pi*pinf2.p_p;
	p_PiP2=pinf1.p_Pi2*pinf2.p_p2;
	ep_PiP=getProductUncert(pinf1.p_Pi,pinf1.ep_Pi,pinf2.p_p,pinf2.ep_p);
      }

    p_KPi=0.0;
    p_KPi2=0.0;
    ep_KPi=0.0;
    if(thrustMethod || (pinf1.boostedMoms[kaonIdx].dot(pinf2.boostedMoms[pionIdx])<0))
      {
      p_KPi=pinf1.p_K*pinf2.p_Pi;
      p_KPi2=pinf1.p_K2*pinf2.p_Pi2;
      ep_KPi=getProductUncert(pinf1.p_K,pinf1.ep_K,pinf2.p_Pi,pinf2.ep_Pi);


      }
    p_KK=0.0;
    p_KK2=0.0;
    ep_KP=0.0;
    if(thrustMethod || (pinf1.boostedMoms[kaonIdx].dot(pinf2.boostedMoms[kaonIdx])<0))
      {
	p_KK=pinf1.p_K*pinf2.p_K;
	p_KK2=pinf1.p_K2*pinf2.p_K2;
	ep_KK=getProductUncert(pinf1.p_K,pinf1.ep_K,pinf2.p_K,pinf2.ep_K);
      }
    p_KP=0.0;
    p_KP2=0.0;
    ep_KP=0.0;
    if(thrustMethod || (pinf1.boostedMoms[kaonIdx].dot(pinf2.boostedMoms[protonIdx])<0))
      {
	p_KP=pinf1.p_K*pinf2.p_p;
	p_KP2=pinf1.p_K2*pinf2.p_p2;
	ep_KP=getProductUncert(pinf1.p_K,pinf1.ep_K,pinf2.p_p,pinf2.ep_p);
      }

    p_PPi=0.0;
    p_PPi2=0.0;
    ep_PPi=0.0;
    if(thrustMethod || (pinf1.boostedMoms[protonIdx].dot(pinf2.boostedMoms[pionIdx])<0))
      {
      p_PPi=pinf1.p_p*pinf2.p_Pi;
      p_PPi2=pinf1.p_p2*pinf2.p_Pi2;
      ep_PPi=getProductUncert(pinf1.p_p,pinf1.ep_p,pinf2.p_Pi,pinf2.ep_Pi);
      }
    p_PK=0.0;
    p_PK2=0.0;
    ep_PK=0.0;
    if(thrustMethod || (pinf1.boostedMoms[protonIdx].dot(pinf2.boostedMoms[kaonIdx])<0))
      {
      p_PK=pinf1.p_p*pinf2.p_K;
      p_PK2=pinf1.p_p2*pinf2.p_K2;
      ep_PK=getProductUncert(pinf1.p_p,pinf1.ep_p,pinf2.p_K,pinf2.ep_K);
      }
    p_PP=0.0;
    p_PP2=0.0;
    ep_PP=0.0;
    if(thrustMethod || (pinf1.boostedMoms[protonIdx].dot(pinf2.boostedMoms[protonIdx])<0))
      {
	p_PP=pinf1.p_p*pinf2.p_p;
	p_PP2=pinf1.p_p2*pinf2.p_p2;
	ep_PP=getProductUncert(pinf1.p_p,pinf1.ep_p,pinf2.p_K,pinf2.ep_K);
      }

    //vector<int> v={pionIdx,kaonIdx,protonIdx};
    //    for(auto& elem:v)
    //    for(vector<int>::iterator it=v.begin();it!=v.end();it++)
      {

      }
    //      kT=firstHadron->p3().perp(secondHadron->p3());

////    kT_PiPi=pinf1.boostedMoms[pionIdx].perp(pinf2.boostedMoms[pionIdx]);
////    kT_PiK=pinf1.boostedMoms[pionIdx].perp(pinf2.boostedMoms[kaonIdx]);
////    kT_PiP=pinf1.boostedMoms[pionIdx].perp(pinf2.boostedMoms[protonIdx]);
////
////
////    kT_KPi=pinf1.boostedMoms[kaonIdx].perp(pinf2.boostedMoms[pionIdx]);
////    kT_KK=pinf1.boostedMoms[kaonIdx].perp(pinf2.boostedMoms[kaonIdx]);
////    kT_KP=pinf1.boostedMoms[kaonIdx].perp(pinf2.boostedMoms[protonIdx]);
////
////    kT_PPi=pinf1.boostedMoms[protonIdx].perp(pinf2.boostedMoms[pionIdx]);
////    kT_PK=pinf1.boostedMoms[protonIdx].perp(pinf2.boostedMoms[kaonIdx]);
////    kT_PP=pinf1.boostedMoms[protonIdx].perp(pinf2.boostedMoms[protonIdx]);
////
    //charlotte switches the definition
      //lower perp higher for debug 
      if(pinf2.boostedMoms[pionIdx].mag() < pinf1.boostedMoms[pionIdx].mag())
	{
	kT_PiPi=pinf2.boostedMoms[pionIdx].perp(pinf1.boostedMoms[pionIdx]);
	}
      else
	{
	  kT_PiPi=pinf1.boostedMoms[pionIdx].perp(pinf2.boostedMoms[pionIdx]);
	}
    kT_PiK=pinf2.boostedMoms[kaonIdx].perp(pinf1.boostedMoms[pionIdx]);
    kT_PiP=pinf2.boostedMoms[protonIdx].perp(pinf1.boostedMoms[pionIdx]);
    //  kT_PiP=pinf1.boostedMoms[pionIdx].perp(pinf2.boostedMoms[protonIdx]);


    kT_KPi=pinf2.boostedMoms[pionIdx].perp(pinf1.boostedMoms[kaonIdx]);
    kT_KK=pinf2.boostedMoms[kaonIdx].perp(pinf1.boostedMoms[kaonIdx]);
        kT_KP=pinf2.boostedMoms[protonIdx].perp(pinf1.boostedMoms[kaonIdx]);
    //    kT_KP=pinf1.boostedMoms[kaonIdx].perp(pinf2.boostedMoms[protonIdx]);

    kT_PPi=pinf2.boostedMoms[pionIdx].perp(pinf1.boostedMoms[protonIdx]);
    kT_PK=pinf2.boostedMoms[kaonIdx].perp(pinf1.boostedMoms[protonIdx]);
    kT_PP=pinf2.boostedMoms[protonIdx].perp(pinf1.boostedMoms[protonIdx]);

    qT_PiPi=getQt(pinf1.boostedLorentzVec[pionIdx],pinf2.boostedLorentzVec[pionIdx]);
    qT_PiK=getQt(pinf1.boostedLorentzVec[pionIdx],pinf2.boostedLorentzVec[kaonIdx]);
    qT_PiP=getQt(pinf1.boostedLorentzVec[pionIdx],pinf2.boostedLorentzVec[protonIdx]);


    qT_KPi=getQt(pinf1.boostedLorentzVec[kaonIdx],pinf2.boostedLorentzVec[pionIdx]);
    qT_KK=getQt(pinf1.boostedLorentzVec[kaonIdx],pinf2.boostedLorentzVec[kaonIdx]);
    qT_KP=getQt(pinf1.boostedLorentzVec[kaonIdx],pinf2.boostedLorentzVec[protonIdx]);

    qT_PPi=getQt(pinf1.boostedLorentzVec[protonIdx],pinf2.boostedLorentzVec[pionIdx]);
    qT_PK=getQt(pinf1.boostedLorentzVec[protonIdx],pinf2.boostedLorentzVec[kaonIdx]);
    qT_PP=getQt(pinf1.boostedLorentzVec[protonIdx],pinf2.boostedLorentzVec[protonIdx]);

    ///////---------

      HepLorentzVector vPhoton=kinematics::firstElectronCM+kinematics::secondElectronCM;
      HepLorentzVector vPhoton2=kinematics::firstElectronCM+kinematics::secondElectronCM;
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
      //now from the getz...
      //      qT=vPhoton.perp(R1Boosted.vect());
      //      cout <<"rsum boosted: "<< RSumBoosted.vect() << " r1 boosted: "<< R1Boosted.vect() <<" r2: " << R2Boosted.vect() <<endl;
      //      ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(firstHadron->userInfo());
      //      ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>(secondHadron->userInfo());

      ;
      //  float getZ(ParticleInfo& pinf1, int massIndex1, ParticleInfo& pinf2, int massIndex2,Hep3Vector& q)
      //Hep3Vector q=vPhoton.vect();
      z1_PiPi=getZ(pinf1,pionIdx,pinf2,pionIdx,vPhoton2);
      z1_PiK=getZ(pinf1,pionIdx,pinf2,kaonIdx,vPhoton2);
      z1_PiP=getZ(pinf1,pionIdx,pinf2,protonIdx,vPhoton2);

      z2_PiPi=getZ(pinf2,pionIdx,pinf1,pionIdx,vPhoton2);
      z2_PiK=getZ(pinf2,kaonIdx,pinf1,pionIdx,vPhoton2);
      z2_PiP=getZ(pinf2,protonIdx,pinf1,pionIdx,vPhoton2);

      z1_KPi=getZ(pinf1,kaonIdx,pinf2,pionIdx,vPhoton2);
      z1_KK=getZ(pinf1,kaonIdx,pinf2,kaonIdx,vPhoton2);
      z1_KP=getZ(pinf1,kaonIdx,pinf2,protonIdx,vPhoton2);

      z2_KPi=getZ(pinf2,pionIdx,pinf1,kaonIdx,vPhoton2);
      z2_KK=getZ(pinf2,kaonIdx,pinf1,kaonIdx,vPhoton2);
      z2_KP=getZ(pinf2,protonIdx,pinf1,kaonIdx,vPhoton2);

      z1_PPi=getZ(pinf1,protonIdx,pinf2,pionIdx,vPhoton2);
      z1_PK=getZ(pinf1,protonIdx,pinf2,kaonIdx,vPhoton2);
      z1_PP=getZ(pinf1,protonIdx,pinf2,protonIdx,vPhoton2);

      z2_PPi=getZ(pinf2,pionIdx,pinf1,protonIdx,vPhoton2);
      z2_PK=getZ(pinf2,kaonIdx,pinf1,protonIdx,vPhoton2);
      z2_PP=getZ(pinf2,protonIdx,pinf1,protonIdx,vPhoton2);

      //set charges/particle types... needed for getZ_Kt
      hadCharge1=getHadCharge(firstHadron->lund());
      hadCharge2=getHadCharge(secondHadron->lund());

      hadPType1=getHadType(firstHadron->lund());
      hadPType2=getHadType(secondHadron->lund());

      hadCharge=getTwoHadCharge(hadCharge1,hadCharge2);
      hadPType=getTwoHadType(hadPType1,hadPType2);

      getZ_Kt(z1,z2,kT,qT);

      z=z1+z2;
      double E1,E2;
      E1=firstHadron->p().e();
      E2=secondHadron->p().e();
      mass=sqrt((E1+E2)*(E1+E2)-(firstHadron->p().vect()+secondHadron->p().vect())*(firstHadron->p().vect()+secondHadron->p().vect()));
      P_h=firstHadron->p().vect()+secondHadron->p().vect();
      PDiff=firstHadron->p().vect()-secondHadron->p().vect();
      diffPhi=PDiff.phi();
      diffTheta=PDiff.theta();
      //      kT=firstHadron->p3().perp(secondHadron->p3());
      //      kT=secondHadron->p3().perp(firstHadron->p3());
      //        cout <<"and now we set the type to : " << hadPType <<endl;
      computePIDWeights();

  }

  //get z according to new definition
  float getZ(ParticleInfo& pinf1, int massIndex1, ParticleInfo& pinf2, int massIndex2,HepLorentzVector& q)
  {
    float momProduct=pinf1.boostedLorentzVec[massIndex1].dot(pinf2.boostedLorentzVec[massIndex2]);
    float firstTerm=(momProduct-kinematics::masses[massIndex1]*kinematics::masses[massIndex1]*kinematics::masses[massIndex2]*kinematics::masses[massIndex2]/momProduct);
   float secondTerm=1.0/(pinf2.boostedLorentzVec[massIndex2].dot(q)-kinematics::masses[massIndex2]*kinematics::masses[massIndex2]*pinf1.boostedLorentzVec[massIndex1].dot(q)/(momProduct));
    float z=firstTerm*secondTerm;
    return z;
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

      cout <<"unknown single hadron "<< llund<<endl;
      return AnaDef::SH_TypeUnknown;
    }

  AnaDef::TwoHadCharge getTwoHadCharge(AnaDef::SingleHadCharge c1, AnaDef::SingleHadCharge c2)
    {
      if(c1==c2)
	return AnaDef::Likesign;
      else
	return AnaDef::Unlikesign;

//      if(AnaDef::SH_ChargeUnknown==c1 || AnaDef::SH_ChargeUnknown==c2)
//	return AnaDef::NA;
//      if((AnaDef::Pos==c1 && AnaDef::Neg==c2)|| (AnaDef::Neg==c1 && AnaDef::Pos==c2))
//	return AnaDef::PN;
//      if(AnaDef::Pos==c1 && AnaDef::Pos==c2)
//	return AnaDef::PP;
//      if(AnaDef::Neg==c1 && AnaDef::Neg==c2)
//	return AnaDef::NN;
//
//      return AnaDef::NA;
    }

  AnaDef::TwoHadPType getTwoHadType(AnaDef::SingleHadType c1, AnaDef::SingleHadType c2)
    {
      if(AnaDef::SH_TypeUnknown==c1 || AnaDef::SH_TypeUnknown==c2)
	{
	  	  cout <<" unknown two had " <<endl;
	return AnaDef::UNKNOWN;
	}

      if(AnaDef::Pion==c1 && AnaDef::Pion==c2)
	return AnaDef::PiPi;
      if(AnaDef::Kaon==c1 && AnaDef::Kaon==c2)
	return AnaDef::KK;
      if(AnaDef::Pion==c1&&AnaDef::Kaon==c2)
	return AnaDef::PiK;
      if(AnaDef::Kaon==c1&& AnaDef::Pion==c2)
	return AnaDef::KPi;
      if(AnaDef::Kaon==c1&& AnaDef::Proton==c2)
	return AnaDef::KP;
      if(AnaDef::Proton==c1&& AnaDef::Kaon==c2)
	return AnaDef::PK;
      if(AnaDef::Pion==c1&& AnaDef::Proton==c2)
	return AnaDef::PiP;
      if(AnaDef::Proton==c1&& AnaDef::Pion==c2)
	return AnaDef::PPi;
      if(AnaDef::Proton==c1&& AnaDef::Proton==c2)
	return AnaDef::PP;

      //      cout <<"no match, unknown" << c1 <<" c2: "<< c2 <<endl;
      return AnaDef::UNKNOWN;
    }

  void getZ_Kt(float &z1, float & z2, float& kT, float& qT)
  {
    switch(hadPType)
      {
      case AnaDef::PiPi:
	z1=z1_PiPi;
	z2=z2_PiPi;
	kT=kT_PiPi;
	qT=qT_PiPi;
	break;

      case AnaDef::PiK:
	z1=z1_PiK;
	z2=z2_PiK;
	kT=kT_PiK;
	qT=qT_PiK;
	break;

      case AnaDef::PiP:
	z1=z1_PiP;
	z2=z2_PiP;
	kT=kT_PiP;
	qT=qT_PiP;
	break;

      case AnaDef::KPi:
	z1=z1_KPi;
	z2=z2_KPi;
	kT=kT_KPi;
	qT=qT_KPi;
	break;
      case AnaDef::KK:
	z1=z1_KK;
	z2=z2_KK;
	kT=kT_KK;
	qT=qT_KK;
	break;

      case AnaDef::KP:
	z1=z1_KP;
	z2=z2_KP;
	kT=kT_KP;
	qT=qT_KP;
	break;


      case AnaDef::PPi:
	z1=z1_PPi;
	z2=z2_PPi;
	kT=kT_PPi;
	qT=qT_PPi;
	break;


      case AnaDef::PK:
	z1=z1_PK;
	z2=z2_PK;
	kT=kT_PK;
	qT=qT_PiK;
	break;

      case AnaDef::PP:
	z1=z1_PP;
	z2=z2_PP;
	kT=kT_PP;
	qT=qT_PP;
	break;

      default:
	//particle is a lepton, let's just use the pion numbers
	z1=z1_PiPi;
	z2=z2_PiPi;
	kT=kT_PiPi;
	qT=qT_PiPi;

	  break;
	
	//	cout <<"unkonwn type" <<hadPType <<endl;
	//	exit(1);
      }

  } 


};
#if defined(BELLE_NAMESPACE)
}
#endif

#endif
