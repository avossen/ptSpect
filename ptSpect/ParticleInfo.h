#ifndef PARTICLEINFO__H
#define PARTICLEINFO__H
#include "particle/ParticleUserInfo.h"
#include "TVector3.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class ParticleInfo:public ParticleUserInfo
{
 public:
  //I don't think that these are actually used..
  double z[5];
  HepLorentzVector boostedLorentzVec[5];
  Hep3Vector boostedMoms[5];
    float labMom;
  //identified as
  int idAs;
  int charge;
  int isWeakDecay;

  float pidProbabilities[5];
  float pidProbabilities2[5];
  float pidUncert[5];
  float pidUncert2[5];
  float p_Pi;
  float p_K;
  float p_p;
  float p_e;
  float p_mu;

  float p_Pi2;
  float p_K2;
  float p_p2;
  float p_e2;
  float p_mu2;

  float ep_Pi;
  float ep_K;
  float ep_p;
  float ep_e;
  float ep_mu;

  float ep_Pi2;
  float ep_K2;
  float ep_p2;
  float ep_e2;
  float ep_mu2;




  double labTheta;
  double cmsTheta;
  double phi;

  double cmsPhi;
  double labPhi;
  double thrustProj;
  //for mc:
  int motherGenId;

 ParticleInfo(){};
 ParticleInfo(const ParticleInfo&x ){}
  ~ParticleInfo(){};
  ParticleInfo* clone() const{return new ParticleInfo(*this);}
  ParticleInfo &operator=(const ParticleInfo& x){}
};


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
