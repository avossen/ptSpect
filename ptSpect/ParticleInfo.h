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

  Hep3Vector boostedMoms[5];
  //identified as
  int idAs;
  int charge;

  float pidProbabilities[5];
  float p_Pi;
  float p_K;
  float p_p;
  float p_e;
  float p_mu;

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
