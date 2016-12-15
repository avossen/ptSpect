#ifndef ANADEFS_H
#define ANADEFS_H

#include <particle/Particle.h>
#include "CLHEP/Vector/LorentzVector.h"
#include <vector>
using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  //follow belle convention
#define pionIdx 2
#define kaonIdx 3
#define protonIdx 4
#define electronIdx 0
#define muonIdx 1


namespace AnaDef
{
  enum TwoHadCharge{Likesign,Unlikesign,NA};//the unknown is for mc, when we
  enum TwoHadPType{PiPi,PiK,PiP,KPi,KK,KP,PPi, PK, PP,UNKNOWN};
  enum ErrorCodes{
    smallRt
  };

  enum SingleHadCharge{Pos, Neg, Neut, SH_ChargeUnknown};
  enum SingleHadType{Pion, Kaon, Proton, Muon,Electron,MuonElectron,SH_TypeUnknown};

enum status
  {
    labSys,
    cmSys
  };

}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
