#ifndef M_EVENT_H
#define  M_EVENT_H
#include <string>
#include <vector>
#include "TChain.h"
#include "ReaderBase.h"
//structure to hold event information

class MEvent:public ReaderBase
{
 public:

  float Thrust;
  float E_miss;
  float thrustThetaCMS;
  float thrustPhiCMS;
  float thrustThetaLab;
  float thrustPhi;
  float thetaEThrust;
  float transProj;
  float longProj;

  int runNr;
  int evtNr;
  int d0Tag;
  int D_Decay;
  int DStarDecay;
  //only for mcAsData
  int D_Decay_mc;
  int DStarDecay_mc;
  int DStarTag;
  int dStarTag;
  float isrPhotonEnergy;

  bool cutEvent;


  const float minThrust;
  const float maxMissingEnergy;
  const float lowerThrustThetaCut;
  const float upperThrustThetaCut;
  const float lowerThrustThetaCutCMS;
  const float upperThrustThetaCutCMS;

  const float thrustThetaCMSMaxProj;
  const bool d0Cut, dStarCut;
  const bool d0CutV2, dStarCutV2;
  const bool DCutMC, DStarCutMC;
  const float isrCut;

  //hard cuts 1.2-1.6
  //1 - 2.14 is 0.7-1.76 in cms
  //acos(0.5e) is 1.05 rad --> in lab 0.73 -.1.7
  //no cuts.. since e.g. cuts on thrust theta together with hadron theta binning biases kT...
  MEvent(TChain* chain, int mMCFlag=mcFlagNone):ReaderBase(mMCFlag),lowerThrustThetaCut(0.0),upperThrustThetaCut(3.5), thrustThetaCMSMaxProj(1.3), lowerThrustThetaCutCMS(0.0), upperThrustThetaCutCMS(1001.8), maxMissingEnergy(3.52), minThrust(0.5),d0Cut(false),d0CutV2(false),dStarCut(false),dStarCutV2(false),isrCut(0),DCutMC(false),DStarCutMC(true)
  {
    myChain=chain;
    if(chain)
      {
	branchPointers.push_back(&thrustThetaCMS);
	branchPointers.push_back(&Thrust);
	if(mMCFlag!=mcFlagWoA)
	  {
	    branchPointers.push_back(&E_miss);
	    branchPointers.push_back(&thrustPhiCMS);
	  }
	branchPointers.push_back(&thrustThetaLab);
	//    branchPointers.push_back(&thrustPhi);
	branchPointers.push_back(&thetaEThrust);
	branchPointersI.push_back(&runNr);
	branchPointersI.push_back(&evtNr);
	if(mMCFlag!=mcFlagWoA)
	  {
	    branchPointers.push_back(&isrPhotonEnergy);
	    branchPointersI.push_back(&d0Tag);
	    branchPointersI.push_back(&dStarTag);
	    branchPointersI.push_back(&D_Decay);
	    branchPointersI.push_back(&DStarDecay);
	    if(mMCFlag==mcAsData)
	      {
		branchPointersI.push_back(&D_Decay_mc);
		branchPointersI.push_back(&DStarDecay_mc);
	      }

	  }
	//the field thrustTheta_mc exists but encodes only the difference, so that leads to all cuts to fail
	//    branchNames.push_back("thrustTheta"+addendum);
	if(mMCFlag==mcFlagWoA)
	  {
	    branchNames.push_back("thrustTheta"+addendum);
	    branchNames.push_back("thrustMag"+addendum);
	  }
	else
	  {
	    branchNames.push_back("thrustTheta");
	    branchNames.push_back("Thrust"+addendum);
	    branchNames.push_back("E_miss"+addendum);
	    branchNames.push_back("thrustPhi");
	  }
	//    branchNames.push_back("thrustThetaLab"+addendum);
	if(mMCFlag==mcFlagWoA)
	  branchNames.push_back("thrustThetaLab"+addendum);////doesn't exist as mc
	else
	  branchNames.push_back("thrustThetaLab");////doesn't exist as mc
	//    branchNames.push_back("thrustPhi"+addendum);
	branchNames.push_back("thetaEThrust"+addendum);


	//same name for all...
	branchNamesI.push_back("runNr");
	branchNamesI.push_back("evtNr");
	if(mMCFlag!=mcFlagWoA)
	  {
	    branchNames.push_back("ISRPhotonEnergy");
	    branchNamesI.push_back("D0Tag");
	    branchNamesI.push_back("DStarTag");
	    branchNamesI.push_back("D_Decay"+addendum);
	    branchNamesI.push_back("DStar_Decay"+addendum);
	    //only makes sense for mcAsData. If it is just mc, the addendum will be _mc anyways
	    if(mMCFlag==mcAsData)
	      {
		branchNamesI.push_back("D_Decay_mc");
		branchNamesI.push_back("DStar_Decay_mc");
	      }
	  }

	doAllBranching();
      }
  }
  //decide if event is good or not
  void afterFill()
  {
    //    cout <<" thrustThetaLab? " << thrustThetaLab <<" cms : " << thrustThetaCMS << " proj: " << fabs(cos(thrustThetaCMS)) <<endl;
    cutEvent=false;
    //institute a cut against too much reconstructed energy...
    if(mMCFlag!=mcFlagWoA)
      {
	if(isrCut && (isrPhotonEnergy> isrCut))
	   cutEvent=true;
	if(d0Cut && (0==d0Tag ))
	  cutEvent=true;
	if(dStarCut && (0==dStarTag))
	  cutEvent=true;
	if(d0CutV2 && (D_Decay<0))
	  cutEvent=true;
	if(dStarCutV2 && (DStarDecay<0))
	  cutEvent=true;
	if(DStarCutMC && (D_Decay_mc<0 && DStarDecay_mc<0))
	  cutEvent=true;


	///check what the reconstruction does
	if(DStarDecay !=2 && D_Decay!=1)
	  cutEvent=true;



    if(E_miss<-1)
      cutEvent=true;
    if(E_miss>maxMissingEnergy)
      cutEvent=true;
      }
    if(Thrust<minThrust)
      cutEvent=true;
   
    if(thrustThetaLab>upperThrustThetaCut)
      cutEvent=true;
    if(thrustThetaLab<lowerThrustThetaCut)
      cutEvent=true;
    if(thrustThetaCMS>upperThrustThetaCutCMS)
      cutEvent=true;
    if(thrustThetaCMS<lowerThrustThetaCutCMS)
      cutEvent=true;

    if(fabs(cos(thrustThetaCMS))>thrustThetaCMSMaxProj)
      cutEvent=true;

    transProj=sin(thetaEThrust)*sin(thetaEThrust)/(1+cos(thetaEThrust)*cos(thetaEThrust));
    longProj=sqrt(1-transProj*transProj);
  }
}; 

#endif
