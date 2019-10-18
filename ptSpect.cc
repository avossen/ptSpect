//#define DEBUG_EVENT 287880//please no output
///--->in AnaDef.h now
//#define DEBUG_EVENT  32950
#define DEBUG_EVENT2 32950

#define pi0Mass 0.1349766
#define etaMass 0.548
//this will be used for the ISR corrections
//const bool onlyGen=false;

//const bool onlyGen=true;
const bool PRINT=false;


#include <iomanip>
#include "TMatrixD.h"
//#include "ptSpect/mc.h"  //one central place to put the define mc
#include "event/BelleEvent.h"
#include "particle/Particle.h"
#include "ptSpect/TwoHadAsymsCommons.h"
#include "particle/utility.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include <TROOT.h>
//#include "LorentzVector.h" //clhep
#include "belle.h"
#include <kid/atc_pid.h>
#include <eid/eid.h>
#include <mdst/Muid_mdst.h>
#include "math.h"
#include "TMath.h"
//for neutral particles:
#include "ip/IpProfile.h"
#include <mdst/findKs.h>
#include <mdst/findLambda.h>
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "benergy/BeamEnergy.h"

#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/khelix2xyz.h"
#include "kfitter/kfitterparticle.h"
//#include "exkvertexfitter... " for more than 1 to 2 decay

#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;


#include MDST_H
#include EVTCLS_H
#define PY_ELECTRON 11
#define PY_MU 13
#define PY_PI 211
#define PY_K 321
#define PY_Pi0 111
#define PY_KS0 310
#define PY_B0 511
#define PY_B 521

//#define SAVE_HISTOS
#define XCHECK

#define D0Mass 1.865
#define D0Width 0.6
#define D0Lower  D0Mass-D0Width
#define D0Upper  D0Mass+D0Width


//#define ThrustMethod
//#define noPID
//#define dataOnlyPID

#include <cmath>
//for thrust etc.
//strange: including this file in front of toolbox: thrust is not found anymore... (maybe "using belle namespace" a problem??)

//#include <mdst/Evtcls_hadron_info.h>
#if defined(BELLE_NAMESPACE)

namespace Belle {

#endif


  Hep3Vector& retSelf(Hep3Vector& vec)
  {
    return vec;
  };
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#include "ptSpect/AnaConsts.h"

#include "ptSpect/ptSpect.h"
#include "ptSpect/HadronPair.h"
#include "ptSpect/ParticleInfo.h"
#include "ptSpect/ParticleInfoMass.h"
#include "ptSpect/DebugHistos.h"
#include "particle/Ptype.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <time.h>
#include <fstream>

//#define W_NEUTRAL

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  using namespace std;
  // Constructor
  ptSpect::ptSpect():onlyGen_(-1.0),smpl_(12345.),cPiPlus("PI+"), cPiNeg("PI-"),cPiZero("PI0"),cKPlus("K+"),cKNeg("K-")
  {
    strcpy(rFileName,"notInitialized.root");
    test=0;
    for(int i=0;i<4;i++)
      {
	//	zVals[i]=0;
      }

    //    histoD0Spect=new TH1D("d0spect","d0spect",1000,0,3.0);
    //    histoDStar=new TH1D("dStarspect","dStarspect",300,1.8,5.0);
    //    histoPiSlowMom=new TH1D("piSlow","piSlow",100,0,3.0);
    //    histoRecDStarSpectToD0Pi=new TH1D("RecdStarSpectD0Pi","RecdStarSpectToD0Pi",300,1.8,2.2);
  }  
  ofstream* pXCheck;
  // Destructor
  ptSpect::~ptSpect(void)
  {

  }
  // initilization
  void ptSpect::init (int *status)
  {
#ifdef XCHECK
    pXCheck=new ofstream("xcheck");
#endif

    zBorders[0]=0.3;
    zBorders[1]=0.5;
    zBorders[2]=0.7;
    zBorders[3]=1.5;

    ptBorders[0]=0.15;
    ptBorders[1]=0.3;
    ptBorders[2]=0.5;
    ptBorders[3]=3.0;

    for(int i=0;i<30;i++)
      {
	numEtas[i]=0;
	numPi0s[i]=0;
      }

    gROOT->SetStyle("Plain");

    thetaPhiLab=new TH2D("thetaPhiLab","thetaPhiLab",100,0,6.3,100,-3.15,3.15);
    thetaPhiCMS=new TH2D("thetaPhiCMS","thetaPhiCMS",100,0,6.3,100,-3.15,3.15);


    const double eler(3.499218);//energies of l, h beam
    const double eher(7.998213);
    const double theta(0.022);
    validRun=true;
    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);
    HepLorentzVector CMBoost2=kinematics::cm;
    CMBoost2.boost(-kinematics::CMBoost);  ///?????->because the sign of the electron vectors is reverted in order to construct a boost vector with a positive sign...
    kinematics::Q=CMBoost2.t();
    //  kinematics::Q=10.52; //die e energien die ich da hab sind on resonance. Dass hier ist aber continuum
    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    kinematics::secondElectronCM.boost(kinematics::CMBoost);
    //    m_histos.setFilenameStart(rFileName);

    srand(time(NULL));

    //declare and read in pid matrices
    pidMatrixPositive=new float***[numMomBins];
    pidMatrixNegative=new float***[numMomBins];
    pidMatrixPositive2=new float***[numMomBins];
    pidMatrixNegative2=new float***[numMomBins];
    pidUncertPositive=new float***[numMomBins];
    pidUncertNegative=new float***[numMomBins];
    pidUncertPositive2=new float***[numMomBins];
    pidUncertNegative2=new float***[numMomBins];


    masses[pionIdx]=0.14;

    masses[protonIdx]=0.938;
    masses[electronIdx]=0.03;
    masses[muonIdx]=0.11;
    masses[kaonIdx]=0.493;

    kinematics::masses[pionIdx]=masses[pionIdx];
    kinematics::masses[protonIdx]=masses[protonIdx];
    kinematics::masses[electronIdx]=masses[electronIdx];
    kinematics::masses[muonIdx]=masses[muonIdx];
    kinematics::masses[kaonIdx]=masses[kaonIdx];

    for(int i=0;i<numMomBins;i++)
      {
	pidMatrixPositive[i]=new float**[numThetaBins];
	pidMatrixNegative[i]=new float**[numThetaBins];
	pidMatrixPositive2[i]=new float**[numThetaBins];
	pidMatrixNegative2[i]=new float**[numThetaBins];
	pidUncertPositive[i]=new float**[numThetaBins];
	pidUncertNegative[i]=new float**[numThetaBins];
	pidUncertPositive2[i]=new float**[numThetaBins];
	pidUncertNegative2[i]=new float**[numThetaBins];


	for(int j=0;j<numThetaBins;j++)
	  {
	    pidMatrixPositive[i][j]=new float*[numPIDs];
	    pidMatrixNegative[i][j]=new float*[numPIDs];
	    pidMatrixPositive2[i][j]=new float*[numPIDs];
	    pidMatrixNegative2[i][j]=new float*[numPIDs];
	    pidUncertPositive[i][j]=new float*[numPIDs];
	    pidUncertNegative[i][j]=new float*[numPIDs];
	    pidUncertPositive2[i][j]=new float*[numPIDs];
	    pidUncertNegative2[i][j]=new float*[numPIDs];

	    for(int k=0;k<numPIDs;k++)
	      {
		pidMatrixPositive[i][j][k]=new float[numPIDs];
		pidMatrixNegative[i][j][k]=new float[numPIDs];
		pidMatrixPositive2[i][j][k]=new float[numPIDs];
		pidMatrixNegative2[i][j][k]=new float[numPIDs];
		pidUncertPositive[i][j][k]=new float[numPIDs];
		pidUncertNegative[i][j][k]=new float[numPIDs];
		pidUncertPositive2[i][j][k]=new float[numPIDs];
		pidUncertNegative2[i][j][k]=new float[numPIDs];


		for(int l=0;l<numPIDs;l++)
		  {
		    //this will be overwritten by Martin's matrices. So set this to zero so we realize that we loaded incorrect values..
		    if(k==l)
		      {
			pidMatrixPositive[i][j][k][l]=0.0;
			pidMatrixNegative[i][j][k][l]=0.0;
			pidMatrixPositive2[i][j][k][l]=0.0;
			pidMatrixNegative2[i][j][k][l]=0.0;

		      }
		  }
	      }
	  }
      }

    loadPIDMatrix();
    //do this after the PID matrix load, so we don't have to deal with changing root directories
    m_file=new TFile(rFileName,"recreate");
  }



  int ptSpect::goodHadronB( ) const 
  {
    // initialize return value
    int b = 0;

    // get manager for evtcls_hadronic_flag
    Evtcls_hadronic_flag_Manager & evtMgr
      = Evtcls_hadronic_flag_Manager::get_manager();

    //   Evtcls_evtcls_flag_Manager & evtMgr2
    //     = Evtcls_evtcls_flag2_Manager::get_manager();

    // get flag for HadronB
    //      hadronic_flag(0) =  10 : old HadronA with R2<0.2
    //                       =  20 : old HadronA with R2>=0.2
    //      hadronic_flag(1) =  10 : new HadronA with R2<0.2
    //                       =  20 : new HadronA with R2>=0.2
    //      hadronic_flag(2) =  10 : HadronB with R2<0.2
    //                       =  20 : HadronB with R2>=0.2
    //      hadronic_flag(3) =  10 : new HadronA with #tracks>=5
    //      hadronic_flag(4) =  10 : HadronB with #tracks>=5
    //      hadronic_flag(5) =  10 : HadronC with R2<0.2
    //                       =  20 : HadronC with R2>=0.2
    Panther_ID id(1);
    Evtcls_hadronic_flag & hadflag = evtMgr(id);

    if ( hadflag.hadronic_flag(2) == 10 ||
	 hadflag.hadronic_flag(2) == 20    ) { b = 1; }

    //hadron J
    //      if ( hadflag.evtcls_flag2(2) >= 0  ) { b = 1; }

    if (b !=1)
      {
	//       printf("bad hadb: %d %d %d %d %d\n",hadflag.hadronic_flag(0),hadflag.hadronic_flag(2),hadflag.hadronic_flag(3),hadflag.hadronic_flag(4),hadflag.hadronic_flag(5));

      }

    // to select tau event
    //  if ( hadflag.hadronic_flag(4) != 0 ) { b += 2; }

    return b;
  }
  // begin_run function
  void ptSpect::begin_run(BelleEvent* evptr, int* status)
  {
      cout <<"debug event : " << DEBUG_EVENT <<endl;
    //    cout <<" in begin run... " <<endl;
    IpProfile::begin_run();
    eid::init_data();

    BeamEnergy::begin_run();
    double eler=BeamEnergy::E_LER();
    double eher=BeamEnergy::E_HER();
    cout <<"got eler: " << eler <<" eher: " << eher<<endl;
    if(eler <3.0 || eher <7.0 || eler > 5.0 || eher > 9.0)
      {
	validRun=false;
	return;
      }
    else
      {
	validRun=true;
      }
    double theta(0.022);
    kinematics::cm=HepLorentzVector(-eher*sin(theta),0.,-eher*cos(theta)+eler,eher+eler);
    kinematics::CMBoost=kinematics::cm.boostVector();
    kinematics::firstElectronCM=HepLorentzVector(eher*sin(theta),0.,eher*cos(theta),eher);
    kinematics::secondElectronCM=HepLorentzVector(0.,0.,-eler,eler);
    HepLorentzVector CMBoost2=kinematics::cm;
    CMBoost2.boost(-kinematics::CMBoost);  ///?????->because the sign of the electron vectors is reverted in order to construct a boost vector with a positive sign...
    kinematics::Q=CMBoost2.t();
    cout <<" Q is : "<<kinematics::Q<< endl;
    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    kinematics::secondElectronCM.boost(kinematics::CMBoost);
    return;
  }

  // hist_def function
  void ptSpect::hist_def()
  {
    Particle p;



    if(rFileName!=0)
      cout <<endl<<":::----- rFileName (handedness): " << rFileName <<endl<<endl;
    else
      cout <<endl <<":::--------No File Name specified (handedness)" <<endl<<endl;

    pTreeSaver=new TreeSaver();
    //    pTreeSaver->setDebugHistos(&m_histos);
    pTreeSaver->addArrayF("z1");
    pTreeSaver->addArrayF("z2");


    //storing entire objects would certainly be better, but we already went down this path...
    //these are the probabilities (weights) for each possible pair
    //these are also too many combinations, could do the same with 2x3
    pTreeSaver->addArrayF("p_PiPi");
    pTreeSaver->addArrayF("p_PiK");
    pTreeSaver->addArrayF("p_PiP");
    pTreeSaver->addArrayF("p_KPi");
    pTreeSaver->addArrayF("p_KK");
    pTreeSaver->addArrayF("p_KP");
    pTreeSaver->addArrayF("p_PPi");
    pTreeSaver->addArrayF("p_PK");
    pTreeSaver->addArrayF("p_PP");


    ///--->the two matrix set from martin
    pTreeSaver->addArrayF("p_PiPi_1");
    pTreeSaver->addArrayF("p_PiK_1");
    pTreeSaver->addArrayF("p_PiP_1");
    pTreeSaver->addArrayF("p_KPi_1");
    pTreeSaver->addArrayF("p_KK_1");
    pTreeSaver->addArrayF("p_KP_1");
    pTreeSaver->addArrayF("p_PPi_1");
    pTreeSaver->addArrayF("p_PK_1");
    pTreeSaver->addArrayF("p_PP_1");

    pTreeSaver->addArrayF("p_PiPi_2");
    pTreeSaver->addArrayF("p_PiK_2");
    pTreeSaver->addArrayF("p_PiP_2");
    pTreeSaver->addArrayF("p_KPi_2");
    pTreeSaver->addArrayF("p_KK_2");
    pTreeSaver->addArrayF("p_KP_2");
    pTreeSaver->addArrayF("p_PPi_2");
    pTreeSaver->addArrayF("p_PK_2");
    pTreeSaver->addArrayF("p_PP_2");

    //uncertainties on the p values
    pTreeSaver->addArrayF("ep_PiPi");
    pTreeSaver->addArrayF("ep_PiK");
    pTreeSaver->addArrayF("ep_PiP");
    pTreeSaver->addArrayF("ep_KPi");
    pTreeSaver->addArrayF("ep_KK");
    pTreeSaver->addArrayF("ep_KP");
    pTreeSaver->addArrayF("ep_PPi");
    pTreeSaver->addArrayF("ep_PK");
    pTreeSaver->addArrayF("ep_PP");

    pTreeSaver->addArrayF("kT_PiPi");
    pTreeSaver->addArrayF("kT_PiK");
    pTreeSaver->addArrayF("kT_PiP");
    pTreeSaver->addArrayF("kT_KPi");
    pTreeSaver->addArrayF("kT_KK");
    pTreeSaver->addArrayF("kT_KP");
    pTreeSaver->addArrayF("kT_PPi");
    pTreeSaver->addArrayF("kT_PK");
    pTreeSaver->addArrayF("kT_PP");

    //dot products
    pTreeSaver->addArrayF("dp_PiPi");
    pTreeSaver->addArrayF("dp_PiK");
    pTreeSaver->addArrayF("dp_PiP");
    pTreeSaver->addArrayF("dp_KPi");
    pTreeSaver->addArrayF("dp_KK");
    pTreeSaver->addArrayF("dp_KP");
    pTreeSaver->addArrayF("dp_PPi");
    pTreeSaver->addArrayF("dp_PK");
    pTreeSaver->addArrayF("dp_PP");


    pTreeSaver->addArrayF("qT_PiPi");
    pTreeSaver->addArrayF("qT_PiK");
    pTreeSaver->addArrayF("qT_PiP");

    pTreeSaver->addArrayF("qT_KPi");
    pTreeSaver->addArrayF("qT_KK");
    pTreeSaver->addArrayF("qT_KP");


    pTreeSaver->addArrayF("qT_PPi");
    pTreeSaver->addArrayF("qT_PK");
    pTreeSaver->addArrayF("qT_PP");



    pTreeSaver->addArrayF("z1_PiPi");
    pTreeSaver->addArrayF("z1_PiK");
    pTreeSaver->addArrayF("z1_PiP");
    pTreeSaver->addArrayF("z1_KPi");
    pTreeSaver->addArrayF("z1_KK");
    pTreeSaver->addArrayF("z1_KP");
    pTreeSaver->addArrayF("z1_PPi");
    pTreeSaver->addArrayF("z1_PK");
    pTreeSaver->addArrayF("z1_PP");

    pTreeSaver->addArrayF("z2_PiPi");
    pTreeSaver->addArrayF("z2_PiK");
    pTreeSaver->addArrayF("z2_PiP");
    pTreeSaver->addArrayF("z2_KPi");
    pTreeSaver->addArrayF("z2_KK");
    pTreeSaver->addArrayF("z2_KP");
    pTreeSaver->addArrayF("z2_PPi");
    pTreeSaver->addArrayF("z2_PK");
    pTreeSaver->addArrayF("z2_PP");



    //theta of the particles...  for mc it is the difference...
    pTreeSaver->addArrayF("labTheta1");
    pTreeSaver->addArrayF("labTheta2");

    pTreeSaver->addArrayF("labPhi1");
    pTreeSaver->addArrayF("labPhi2");

    pTreeSaver->addArrayF("cmsTheta1");
    pTreeSaver->addArrayF("cmsTheta2");

    pTreeSaver->addArrayF("cmsPhi1");
    pTreeSaver->addArrayF("cmsPhi2");


    pTreeSaver->addArrayF("thrustProj1");
    pTreeSaver->addArrayF("thrustProj2");
    pTreeSaver->addArrayF("kT");
    pTreeSaver->addArrayF("dp");
    pTreeSaver->addArrayF("HadDiffTheta");
    pTreeSaver->addArrayF("HadDiffPhi");
    pTreeSaver->addArrayF("qT");



    //important that first all arrays are defined, F, than I 
#ifdef MC
    pTreeSaver->addArrayF("z1_mc");
    pTreeSaver->addArrayF("z2_mc");

    pTreeSaver->addArrayF("labTheta1_mc");
    pTreeSaver->addArrayF("labTheta2_mc");
    pTreeSaver->addArrayF("labPhi1_mc");
    pTreeSaver->addArrayF("labPhi2_mc");

    pTreeSaver->addArrayF("cmsTheta1_mc");
    pTreeSaver->addArrayF("cmsTheta2_mc");
    pTreeSaver->addArrayF("cmsPhi1_mc");
    pTreeSaver->addArrayF("cmsPhi2_mc");


    pTreeSaver->addArrayF("thrustProj1_mc");
    pTreeSaver->addArrayF("thrustProj2_mc");
    pTreeSaver->addArrayF("kT_mc");
    //dot product
    pTreeSaver->addArrayF("dp_mc");
    pTreeSaver->addArrayF("HadDiffTheta_mc");
    pTreeSaver->addArrayF("HadDiffPhi_mc");
    pTreeSaver->addArrayF("qT_mc");

#endif
    //!!!! this charge type is in principle worthless, look at the charges of the hadrons!!-->change this to likesign and unlikesign
    pTreeSaver->addArrayI("chargeType");
    pTreeSaver->addArrayI("particleType");
    pTreeSaver->addArrayI("chargeType1");
    pTreeSaver->addArrayI("particleType1");
    pTreeSaver->addArrayI("chargeType2");
    pTreeSaver->addArrayI("particleType2");
#ifdef MC
    pTreeSaver->addArrayI("chargeType_mc");
    pTreeSaver->addArrayI("particleType_mc");
    pTreeSaver->addArrayI("chargeType1_mc");
    pTreeSaver->addArrayI("particleType1_mc");
    pTreeSaver->addArrayI("chargeType2_mc");
    pTreeSaver->addArrayI("particleType2_mc");

#endif
    pTreeSaver->addFieldF("Thrust");
    pTreeSaver->addFieldF("E_miss");
    pTreeSaver->addFieldF("thrustTheta"); //theta of thrust
    pTreeSaver->addFieldF("thrustPhi"); //phi of thrust
    pTreeSaver->addFieldF("thrustThetaLab");
    pTreeSaver->addFieldF("thetaEThrust"); //angle between e+-e- axis and thrust -> for correction factor


    ///add


#ifdef MC
    pTreeSaver->addFieldF("Thrust_mc");
    pTreeSaver->addFieldF("E_miss_mc");
    pTreeSaver->addFieldF("thrustTheta_mc"); //angle between thrust rec thrust mc
    //  pTreeSaver->addFieldF("thrustPhi_mc");  //<--- not needed because we have phi of orginal and diff phi...
    pTreeSaver->addFieldF("diffThetaThrust_mc");
    pTreeSaver->addFieldF("diffPhiThrust_mc");

    pTreeSaver->addFieldF("VP_Energy"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PX"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PY"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PZ"); //energy of virtual photon
    pTreeSaver->addFieldF("quarkAngle");
    pTreeSaver->addFieldF("thetaEThrust_mc");
    pTreeSaver->addFieldF("ISRPhotonEnergy");
    pTreeSaver->addFieldI("numQuarks");
    pTreeSaver->addFieldI("D_Decay_mc");
    pTreeSaver->addFieldI("DStar_Decay_mc");


    //the fields for the mc w/o acceptance
#endif
    pTreeSaver->addFieldI("runNr");
    pTreeSaver->addFieldI("evtNr");


    pTreeSaver->addFieldI("D0Tag");
    pTreeSaver->addFieldI("DStarTag");

    pTreeSaver->addFieldI("D_Decay");
    pTreeSaver->addFieldI("DStar_Decay");


    //doesn't work anymore with arrays...
    //  pTreeSaver->createNTuple(tm);
#ifdef MC

    //  pTreeSaver->addArrayPi0AsymmetryF("realPi0_gammaAsymmetry");
#endif
  }

  // event function
  void ptSpect::event(BelleEvent* evptr, int* status)
  {
    //        cout <<"debug event : " << DEBUG_EVENT <<endl;
    bool onlyGen=false;
    if(onlyGen_>0.0)
      {
	onlyGen=true;
      }


    //    cout <<"value of onlygen: "<< onlyGen_ << ", onlyGen: " << onlyGen <<endl;

    bool eventCut=false;
    int evtNr;
    int runNr;
    /////for xcheck

    if(Belle_event_Manager::get_manager().begin()==Belle_event_Manager::get_manager().end())
      return;

    evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    runNr=Belle_event_Manager::get_manager().begin()->RunNo();

    kinematics::runNr=runNr;
    kinematics::evtNr=evtNr;

    //      cout <<"in event " <<endl;
    if(onlyGen)
      {
	eventCut=true;
	//	cout <<"save gen info" <<endl;
	pTreeSaver->saveGenInfo(v_allParticles, eventCut);
	if(eventCut)
	  {
	    exitEvent();
	    return;
	  }
      }
    const double m_pi0=0.1349766;
    vector<float> v_drH1;
    vector<float> v_drH2;
    if(!validRun)
      {
	eventCut=true;
	//	return;
      }
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"looking at debug event " <<endl;
      }
    //    (*pXCheck) << endl<<  "processing event nr " << evtNr <<endl;
    //        cout <<endl<<"--> run = " << runNr <<" evtNr = "  <<evtNr <<endl;
    //    cout <<" --> exp = " << " run = " << runNr << " event = " << evtNr <<endl;
    if(!IpProfile::usable())
      {
	//	cout <<"ip profile not usable" <<endl;
	eventCut=true;
	//	return;
      }
    if(!(test%1000))
      {
	//	      cout << "evt " <<test <<endl;
	//	        cout << "nr " <<evtNr <<endl;
      }
    test++;

    //#ifndef MC
    if(!goodHadronB())
      {
	//	cout <<"no good hadronb " <<endl;
	eventCut=true;

	//      return;
      }
    //#endif
    float visEnergy=0;
    float visEnergyOnFile=0;
    vector<float> v_pLab;
    vector<float> v_pCMS;
    char ptypName[200];
    int lundPC=-1;
    vector<float> v_z;
    vector<float> v_phi;
    vector<float> v_theta;
    vector<float> v_qt;

    vector<float> v_pxCMS;
    vector<float> v_pyCMS;
    vector<float> v_pzCMS;

    vector<float> v_pxLab;
    vector<float> v_pyLab;
    vector<float> v_pzLab;

    vector<float> v_pxG;
    vector<float> v_pyG;
    vector<float> v_pzG;
    vector<float> v_pGLab;
    vector<Hep3Vector> allParticlesBoosted;  
    vector<Hep3Vector> boostedPi0s;
    vector<Hep3Vector> boostedEtas;
    vector<int> allPB_particleClass;
    vector<int> allPB_particleCharge;

    vector<float> nonBoostedE;
    vector<Hep3Vector> allParticlesNonBoosted;

    vector<float> allPB_E; //energy for the three vec above...
    Mdst_charged_Manager& mdst_chr_Mgr=Mdst_charged_Manager::get_manager();
    //  Mdst_klong_Manager& mdst_klong_Mgr=Mdst_klong_Manager::get_manager();
    Mdst_trk_Manager& mdst_trk_Mgr=Mdst_trk_Manager::get_manager();
    atc_pid selKPi(3,1,5,3,2);  //K/pi separation
    atc_pid selPK(3,1,5,4,3); //proton kaon separation
    atc_pid selPiP(3,1,5,2,4); //pion proton separation
    atc_pid selKP(3,1,5,3,4);

    int itmp=0;
    int iChTrks=0;//num of charged Tracks
    int iChTrkPos=0;
    int iChTrkNeg=0;
    //   cout <<"chargedTracks: " << mdst_chr_Mgr.size() <<" num Trk: " << mdst_trk_Mgr.size() <<" klong: " << mdst_klong_Mgr.size() <<endl;

    //    cout <<"there are " << mdst_chr_Mgr.size() << " charge tracks in mdst_chr " <<endl;
    //
    //saves generator level pairs in MC (doesn't do anything for no MC)



    for(Mdst_charged_Manager::iterator chr_it=mdst_chr_Mgr.begin();chr_it!=mdst_chr_Mgr.end();chr_it++)
      {
	Hep3Vector tmpv3(chr_it->p(0),chr_it->p(1),chr_it->p(2));
	if(DEBUG_EVENT==evtNr)
	  {
	    cout <<"very first looking at track with theta: "<< tmpv3.theta() <<" lab p : "<< tmpv3.mag()<<endl;
	  }

	if(evtNr==DEBUG_EVENT)
	  {
	    atc_pid selPID2(3,1,5,2,3);
	    double pid_K = selPID2.prob(*chr_it);
	    printf("pid_K %f \n",pid_K);
	  }
	if(!enoughSVDHits(chr_it))
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		cout <<" cut track with naive momentum: "<< (*chr_it).p(0)<<", " << (*chr_it).p(1) <<", " << (*chr_it).p(2)<<endl;
	      }
	    continue;
	  }
	bool positivId=false;
	double m_mass=m_pi;
	int massHyp=pionIdx;
	if(evtNr==DEBUG_EVENT)
	  {
	    cout <<"setting massHyp to " << massHyp <<endl;
	  }


	double m_theta=0;
	double m_phi=0;
	double m_qt=0;
	double m_z[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double charge=(*chr_it).charge();
	///we take all (and the name doesn't really matter, so make this pion default)
	//	strcpy(ptypName,"unknown");
	if(charge>0)
	  strcpy(ptypName,"PI+");
	else
	  strcpy(ptypName,"PI-");
	lundPC=-1;


	//      HepLorentzVector hepvec;
	//immer daran denken in die richtung -boostvector zu boosten! ... nein ist anscheinend schon in der boost vector definition drin...
	eid sel_e(*chr_it);
	double mu_id=0;
	Muid_mdst muID(*chr_it);
	if(muID.Chi_2()>0)
	  mu_id=muID.Muon_likelihood();
	double atcKPi=selKPi.prob(*chr_it);
	double atcKP=selKP.prob(*chr_it);
	double atcPiP=selPiP.prob(*chr_it);
	float e_cut=0.85;
	float mu_cut=0.9;
	//	double e_id=sel_e.prob(0,-1,0);
	double e_id=sel_e.prob(3,-1,5);
	//	cout <<"atcKP " << atcKP <<" atcPiP : "<< atcPiP <<endl;
	if(DEBUG_EVENT==evtNr)
	  {
	    cout <<"pids "<< e_id <<" mu_id: "<< mu_id <<" atcKPi: "<< atcKPi << " atcKP: "<< atcKP <<" atcPiP: "<< atcPiP <<" px: "<< (*chr_it).p(0)<<endl;
	Hep3Vector tmpv2(chr_it->p(0),chr_it->p(1),chr_it->p(2));
	if(DEBUG_EVENT==evtNr)
	  {
	    cout <<"first looking at track with theta: "<< tmpv2.theta() <<" lab p : "<< tmpv2.mag()<<endl;
	  }

	    //cout <<"pid kpi: " << atcKPiAlt <<" pid KP: " << atcKPAlt << " e_id: " << e_id << " mu_id: " << mu_id <<endl;
	  }
	bool isLepton=false;
	bool isPionKaon=false;
	bool isProton=false;
	//	cout <<"e id: "<< e_id <<" mu id: "<< mu_id << endl;
	if(e_id>e_cut&& mu_id<0.9 && e_id<=1.0 &&mu_id >=-999)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		cout <<"is electron "<< e_id <<" mu_id: "<< mu_id <<" atcKPi: "<< atcKPi << " atcKP: "<< atcKP <<" atcPiP: "<< atcPiP <<endl;
	      }
	    m_mass=m_e;
	    positivId=true;
	    massHyp=electronIdx;
	    isLepton=true;
	    //	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	    //	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	    if(charge>0)
	      {
		strcpy(ptypName,"E+");
		lundPC=-11;
	      }
	    else
	      {
		strcpy(ptypName,"E-");
		lundPC=11;
	      }

	  }
	//used to be mu_id > e_cut 
	if(mu_id>mu_cut && e_id<e_cut && mu_id<=1.0&&e_id>=0.0)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
	       	cout <<"is muon" <<endl;
	      }
	    //	    m_histos.hPidMuPi->Fill(chr_it->trk().pid_mu(),chr_it->trk().pid_pi());
	    m_mass=m_muon;
	    positivId=true;
	    massHyp=muonIdx;
	    isLepton=true;
	    if(charge>0)
	      {
		strcpy(ptypName,"MU+");
		lundPC=13;
	      }
	    else
	      {
		strcpy(ptypName,"MU-");
		lundPC=-13;
	      }

	    //	    m_histos.hPidMu->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }
	if(mu_id>0.9&& e_id>e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		cout <<"is electron" <<endl;
	      }
	    //	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	    m_mass=m_e;
	    positivId=true;
	    //massHyp is already set to zero above 
	    isLepton=true;
	    //	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	    if(charge>0)
	      {
		strcpy(ptypName,"E+");
		lundPC=-11;
	      }
	    else
	      {
		strcpy(ptypName,"E-");
		lundPC=11;
	      }

	  }


	
	//	if(!isLepton)// && e_id>=0.0 && e_id<0.85 && mu_id>=-999 && mu_id<0.9 )
	if(!isLepton && e_id>=0.0 && e_id<0.85 && mu_id>=-999 && mu_id<0.9 )
	  {
	    if(atcKPi>0.6 && atcKPi <= 1.0 && atcKP >0.2 && atcKP <= 1.0 && atcPiP<=1.0 && atcPiP>=0.0) //kaon
	      {
		m_mass=m_k;
		massHyp=kaonIdx;
		isPionKaon=true;
		positivId=true;
		if(charge>0)
		  {
		    strcpy(ptypName,"K+");
		    lundPC=321;
		  }
		else
		  {
		    strcpy(ptypName,"K-");
		    lundPC=-321;
		  }
		if(DEBUG_EVENT==evtNr)
		  {
		    cout <<"is kaon "<< e_id <<" mu_id: "<< mu_id <<" atcKPi: "<< atcKPi << " atcKP: "<< atcKP <<" atcPiP: "<< atcPiP <<endl;
		    //cout <<"pid kpi: " << atcKPiAlt <<" pid KP: " << atcKPAlt << " e_id: " << e_id << " mu_id: " << mu_id <<endl;
		  }

		//		  m_histos.hPidK->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
	      }
	    else
	      {
		if(atcKP<0.2 && atcKP <=1.0 && atcPiP<0.2 && atcPiP<=1.0)
		  {
		    m_mass=m_pr;
		    massHyp=protonIdx;
		    //			m_histos.hPidPr->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_p());
		    //			m_histos.hPidPrPi->Fill(chr_it->trk().pid_p(),chr_it->trk().pid_pi());

		    isProton=true;
		    positivId=true;
		    if(charge>0)
		      {
			strcpy(ptypName,"P+");
			lundPC=2212;
		      }
		    else
		      {
			//from decay dec, still negative charge...
			strcpy(ptypName,"AP+");
			lundPC=-2212;
		      }
		  }
		else  //pion
		  {
		    //default mass assignment if nothing is found is pion anywasy....
		    //			if(atcKPi<0.3)
		    if(atcKPi<0.6 && atcPiP>=0.2 && atcKPi>=0.0 && atcPiP<=1.0)
		      {
			massHyp=pionIdx;
			m_mass=m_pi;
			positivId=true;
			if(evtNr==DEBUG_EVENT)
			  {
			    cout <<"setting massHyp positively to " << massHyp <<endl;
			  }



			isPionKaon=true;
			if(charge>0)
			  {
			    strcpy(ptypName,"PI+");
			    lundPC=211;
			  }
			else
			  {
			    strcpy(ptypName,"PI-");
			    lundPC=-211;
			  }
			/*		      m_histos.hPidPi->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
					      m_histos.hPidPiMu->Filpl(chr_it->trk().pid_pi(),chr_it->trk().pid_mu());
					      m_histos.hPidPiE->Fill(chr_it->trk().pid_pi(),chr_it->trk().pid_e());
					      m_histos.hPidPiPr->Fill(chr_it->trk().pid_pi(),chr_it->trk().pid_p());*/

		      }
		  }
	      }
	  }
	//new...
	if(!positivId)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		cout <<"no pos id " <<endl;
	      }
	    continue;
	  }
	Hep3Vector tmpv(chr_it->p(0),chr_it->p(1),chr_it->p(2));
	if(DEBUG_EVENT==evtNr)
	  {
	    cout <<"looking at track with theta: "<< tmpv.theta() <<" lab p : "<< tmpv.mag()<<" px: "<< chr_it->p(0) <<" ptyp: "<< ptypName<< endl;
	  }
	//	(*pXCheck) <<"looking at charged " << ptypName <<"  with lab theta: "<< tmpv.theta() << " costheta: " << tmpv.cosTheta() << " and mom: "<< tmpv.mag()<<endl;

	//	cout <<"massHyp: "<< massHyp <<endl;
	double dr, dz, refitPx, refitPy, refitPz;
	getDrDz(chr_it, massHyp,dr,dz, refitPx, refitPy, refitPz);
	v_vertexR.push_back(dr);
	v_vertexZ.push_back(dz);
	//	Particle* tmpP=new Particle(*chr_it,string(ptypName));
	//	(*pXCheck).precision(3);
	//	(*pXCheck) << std::setw(7);
	//	(*pXCheck)<< "lab refit " << string(ptypName) << " massHyp: "<< massHyp <<", x"  << refitPx << " y: " << refitPy<< " z: " << refitPz << endl;
	//	(*pXCheck)<< "lab non-refit rho "<<tmpP->p().rho()<< "x :"  << tmpP->p().x() << " y: " <<  tmpP->p().y() << " z: " << tmpP->p().z() << " dr: "<< dr <<" dz: "<< dz << endl;
	//	(*pXCheck) <<" mass: "<< masses[massHyp] <<endl;

	Hep3Vector h3Vect(refitPx,refitPy,refitPz);
	////
	float labMom=h3Vect.mag();
	//	cout <<"track with " << labMom << " momentum " <<endl;
	///
	//      cout <<"looking at " <<(*chr_it).p(0) <<" " << (*chr_it).p(1) <<" " << (*chr_it).p(2) <<endl;
	if ( fabs(dr) > cuts::vertexR )//slides from kibayashi 
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		cout <<"dr cut: " << fabs(dr) <<endl;
		cout <<" cut track due to vertex.. r: " << fabs(dr) <<" px lab of track:  " <<(*chr_it).p(0)<< endl;
	      }

	    continue;
	  }
	if ( fabs(dz) > cuts::vertexZ ) 
	  {

	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		cout <<"dz cut: " << fabs(dz) <<endl;
		cout <<" cut track due to vertex.. z" <<endl;
	      }
	   
	    continue;//used to be 4
	  }
	if(charge>0)
	  iChTrkPos++;
	else
	  iChTrkNeg++;
	//	cout <<"compare " << (*chr_it).p(0) << ", " << (*chr_it).p(1) <<", " << (*chr_it).p(2) <<endl;
	//		cout <<" to: "<< refitPx << " " << refitPy <<" " << refitPz <<endl;
	
	//		Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
	double E[5];
	for(int i=0;i<5;i++)
	  {
	    E[i]=sqrt(masses[i]*masses[i]+h3Vect.mag2());
	  }
	HepLorentzVector nonBoostedVec[5];
	HepLorentzVector boostedVec[5];
	for(int i=0;i<5;i++)
	  {
	    nonBoostedVec[i]=HepLorentzVector(h3Vect,E[i]);
	    boostedVec[i]=HepLorentzVector(h3Vect,E[i]);
	    if(i==massHyp)
	      {
		//(*pXCheck) <<" before boost px: "<<    boostedVec[i].px() <<", y: "<< boostedVec[i].py() <<"pz: "<< boostedVec[i].pz() <<" e: "<< boostedVec[i].e() <<endl;
	      }
	    boostedVec[i].boost(kinematics::CMBoost);
	    if(i==massHyp)
	      {
		//		(*pXCheck) <<" after boost px: "<<    boostedVec[i].px() <<", y: "<< boostedVec[i].py() <<"pz: "<< boostedVec[i].pz() <<" e: "<< boostedVec[i].e() <<" naive z: "<< boostedVec[i].vect().mag()/5.25 <<endl;
	      }
	  }



	if(h3Vect.perp()<cuts::minPtThrust)
	  {

	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		cout <<"removing pt=: " << h3Vect.perp() <<endl;
	      }
	    //	      cout <<"removing pt=: " << h3Vect.perp() <<endl;
	    //	    (*pXCheck) <<" cut due to min pt " << endl;
	    continue;
	  }
	for(int i=0;i<5;i++)
	  {
	    m_z[i]=2*boostedVec[i].e()/kinematics::Q;
	  }
	iChTrks++;
	if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	  {
	    cout <<"charged z: " << m_z<<endl;
	  }
	//just use min p cut
	if(m_z[massHyp]<cuts::minZThrust)
	  {
	    //	    continue;
	  }
	if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	  {
	    //	    cout <<"adding charged track: " << boostedVec.x() <<" y: " << boostedVec.y() << " z: " << boostedVec.z() << " e: " << boostedVec.e() <<endl;
	  }
	allParticlesBoosted.push_back(boostedVec[massHyp].vect());



	//	pXCheck


	allPB_particleClass.push_back(massHyp);
	//disasterous code...
	allPB_particleCharge.push_back(charge);

	nonBoostedE.push_back(E[massHyp]);

	//	cout <<"add to allparticles boosted " <<endl;
       	allParticlesNonBoosted.push_back(h3Vect);
	allPB_E.push_back(boostedVec[massHyp].e());
	if(DEBUG_EVENT==evtNr)
	  {
	    //	  (*pXCheck).precision(3);
	    //	  (*pXCheck)<<boostedVec.vect().x() << " " <<boostedVec.vect().y() << " " << boostedVec.vect().z() << " " << m_mass <<" charged " <<endl;
	  }
	visEnergy+=boostedVec[massHyp].e();

	//this is now done by the new and better functions...
	///--- temp disables
	/////---       	findDStar(allParticlesBoosted, allPB_particleClass, allPB_particleCharge);

	Particle* pNonBoosted=0;
	if(isLepton){

	}
	else{
	  if(fabs(lundPC)==PY_PI || fabs(lundPC)==PY_K)
	    {
	      pNonBoosted=new Particle(*chr_it,string(ptypName));
	      pNonBoosted->momentum().momentum(nonBoostedVec[massHyp]);
	      bool found=false;
	      if(fabs(lundPC)==PY_PI)
		{
		  chargedPiCandidates.push_back(pNonBoosted);
		}
	      if(fabs(lundPC)==PY_K)
		{
		  chargedKCandidates.push_back(pNonBoosted);
		}
	    }
	}

	//////
	////// To use the PID unfolding we also have to save leptons and protons
	/////


	//just use minP (next line)
	if(m_z[massHyp]<cuts::minZ)
	  {
	    //	    cout <<"didn't pass min z...: "<< m_z <<", energy: " << boostedVec.e()<<endl;
	    //	    continue;
	  }
	if(labMom<cuts::minPLab || labMom>cuts::maxPLab) 
	  {

	    if(DEBUG_EVENT==evtNr)
	      {
		cout <<"cut on plabl " << labMom<<endl;
	      }
	    //	    (*pXCheck) <<" cut due to min plab " << endl;
	    continue;
	  }
	if(cos(h3Vect.theta())<cuts::minCosTheta||cos(h3Vect.theta())>cuts::maxCosTheta)
	  {
	    //	    (*pXCheck) <<" cut duecos theata " << endl;
	    if(DEBUG_EVENT==evtNr)
	      {
		cout << "CUT cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	      }
	    //	    cout << "CUT cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	    continue;
	  }
      
	if(DEBUG_EVENT==evtNr)
	  {
	    cout << "cos theta: " << cos(h3Vect.theta()) <<" charge: " << charge <<" massHyp: "<< massHyp <<endl;
	  }

	Particle* p=new Particle(*chr_it,string(ptypName));
      
	//has to be in parantheses
       	p->userInfo(ParticleInfo());
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(p->userInfo());
      
	pinf.idAs=massHyp;
	pinf.charge=charge;
       	pinf.labMom=labMom;
	for(int i=0;i<5;i++)
	  {
	    pinf.z[i]=m_z[i];
	    pinf.boostedMoms[i]=boostedVec[i].vect();
	    pinf.boostedLorentzVec[i]=boostedVec[i];
	  }
	////need to set the PID matrices
	pinf.labTheta=h3Vect.theta();

	if(setHadronPIDProbs(&pinf, labMom)<0)
	  {
	    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
	      {
		cout <<"pid exit" <<endl;
	      }
	    delete p;
	    //	    (*pXCheck) <<" problem with pid , exit event " << endl;
	    //	    cout <<"exit event" <<endl;
	    exitEvent();
	    return;
	  }

	pinf.cmsTheta=boostedVec[massHyp].theta();
	//	cout <<"theta lab:" << h3Vect.theta() <<"cms: "<< boostedVec.theta()<<endl;
	//	cout <<"theta phi:" << h3Vect.phi() <<"cms: "<< boostedVec.phi()<<endl;

	pinf.labPhi=h3Vect.phi();
	pinf.cmsPhi=boostedVec[massHyp].phi();
	Ptype& m_pt=p->pType();
	//is it ok, to leave the default error matrix?
	p->momentum().momentum(boostedVec[massHyp]);

	//	cout <<"add to all particles for comp " <<endl;
	//only hadrons, even if unidentified trust that leptons can be separated
	//	if(!isLepton)  --> changed since we apply PID unfolding
	v_allParticles.push_back(p);

      }
    /////find memory leaks...
    //    exitEvent();
    //    return;
    ////do this after we collected the reconstructed particles, so we know what we don't have to save
    pTreeSaver->saveGenInfo(v_allParticles, eventCut);
    if(eventCut)
      {
	exitEvent();
	return;
      }

    //        cout <<"done with charged particles " << endl;
    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();

    //////---->just for the collins charm stuff

    Mdst_pi0_Manager &pi0_mgr=Mdst_pi0_Manager::get_manager();
    for(std::vector<Mdst_pi0>::const_iterator i =pi0_mgr.begin();i!=pi0_mgr.end();i++)
      {
	const Mdst_pi0& pi0=*i;
	int id =(int)pi0.get_ID();
	//	cout <<"pi0 children: " << pi0.nChildren()<<endl;

	double px=pi0.px();
	double py=pi0.py();
	double pz=pi0.pz();


	Mdst_ecl_aux &aux1 =eclaux_mgr(Panther_ID(pi0.gamma(0).ecl().get_ID()));
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25_1 =aux1.e9oe25();
	Mdst_ecl_aux &aux2 =eclaux_mgr(Panther_ID(pi0.gamma(1).ecl().get_ID()));
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25_2 =aux2.e9oe25();
	double mass=pi0.mass(); //mass before fitting ???
	///	if(mass>0.15 || mass<0.12)
	//almost the same as above...
	float pLab=sqrt(px*px+py*py+pz*pz);
	//      cout <<"pi0mass: "<< mass <<endl;≈

	float g1Energy= sqrt(pi0.gamma(0).px()*pi0.gamma(0).px()+pi0.gamma(0).py()*pi0.gamma(0).py()+pi0.gamma(0).pz()*pi0.gamma(0).pz());
	float g2Energy= sqrt(pi0.gamma(1).px()*pi0.gamma(1).px()+pi0.gamma(1).py()*pi0.gamma(1).py()+pi0.gamma(1).pz()*pi0.gamma(1).pz());
	//	cout <<"pi0 gamma1: "<< g1Energy <<" gamma2: "<< g2Energy <<endl;
       	if(g1Energy < 0.1 || g2Energy < 0.1)
	  continue;


	Hep3Vector h3Vect(px,py,pz);
	float E=sqrt(mass*mass+h3Vect.mag2());
	HepLorentzVector boostedVec(h3Vect,E);
	boostedVec.boost(kinematics::CMBoost);


	//	cout <<"testing diphoton mass: "<< mass <<endl;
	if(fabs(mass-pi0Mass)<0.03)
	  {
	    //pi0
	    boostedPi0s.push_back(boostedVec);
	  }
	if(fabs(mass-etaMass)<0.03)
	  {
	    //eta
	    boostedEtas.push_back(boostedVec);
	  }
      }


    /////-----


    for(std::vector<Mdst_gamma>::const_iterator i =gamma_mgr.begin();i!=gamma_mgr.end();i++)
      {
	float energy1=sqrt(i->px()*i->px()+i->py()*i->py()+i->pz()*i->pz());
	if(energy1<0.3)
	  continue;
	HepLorentzVector h1(i->px(),i->py(),i->pz(),energy1);
	for(std::vector<Mdst_gamma>::const_iterator j=i+1;j!=gamma_mgr.end();j++)
	  {
	    float energy2=sqrt(j->px()*j->px()+j->py()*j->py()+j->pz()*j->pz());
	    if(energy2<0.3)
	      continue;
	    HepLorentzVector h2(j->px(),j->py(),j->pz(),energy2);
	    HepLorentzVector etaCandidate=h1+h2;
	    etaCandidate.boost(kinematics::CMBoost);
	    float mass=etaCandidate.mag();
	    if(fabs(mass-etaMass)<0.04)
	      {
		//eta
		boostedEtas.push_back(etaCandidate);
	      }


	  }
      }


    //    cout <<"do gammas " <<endl;
    int gammaCount=0;
    for(std::vector<Mdst_gamma>::const_iterator i =gamma_mgr.begin();i!=gamma_mgr.end();i++)
      {
	Hep3Vector h(i->px(),i->py(),i->pz());
	const Mdst_gamma& gam=*i;
	int id=(int)gam.get_ID();
	double px=gam.px();
	double py=gam.py();
	double pz=gam.pz();
	//does not make sensee because if we only look in the central region for the
	//computation of the thrust axis, that would change the axis, same for charged
	///there I take all for the thrust axis and the ones in the central region for
	//the asymmetry
	Mdst_ecl_aux &aux =eclaux_mgr(Panther_ID(gam.ecl().get_ID()));
	if(gam.ecl().quality()!=0)
	  {
	    //	    cout <<"loosing photon due to  quality " <<endl;
	    continue;
	  }
	//ratio of energy in 3x3 cluster compared to 5x5 cluster in emcal
	double e9oe25 =aux.e9oe25();
	double gammaE=sqrt(px*px+py*py+pz*pz);
	HepLorentzVector boostedVec(px,py,pz,gammaE);
	Hep3Vector photVec(px,py,pz);
	boostedVec.boost(kinematics::CMBoost);
	v_gammaE.push_back(gammaE);

	//barrel energy cut is the lowest
	if(gammaE<cuts::minGammaEBarrel)
	  {
	    //	    cout <<" loosing photon in barrel due to energy cut " << gammaE<<endl;
	    continue;
	  }
	float photTheta= photVec.theta();
	photTheta=180*(photTheta/TMath::Pi());
	if(photTheta >cuts::barrelThetaMax) 
	  {
	    if(gammaE<cuts::minGammaEEndcapBkwd)
	      {
		//	    cout <<" loosing photon in forwrad endcap due to energy cut " << gammaE<<endl;
		continue;
	      }
	  }
	if(photTheta< cuts::barrelThetaMin)
	  {
	    if(gammaE<cuts::minGammaEEndcapFwd)
	      {
		//	    cout <<" loosing photon in backward endcap due to energy cut " << gammaE<<endl;
		continue;
	      }
	  }
	allParticlesNonBoosted.push_back(photVec);
	//photon
	allPB_particleClass.push_back(-1);
	allParticlesBoosted.push_back(boostedVec.vect());
	nonBoostedE.push_back(gammaE);

	allPB_E.push_back(boostedVec.e());
	visEnergy+=boostedVec.e();
	gammaCount++;
      }
    //    cout <<" number of charged tracks in the event: pos - " << iChTrkPos <<" neg - " << iChTrkNeg <<endl;
    //    cout <<"number of photons in the event: " << gammaCount <<endl;

    //    cout <<allParticlesBoosted.size() <<" particles in jet and thrust computation " << endl;
    //    if(allParticlesBoosted.size()!=nonBoostedE.size())
    //      cout <<"e not the same size! " <<endl;
    for(int i=0;i<allParticlesNonBoosted.size();i++)
      {
	Hep3Vector& vec=allParticlesNonBoosted[i];
	//	  cout <<"------> in lab frame Px, Py, Pz, E: " << vec.x()<<" "<< vec.y() <<" " << vec.z() <<" " << nonBoostedE[i]<<endl;
      }
    for(int i=0;i<allParticlesBoosted.size();i++)
      {
	//	  cout <<"Px: " << allParticlesBoosted[i].x()<<" Py: " << allParticlesBoosted[i].y()<<" "<< allParticlesBoosted[i].z()  <<endl;
	Hep3Vector& vec=allParticlesBoosted[i];
	//	  cout <<"------> in CM frame Px, Py, Pz, E: " << vec.x()<<" "<< vec.y() <<" " << vec.z() <<" " << allPB_E[i]<< " mass_hyp: " << allPB_particleClass[i] <<endl;
      }

    ////-----> look for Ds

    //find possible Ks..
    Mdst_vee2_Manager &vee2_m=Mdst_vee2_Manager::get_manager();
    for(vector<Mdst_vee2>::iterator vee_it=vee2_m.begin();vee_it!=vee2_m.end();vee_it++)
      {
	//not a K short...
	if(vee_it->kind()!=1)
	  continue;

	/// also quality checks?
	//second parameter keeps relation
	Particle* p=new Particle(*vee_it,true);
	FindKs findks;

	findks.candidates(*vee_it,IpProfile::position());
	if(findks.goodKs())
	  {
	    KsCandidates.push_back(p);
	  }
	else
	  {
	    delete p;
	  }

      }


    //    cout <<"do Ds " <<endl;
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    ///got all candidates, now reconstruct Ds...
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    if(m_mc)
      {
	/////	recD0MC();
	/////	recDStarMC();
      }



    /////    reconstructD0();
    /////    //    if(chargedDCandidates.size()>0 || D0Candidates.size()>0)
    /////    //      cout <<"we have "  << D0Candidates.size() <<" D0s before fit" <<endl;
    /////    vector<Particle*> tempParticles;
    /////    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
    /////      {
    /////	double confLevel;
    /////	//		cout <<"mass before d0 mit: "<< (*itD)->p().mag()<<endl;
    /////	///	if(false)
    /////	if(!doKmVtxFit2(*(*itD),  confLevel,0))
    /////	  {
    /////	    //	    	    cout <<"removing D0 due to bad fit " << endl;
    /////	    delete *itD;
    /////	  }
    /////	else
    /////	  {
    /////	    //	    cout <<"D0 fit is good " << endl;
    /////	    tempParticles.push_back(*itD);
    /////	    //	    	    cout <<"mass  d0 mit: "<< (*itD)->p().mag()<<endl;
    /////	  }
    /////	
    /////      }
    /////    D0Candidates.clear();
    /////    
    /////    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
    /////      {
    /////	D0Candidates.push_back(*itD);
    /////      }
    /////    tempParticles.clear();
    /////
    /////    reconstructChargedD();
    /////    //    if(chargedDCandidates.size()>0)
    /////      //    cout <<"we have " << chargedDCandidates.size() <<" before fit " <<endl;
    /////    for(vector<Particle*>::iterator itD=chargedDCandidates.begin();itD!=chargedDCandidates.end();itD++)
    /////      {
    /////	double confLevel;
    /////	//	if(false)
    /////    	if(!doKmVtxFit2(*(*itD),  confLevel,0))
    /////	  {
    /////	    //	    printD();
    /////	    //	    	    cout <<" no good charged D " <<endl;
    /////	    //	    cout <<"removing charged D due to bad fit " << endl;
    /////	    delete *itD;
    /////	  }
    /////	else
    /////	  {
    /////	    //    cout <<"good fit for charged D " << endl;
    /////	    //	    cout <<" good charged D " << endl;
    /////	    tempParticles.push_back(*itD);
    /////	  }
    /////      }
    /////    chargedDCandidates.clear();
    /////
    /////    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
    /////      {
    /////	chargedDCandidates.push_back(*itD);
    /////      }
    /////    tempParticles.clear();
    /////    //    if(chargedDCandidates.size()>0 || D0Candidates.size()>0)
    /////    //      cout <<"we have " << chargedDCandidates.size() <<" charged and " << D0Candidates.size() <<" D0s after fit" <<endl;
    /////    reconstructDStar();
    /////    for(vector<Particle*>::iterator itD=DStarCandidates.begin();itD!=DStarCandidates.end();itD++)
    /////      {
    /////	double confLevel;
    /////	//	cout <<" d star before: "<<	(*itD)->p().mag() <<endl;
    /////	//in principle mass-vertex constrained fit here (not just mass...)
    /////	if(!doKmVtxFit2(*(*itD),  confLevel,0))
    /////	  {
    /////	    //	    	    cout <<" no good dstar .. " <<endl;
    /////	    delete *itD;
    /////	  }
    /////	else
    /////	  {
    /////	    //	          if(kinematics::DStarDecay==2){
    /////		    //		    cout <<"found good dstar in correct decay " <<endl;
    /////	    //		  }
    /////	    //	    	      cout <<"found good dstar .." <<endl;
    /////	    tempParticles.push_back(*itD);
    /////	    //	cout <<" d star after: "<<	(*itD)->p().mag() <<endl;
    /////	  }
    /////
    /////      }
    /////
    /////    DStarCandidates.clear();
    /////    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
    /////      {
    /////	DStarCandidates.push_back(*itD);
    /////      }
    /////    tempParticles.clear();
    /////
    /////      /////----> end look for Ds

    //    if(D0Candidates.size()>0|| chargedDCandidates.size()>0 )
    //      {
    //		kinematics::D0Tag=1;
    //      }
    //    else
    //      {
    //	kinematics::D0Tag=0;
    //      }
    //    if(DStarCandidates.size()>0)
    //      {
    //	kinematics::DStarTag=1;
    //      }
    //    else
    //      {
    //	kinematics::DStarTag=0;
    //      }
    //
    //    cout <<"doing thrust "<<endl;
    Thrust t=thrustall(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    //    Thrust t2=thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);


    Thrust labThrust=thrustall(allParticlesNonBoosted.begin(),allParticlesNonBoosted.end(),retSelf);
    ///jet computations:

    //     cout <<"got  lab thrust " <<endl;
    //            cout <<"lab thrust theta: " << labThrust.axis.theta() <<endl;
    //   	cout <<"thrust cms theta: "<< t.axis.theta() <<" phi: "<< t.axis.phi() << " mag: "<< t.thru<<endl;
    //     cout <<"lab thrust phi: " << labThrust.axis.phi() <<endl;

    //    cout <<" Thrust cos theta in lab system: " << cos(labThrust.axis.theta())<<endl;
    //    cout <<" Thrust cos theta in CMS system: " << cos(t.axis.theta())<<endl;
    //    cout <<"thrust px: " << t.axis.x() << " py: "<< t.axis.y()<<" pz: " << t.axis.z()<<endl;

    //

    if(cos(labThrust.axis.theta()) > cuts::maxLabThrustCosTheta || cos(labThrust.axis.theta())<cuts::minLabThrustCosTheta)
      {
	//	cout <<" cut on lab axis " <<endl;
	//	    (*pXCheck) <<" exit event due to lab thrust " << endl;
	if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
	  {
	    cout <<"lab thrust exit" <<endl;
	  }
	//we don't cut on thrust values anymore...
	//	exitEvent();
	//	return;
      }

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      {
	//	cout <<"particles in thrust computation" <<endl;
	for(int i=0;i<allParticlesBoosted.size();i++)
	  {
	    Hep3Vector l=allParticlesBoosted[i];
	    //	    cout <<"x: " << l.x() << " y: " << l.y() <<" z: " << l.z() << endl;
	  }
      }


    //  kinematics::thrustDirCM=thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    kinematics::thrustDirCM=t.axis;


    //        if(rand() % 100 <50)
    if(kinematics::thrustDirCM.z()<0)
      {
	kinematics::thrustDirCM.setZ((-1)*kinematics::thrustDirCM.z());
	kinematics::thrustDirCM.setY((-1)*kinematics::thrustDirCM.y());
	kinematics::thrustDirCM.setX((-1)*kinematics::thrustDirCM.x());
	kinematics::thrustZReverted=true;
      }
    else
      {
	kinematics::thrustZReverted=false;
      }

    kinematics::thrustDirLab=labThrust.axis;
    kinematics::thrustMag=t.thru;

    Hep3Vector tDiff=t.axis-thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    visEnergyOnFile=hadronInfo_mgr.begin()->Evis();
    //  (*pXCheck)<<"diff: " << tDiff.x() << " " << tDiff.y() << " " << tDiff.z() <<endl;
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      {	//  cout <<"boost: " << kinematics::CMBoost.x() << " y: " << kinematics::CMBoost.y() << " " << kinematics::CMBoost.z() <<endl;
	//      cout <<"thrustDirCM: " << kinematics::thrustMag << " thrustDir z: " << abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()<<" vorher: " << kinematics::thrustDirCM.z()<< " visE: " << visEnergy<<" on file: " << visEnergyOnFile << " ichTrak: " << iChTrks <<" tdx: " << kinematics::thrustDirCM.x() << " tdy: " << kinematics::thrustDirCM.y() <<endl;
      }

    kinematics::E_miss=kinematics::Q-visEnergyOnFile;
    //    if(kinematics::thrustMag<cuts::minThrust || abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ|| visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
    if( visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
      {
	bool foundReason=false;
	if(visEnergyOnFile < cuts::minVisEnergy)
	  {
	    //	    (*pXCheck) <<" exit event due to min vis e " << endl;
	    //	    	    	  cout <<"---------------------------"<<endl<<"cut on vis energy: " << cuts::minVisEnergy <<" found only: " <<visEnergyOnFile <<" GeV" <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
	      {
		cout <<"exit due to minvis energy" <<visEnergyOnFile <<endl;
	      }
	    
	    foundReason=true;
	  }
	if(abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ)
	  {
	    //	    (*pXCheck) <<" exit event due thrustz " << endl;
	    //	    	    cout <<"------event nr: "<< kinematics::evtNr <<" ---------------------"<<endl<<"cut on thrust z projection in CM: " << cuts::maxThrustZ <<" found  " <<kinematics::thrustDirCM.z()/kinematics::thrustDirCM.mag() <<endl;
	    //	  cout <<"-------------------------------"<<endl;


	    foundReason=true;
	  }
	if(kinematics::thrustMag<cuts::minThrust)
	  {
	    //	    (*pXCheck) <<" exit event due to minthrust " << endl;
	    //	    cout <<"---------------------------"<<endl<<"cut magnitude: " << cuts::minThrust <<" found  " <<kinematics::thrustMag <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	    foundReason=true;
	  }
	if( iChTrks < cuts::minNTracks)
	  {
	    
	    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
	      {
		cout <<"exit due to ntracks" << iChTrks<<endl;
	      }
	    //	    (*pXCheck) <<" exit event due to min ntracks " << endl;
	    //	        cout <<"-----------------------"<<endl <<" cut on min tracks, need; "<< cuts::minNTracks <<" have: " << iChTrks <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	    foundReason=true;
	  }

	if(!foundReason)
	  {
	    //	  cout <<"exiting event due to some cut" <<endl;
	  }
	if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
	  {
	    cout <<"exit due to some other reason" <<endl;
	  }

	exitEvent();
	return;

      }
    //    cout <<"passed event cut " <<endl;
    for(int i=0;i<allParticlesBoosted.size();i++)
      {
	float ltheta=AuxFunc::getTheta(kinematics::thrustDirCM,allParticlesBoosted[i]);
	float m_z=2*allPB_E[i]/kinematics::Q;
	//	m_histos.hEFlowFromThrust->Fill(ltheta,m_z);
      }
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"evt good" <<endl;
      }
    //  cout <<"thrust cms: " << kinematics::thrustDirCM.theta() <<endl;

    //  cout <<"thrust theta, phi: " << kinematics::thrustDirCM.theta() <<" / " << kinematics::thrustDirCM.phi()<<endl;
    //  cout <<"thrust cms after possible flip: " << kinematics::thrustDirCM.theta() <<endl;
    kinematics::thetaEThrust=kinematics::thrustDirCM.angle(kinematics::firstElectronCM.vect());
    //  cout <<"theta - e beam angle: " << kinematics::thetaEThrust <<" thrust theta: "  <<kinematics::thrustDirCM.theta()<<endl;
    //    cout <<"pi0"<<endl;

    int sB=v_allParticles.size();

    //////////

    setParticleProperties();
#ifdef ThrustMethod
    findHadronPairsThrust();
#else
    findHadronPairs();
#endif
    Hep3Vector axis=kinematics::thrustDirCM;
    ////------ stuff for hairongs analysis
    //     cout <<"found " << boostedPi0s.size() <<" pi0s and " << boostedEtas.size() <<" etas " <<endl;
    for(int i=0;i<boostedPi0s.size();i++)
      {
	Hep3Vector& vec=boostedPi0s[i];
	float thrustProj=axis.dot(vec);
	float pi0E=sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z()+pi0Mass*pi0Mass);
	float pi0Z=2*pi0E/10.5;
	if(pi0Z<0.2)
	  continue;
	float pi0Pt=axis.perp(vec);
	//	  cout <<"pi0Z: "<< pi0Z<<" pt: " << pi0Pt <<endl;
	if(pi0Pt<1.5)
	  {
	    int zBin=getBin(zBorders,4,pi0Z);
	    int ptBin=getBin(ptBorders,4,pi0Pt);
	    //  cout <<"zBin: " << zBin <<" ptBin: "<< ptBin <<endl;
	    numPi0s[zBin*4+ptBin]++;
	  }
	//     pinf.thrustProj=axis.dot((*it)->p().vect())/(axis.mag()*(*it)->p().vect().mag());
      }

    for(int i=0;i<boostedEtas.size();i++)
      {
	Hep3Vector& vec=boostedEtas[i];
	float thrustProj=axis.dot(vec);
	float etaE=sqrt(vec.x()*vec.x()+vec.y()*vec.y()+vec.z()*vec.z()+pi0Mass*pi0Mass);
	float etaZ=2*etaE/10.5;
	if(etaZ<0.3)
	  continue;
	float etaPt=axis.perp(vec);
	if(etaPt<1.5)
	  {
	    int zBin=getBin(zBorders,4,etaZ);
	    int ptBin=getBin(ptBorders,4,etaPt);
	    numEtas[zBin*4+ptBin]++;
	  }
	//     pinf.thrustProj=axis.dot((*it)->p().vect())/(axis.mag()*(*it)->p().vect().mag());
      }




    ////----




#ifdef SAVE_HISTOS
    saveHistos(allParticlesBoosted, allParticlesNonBoosted);
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      {
	cout <<"about to histos" <<endl;}

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2){
      cout <<"done saving histos" <<endl;}
#endif
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2){
      cout <<"before save tree" <<endl;}
    //  cout <<"vor tree" <<endl;



#ifdef XCHECK
    int qCounter=-1;
    (*pXCheck) << std::setw(4);
    (*pXCheck) <<" runNr: "<<runNr <<" eventNr: "<< evtNr;
    (*pXCheck) <<" thrust Mag: " <<kinematics::thrustMag;
    (*pXCheck) <<" thrust dir (CM) theta: "<<     kinematics::thrustDirCM.theta();
    (*pXCheck) <<" phi "<<     kinematics::thrustDirCM.phi()<<endl;

    //    (*pXCheck) <<" thrust dir (Lab) theta: "<<     kinematics::thrustDirLab.theta();
    //    (*pXCheck) <<" phi "<<     kinematics::thrustDirLab.phi()<<endl;

    for(vector<HadronPair*>::iterator it=v_hadronPairs.begin();it!=v_hadronPairs.end();it++)
      {
	qCounter++;
	HadronPair* hp=(*it);
	//      if(evtNr==1)
	{
	  ParticleInfo& pinf1=dynamic_cast<ParticleInfo&>(hp->firstHadron->userInfo());
	  ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>(hp->secondHadron->userInfo());

	  //momentum ordering for easier comparison
	  if(pinf1.labMom<pinf2.labMom)
	    {
	      pinf1=dynamic_cast<ParticleInfo&>(hp->secondHadron->userInfo());
	      pinf2=dynamic_cast<ParticleInfo&>(hp->firstHadron->userInfo());
	    }

	  //	(*pXCheck) <<" two had type; "<< hp->hadPType<<endl;
	  //seems to be in cms...
	  //	(*pXCheck) << "other mom: "<<  hp->firstHadron->p().vect().mag() <<" or: "<< hp->secondHadron->p().vect().mag()<<endl;
	  (*pXCheck) <<"plab h1 "<< pinf1.labMom << " h2 " << pinf2.labMom;
	  (*pXCheck) << " costheta h1 " << cos(pinf1.labTheta)<<"  h2 " << cos(pinf2.labTheta)<<endl;

	  (*pXCheck) <<"i 0 pid h1 data "<< pinf1.p_Pi <<" mc " << pinf1.p_Pi2;
	  (*pXCheck) <<" pid h2 data "<< pinf2.p_Pi <<" mc " << pinf2.p_Pi2<<endl;

	  (*pXCheck) <<"i 1 pid h1 data "<< pinf1.p_K <<" mc " << pinf1.p_K2;
	  (*pXCheck) <<" pid h2 data"<< pinf2.p_K <<" mc " << pinf2.p_K2<<endl;

	  (*pXCheck) <<"i 2 pid h1 data "<< pinf1.p_p <<" mc " << pinf1.p_p2;
	  (*pXCheck) <<" pid h2 data "<< pinf2.p_p <<" mc " << pinf2.p_p2<<endl;


	  (*pXCheck) << "i 0 z h1 " << hp->z1_PiPi << " z h2 "<< hp->z2_PiPi <<" kT " << hp->kT_PiPi << " qT " << hp->qT_PiPi <<endl;
	  (*pXCheck) << "i 1 z h1 " << hp->z1_PiK << " z h2 "<< hp->z2_PiK <<" kT " << hp->kT_PiK << " qT " << hp->qT_PiK <<endl;
	  (*pXCheck) << "i 2 z h1 " << hp->z1_PiP << " z h2 "<< hp->z2_PiP <<" kT " << hp->kT_PiP << " qT " << hp->qT_PiP <<endl;

	  (*pXCheck) << "i 3 z h1 " << hp->z1_KPi << " z h2 "<< hp->z2_KPi <<" kT " << hp->kT_KPi << " qT " << hp->qT_KPi <<endl;
	  (*pXCheck) << "i 4 z h1 " << hp->z1_KK << " z h2 "<< hp->z2_KK <<" kT " << hp->kT_KK << " qT " << hp->qT_KK <<endl;
	  (*pXCheck) << "i 5 z h1 " << hp->z1_KP << " z h2 "<< hp->z2_KP <<" kT " << hp->kT_KP << " qT " << hp->qT_KP <<endl;

	  (*pXCheck) << "i 6 z h1 " << hp->z1_PPi << " z h2 "<< hp->z2_PPi <<" kT " << hp->kT_PPi << " qT " << hp->qT_PPi <<endl;
	  (*pXCheck) << "i 7 z h1 " << hp->z1_PK << " z h2 "<< hp->z2_PK <<" kT " << hp->kT_PK << " qT " << hp->qT_PK <<endl;
	  (*pXCheck) << "i 8 z h1 " << hp->z1_PP << " z h2 "<< hp->z2_PP <<" kT " << hp->kT_PP << " qT " << hp->qT_PP <<endl;

	  //(*pXCheck)<< "qT " <<hp->qT <<endl;
	  //	(*pXCheck) << "z1: " << hp->z1 <<" z2: " << hp->z2 <<" kT: " << hp->kT <<  " pid1 " << hp->firstHadron->lund() <<" pid2 " << hp->firstHadron->lund() <<endl;

	  //cos(pinf1.labTheta) <<" cos labTheta2: "<< cos(pinf2.labTheta) <<" z1: "<< hp->z1 << " z2: "<< hp->z2 << " kT: "<< hp->kT <<", qT: "<< hp->qT << " had type1 "<< hp->hadPType1 <<" second had type: "<< hp->hadPType2  <<" first charge: "<< hp->hadCharge1 << " second charge: "<< hp->hadCharge2 <<endl;
	  //	(*pXCheck) << "naive z1: "<< pinf1.boostedMoms[pinf1.idAs].mag()/5.25 <<" naive z2: "<< pinf2.boostedMoms[pinf2.idAs].mag()/5.25 <<endl;
	  //	(*pXCheck) << " lab Theta1 "<< pinf1.labTheta <<" cosLabTheta: "<< cos(pinf1.labTheta) <<" labTheta2: "<< pinf2.labTheta <<" cosLt2: ";
	  //	(*pXCheck) << cos(pinf2.labTheta)<<endl;
	  //	(*pXCheck) << " cms Theta1 "<< pinf1.cmsTheta <<" cosCmsTheta: "<< cos(pinf1.cmsTheta) <<" cmsTheta2: "<< pinf2.cmsTheta <<" cosLt2: ";
	  //	(*pXCheck) << cos(pinf2.cmsTheta)<<endl;
	}
	//      (*pXCheck)<<setprecision(4);
      }
    (*pXCheck) <<endl;
#endif
    saveTree();

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2){
      cout <<"dones saving tree" <<endl;}

    cleanUp();
    //      cout <<"aft cleaning"<<endl;
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2){
      cout <<"cleaning"<<endl;}
   
  }

  bool ptSpect::enoughSVDHits(Mdst_charged_Manager::iterator chr_it)
  {
    Mdst_trk& mdsttrk = chr_it->trk();
    Mdst_trk_fit& trk_fit=mdsttrk.mhyp(2);
    if(!mdsttrk)
      return false;
    //    HepPoint3D pivot(trk_fit.pivot_x(),trk_fit.pivot_y(),trk_fit.pivot_z());
    if(trk_fit.nhits(3)+trk_fit.nhits(4)<3)
      return false;

    return true;

  }


  void ptSpect::getDrDz(Mdst_charged_Manager::iterator chr_it, int masshyp, double& dr, double& dz, double& refitPx, double& refitPy, double& refitPz)
  {
    Mdst_trk& mdsttrk = chr_it->trk();
    Mdst_trk_fit &mdsttrkfit=mdsttrk.mhyp(masshyp);
    HepPoint3D pivot(mdsttrkfit.pivot_x(),mdsttrkfit.pivot_y(),mdsttrkfit.pivot_z());

    HepVector a( 5, 0 );
    a[0] = mdsttrkfit.helix( 0 ); // helix parameters defined at the pivot
    a[1] = mdsttrkfit.helix( 1 );
    a[2] = mdsttrkfit.helix( 2 );
    a[3] = mdsttrkfit.helix( 3 );
    a[4] = mdsttrkfit.helix( 4 );

    HepSymMatrix Ea( 5, 0 );
    Ea[0][0] = mdsttrkfit.error( 0 );
    Ea[1][0] = mdsttrkfit.error( 1 );
    Ea[1][1] = mdsttrkfit.error( 2 );
    Ea[2][0] = mdsttrkfit.error( 3 );
    Ea[2][1] = mdsttrkfit.error( 4 );
    Ea[2][2] = mdsttrkfit.error( 5 );
    Ea[3][0] = mdsttrkfit.error( 6 );
    Ea[3][1] = mdsttrkfit.error( 7 );
    Ea[3][2] = mdsttrkfit.error( 8 );
    Ea[3][3] = mdsttrkfit.error( 9 );
    Ea[4][0] = mdsttrkfit.error( 10 );
    Ea[4][1] = mdsttrkfit.error( 11 );
    Ea[4][2] = mdsttrkfit.error( 12 );
    Ea[4][3] = mdsttrkfit.error( 13 );
    Ea[4][4] = mdsttrkfit.error( 14 );


    Helix helix( pivot, a, Ea );
    helix.pivot( IpProfile::position(1));
    //    cout <<"helix momentum: "<< helix.momentum()<<endl;
    refitPx=helix.momentum().x();
    refitPy=helix.momentum().y();
    refitPz=helix.momentum().z();
    //    HepLorentzVector boostedVec(helix.momentum(),sqrt(helix.momentum().mag2()+m_mass*m_mass));
    dr  = helix.dr();
    dz  = helix.dz();
  }

  void ptSpect::exitEvent()
  {
#ifdef SAVE_HISTOS
    // saveHistos(allParticlesBoosted, allParticlesNonBoosted);
#endif
    cleanUp();
  }
  // begin_run function

  void ptSpect::end_run(BelleEvent* evptr, int* status)
  {
    std::cout << "ptSpect's end_run function" << std::endl;
    //
  }

  void ptSpect::saveTree()
  {


    if(DEBUG_EVENT==kinematics::evtNr)
      {
	cout <<" saving " << v_hadronPairs.size() <<" pairs" <<endl;
      }
    pTreeSaver->fillWPairData(v_hadronPairs,m_evtInfo);
  }
  void ptSpect::saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted)
  {

    for(int i=0;i<v_vertexR.size();i++)
      {
	//	m_histos.hVertexR->Fill(v_vertexR[i]);
	//	m_histos.hVertexZ->Fill(v_vertexZ[i]);
      }
    for(int i=0;i<v_pi0GammaE.size();i++)
      {
	//	m_histos.hPi0GammaE->Fill(v_pi0GammaE[i]);
      }
    for(int i=0;i<v_gammaE.size();i++)
      {
	//	m_histos.hGammaE->Fill(v_gammaE[i]);
      }
    for(int i=0;i<v_asyms.size();i++)
      {
	//	m_histos.hGammaAsym->Fill(v_asyms[i]);
      }
    //for comparison with the thrust from the manager
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    //    m_histos.hThrustOnFile->Fill(hadronInfo_mgr.begin()->Thrust());



    for(int i=0;i<v_allParticles.size();i++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_allParticles[i]->userInfo());
	//      cout <<"sveing phi: " << pinf.phi <<", theta: " << pinf.theta<<endl;
	if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	  continue;

	//	m_histos.hz->Fill(pinf.z[pinf.idAs]);
	//	m_histos.hTheta->Fill(pinf.cmsTheta);
	//	m_histos.hCosTheta->Fill(cos(pinf.cmsTheta));
	//	m_histos.hThrustProj->Fill(pinf.thrustProj);



	if(v_allParticles[i]->charge()>0)
	  {

	    //	    m_histos.hThetaPos->Fill(pinf.cmsTheta);
	    //	    m_histos.hZPos->Fill(pinf.z[pinf.idAs]);
	  }
	else
	  {
	    if(v_allParticles[i]->charge()==0)
	      {
		//		m_histos.hPi0Mass->Fill(dynamic_cast<ParticleInfoMass&>(pinf).mass);

		//		m_histos.hThetaNeut->Fill(pinf.cmsTheta);
		//		m_histos.hZNeut->Fill(pinf.z[pinf.idAs]);
	      }
	    else
	      {

		//		m_histos.hThetaNeg->Fill(pinf.cmsTheta);
		//		m_histos.hZNeg->Fill(pinf.z[pinf.idAs]);
	      }
	  }
	//	m_tup->dumpData();
      }



    //    m_histos.hThrust->Fill(kinematics::thrustMag);
    //    m_histos.hThrustX->Fill(kinematics::thrustDirCM.x());
    //    m_histos.hThrustY->Fill(kinematics::thrustDirCM.y());
    //    m_histos.hThrustZ->Fill(kinematics::thrustDirCM.z());
  }

  bool ptSpect::findJetMatch(Hep3Vector& vec, vector<PseudoJet>* const1, vector<PseudoJet>* const2)
  {
    for(unsigned j=0;j<const1->size();j++)            {
      if((fabs((*const1)[j].px()-vec.x()) <0.001)&& (fabs((*const1)[j].py()-vec.y()) <0.001) && (fabs((*const1)[j].pz()-vec.z()) <0.001))		{
	return true;
      }
    }

    for(unsigned j=0;j<const2->size();j++)            {
      if((fabs((*const2)[j].px()-vec.x()) <0.001)&& (fabs((*const2)[j].py()-vec.y()) <0.001) && (fabs((*const2)[j].pz()-vec.z()) <0.001))		{
	return true;
      }
    }

    return false;
  }

  void ptSpect::setParticleProperties()
  {
    int  evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    for(vector<Particle*>::const_iterator it=v_allParticles.begin();it!=v_allParticles.end();it++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	bool firstHemi=false;
	Hep3Vector axis=kinematics::thrustDirCM;
	if(kinematics::thrustDirCM.dot((*it)->p().vect())>0)
	  {
	    firstHemi=true;
	  }

	float phi=getPhi(axis,**it);
	//      cout <<"phi is: " << phi<< endl;
	float theta=getTheta(axis,**it);

	pinf.phi=phi;
	///leave theta alone!!! Leave as lab theta...
	///      pinf.cmsTheta=theta;

	//      pinf.thxrustProj=kinematics::thrustDirCM.dot((*it)->p().vect())/(kinematics::thrustDirCM.mag()*(*it)->p().vect().mag());
	pinf.thrustProj=axis.dot((*it)->p().vect())/(axis.mag()*(*it)->p().vect().mag());
	if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	  {
	    //    if(evtNr==DEBUG_EVENT)
	    {
	      //		cout <<"cut min thrustProj ("<<pinf.thrustProj << " required: " << cuts::minThrustProj<<endl;     
	    }
	    //	  cout <<"cut with " <<fabs(pinf.thrustProj)<<endl;
	    continue;
	  }

	/////check if this matches....
	Hep3Vector tmpVec=(*it)->p().vect();

    
	//      if(pinf.thrustProj>0)//one hemisphere
	if(firstHemi)
	  {
	    if((*it)->pType().charge()!=0)
	      {
		//		(*pXCheck )<< "putting particle with cos lab theta: "<< cos(pinf.labTheta) <<" in first hemi, naive z: " << tmpVec.mag()/5.25<<endl;
		if(DEBUG_EVENT==evtNr)
		  {
		    cout <<"pushing in first hemi p: "<< pinf.labMom<<endl;
		  }

		v_firstHemi.push_back(*it);
	      }

	  }
	else//particle is in the other hemisphee
	  {
	    if((*it)->pType().charge()!=0)
	      {
		//		(*pXCheck )<< "putting particle with cos lab theta: "<< cos(pinf.labTheta) <<" in second hemi, naive z: " << tmpVec.mag()/5.25<< endl;
		//rather pointers...
		if(DEBUG_EVENT==evtNr)
		  {
		    cout <<"pushing in second hemi p: "<< pinf.labMom<<endl;
		  }

		v_secondHemi.push_back(*it);
	      }

	  }

      }
  }

  //don't use thrust to make the pairs
  //we have the hadrons in two sets (first hemi, second hemi) which are meaningless for this application
  // so we can either put them in the same and then do all combinations 
  //it should be the same to do all combinations of the first with the second and then the first with the first and the second with the second
  //probably (hopefully)  get teh similar result by just running over v_allParticles
  void ptSpect::findHadronPairs()
  {
    float dotProduct[25];

    if(DEBUG_EVENT==kinematics::evtNr)
      {
	cout <<"combining hadrons: "<< v_firstHemi.size() <<" second: "<< v_secondHemi.size() <<endl;
      }
    

    for(vector<Particle*>::const_iterator it=v_firstHemi.begin();it!=v_firstHemi.end();it++)
      {
 	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_secondHemi.begin();it2!=v_secondHemi.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	    //not back to back according to thrustless definition (note that being in first and second hemisphere is not necessarily enough
	    //	    if((*it2)->sp().vect().dot((*it)->p().vect())>0)
	    bool acceptPair=false;
	    for(int i=0;i<5;i++)
	      {
		for(int j=0;j<5;j++)
		  {
		    dotProduct[i*5+j]=pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]);
		    //		      cout <<"dot prodcut (1) " <<pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]) <<endl;
		    if(pinf.boostedMoms[i].dot(pinf2.boostedMoms[j])<0)
		      {
			acceptPair=true;

		      }

		    if(DEBUG_EVENT==kinematics::evtNr)
		      {
			cout <<"looking to combine first/second p: "<< pinf.labMom <<" and " << pinf2.labMom <<endl;
			cout << " dot product: "<< pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]) <<endl;
			cout <<" accept: ("<<i <<" ) " << acceptPair <<endl;
		      }

		  }
	      }

	    if(!acceptPair)
	      continue;
	    //now unknowns...
	    HadronPair* hp=new HadronPair();
	    //	    HadronPair* hp2=new HadronPair();
	    //	    hp2->secondRun=true;
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->setDotProducts(dotProduct);
	    //	    hp2->firstHadron=*it2;
	    //	    hp2->secondHadron=*it;
	    
	    //	    hp->hadCharge=AnaDef::PN; -->let this be set automatically
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType()); //meaningless, since we know deal in probabilities
	    //	    hp2->hadPType=AuxFunc::getPType((*it2)->pType(),(*it)->pType()); //meaningless, since we know deal in probabilities
	    //	    cout <<"setting ptype in data: "<< hp->hadPType<<endl;
	    //	  cout <<"R1: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);



	    hp->compute();
	    //	    hp2->compute();
	    if(kinematics::evtNr==DEBUG_EVENT || kinematics::evtNr==DEBUG_EVENT2)
	      {
		cout <<"evt: " << kinematics::evtNr <<" putting hadron pair" <<endl;
		cout <<"prop pion 1: "<< hp->p_PiPi << " first: "<< pinf.p_Pi <<" second: "<< pinf2.p_Pi <<endl;
	      }
	    v_hadronPairs.push_back(hp);

	    //decided not put the other combination in anymore...
	    //	    v_hadronPairs.push_back(hp2);
	  }
      }
    //  cout <<v_hadronPairs.size() <<" pairs after first with second" <<endl;
    for(vector<Particle*>::const_iterator it=v_firstHemi.begin();it!=v_firstHemi.end();it++)
      {
 	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=it+1;it2!=v_firstHemi.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());

	    //	    if((*it2)->p().vect().dot((*it)->p().vect())>0)
	    bool acceptPair=false;
	    for(int i=0;i<5;i++)
	      {
		for(int j=0;j<5;j++)
		  {
		    dotProduct[i*5+j]=pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]);
		    //		      cout <<"dot prodcut (2) " <<pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]) <<endl;
		    if(DEBUG_EVENT==kinematics::evtNr)
		      {
			cout <<"looking to combine p: "<< pinf.labMom <<" and " << pinf2.labMom <<endl;
			cout << " dot product: "<< pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]) <<endl;
			
		      }
		    if(pinf.boostedMoms[i].dot(pinf2.boostedMoms[j])<0)
		      acceptPair=true;
		  }
	      }

	    if(!acceptPair)
	      continue;
	    //at least one mass hypothesis is fine
	    //now unknowns...
	    HadronPair* hp=new HadronPair();
	    //		HadronPair* hp2=new HadronPair();
	    //		hp2->secondRun=true;

	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->setDotProducts(dotProduct);

	    //		hp2->firstHadron=*it2;
	    //		hp2->secondHadron=*it;
		
	    //	    hp->hadCharge=AnaDef::PN; -->let this be set automatically
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType()); //meaningless, since we know deal in probabilities
	    //		hp2->hadPType=AuxFunc::getPType((*it2)->pType(),(*it)->pType()); //meaningless, since we know deal in probabilities
	    //	    cout <<"setting ptype in data: "<< hp->hadPType<<endl;
	    //	  cout <<"R1: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->compute();
	    //		hp2->compute();
	    if(kinematics::evtNr==DEBUG_EVENT || kinematics::evtNr==DEBUG_EVENT2)
	      {
		cout <<"prop pion 1: "<< hp->p_PiPi << " first: "<< pinf.p_Pi <<" second: "<< pinf2.p_Pi <<endl;
		cout <<"putting hadron pair" <<endl;
	      }

	    v_hadronPairs.push_back(hp);
	    //		v_hadronPairs.push_back(hp2);
	  }
      }  
    //  cout <<v_hadronPairs.size() <<" pairs after first with first" <<endl;
    for(vector<Particle*>::const_iterator it=v_secondHemi.begin();it!=v_secondHemi.end();it++)
      {
 	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=it+1;it2!=v_secondHemi.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());

	    //	    if((*it2)->p().vect().dot((*it)->p().vect())>0)
	    bool acceptPair=false;
	    for(int i=0;i<5;i++)
	      {
		for(int j=0;j<5;j++)
		  {
		    dotProduct[i*5+j]=pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]);
		    //		      cout <<"dot prodcut (3) " <<pinf.boostedMoms[i].dot(pinf2.boostedMoms[j]) <<endl;
		    if(pinf.boostedMoms[i].dot(pinf2.boostedMoms[j])<0)
		      acceptPair=true;
		  }
	      }

	    if(!acceptPair)
	      continue;

	    //now unknowns...
	    HadronPair* hp=new HadronPair();
	    //		HadronPair* hp2=new HadronPair();
	    //		hp2->secondRun=true;

	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->setDotProducts(dotProduct);

	    //		hp2->firstHadron=*it2;
	    //		hp2->secondHadron=*it;
		
	    //	    hp->hadCharge=AnaDef::PN; -->let this be set automatically
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType()); //meaningless, since we know deal in probabilities
	    //		hp2->hadPType=AuxFunc::getPType((*it2)->pType(),(*it)->pType()); //meaningless, since we know deal in probabilities
	    //	    cout <<"setting ptype in data: "<< hp->hadPType<<endl;
	    //	  cout <<"R1: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->compute();
	    //		hp2->compute();
	    if(kinematics::evtNr==DEBUG_EVENT || kinematics::evtNr==DEBUG_EVENT2)
	      {
		cout <<"prop pion 1: "<< hp->p_PiPi << " first: "<< pinf.p_Pi <<" second: "<< pinf2.p_Pi <<endl;
		cout <<"putting hadron pair" <<endl;
	      }

	    v_hadronPairs.push_back(hp);
	    //		v_hadronPairs.push_back(hp2);
	  }
      }  
    //  cout <<v_hadronPairs.size() <<" pairs after second with second" <<endl;

  }

  void ptSpect::findHadronPairsThrust()
  {
    int  evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();

    for(vector<Particle*>::const_iterator it=v_firstHemi.begin();it!=v_firstHemi.end();it++)
      {
	//don't continue if it is not pion or kaon, update: no kaons, save space
	if(!((*it)->pType()==cPiPlus))// || (*it)->pType()==cKPlus))
	  {
	    if(evtNr==DEBUG_EVENT)
	      {
		cout <<"one pos is not pi/k" <<endl;
	      }
	    //	    continue;
	  }
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	for(vector<Particle*>::const_iterator it2=v_secondHemi.begin();it2!=v_secondHemi.end();it2++)
	  {
	    ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());

	    //now unknowns...
	    HadronPair* hp=new HadronPair();
	    //	    HadronPair* hp2=new HadronPair();
	    hp->thrustMethod=true;
	    //	    hp2->thrustMethod=true;
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;

	    //	    hp2->firstHadron=*it2;
	    //	    hp2->secondHadron=*it;


	    //	    hp->hadCharge=AnaDef::PN; -->let this be set automatically
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType()); //meaningless, since we know deal in probabilities
	    //	    hp2->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType()); //meaningless, since we know deal in probabilities
	    //	    cout <<"setting ptype in data: "<< hp->hadPType<<endl;
	    //	  cout <<"R1: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->compute();
	    //	    hp2->compute();
	    v_hadronPairs.push_back(hp);
	    //	    v_hadronPairs.push_back(hp2);
	  }
      }

  }


  void ptSpect::cleanUp()
  {
    v_vertexR.clear();
    v_vertexZ.clear();
    v_pi0GammaE.clear();
    v_gammaE.clear();
    v_asyms.clear();
    //is a local variable of event() anyways
    //    allParticlesBoosted.clear();


    for(int i=0;i<D0Candidates.size();i++){
      delete D0Candidates[i];
    }
    D0Candidates.clear();

    //    cout <<" done D0"<<D0Candidates.size()<<endl;

    //    cout <<"  DStar size " <<DStarCandidates.size()<< endl;
    for(int i=0;i<DStarCandidates.size();i++){
      delete DStarCandidates[i];
    }
    DStarCandidates.clear();
    //    cout <<" done DStar " <<DStarCandidates.size()<< endl;

    //    cout <<" charged D "<<chargedKCandidates.size()<<endl;
    for(int i=0;i<chargedDCandidates.size();i++){

      delete chargedDCandidates[i];
    }
    chargedDCandidates.clear();
    //    cout <<" done charged D "<<chargedDCandidates.size()<<endl;
    //    cout <<" charged K "<<chargedKCandidates.size()<<endl;

    //don't delete, since they are already deleted as part of v_allParticles
    //needs to be deleted, is done later, but if we already clear her, it won't get deleted...
    //    chargedKCandidates.clear();
    //    cout <<" done charged K"<< chargedKCandidates.size()<<endl;


    //    cout <<"pi0s: "<< pi0Candidates.size()<<endl;
    for(int i=0;i<pi0Candidates.size();i++){
      delete pi0Candidates[i];
    }
    pi0Candidates.clear();
    //    cout <<" done pi0"<<endl;


    //    cout <<"Ks candidates: " << KsCandidates.size()<<endl;
    for(int i=0;i<KsCandidates.size();i++){
      delete KsCandidates[i];
    }
    KsCandidates.clear();


    //    cout <<" done Ks"<<endl;
    //don't delete, since they are already deleted as part of 'all particles'
    //    chargedPiCandidates.clear();
    //    cout <<" done charged Pi"<<endl;


    //      cout <<"cleaning up.."<<endl;
    for(vector<HadronPair*>::iterator it=v_hadronPairs.begin();it!=v_hadronPairs.end();it++)
      {
	delete *it;
      }


    //      cout <<"c1"<<endl;
    v_hadronPairs.clear();
    v_firstHemi.clear();
    v_secondHemi.clear();

    for(int i=0;i<chargedPiCandidates.size();i++)
      {
	delete chargedPiCandidates[i];
      }
    chargedPiCandidates.clear();

    for(int i=0;i<chargedKCandidates.size();i++)
      {
	delete chargedKCandidates[i];
      }
    chargedKCandidates.clear();



    for(int i=0;i<v_allParticles.size();i++)
      {
	//	  ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_allParticles[i]->userInfo());
	//	  cout <<"del pin"<<endl;
	//	  delete &pinf; //leads to crash in delete particle... maybe the desctructor of particle also deletes this
	//	  cout <<"aft del pin" <<endl;
	delete v_allParticles[i];
	//	  cout <<"aft del part" <<endl;

      }
    v_allParticles.clear();
  }

  // begin_run function
  void ptSpect::term()
  {
    for(int i=0;i<4;i++)
      {
	for(int j=0;j<4;j++)
	  {
	    //	    cout << "zBin " << i <<" ptBin: "<< j << " numPi0s: "<<numPi0s[i*4+j]<<" numEtas: "<< numEtas[i*4+j]<<endl;
	  }

      }

    pTreeSaver->finalize();

    //    histoD0Spect->Write();
    //    histoDStar->Write();
    //    histoPiSlowMom->Write();
    //    histoRecDStarSpectToD0Pi->Write();

    //    thetaPhiLab->Write();
    //    thetaPhiCMS->Write();


        cout <<"writing file.." <<endl;
    m_file->Write();
    cout <<"closing file.." <<endl;
    m_file->Close();
    cout <<"on to histo " <<endl;

    cout <<"numKPlus FH: " << numKPlusFH <<" KMinusFH: " << numKMinusFH <<" numPiPlusFH " << numPiPlusFH <<" piMinusFH: " << numPiMinusFH <<endl;
    cout <<"numKPlus SH: " << numKPlusSH <<" KMinusSH: " << numKMinusSH <<" numPiPlusSH " << numPiPlusSH <<" piMinusSH: " << numPiMinusSH <<endl;
    cout <<"numPi0 FH: " << numPi0FH <<" Pi0 SH: " << numPi0SH <<endl;
    cout <<"numKPiFH; " << numKPiFH << " numPiKFH: " << numPiKFH <<endl;
    cout <<"numKPiSH; " << numKPiSH << " numPiKSH: " << numPiKSH <<endl;
    cout <<"numKPQ: " << numKPQ << ", numPKQ: " << numPKQ <<endl;
    std::cout << "ptSpect's term function" << std::endl;

#ifdef SAVE_HISTOS

    TFile histoFile("myHistos.root","recreate");
    thetaPhiCMS->Write();
    m_histos.hEFlowNorm->Write();


    m_histos.hEFlowMC->Write();

    TCanvas c3;
    m_histos.hHPairMassMC->Write();


    m_histos.hHPairMassMCBin0->Write();


    m_histos.hHPairMassMCBin1->Write();

    m_histos.hHPairMassMCBin2->Write();

    m_histos.hHPairMassMCBin3->Write();

    m_histos.hEFlowFromThrust->Write();

    m_histos.hEFlowFromThrust->Scale(1/(float)m_histos.hEFlowMC->GetMaximum());
    m_histos.hEFlowMC->Scale(1/(float)m_histos.hEFlowMC->GetMaximum());

    m_histos.hEFlowMC->Add(m_histos.hEFlowFromThrust,-1);

    m_histos.hEFlowMC->Write();

    m_histos.hEFlowMC->Add(m_histos.hEFlowFromThrust,1);
    m_histos.hEFlowMC->Divide(m_histos.hEFlowFromThrust);
    m_histos.hEFlowMC->Write();

    m_histos.hEFlowNorm->Write();
    histoFile.Write();

    for(int i =0;i<m_histos.v_histos.size();i++)
      {
	m_histos.v_histos[i]->Write();
      }
    for(int i =0;i<m_histos.v_histos2D.size();i++)
      {
	m_histos.v_histos2D[i]->Write();
      }
    for(int i =0;i<m_histos.v_histosI.size();i++)
      {
	m_histos.v_histosI[i]->Write();
      }

    cout <<"numEvents; " << test<<endl;
#endif

  }


  void ptSpect::findDStar(vector<Hep3Vector>& allPB, vector<int>& allPB_Class, vector<int>& allPB_Charge)
  {
    const float D0mass=1.865;
    const float DStarMass=2.010;
    const float pionMass=0.139570;
    const float kaonMass=0.493667;

    kinematics::D0Tag=0;
    kinematics::DStarTag=0;
    kinematics::DDecay=-1;
    kinematics::DStarDecay=-1;

    if(allPB.size()!=allPB_Class.size() || allPB_Class.size()!=allPB_Charge.size())
      {
	cout<<"findstar size problem " << endl;
	exit(1);
      }
    for(int i=0;i<allPB.size();i++)
      {
	//find kaon- (should also work with k+...)
	if(allPB_Class[i]==3 && allPB_Charge[i]<0)
	  { 
	    if(allPB[i].mag()<0.1)  
	      continue;
	    double m_mass=kaonMass;
	    double EK=sqrt(m_mass*m_mass+allPB[i].mag2());

	    //find pion to form D0 (no vertex fit yet...)
	    for(int j=0;j<allPB.size();j++)
	      {
		if(allPB[j].mag()<0.1)
		  continue;

		//found pion, see if we can form D0 mass...
		if(allPB_Class[j]==2 && allPB_Charge[j]>0)
		  {
		    m_mass=pionMass;
		    float EP1=sqrt(m_mass*m_mass+allPB[j].mag2());
		    float d0candidateMass=sqrt((EK+EP1)*(EK+EP1)-(allPB[i]+allPB[j]).mag2());
		    //		    histoD0Spect->Fill(d0candidateMass);
		    if(fabs(d0candidateMass-D0mass)<0.1)
		      {
			kinematics::D0Tag=1;
			kinematics::DDecay=1;

			Hep3Vector d0candidateMom=(allPB[i]+allPB[j]);
			float d0candidateE=EP1+EK;

			for(int k=0;k<allPB.size();k++)
			  {
			    if(k==j)
			      continue;
			    if(allPB_Class[k]==2 && allPB_Charge[k]>0)
			      {
				if(allPB[k].perp()<0.1)
				  continue;
				float EP2=sqrt(m_mass*m_mass+allPB[k].mag2());
				float dStarcandidateMass=sqrt((d0candidateE+EP2)*(d0candidateE+EP2)-(d0candidateMom+allPB[k]).mag2());
				//martins cut... 
				Hep3Vector dStarCandMom=d0candidateMom+allPB[k];
				if(dStarCandMom.mag()>2.0 && dStarCandMom.mag()<4.9)
				  {
				    //				    histoDStar->Fill(dStarcandidateMass);
				    if(fabs(dStarcandidateMass-DStarMass)<0.05)
				      {
					kinematics::DStarTag=1;
					//					histoPiSlowMom->Fill(allPB[k].mag());
				      }

				  }
				
			      }

			  }

		      }
		  }


	      }

	  }
      }


  }



  //reconstruct D0 from pion, kaon, ks..
  void ptSpect::reconstructD0()
  {
    double m_D0=1.86484;
    double m_d0mass_max=m_D0+0.015;
    double m_d0mass_min=m_D0-0.015;


    //the mass range for which we do a mass/vertex constrained fit
    double m_d0mass_maxLoose=m_D0+3*0.015;
    double m_d0mass_minLoose=m_D0-3*0.015;



    //nominal mass: 1.86484, use +- 0.015, except for K-pi+pi0 (0.025)
    //D-->Kpi
    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
      {
	for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	  {
	    Particle& pion= *(*itP);
	    Particle& kaon= *(*itK);
	    if(kaon.charge()+pion.charge()!=0) continue;
	    HepLorentzVector p_d0=pion.p()+kaon.p();
	    double m=p_d0.mag();




	    //	if(!doKmVtxFit2(*(*itD),  confLevel,0))


	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
	    Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));

	    d0->relation().append(kaon);
	    d0->relation().append(pion);
	    if(m_mc)
	      {
		const Gen_hepevt &h_kaon=kaon.genHepevt();
		const Gen_hepevt &h_pion=pion.genHepevt();
		if(h_kaon && h_pion && h_kaon.mother() && h_pion.mother() && h_kaon.mother().get_ID()==h_pion.mother().get_ID()){
		  d0->relation().genHepevt(h_kaon.mother());
		}
	      } 
	    D0Candidates.push_back(d0);
	    kinematics::DDecay=1;

	  }
      }

    ////--->D to K-pi+pi0
    //uses different mass cut:
    m_d0mass_max=m_D0+0.025;
    m_d0mass_min=m_D0-0.025;
    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
      {
	for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	  {
	    for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
	      {
		Particle& pi0=*(*itPi0);
		Particle& pion= *(*itP);
		Particle& kaon= *(*itK);
		if(kaon.charge()+pion.charge()!=0) continue;
		HepLorentzVector p_d0=pion.p()+kaon.p()+pi0.p();
		double m=p_d0.mag();
		//	    cout <<"2 filling with " << m <<endl;
		if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;

		Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));

		d0->relation().append(kaon);
		d0->relation().append(pion);
		d0->relation().append(pi0);
		if(m_mc)
		  {
		    const Gen_hepevt &h_kaon=kaon.genHepevt();
		    const Gen_hepevt &h_pion=pion.genHepevt();
		    const Gen_hepevt &h_pi0=pi0.genHepevt();
		    if(h_pi0&&h_kaon && h_pion && h_kaon.mother() && h_pion.mother() && h_kaon.mother().get_ID()==h_pion.mother().get_ID() && h_pi0.mother() && h_pi0.mother()==h_kaon.mother()){
		      d0->relation().genHepevt(h_kaon.mother());
		    }
		  }
		D0Candidates.push_back(d0);
		kinematics::DDecay=2;
	      }
	  }
      }
    //and set the mass cuts back..
    m_d0mass_max=m_D0+0.015;
    m_d0mass_min=m_D0-0.015;
    ///D->K-pi+pi+pi-
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=(itP1+1);itP2!=chargedPiCandidates.end();itP2++)  
	  {
	    for(vector<Particle*>::iterator itP3=(itP2+1);itP3!=chargedPiCandidates.end();itP3++)  
	      {
		for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
		  {
		    Particle& pion1= *(*itP1);
		    Particle& pion2= *(*itP2);
		    Particle& pion3= *(*itP3);
		    Particle& kaon= *(*itK);
		    if(kaon.charge()+pion1.charge()+pion2.charge()+pion3.charge()!=0) continue;
		    HepLorentzVector p_d0=kaon.p()+pion1.p()+pion2.p()+pion3.p();
		    double m=p_d0.mag();

		    //	    cout <<"3 filling with " << m <<endl;


		    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;

		    Particle* d0 =new Particle(p_d0,Ptype(kaon.charge()<0 ? "D0" : "D0B"));

		    d0->relation().append(kaon);
		    d0->relation().append(pion1);
		    d0->relation().append(pion2);
		    d0->relation().append(pion3);
		    if(m_mc)
		      {
			const Gen_hepevt &h_kaon=kaon.genHepevt();
			const Gen_hepevt &h_pion1=pion1.genHepevt();
			const Gen_hepevt &h_pion2=pion2.genHepevt();
			const Gen_hepevt &h_pion3=pion3.genHepevt();
		    
			if(h_kaon && h_pion1 && h_pion2&&h_pion3){
			  if(h_kaon.mother() && h_pion1.mother() &&h_pion2.mother() && h_pion3.mother()){
			    if( h_kaon.mother().get_ID()==h_pion1.mother().get_ID() && h_pion1.mother()==h_pion2.mother() && h_pion2.mother()==h_pion2.mother()){
			      d0->relation().genHepevt(h_kaon.mother());
			    }
			  }
			}
		      }
		    D0Candidates.push_back(d0);
		    kinematics::DDecay=3;
		  }
	      }
	  }
      }
 
    ////D-->Ks pi+pi-
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=itP1+1;itP2!=chargedPiCandidates.end();itP2++)
	  {
	    for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
	      {
		Particle& pion1=*(*itP1);
		Particle& pion2= *(*itP2);
		Particle& Ks= *(*itKs);
		if(pion1.charge()+pion2.charge()!=0) continue;
		HepLorentzVector p_d0=pion1.p()+pion2.p()+Ks.p();
		double m=p_d0.mag();

		//	    cout <<"4 filling with " << m <<endl;


		if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));

		Particle* d0 =new Particle(p_d0,Ptype("D0"));

		d0->relation().append(pion1);
		d0->relation().append(pion2);
		d0->relation().append(Ks);
		if(m_mc)
		  {
		    const Gen_hepevt &h_pion1=pion1.genHepevt();
		    const Gen_hepevt &h_pion2=pion2.genHepevt();
		    const Gen_hepevt &h_Ks=Ks.genHepevt();
		    if(h_pion1&&h_pion2 && h_Ks&& h_pion1.mother() && h_pion2.mother() && h_Ks.mother() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_pion1.mother()==h_Ks.mother()){
		      d0->relation().genHepevt(h_pion1.mother());
		    }
		    //		printD();
		  }
		D0Candidates.push_back(d0);
		kinematics::DDecay=4;


	      }
	  }
      }

    ////D-->K+K-
    for(vector<Particle*>::iterator itK1=chargedKCandidates.begin();itK1!=chargedKCandidates.end();itK1++)
      {
	for(vector<Particle*>::iterator itK2=itK1+1;itK2!=chargedKCandidates.end();itK2++)
	  {
	    Particle& kaon1= *(*itK1);
	    Particle& kaon2= *(*itK2);
	    if(kaon1.charge()+kaon2.charge()!=0) continue;
	    HepLorentzVector p_d0=kaon1.p()+kaon2.p();
	    double m=p_d0.mag();


	    //	    cout <<"found k/k combination, filling with m: "<< m <<endl;
	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;

	    Particle* d0 =new Particle(p_d0,Ptype("D0"));

	    d0->relation().append(kaon1);
	    d0->relation().append(kaon2);
	    if(m_mc)
	      {
		const Gen_hepevt &h_kaon1=kaon1.genHepevt();
		const Gen_hepevt &h_kaon2=kaon2.genHepevt();
		if(h_kaon1 && h_kaon2 && h_kaon1.mother() && h_kaon2.mother() && h_kaon1.mother().get_ID()==h_kaon2.mother().get_ID()){
		  d0->relation().genHepevt(h_kaon1.mother());
		}
	      } 
	    D0Candidates.push_back(d0);
	    kinematics::DDecay=5;
	  }
      }

    //D-->Kspi0
    for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
      {
	for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
	  {
	    Particle& pi0= *(*itPi0);
	    Particle& Ks= *(*itKs);
	    HepLorentzVector p_d0=pi0.p()+Ks.p();
	    double m=p_d0.mag();



	    if(m>m_d0mass_max || m < m_d0mass_min ||isnan(m)) continue;

	    Particle* d0 =new Particle(p_d0,Ptype( "D0"));

	    d0->relation().append(pi0);
	    d0->relation().append(Ks);
	    if(m_mc)
	      {
		//		printD();
		const Gen_hepevt &h_pi0=pi0.genHepevt();
		const Gen_hepevt &h_Ks=Ks.genHepevt();
		if(h_pi0 && h_Ks && h_pi0.mother() && h_Ks.mother() && h_pi0.mother().get_ID()==h_Ks.mother().get_ID()){
		  d0->relation().genHepevt(h_pi0.mother());
		}
	      } 
	    D0Candidates.push_back(d0);
	    kinematics::DDecay=6;
	  }
      }
  }

  void ptSpect::reconstructChargedD()
  {
    double m_DPlus=1.86962;
    double m_dPlusmass_max=m_DPlus+0.015;
    double m_dPlusmass_min=m_DPlus-0.015;
    ////D-->KsPi
    for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
      {
	for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	  {
	    Particle& Ks= *(*itKs);
	    Particle& pion= *(*itP);
	    HepLorentzVector p_dPlus=Ks.p()+pion.p();
	    double m=p_dPlus.mag();

	    if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;

	    Particle* dPlus =new Particle(p_dPlus,Ptype(pion.charge() > 0 ? "D+" : "D-"));

	    dPlus->relation().append(Ks);
	    dPlus->relation().append(pion);
	    if(m_mc)
	      {
		//		printD();
		const Gen_hepevt &h_Ks=Ks.genHepevt();
		const Gen_hepevt &h_pion=pion.genHepevt();
		if(h_Ks && h_pion && h_Ks.mother() && h_pion.mother() && h_Ks.mother().get_ID()==h_pion.mother().get_ID()){
		  dPlus->relation().genHepevt(h_Ks.mother());
		}
	      } 
	    chargedDCandidates.push_back(dPlus);
	    kinematics::DDecay=7;
	  }
      }
    

    ///D->Kspi+pi+pi-
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=(itP1+1);itP2!=chargedPiCandidates.end();itP2++)  
	  {
	    for(vector<Particle*>::iterator itP3=(itP2+1);itP3!=chargedPiCandidates.end();itP3++)  
	      {
		for(vector<Particle*>::iterator itKs=KsCandidates.begin();itKs!=KsCandidates.end();itKs++)
		  {
		    Particle& pion1= *(*itP1);
		    Particle& pion2= *(*itP2);
		    Particle& pion3= *(*itP3);
		    Particle& Ks= *(*itKs);
		    double charge =pion1.charge()+pion2.charge()+pion3.charge();
		    if(fabs(pion1.charge()+pion2.charge()+pion3.charge())>1) continue;
		    HepLorentzVector p_dPlus=Ks.p()+pion1.p()+pion2.p()+pion3.p();
		    double m=p_dPlus.mag();


		    if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;

		    Particle* dPlus =new Particle(p_dPlus,Ptype(charge >0 ? "D+" : "D-"));

		    dPlus->relation().append(Ks);
		    dPlus->relation().append(pion1);
		    dPlus->relation().append(pion2);
		    dPlus->relation().append(pion3);
		    if(m_mc)
		      {
			//			printD();
			const Gen_hepevt &h_Ks=Ks.genHepevt();
			const Gen_hepevt &h_pion1=pion1.genHepevt();
			const Gen_hepevt &h_pion2=pion2.genHepevt();
			const Gen_hepevt &h_pion3=pion3.genHepevt();
		    
			if(h_Ks && h_pion1 && h_pion2&&h_pion3){
			  if(h_Ks.mother() && h_pion1.mother() &&h_pion2.mother() && h_pion3.mother()){
			    if( h_Ks.mother().get_ID()==h_pion1.mother().get_ID() && h_pion1.mother()==h_pion2.mother() && h_pion2.mother()==h_pion2.mother()){
			      dPlus->relation().genHepevt(h_pion1.mother());
			    }
			  }
			}
		      }
		    chargedDCandidates.push_back(dPlus);
		    kinematics::DDecay=8;
		  }
	      }
	  }
      }


    //D-->K-P+P+
    for(vector<Particle*>::iterator itP1=chargedPiCandidates.begin();itP1!=chargedPiCandidates.end();itP1++)
      {
	for(vector<Particle*>::iterator itP2=itP1+1;itP2!=chargedPiCandidates.end();itP2++)
	  {
	    for(vector<Particle*>::iterator itK=chargedKCandidates.begin();itK!=chargedKCandidates.end();itK++)
	      {
		Particle& pion1=*(*itP1);
		Particle& pion2= *(*itP2);
		Particle& kaon= *(*itK);
		
		if(pion1.charge()!=pion2.charge()) continue;
		if(kaon.charge()==pion1.charge()) continue;
		double charge = pion1.charge()+pion2.charge()+kaon.charge();
		HepLorentzVector p_dPlus=pion1.p()+pion2.p()+kaon.p();
		double m=p_dPlus.mag();


		if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));

		Particle* dPlus =new Particle(p_dPlus,Ptype(charge > 0 ? "D+": "D-"));

		dPlus->relation().append(pion1);
		dPlus->relation().append(pion2);
		dPlus->relation().append(kaon);
		if(m_mc)
		  {
		    const Gen_hepevt &h_pion1=pion1.genHepevt();
		    const Gen_hepevt &h_pion2=pion2.genHepevt();
		    const Gen_hepevt &h_kaon=kaon.genHepevt();
		    if(h_pion1&&h_pion2 && h_kaon&& h_pion1.mother() && h_pion2.mother() && h_kaon.mother() && h_pion1.mother().get_ID()==h_pion2.mother().get_ID() && h_pion1.mother()==h_kaon.mother()){
		      dPlus->relation().genHepevt(h_pion1.mother());
		    }
		  }
		chargedDCandidates.push_back(dPlus);
		kinematics::DDecay=9;
	      }
	  }
      }

    //D-->K+K-pi+
    for(vector<Particle*>::iterator itK1=chargedKCandidates.begin();itK1!=chargedKCandidates.end();itK1++)
      {
	for(vector<Particle*>::iterator itK2=itK1+1;itK2!=chargedKCandidates.end();itK2++)
	  {
	    for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	      {
		Particle& kaon1=*(*itK1);
		Particle& kaon2= *(*itK2);
		Particle& pion= *(*itP);
		
		if(kaon1.charge()==kaon2.charge()) continue;
		double charge = kaon1.charge()+kaon2.charge()+pion.charge();
		HepLorentzVector p_dPlus=kaon1.p()+kaon2.p()+pion.p();
		double m=p_dPlus.mag();



		if(m>m_dPlusmass_max || m < m_dPlusmass_min ||isnan(m)) continue;
		//		Particle* d0 =new Particle(p_d0,Ptype(k/.charge()<0 ? "D0" : "D0B"));

		Particle* dPlus =new Particle(p_dPlus,Ptype(charge > 0 ? "D+": "D-"));

		dPlus->relation().append(kaon1);
		dPlus->relation().append(kaon2);
		dPlus->relation().append(pion);
		if(m_mc)
		  {
		    const Gen_hepevt &h_kaon1=kaon1.genHepevt();
		    const Gen_hepevt &h_kaon2=kaon2.genHepevt();
		    const Gen_hepevt &h_pion=pion.genHepevt();
		    if(h_kaon1&&h_kaon2 && h_pion&& h_kaon1.mother() && h_kaon2.mother() && h_pion.mother() && h_kaon1.mother().get_ID()==h_kaon2.mother().get_ID() && h_kaon1.mother()==h_pion.mother()){
		      dPlus->relation().genHepevt(h_kaon1.mother());
		    }
		  }
		chargedDCandidates.push_back(dPlus);
		kinematics::DDecay=10;
	      }
	  }
      }



  }

  void ptSpect::reconstructDStar()
  {
    double m_DStarPlus=2.01027;
    double m_DStar0=2.00697;
    double m_D0=1.86484;
    double m_DPlus=1.86962;


    //Dmitry does not have max and min on DStar, just the difference (probably due to the uncertainty on the DStar)
    double m_dStarPlusmass_max=m_DStarPlus+10.015;
    double m_dStarPlusmass_min=m_DStarPlus-10.015;

    double m_dStar0mass_max=m_DStar0+10.015;
    double m_dStar0mass_min=m_DStar0-10.015;

    double max_massDifference=0.003;

    //    cout <<" combining " << chargedDCandidates.size() <<" charged Ds with " << pi0Candidates.size() <<" pi0s"<<endl;
    for(vector<Particle*>::iterator itD=chargedDCandidates.begin();itD!=chargedDCandidates.end();itD++)
      {
	//	break;
	for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
	  {
	    Particle& D= *(*itD);
	    Particle& pi0= *(*itPi0);

	    bool doubleUse=false;
	    //make sure that pi0 is not child of D
	    for(int i =0;i<D.nChildren();i++)
	      {
		if(D.child(i).relation().isIdenticalWith(pi0.relation()))
		  {
		    //		    cout <<"found double use 1 .." <<endl;
		    doubleUse=true;
		    break;
		  }
	      }
	    if(doubleUse)
	      continue;


	    HepLorentzVector p_dStar=D.p()+pi0.p();
	    double m=p_dStar.mag();

	    if(m>m_dStarPlusmass_max || m < m_dStarPlusmass_min ||isnan(m)) continue;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStarPlus-m_DPlus) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStarPlus-m_DPlus)) <<endl;

	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_DPlus)) > max_massDifference) continue;
	    //	    cout <<"done " <<endl;

	    Particle* dStar =new Particle(p_dStar,Ptype(D.charge()>0 ? "D*+" : "D*-"));

	    dStar->relation().append(D);
	    dStar->relation().append(pi0);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Gen_hepevt &h_pi0=pi0.genHepevt();
		if(h_D && h_pi0 && h_D.mother() && h_pi0.mother() && h_D.mother().get_ID()==h_pi0.mother().get_ID()){
		  dStar->relation().genHepevt(h_D.mother());
		}
	      } 
	    DStarCandidates.push_back(dStar);
	    kinematics::DStarDecay=1;
	  }
      }


    //    cout <<" combining " << D0Candidates.size() <<"  D0s with " << chargedPiCandidates.size() <<" charged "<<endl;
    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	for(vector<Particle*>::iterator itP=chargedPiCandidates.begin();itP!=chargedPiCandidates.end();itP++)
	  {
	    /////	    cout <<"looking at D0 - pi combination- "<<kinematics::evtNr <<endl;
	    Particle& D= *(*itD);
	    Particle& pion= *(*itP);

	    ////	    cout <<"D: "<< D.p().px()<<" y: "<< D.p().py() << " z: "<< D.p().pz()<<endl;
	    /////	    cout <<"pion: "<< pion.p().px()<<" y: "<< pion.p().py() << " z: "<< pion.p().pz()<<endl;

	    bool doubleUse=false;
	    //make sure that pion is not child of D
	    for(int i =0;i<D.nChildren();i++)
	      {
		if(D.child(i).relation().isIdenticalWith(pion.relation()))
		  {
		    //		    cout <<"found double use 2 .." <<endl;
		    doubleUse=true;
		    break;
		  }
	      }
	    if(doubleUse)
	      continue;

	    //	    cout <<"no double use" <<endl;
	    HepLorentzVector p_dStar=D.p()+pion.p();
	    double m=p_dStar.mag();
	    //	    histoRecDStarSpectToD0Pi->Fill(m);

	    ///	    cout <<"dstar cand mass: "<< m <<endl;
	    if(m>m_dStarPlusmass_max || m < m_dStarPlusmass_min ||isnan(m)) continue;
	    //	   	    cout <<"m -D: "<< m-D.p().mag() <<endl;
	    //	       cout <<"looking at dstar, mass diff: " <<(m_DStarPlus-m_D0) <<" vs : " << (m-D.p().mag());
	    //	        cout <<" gives: " << fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) <<endl;
	    //		cout <<"individual: m: "<< m <<" D mass: "<< D.p().mag() <<endl;

	    //	    cout<<"rec D: "<<D.pType().name() <<" (" << D.pType().lund()<<"), "<<" |p|: " <<D.p().rho()<< " ("<<D.p().px()<<", " << D.p().py()<< ", " << D.p().pz() <<", "<< D.p().t()<<")" <<endl;

	    //	    for(int i =0;i<D.nChildren();i++)
	    {
	      //		cout<<"rec D: "<<D.child(i).pType().name() <<" (" << D.child(i).pType().lund()<<"), "<<" |p|: " <<D.child(i).p().rho()<< " ("<<D.child(i).p().px()<<", " << D.child(i).p().py()<< ", " << D.child(i).p().pz() <<", "<< D.child(i).p().t()<<")" <<endl;
	    }

	    //	    cout <<"rec pion: "<<pion.pType().name() <<" (" << pion.pType().lund()<<"), "<<" |p|: " <<pion.p().rho()<< " ("<<pion.p().px()<<", " << pion.p().py()<< ", " << pion.p().pz() <<", "<< pion.p().t()<<")" <<endl;
	    //	    cout <<"momentum of reconstructed D: "<< D.p().px() <<", " << D.p().py()<<", " << D.p().pz()<<endl;
	    //	    cout <<"momentum of reconstructed pion: "<< pion.p().px() <<", " << pion.p().py()<<", " << pion.p().pz()<<endl;
	    //	    HepLorentzVector p_dStar=D.p()+pi0.p();

	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) <0.1)
	      {
		//			cout <<"printing D.."<<endl;
		//				printD(true);
	      }
	    //	    cout <<"checking for mass difference: "<< fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) <<endl;
	    if(fabs(m-D.p().mag()-(m_DStarPlus-m_D0)) > max_massDifference) continue;
	    ///	    cout <<"done" <<endl;

	    Particle* dStar =new Particle(p_dStar,Ptype(pion.charge()>0 ? "D*+" : "D*-"));

	    dStar->relation().append(D);
	    dStar->relation().append(pion);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Gen_hepevt &h_pion=pion.genHepevt();
		if(h_D && h_pion && h_D.mother() && h_pion.mother() && h_D.mother().get_ID()==h_pion.mother().get_ID()){
		  dStar->relation().genHepevt(h_D.mother());
		}
	      } 
	    DStarCandidates.push_back(dStar);


	    kinematics::DStarDecay=2;
	  }
      }

    ///D*0 from D0pi0
    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	//	break;
	for(vector<Particle*>::iterator itPi0=pi0Candidates.begin();itPi0!=pi0Candidates.end();itPi0++)
	  {
	    Particle& D= *(*itD);
	    Particle& pi0= *(*itPi0);

	    bool doubleUse=false;
	    //make sure that pi0 is not child of D
	    for(int i =0;i<D.nChildren();i++)
	      {
		if(D.child(i).relation().isIdenticalWith(pi0.relation()))
		  {
		    //		    cout <<"found double use 3 .." <<endl;
		    doubleUse=true;
		    break;
		  }
	      }
	    if(doubleUse)
	      continue;


	    HepLorentzVector p_dStar=D.p()+pi0.p();
	    double m=p_dStar.mag();


	    if(m>m_dStar0mass_max || m < m_dStar0mass_min ||isnan(m)) continue;
	    //	    cout <<"looking at dstar, mass diff: " <<(m_DStar0-m_D0) <<" vs : " << (m-D.p().mag());
	    //	    cout <<" gives: " << fabs(m-D.p().mag()-(m_DStar0-m_D0)) <<endl;
	    if(fabs(m-D.p().mag()-(m_DStar0-m_D0)) > max_massDifference) continue;
	    //	    cout <<" done " <<endl;

	    Particle* dStar =new Particle(p_dStar,Ptype("D*0"));

	    dStar->relation().append(D);
	    dStar->relation().append(pi0);
	    if(m_mc)
	      {
		const Gen_hepevt &h_D=D.genHepevt();
		const Gen_hepevt &h_pi0=pi0.genHepevt();
		if(h_D && h_pi0 && h_D.mother() && h_pi0.mother() && h_D.mother().get_ID()==h_pi0.mother().get_ID()){
		  dStar->relation().genHepevt(h_D.mother());
		}
	      } 
	    DStarCandidates.push_back(dStar);
	    kinematics::DStarDecay=3;
	  }
      }





  }


  //this seems to have some impact...
  unsigned ptSpect::doKmFit(Particle &p, double& confLevel, int debug, double mass)
  {
    //    return true;
    kmassfitter km;
    if(mass!=0)
      {
	//	cout <<" using mass; "<< mass <<endl;
	km.invariantMass(mass);
      }
    //    km.invariantMass(mass==0 ? p.pType().mass(): mass);
    else
      km.invariantMass(mass==0 ? p.pType().mass(): mass);
    for(unsigned j=0;j<p.relation().nChildren();j++)
      {
	Particle child=p.relation().child(j);
	km.addTrack(child.momentum().p(),child.momentum().x(),child.momentum().dpx(),child.pType().charge(),child.pType().mass());
      }
    km.notDecayPoint();
    unsigned err = km.fit();
    if(err){
      //           cout <<"Err in kmassvertexfitter: "<< err <<endl;
      return 0;
    }
    //    else{cout <<"fit was ok.." <<endl;}
    confLevel=km.cl();
    return makeMother(km,p);
  }

  //dmitries code..
  unsigned ptSpect::doKmVtxFit2(Particle &p, double& confLevel, int debug, double mass)
  {
    kmassvertexfitter kmvfitter;
    kmvfitter.invariantMass(mass==0 ? p.pType().mass() : mass);
    for(unsigned i=0; i<p.nChildren(); ++i)
      addTrack2fit(kmvfitter, p.child(i));

    if(p.nChildren()>2)
      {
	
      }
    //this is in the example, but probably old interface?
    //    kmvfitter.vertex(IpProfile::position());
    //    kmvfitter.errVertex(IpProfile::position_err());
    kmvfitter.initialVertex(IpProfile::position());
    //no error  
    if(!kmvfitter.fit()) {
      makeMother(kmvfitter, p);
      p.momentum().vertex(kmvfitter.vertex(),kmvfitter.errVertex());
      confLevel=kmvfitter.cl();
      return true;
    }
    return false;
  }
  //this seems to have some impact...
  unsigned ptSpect::doKmVtxFit(Particle &p, double& confLevel, int debug)
  {
    //first get vertex:
    kvertexfitter vtxFit;
    vtxFit.initialVertex(p.x());
    for(unsigned j=0;j<p.relation().nChildren();j++)
      {
	//	Particle child=p.relation().child(j);
	//	vtxFit.addTrackToFit(child.momentum().p(),child.momentum().x(),child.momentum().dpx(),child.pType().charge(),child.pType().mass());
      }

    //    return true;
    kmassvertexfitter km;
    km.initialVertex(p.x());
    km.invariantMass(p.pType().mass());
    for(unsigned j=0;j<p.relation().nChildren();j++)
      {
	Particle child=p.relation().child(j);
	km.addTrack(child.momentum().p(),child.momentum().x(),child.momentum().dpx(),child.pType().charge(),child.pType().mass());
      }
    //    km.notDecayPoint();
    unsigned err = km.fit();
    if(err){
      //           cout <<"Err in kmassvertexfitter: "<< err <<endl;
      return 0;
    }
    //    else{cout <<"fit was ok.." <<endl;}
    confLevel=km.cl();
    return makeMother(km,p);
  }


  genhep_vec* ptSpect::getDaughters(const Gen_hepevt &mother)
  {
    /* get a vector with all daughters as Gen_hepevt* */

    Gen_hepevt_Manager& gen_hepevt_mgr = Gen_hepevt_Manager::get_manager();

    int n_children = mother.daLast() - mother.daFirst() + 1;

    genhep_vec *children = new genhep_vec();

    for(int i=0; i<n_children; i++) {

      Panther_ID ID0(mother.daFirst()+i);
      if(ID0==0)
	{
	  if(PRINT)
	    cout <<"wrong!!!" <<endl;
	  break;
	}
      Gen_hepevt& temp = gen_hepevt_mgr(ID0);

      if (temp) 
	{
	  children->push_back(&temp);
	}
    }

    return children;
  }


  bool ptSpect::recursivePrint(const Gen_hepevt gen_it, string s)
  {
    genhep_vec* daughters=getDaughters(gen_it);
    int lund=fabs(gen_it.idhep());

    if(lund==911|| lund>9000000)
      return false;
    Particle p(gen_it);
    if(lund== 100423 || lund ==100421 ||lund==100411 || lund==100413 )
      cout <<s << lund << endl;
    else
      cout <<s<<p.pType().name() <<" (" << p.pType().lund()<<"), "<<" |p|: " <<p.p().rho()<< " ("<<p.p().px()<<", " << p.p().py()<< ", " << p.p().pz() <<", "<< p.p().t()<<")" <<endl;
    if(daughters->size()<=0)
      {
	//	cout <<"has " << p.nChildren()<<" children " <<endl;
	delete daughters;
	return true;
      }
    //don't print eventual decay products of kaons etc...
    if(lund==211 || lund==321 || lund==13 || lund==111)
      return true;
    
    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
      {
	recursivePrint(**it,s+("-->"));
      }
    delete daughters;
    //    cout <<s<<endl;
  }




  void ptSpect::printD(bool star)
  {

    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
    int neutralD=411;
    int charged=421;
    if(star)
      {
	neutralD=413;
	charged=423;
      }
    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {

	if(fabs(gen_it->idhep())==neutralD || fabs(gen_it->idhep())==charged)
	  {
	    recursivePrint(*gen_it,"");
	  }
      }
  }


  bool ptSpect::recD0MC()
  {


    return true;
  }
  bool ptSpect::recDStarMC()
  {
    kinematics::DStarDecayMC=-1;
    kinematics::DDecayMC=-1;
    Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();

    for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
      {
	if(fabs(gen_it->idhep())==413 || fabs(gen_it->idhep())==423)
	  {
	    int numK=0;
	    int numPi=0;
	    int numDaughters=0;
	    int numGDaughters=0;

	    genhep_vec* daughters=getDaughters(*gen_it);
	    numDaughters=daughters->size();
	    if(numDaughters<=0)
	      {
		delete daughters;
		return false;
	      }

	    for(genhep_vec::iterator it=daughters->begin();it!=daughters->end();it++)
	      {		    
		if(fabs((*it)->idhep())==211)
		  {
		    numPi++;
		  }
		if(fabs((*it)->idhep())==321)
		  {
		    numK++;
		  }
		if((fabs((*it)->idhep())==411) || (fabs((*it)->idhep())==421))
		  {
		    genhep_vec* grandDaughters=getDaughters(**it);
		    numGDaughters=grandDaughters->size();
		    for(genhep_vec::iterator it2=grandDaughters->begin();it2!=grandDaughters->end();it2++)
		      {
			if(fabs((*it2)->idhep())==211)
			  {
			    numPi++;
			  }
			if(fabs((*it2)->idhep())==321)
			  {
			    numK++;
			  }
		      }
		    delete grandDaughters;
		       
		  }

	      }
	    //delete daughters of D*
	    delete daughters;	
	    //don't print eventual decay products of kaons etc...
		
	    if(numK==1 && numPi==2 && numDaughters==2 && numGDaughters==2)
	      {
		//ddecay ==1 is pi/k ,dstardecay 2 is D0/pi
		//		    recursivePrint(*gen_it,"");
		kinematics::DDecayMC=1;
		kinematics::DStarDecayMC=2;
		//		    cout <<"found d in mc "<< kinematics::evtNr<<endl;
		//		    printD(true);
		return true;
	      }
	    else
	      {
		//start new
		numDaughters=0;
		numGDaughters=0;
		numK=0;
		numPi=0;
	      }

	  }

      }
    return false;
  }

  int ptSpect::setHadronPIDProbs(ParticleInfo* info, float mom)
  {
    //    cout <<"pb+1 " << pb+1 <<" phb+1: "<< thb+1 <<endl;
    //    cout <<"getting mombin for mom: " << mom <<endl;
    if(kinematics::evtNr==DEBUG_EVENT)
      {
	cout <<"setting pid probs for p: "<< mom<<endl;
	cout <<"id as : "<< info->idAs<<endl;
      }
    int momBin=getBin(plabb,pb+1,mom);
    //    cout <<"getting mombin for theta: " << info->labTheta <<endl;
    int thetaBin=getBin(costhetab,thb+1,cos(info->labTheta));

    //lowest bin is the border..
    thetaBin--;
    momBin--;
    //     cout <<" plab: " << mom<<" bin: " << momBin <<" cosTheta: "<< cos(info->labTheta) <<"/bin " << thetaBin << " charge: " << info->charge <<endl;
    //probability for each mass hyposthesis
    //    cout <<"mom: " << mom << " bin: " << momBin <<" that: " << info->labTheta <<" bin: "<< thetaBin <<endl;
    int idAs=info->idAs;
    //    cout <<" set hadron pid probs idAs: " << idAs <<" thetaBin: "<< thetaBin << " momBin: "<< momBin <<endl;

    int momBinUp=momBin+1;
    int momBinDown=momBin-1;
    int thetaBinUp=thetaBin+1;
    int thetaBinDown=thetaBin-1;
    if(momBinDown<0)
      momBinDown=0;
    if(thetaBinDown<0)
      thetaBinDown=0;
    if(momBinUp>=numMomBins)
      momBinUp=numMomBins-1;
    if(thetaBinUp>=numThetaBins)
      thetaBinUp=numThetaBins-1;


    //    cout <<" looking at pi/pi hypo: "<<   pidMatrixPositive[momBin][thetaBin][2][2] <<" or: "<<  pidMatrixNegative[momBin][thetaBin][2][2]<<endl;


    //    cout <<" mom/theta up looking at pi/pi hypo: "<<   pidMatrixPositive[momBinUp][thetaBinUp][2][2] <<" or: "<<  pidMatrixNegative[momBinUp][thetaBinUp][2][2]<<endl;

    //    cout <<" mom/theta down looking at pi/pi hypo: "<<   pidMatrixPositive[momBinDown][thetaBinDown][2][2] <<" or: "<<  pidMatrixNegative[momBinDown][thetaBinDown][2][2]<<endl;




    for(int hypo=0;hypo<5;hypo++)
      {
	//default
	info->pidProbabilities[hypo]=0.0;
	info->pidProbabilities2[hypo]=0.0;
	info->pidUncert[hypo]=0.0;
	info->pidUncert2[hypo]=0.0;
	//	cout <<"momBin: " << momBin <<" thetaBin: "<< thetaBin <<" hypo: " << hypo <<  " idAs: "<< idAs << " charge: "<< info->charge <<endl;
	if(info->charge > 0)
	  {
	    //we don't have matrices for everything
	    if((thetaBin>=0 && momBin>=0) && mom >0.5){


	      if(kinematics::evtNr==DEBUG_EVENT)
		{
		  cout <<"theta and mom bin in range " <<endl;
		  cout <<"setting prop to : "<< pidMatrixPositive[momBin][thetaBin][hypo][idAs] << " for hypo " << hypo <<endl;
		}
	    
	      info->pidProbabilities[hypo]=pidMatrixPositive[momBin][thetaBin][hypo][idAs];
	      info->pidProbabilities2[hypo]=pidMatrixPositive2[momBin][thetaBin][hypo][idAs];
	      info->pidUncert[hypo]=pidUncertPositive[momBin][thetaBin][hypo][idAs];
	      info->pidUncert2[hypo]=pidUncertPositive2[momBin][thetaBin][hypo][idAs];
	      //	           cout <<"hypothesis: " << hypo <<" weight: "<< info->pidProbabilities[hypo] <<endl;
	      if(fabs(info->pidProbabilities[hypo])>10000)
		{
		  return -1;
		  //info->pidProbabilities[
		  //		return -1;
		}
	      //	      cout <<"pid probability for hypo " << hypo<< " is: "<< pidMatrixPositive[momBin][thetaBin][hypo][idAs] <<endl;
	    }
	    else
	      {
		if(idAs==hypo)
		  {
#ifdef noPID
		  
		    //even with noPID reject events outside acceptance
		    info->pidProbabilities[hypo]=0.0;
		    info->pidProbabilities2[hypo]=0.0;
#else
		    info->pidProbabilities[hypo]=0.0;
		    info->pidProbabilities2[hypo]=0.0;
#endif
		  }
		else
		  {
		    info->pidProbabilities[hypo]=0.0;
		    info->pidProbabilities2[hypo]=0.0;
		  }
	      }

	  }
	else
	  {
	    if((thetaBin>=0 && momBin>=0) && mom >0.5)
	      {
		//	      cout <<"pid neg probability for hypo " << hypo<< " is: "<< pidMatrixPositive[momBin][thetaBin][hypo][idAs] <<endl;
		info->pidProbabilities[hypo]=pidMatrixNegative[momBin][thetaBin][hypo][idAs];	
		info->pidProbabilities2[hypo]=pidMatrixNegative2[momBin][thetaBin][hypo][idAs];	
		info->pidUncert[hypo]=pidUncertNegative[momBin][thetaBin][hypo][idAs];
		info->pidUncert2[hypo]=pidUncertNegative2[momBin][thetaBin][hypo][idAs];
		//	      cout <<"hypothesis: " << hypo <<" weight: "<< info->pidProbabilities[hypo] <<endl;    
		if(fabs(info->pidProbabilities[hypo])>10000)
		  return -1;
	      }
	    else
	      {
		if(idAs==hypo)
		  {
#ifdef noPID

		    //outside acceptance
		    info->pidProbabilities[hypo]=0.0;
		    info->pidProbabilities2[hypo]=0.0;
#else
		    info->pidProbabilities[hypo]=0.0;
		    info->pidProbabilities2[hypo]=0.0;
#endif
		  }
		else
		  {
		    info->pidProbabilities[hypo]=0.0;
		    info->pidProbabilities2[hypo]=0.0;
		  }
	      }
	  }
      }
    info->p_Pi=info->pidProbabilities[pionIdx];
    info->p_K=info->pidProbabilities[kaonIdx];
    info->p_p=info->pidProbabilities[protonIdx];


    info->p_e=info->pidProbabilities[electronIdx];
    info->p_mu=info->pidProbabilities[muonIdx];


    info->p_Pi2=info->pidProbabilities2[pionIdx];
    info->p_K2=info->pidProbabilities2[kaonIdx];
    info->p_p2=info->pidProbabilities2[protonIdx];


    info->p_e2=info->pidProbabilities2[electronIdx];
    info->p_mu2=info->pidProbabilities2[muonIdx];


    info->ep_Pi=info->pidUncert[pionIdx];
    info->ep_K=info->pidUncert[kaonIdx];
    info->ep_p=info->pidUncert[protonIdx];


    info->ep_e=info->pidUncert[electronIdx];
    info->ep_mu=info->pidUncert[muonIdx];

    info->ep_Pi2=info->pidUncert2[pionIdx];
    info->ep_K2=info->pidUncert2[kaonIdx];
    info->ep_p2=info->pidUncert2[protonIdx];


    info->ep_e2=info->pidUncert2[electronIdx];
    info->ep_mu2=info->pidUncert2[muonIdx];





    return 1;
    //    cout <<"p_Pi " << info->p_Pi <<" p_K: " << info->p_K <<" p_p: "<< info->p_p << " p_e " << info->p_e <<" p_mu : " << info->p_mu <<endl;
  }

  //code from FRancesca
  void  ptSpect::loadPIDMatrix()
  {
    //   TFile* fpid = new TFile("newpid.root","read");

    TFile* fpid = new TFile("~vossen/myProjects/ptSpect/invertedpidmatrices_setb061810I_inv2_realdataalways.root","read");
    TFile* fpid2 = new TFile("~vossen/myProjects/ptSpect/invertedpidmatrices_setb061810I_inv1_MConlyatlooseends.root","read");

    char matrix_name[300];
    char uncert_minus_name[300];
    char uncert_plus_name[300];
    if (fpid->IsZombie() || fpid2->IsZombie()) {
      printf("File code.root does not exist.\n");
      return;
    }
  
    // reading all matrices together, as the code crashes if I try to open the file too many times...
    for (Int_t u = 0; u < pb; u++)
      for (Int_t v = 0; v < thb; v++)
	for (Int_t w = 0; w < 2; w++){
	  //w is charge ( 0  negative, 1 positive);
	  sprintf(matrix_name,"invanalyticmatrix_u%d_v%d_w%d",u,v,w); // invanalyticmatrix_u16_v0_w0 does not exists!!! 
	  sprintf(uncert_plus_name,"invuncertrejplmatrix_u%d_v%d_w%d",u,v,w); 
	  sprintf(uncert_minus_name,"invuncertrejmimatrix_u%d_v%d_w%d",u,v,w); 
	  //u16 starts from v3 (invanalyticmatrix_u16_v3_w0)
	  // corresponds to z > 1 (smearing?) additional uncertainties would be 1-2
	  // up to 9 July 2015! int k= u*8*2+v*2+w; using 8 instead of 9!
	  //    int k= u*thb*2+v*2+w;
	  if(u<16 || v>2 ){
	    //       cout <<" u:  " << u << " v: " << v << " w: " << w <<endl;
	    TMatrixD mat = *(TMatrixD*)fpid->Get(matrix_name);
	    TMatrixD mat2 = *(TMatrixD*)fpid2->Get(matrix_name);
	    TMatrixD matUncertPos = *(TMatrixD*)fpid->Get(uncert_plus_name);
	    TMatrixD matUncertNeg = *(TMatrixD*)fpid->Get(uncert_minus_name);

	    TMatrixD matUncertPos2 = *(TMatrixD*)fpid2->Get(uncert_plus_name);
	    TMatrixD matUncertNeg2 = *(TMatrixD*)fpid2->Get(uncert_minus_name);

	    for (Int_t i = 0; i <= 4; i++)
	      for (Int_t j = 0; j <= 4; j++){

		//	   cout <<"mat: " << mat(i,j) <<" mat2: "<< mat2(i,j) <<endl;
		//	   cout <<"pos uncert: " << matUncertPos(i,j) <<" neg: "<< matUncertNeg(i,j) <<endl;
		//	   	   cout <<"pos uncert2: " << matUncertPos2(i,j) <<" neg: "<< matUncertNeg2(i,j) <<endl;
		float symUncert=matUncertPos(i,j);
		float symUncert2=matUncertPos2(i,j);
		if(isnan(matUncertPos(i,j)))
		  {
		    cout <<" matuncert pos nan"<<endl;
		    //to force using the neg uncert
		    symUncert=0;
		  }
		if(isnan(matUncertNeg(i,j)))
		  {
		    cout <<" matuncert neg nan"<<endl;
		  }
		else
		  {
		    if(fabs(matUncertNeg(i,j))>symUncert)
		      symUncert=fabs(matUncertNeg(i,j));
		  }
		if(isnan(matUncertNeg2(i,j)))
		  {
		    cout <<" matuncert neg2 nan"<<endl;
		  }
		else
		  {
		    if(fabs(matUncertNeg2(i,j))>symUncert2)
		      symUncert2=fabs(matUncertNeg2(i,j));
		  }

		//not sure if that makes sense here since we do the averaging later on again
		//		if(fabs(symUncert2)<100.0)
		//		  symUncert=sqrt((symUncert*symUncert+symUncert2*symUncert2))/2;
		if(isnan(symUncert))
		  cout<<"resulting uncertainty still nan!!!!!!!!"<<endl<<endl;


		//now we got the symmetric uncertainties

		//at this point symUncert is the data uncertainty if there is a problem with MC, otherwise the geometric mean
		//

		if(w==1)
		  {
#ifdef noPID
		    if(i==j)
		      {
			pidMatrixPositive[u][v][i][j]=1;
			pidMatrixPositive2[u][v][i][j]=1;
			pidUncertPositive[u][v][i][j]=0.0;
			pidUncertPositive2[u][v][i][j]=0.0;
		      }
		    else
		      {
			pidMatrixPositive[u][v][i][j]=0;
			pidMatrixPositive2[u][v][i][j]=0;
			pidUncertPositive[u][v][i][j]=0.0;
			pidUncertPositive2[u][v][i][j]=0.0;
		      }
#else
		    pidMatrixPositive[u][v][i][j]=mat(i,j);
		    pidMatrixPositive2[u][v][i][j]=mat2(i,j);
		    //			    cout <<"loading " << mat(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" positive " <<endl;
		    //		    cout <<"loading 2" << mat2(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" positive " <<endl;

		    if(fabs(pidMatrixPositive2[u][v][i][j])>100.0 || fabs(symUncert2)>100.0)
		      {
			symUncert2=symUncert;
			pidMatrixPositive2[u][v][i][j]=pidMatrixPositive[u][v][i][j];
			//the uncert asseignment doesnt' matter, since we use the average between the two anyways
			//			symUncert2=symUncert;
		      }

		    if(isnan(pidMatrixPositive2[u][v][i][j]) || isnan(matUncertPos2(i,j)) || isnan(matUncertNeg2(i,j)))
		      {
			pidMatrixPositive2[u][v][i][j]=pidMatrixPositive[u][v][i][j];
			symUncert2=symUncert;
		      }
		    //		    symUncert=sqrt(symUncert*symUncert+symUncert2*symUncert2)/2;
		    pidUncertPositive[u][v][i][j]=symUncert;
		    pidUncertPositive2[u][v][i][j]=symUncert2;
		    if(isnan(symUncert))
		      {
			cout <<"poblem with pos symUncert!!!" <<endl;
		      }
		    if(symUncert>100.0)
		      {
			cout <<" uncert > 100 " <<endl;
		      }	   

#endif
		    //	     cout <<"loading positive "<< i <<", " << j << " "  << mat(i,j) <<endl;




		  }
		else
		  {
#ifdef noPID
		    if(i==j)
		      {
			pidMatrixNegative[u][v][i][j]=1;
			pidMatrixNegative2[u][v][i][j]=1;
			pidUncertNegative[u][v][i][j]=0.0;
			pidUncertNegative2[u][v][i][j]=0.0;
		      }
		    else
		      {
			pidMatrixNegative[u][v][i][j]=0;
			pidMatrixNegative2[u][v][i][j]=0;
			pidUncertNegative[u][v][i][j]=0.0;
			pidUncertNegative2[u][v][i][j]=0.0;
		      }
#else
		    //		    cout <<"loading " << mat(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" negative " <<endl;
		    //		    cout <<"symUncert: "<< symUncert <<endl;
		    //		    cout <<"loading 2 " << mat2(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" negative " <<endl;
		    pidMatrixNegative[u][v][i][j]=mat(i,j);
		    pidMatrixNegative2[u][v][i][j]=mat2(i,j);
		    if(fabs(pidMatrixNegative2[u][v][i][j])>100.0 || fabs(symUncert2)>100.0) 
		      {
			pidMatrixNegative2[u][v][i][j]=pidMatrixNegative[u][v][i][j];
			//			symUncert2=symUncert;
		      }
		    if(isnan(pidMatrixNegative2[u][v][i][j]) || isnan(matUncertPos2(i,j)) || isnan(matUncertNeg2(i,j)))
		      {
			pidMatrixNegative2[u][v][i][j]=pidMatrixNegative[u][v][i][j];
			symUncert2=symUncert;
		      }

		    //		    symUncert=sqrt(symUncert*symUncert+symUncert2*symUncert2)/2;

		    if(isnan(symUncert))
		      cout <<"problem with neg symUncert!!!" <<endl;
		    pidUncertNegative[u][v][i][j]=symUncert;
		    pidUncertNegative2[u][v][i][j]=symUncert2;

		  
#endif

		    //	     cout <<"loading negative " << i << ", " << j << " "  << mat(i,j) <<endl;




		  }

		//        matrix[k][i][j] = mat(i,j);
		//if(k==72) cout << i << " " << j <<" " <<  mat(i,j) << endl;
	      }
	  }
	}
    fpid->Close();
    cout << std::fixed;
    cout << setprecision(4);
    for(Int_t w=0;w<2;w++)
      {
	for (Int_t u = 0; u < pb; u++)
	  {

	    for (Int_t v = 0; v < thb; v++)
	      {
		cout <<"momentum bin " << u;
		cout <<", cosTheta bin " <<v<<endl;
		for (Int_t i = 0; i <= 4; i++)
		  {
		    for (Int_t j = 0; j <= 4; j++)
		      {
			if(w==0)
			  cout <<"data: " << pidMatrixNegative[u][v][i][j] <<" mc: " << pidMatrixNegative2[u][v][i][j] <<" ";
			else
			  cout <<"data: " <<pidMatrixPositive[u][v][i][j]<< " mc: "<< pidMatrixPositive2[u][v][i][j] <<" ";
			//			pidUncertNegative[u][v][i][j];
		      }
		    cout <<endl;
		  }
		cout <<"symmetrized uncertainties  "<<endl;
		for (Int_t i = 0; i <= 4; i++)
		  {
		    for (Int_t j = 0; j <= 4; j++)
		      {
			if(w==0)
			  cout <<"data: " << pidUncertNegative[u][v][i][j] <<" mc: "<< pidUncertNegative2[u][v][i][j]<< " " ;
			else
			  cout <<"data: " << pidUncertPositive[u][v][i][j] <<" mc: "<< pidUncertPositive2[u][v][i][j]<<" ";
			//			pidUncertNegative[u][v][i][j];
		      }
		    cout <<endl;
		  }

		cout <<endl;
	      }
	  }
      }

  }
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
