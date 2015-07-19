#include <iomanip>
#include "ptSpect/mc.h"  //one central place to put the define mc
#include "event/BelleEvent.h"
#include "particle/Particle.h"
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

//#define SAVE_HISTOS
#define XCHECK


#define D0Mass 1.865
#define D0Width 0.6
#define D0Lower  D0Mass-D0Width
#define D0Upper  D0Mass+D0Width

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

//#define DEBUG_EVENT 287880//please no output
#define DEBUG_EVENT 1//please no output
//#define DEBUG_EVENT 15859
#define DEBUG_EVENT2 -287880
#include "ptSpect/AnaConsts.h"
#include "ptSpect/HadronPair.h"
#include "ptSpect/ptSpect.h"
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
  ptSpect::ptSpect():smpl_(12345.),cPiPlus("PI+"), cPiNeg("PI-"),cPiZero("PI0"),cKPlus("K+"),cKNeg("K-")
  {
    strcpy(rFileName,"notInitialized.root");
    test=0;
    for(int i=0;i<4;i++)
      {
	//	zVals[i]=0;
      }

    histoD0Spect=new TH1D("d0spect","d0spect",1000,0,3.0);
    histoDStar=new TH1D("dStarspect","dStarspect",300,1.8,5.0);
    histoPiSlowMom=new TH1D("piSlow","piSlow",100,0,3.0);

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
    gROOT->SetStyle("Plain");
    m_file=new TFile(rFileName,"recreate");
    thetaPhiLab=new TH2D("thetaPhiLab","thetaPhiLab",100,0,6.3,100,-3.15,3.15);
    thetaPhiCMS=new TH2D("thetaPhiCMS","thetaPhiCMS",100,0,6.3,100,-3.15,3.15);

    jetThrustDiff=new TH1D("jetThrustDiff","jetThrustDiff",100,0,1.0);
    jetJetDiff=new TH1D("jetJetDiff","jetJetDiff",100,-1.0,1.0);
    jetEnergy=new TH1D("jetEnergy","jetEnergy",100,0,6);
    numJets=new TH1D("numJets","numJets",10,0,10);

    numPartInJet=new TH1D("numPartInJet","numPartInJet",20,0,20);
    partEnergyInJet=new TH1D("partEnergyInJet","partEnergyInJet",100,0,6);


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
    m_histos.setFilenameStart(rFileName);

    srand(time(NULL));

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
    IpProfile::begin_run();
    eid::init_data();

    BeamEnergy::begin_run();
    double eler=BeamEnergy::E_LER();
    double eher=BeamEnergy::E_HER();

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
    kinematics::firstElectronCM.boost(kinematics::CMBoost);
    kinematics::secondElectronCM.boost(kinematics::CMBoost);
    //    cout <<"end beginrun " <<endl;
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
    pTreeSaver->setDebugHistos(&m_histos);
    pTreeSaver->addArrayF("z1");
    pTreeSaver->addArrayF("z2");

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
    pTreeSaver->addArrayF("HadDiffTheta_mc");
    pTreeSaver->addArrayF("HadDiffPhi_mc");

#endif
    //!!!! this charge type is in principle worthless, look at the charges of the hadrons!!
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
    pTreeSaver->addFieldF("fluffiness1");
    pTreeSaver->addFieldF("fluffiness2");
    pTreeSaver->addFieldF("jetE1");
    pTreeSaver->addFieldF("jetE2");

    pTreeSaver->addFieldF("jet1Phi");
    pTreeSaver->addFieldF("jet2Phi");
    pTreeSaver->addFieldF("jet1Theta");
    pTreeSaver->addFieldF("jet2Theta");

#ifdef MC
    pTreeSaver->addFieldF("Thrust_mc");
    pTreeSaver->addFieldF("E_miss_mc");
    pTreeSaver->addFieldF("thrustTheta_mc"); //angle between thrust rec thrust mc
    //  pTreeSaver->addFieldF("thrustPhi_mc");  //<--- not needed because we have phi of orginal and diff phi...
    pTreeSaver->addFieldF("diffThetaThrust_mc");
    pTreeSaver->addFieldF("diffPhiThrust_mc");

    pTreeSaver->addFieldF("fluffiness1_mc");
    pTreeSaver->addFieldF("fluffiness2_mc");
    pTreeSaver->addFieldF("jetE1_mc");
    pTreeSaver->addFieldF("jetE2_mc");

    pTreeSaver->addFieldF("jet1Phi_mc");
    pTreeSaver->addFieldF("jet2Phi_mc");
    pTreeSaver->addFieldF("jet1Theta_mc");
    pTreeSaver->addFieldF("jet2Theta_mc");


    pTreeSaver->addFieldF("VP_Energy"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PX"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PY"); //energy of virtual photon
    pTreeSaver->addFieldF("VP_PZ"); //energy of virtual photon
    pTreeSaver->addFieldF("quarkAngle");
    pTreeSaver->addFieldF("thetaEThrust_mc");
    pTreeSaver->addFieldF("ISRPhotonEnergy");
    pTreeSaver->addFieldI("numQuarks");


    //the fields for the mc w/o acceptance
#endif
    pTreeSaver->addFieldI("runNr");
    pTreeSaver->addFieldI("evtNr");

    pTreeSaver->addFieldI("jetNumPart1");
    pTreeSaver->addFieldI("jetNumPart2");

    pTreeSaver->addFieldI("D0Tag");
    pTreeSaver->addFieldI("DStarTag");



    //doesn't work anymore with arrays...
    //  pTreeSaver->createNTuple(tm);
#ifdef MC

    //  pTreeSaver->addArrayPi0AsymmetryF("realPi0_gammaAsymmetry");
#endif
  }

  // event function
  void ptSpect::event(BelleEvent* evptr, int* status)
  {
    vector<float> v_drH1;
    vector<float> v_drH2;
    if(!validRun)
      return;
    int evtNr;
    int runNr;
    /////for xcheck

    evtNr=Belle_event_Manager::get_manager().begin()->EvtNo();
    runNr=Belle_event_Manager::get_manager().begin()->RunNo();

    kinematics::runNr=runNr;
    kinematics::evtNr=evtNr;

    //    cout <<"--> run = " << runNr <<" evtNr = "  <<evtNr <<endl;
    //    cout <<" --> exp = " << " run = " << runNr << " event = " << evtNr <<endl;
    if(!IpProfile::usable())
      return;
    if(!(test%1000))
      {
	//	      cout << "evt " <<test <<endl;
	//	        cout << "nr " <<evtNr <<endl;
      }
    test++;

    //#ifndef MC
    if(!goodHadronB())
      return;
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
    vector<int> allPB_particleClass;
    vector<int> allPB_particleCharge;

    vector<float> nonBoostedE;
    vector<Hep3Vector> allParticlesNonBoosted;


    vector<PseudoJet> fjParticles;

    vector<float> allPB_E; //energy for the three vec above...
    Mdst_charged_Manager& mdst_chr_Mgr=Mdst_charged_Manager::get_manager();
    //  Mdst_klong_Manager& mdst_klong_Mgr=Mdst_klong_Manager::get_manager();
    Mdst_trk_Manager& mdst_trk_Mgr=Mdst_trk_Manager::get_manager();
    atc_pid selKPi(3,1,5,3,2);  //K/pi separation
    atc_pid selPK(3,1,5,4,3); //proton kaon separation
    atc_pid selPiP(3,1,5,2,3); //pion proton separation
    atc_pid selKP(3,1,5,3,4);

    int itmp=0;
    int iChTrks=0;//num of charged Tracks
    int iChTrkPos=0;
    int iChTrkNeg=0;
    //   cout <<"chargedTracks: " << mdst_chr_Mgr.size() <<" num Trk: " << mdst_trk_Mgr.size() <<" klong: " << mdst_klong_Mgr.size() <<endl;

    //    cout <<"there are " << mdst_chr_Mgr.size() << " charge tracks in mdst_chr " <<endl;
    for(Mdst_charged_Manager::iterator chr_it=mdst_chr_Mgr.begin();chr_it!=mdst_chr_Mgr.end();chr_it++)
      {
	if(!enoughSVDHits(chr_it))
	  continue;
	double m_mass=m_pi;
	int massHyp=2;
	double m_theta=0;
	double m_phi=0;
	double m_qt=0;
	double m_z=0;
	strcpy(ptypName,"unknown");
	double charge=(*chr_it).charge();

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
	double e_id=sel_e.prob(0,-1,0);

	if(DEBUG_EVENT==evtNr)
	  {
	    //	    cout <<"pid kpi: " << atcKPiAlt <<" pid KP: " << atcKPAlt << " e_id: " << e_id << " mu_id: " << mu_id <<endl;
	  }


	bool isLepton=false;
	bool isPionKaon=false;
	bool isProton=false;
      
	if(e_id>e_cut&& mu_id<0.9)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//	      cout <<"is electron" <<endl;
	      }
	    m_mass=m_e;
	    massHyp=0;
	    isLepton=true;
	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	  }
	//used to be mu_id > e_cut 
	if(mu_id>mu_cut && e_id<e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//		cout <<"is muon" <<endl;
	      }
	    m_histos.hPidMuPi->Fill(chr_it->trk().pid_mu(),chr_it->trk().pid_pi());
	    m_mass=m_muon;
	    massHyp=1;
	    isLepton=true;
	    m_histos.hPidMu->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }
	if(mu_id>0.9&& e_id>e_cut)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//	      cout <<"is electron" <<endl;
	      }
	    m_histos.hPidEPi->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_pi());
	    m_mass=m_e;
	    isLepton=true;
	    m_histos.hPidE->Fill(chr_it->trk().pid_e(),chr_it->trk().pid_mu());
	  }


	if(!isLepton)
	  {
	    if(atcKPi>0.6 && atcKPi < 1.0) //kaon
		{
		  m_mass=m_k;
		  massHyp=3;
		  isPionKaon=true;
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
		  m_histos.hPidK->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_pi());
		}
	      else
		{
		  if(atcKP<0.2 && atcKP <1.0 && atcPiP<0.2 && atcPiP<1.0)
		      {
			m_mass=m_pr;
			massHyp=4;
			m_histos.hPidPr->Fill(chr_it->trk().pid_K(),chr_it->trk().pid_p());
			m_histos.hPidPrPi->Fill(chr_it->trk().pid_p(),chr_it->trk().pid_pi());

			isProton=true;
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
			if(atcKPi<0.6)
			    {
			      massHyp=2;
			      m_mass=m_pi;
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



	double dr, dz, refitPx, refitPy, refitPz;
	getDrDz(chr_it, massHyp,dr,dz, refitPx, refitPy, refitPz);
	v_vertexR.push_back(dr);
	v_vertexZ.push_back(dz);

	///
	//      cout <<"looking at " <<(*chr_it).p(0) <<" " << (*chr_it).p(1) <<" " << (*chr_it).p(2) <<endl;
	if ( fabs(dr) > cuts::vertexR )//slides from kibayashi 
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//	      cout <<"dr cut: " << fabs(dr) <<endl;
	      }
	    //	    cout <<" cut track due to vertex.. r: " << fabs(dr) <<" px lab of track:  " <<(*chr_it).p(0)<< endl;
	    continue;
	  }
	if ( fabs(dz) > cuts::vertexZ ) 
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//cout <<"dz cut: " << fabs(dz) <<endl;
	      }
	    //	    cout <<" cut track due to vertex.. z" <<endl;
	    continue;//used to be 4
	  }
	if(charge>0)
	  iChTrkPos++;
	else
	  iChTrkNeg++;
	//	cout <<"compare " << (*chr_it).p(0) << ", " << (*chr_it).p(1) <<", " << (*chr_it).p(2) <<endl;
	//	cout <<" to: "<< refitPx << " " << refitPy <<" " << refitPz <<endl;
	//	Hep3Vector h3Vect((*chr_it).p(0),(*chr_it).p(1),(*chr_it).p(2));
	Hep3Vector h3Vect(refitPx,refitPy,refitPz);
	////

	double E=sqrt(m_mass*m_mass+h3Vect.mag2());
	HepLorentzVector boostedVec(h3Vect,E);
	boostedVec.boost(kinematics::CMBoost);




	if(h3Vect.perp()<cuts::minPtThrust)
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//	      cout <<"removing pt=: " << h3Vect.perp() <<endl;
	      }
	    //	      cout <<"removing pt=: " << h3Vect.perp() <<endl;
	    continue;
	  }
	m_z=2*boostedVec.e()/kinematics::Q;
	iChTrks++;
	if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	  {
	    //	  cout <<"charged z: " << m_z<<endl;
	  }


	if(m_z<cuts::minZThrust)
	  continue;
	if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	  {
	    //	  cout <<"adding charged track: " << boostedVec.x() <<" y: " << boostedVec.y() << " z: " << boostedVec.z() << " e: " << boostedVec.e() <<endl;
	  }
	allParticlesBoosted.push_back(boostedVec.vect());
	allPB_particleClass.push_back(massHyp);
	//disasterous code...
	allPB_particleCharge.push_back(charge);

	nonBoostedE.push_back(E);
	fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.pz(),boostedVec.e()));
	//	cout <<"add to allparticles boosted " <<endl;
	allParticlesNonBoosted.push_back(h3Vect);
	allPB_E.push_back(boostedVec.e());
	if(DEBUG_EVENT==evtNr)
	  {
	    //	  (*pXCheck).precision(3);
	    //	  (*pXCheck)<<boostedVec.vect().x() << " " <<boostedVec.vect().y() << " " << boostedVec.vect().z() << " " << m_mass <<" charged " <<endl;
	  }
	visEnergy+=boostedVec.e();
	findDStar(allParticlesBoosted, allPB_particleClass, allPB_particleCharge);

	if(isLepton)
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//	      cout <<"is lepton z: " << m_z <<endl;
	      }
	    //	    cout <<"is lepton z: " << m_z <<endl;
	    //	    continue;
	  }
	if(!(isPionKaon || isProton))
	  {
	    if(DEBUG_EVENT==evtNr|| DEBUG_EVENT2==evtNr)
	      {
		//cout <<"is not pion/kaon: " << m_z <<endl;
	      }
	    //	    cout <<"is not pion/kaon: " << m_z <<endl;
	    continue;
	  }

	if(DEBUG_EVENT==evtNr)
	  {
	    //	    cout << "z: " << m_z <<"charge: " << charge <<endl;
	  }



	if(m_z<cuts::minZ)
	  {
	    //	    cout <<"didn't pass min z...: "<< m_z <<", energy: " << boostedVec.e()<<endl;
	    continue;
	  }
	if(cos(h3Vect.theta())<cuts::minCosTheta||cos(h3Vect.theta())>cuts::maxCosTheta)
	  {
	    if(DEBUG_EVENT==evtNr)
	      {
		//		cout << "CUT cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	      }
	    //	    cout << "CUT cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	    continue;
	  }
      
	if(DEBUG_EVENT==evtNr)
	  {
	    // cout << "cos theta: " << cos(h3Vect.theta()) <<"charge: " << charge <<endl;
	  }




	Particle* p=new Particle(*chr_it,string(ptypName));




	//has to be in parantheses
	p->userInfo(*(new ParticleInfo()));
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(p->userInfo());
	pinf.z=m_z;
	pinf.labTheta=h3Vect.theta();
	pinf.cmsTheta=boostedVec.theta();
	//	cout <<"theta lab:" << h3Vect.theta() <<"cms: "<< boostedVec.theta()<<endl;
	//	cout <<"theta phi:" << h3Vect.phi() <<"cms: "<< boostedVec.phi()<<endl;

	pinf.labPhi=h3Vect.phi();
	pinf.cmsPhi=boostedVec.phi();
	Ptype& m_pt=p->pType();
	//is it ok, to leave the default error matrix?
	p->momentum().momentum(boostedVec);
	//	cout <<"add to all particles for comp " <<endl;
	//only hadrons, even if unidentified trust that leptons can be separated
	if(!isLepton)
	  v_allParticles.push_back(p);
      }
    Mdst_gamma_Manager& gamma_mgr=Mdst_gamma_Manager::get_manager();
    Mdst_ecl_aux_Manager& eclaux_mgr = Mdst_ecl_aux_Manager::get_manager();


    //////--->inserted pi0

    Mdst_pi0_Manager &pi0_mgr=Mdst_pi0_Manager::get_manager();
    for(std::vector<Mdst_pi0>::const_iterator i =pi0_mgr.begin();i!=pi0_mgr.end();i++)
      {
	const Mdst_pi0& pi0=*i;
	int id =(int)pi0.get_ID();

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
	if(mass>0.15 || mass<0.12)
	  continue;
	float pLab=sqrt(px*px+py*py+pz*pz);
	//      cout <<"pi0mass: "<< mass <<endl;≈

	float g1Energy= sqrt(pi0.gamma(0).px()*pi0.gamma(0).px()+pi0.gamma(0).py()*pi0.gamma(0).py()+pi0.gamma(0).pz()*pi0.gamma(0).pz());
	float g2Energy= sqrt(pi0.gamma(1).px()*pi0.gamma(1).px()+pi0.gamma(1).py()*pi0.gamma(1).py()+pi0.gamma(1).pz()*pi0.gamma(1).pz());

	if(g1Energy < 0.05 || g2Energy < 0.05)
	  continue;
	Particle* p=new Particle(pi0);
	double confLevel;

	HepPoint3D pi0DecPoint;
	HepSymMatrix pi0ErrMatrix;

	setGammaError(p->child(0),IpProfile::position(), IpProfile::position_err_b_life_smeared());
	setGammaError(p->child(1),IpProfile::position(), IpProfile::position_err_b_life_smeared());

	if(!doKmFit(*p,  confLevel,0,m_pi0))
	  {
	    continue;
	  }

	p->userInfo(*(new ParticleInfoMass()));
	ParticleInfoMass& pinf=dynamic_cast<ParticleInfoMass&>(p->userInfo());
	pinf.gammaE1=g1Energy;
	pinf.gammaE2=g2Energy;
	pinf.e9oe25_1=e9oe25_1;
	pinf.e9oe25_2=e9oe25_2;
	pi0Candidates.push_back(p);

      }


    /////---->done



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
	fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.pz(),boostedVec.e()));
	allPB_E.push_back(boostedVec.e());
	visEnergy+=boostedVec.e();
	gammaCount++;
      }
    //    cout <<" number of charged tracks in the event: pos - " << iChTrkPos <<" neg - " << iChTrkNeg <<endl;
    //    cout <<"number of photons in the event: " << gammaCount <<endl;

    //    cout <<allParticlesBoosted.size() <<" particles in jet and thrust computation " << endl;
    if(allParticlesBoosted.size()!=nonBoostedE.size())
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



    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    ///got all candidates, now reconstruct Ds...
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    reconstructD0();
    vector<Particle*> tempParticles;
    for(vector<Particle*>::iterator itD=D0Candidates.begin();itD!=D0Candidates.end();itD++)
      {
	double confLevel;
	//		cout <<"mass before d0 mit: "<< (*itD)->p().mag()<<endl;
	if(!doKmVtxFit2(*(*itD),  confLevel,0))
	  {
	    delete *itD;
	  }
	else
	  {
	    tempParticles.push_back(*itD);
	    //	    cout <<"mass  d0 mit: "<< (*itD)->p().mag()<<endl;
	  }

      }
    D0Candidates.clear();
    
    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
      {
	D0Candidates.push_back(*itD);
      }
    tempParticles.clear();

    reconstructChargedD();
    for(vector<Particle*>::iterator itD=chargedDCandidates.begin();itD!=chargedDCandidates.end();itD++)
      {
	double confLevel;
	if(!doKmVtxFit2(*(*itD),  confLevel,0))
	  {
	    //	    cout <<" no good charged D " <<endl;
	    delete *itD;
	  }
	else
	  {
	    //	    cout <<" good charged D " << endl;
	    tempParticles.push_back(*itD);
	  }
      }
    chargedDCandidates.clear();

    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
      {
	chargedDCandidates.push_back(*itD);
      }
    tempParticles.clear();

    reconstructDStar();
    for(vector<Particle*>::iterator itD=DStarCandidates.begin();itD!=DStarCandidates.end();itD++)
      {
	double confLevel;
	//	cout <<" d star before: "<<	(*itD)->p().mag() <<endl;
	//in principle mass-vertex constrained fit here (not just mass...)
	if(!doKmVtxFit2(*(*itD),  confLevel,0))
	  {
	    //	    cout <<" no good dstar .. " <<endl;
	    delete *itD;
	  }
	else
	  {
	    //	      cout <<"found good dstar .." <<endl;
	    tempParticles.push_back(*itD);
	    //	cout <<" d star after: "<<	(*itD)->p().mag() <<endl;
	  }

      }

    DStarCandidates.clear();
    for(vector<Particle*>::iterator itD=tempParticles.begin();itD!=tempParticles.end();itD++)
      {
	DStarCandidates.push_back(*itD);
      }
    tempParticles.clear();











      /////----> end look for Ds











    Thrust t=thrustall(allParticlesBoosted.begin(),allParticlesBoosted.end(),retSelf);
    Thrust labThrust=thrustall(allParticlesNonBoosted.begin(),allParticlesNonBoosted.end(),retSelf);
    ///jet computations:

    //  cout <<"got  lab thrust " <<endl;
    ///    cout <<"lab thrust theta: " << labThrust.axis.theta() <<endl;
    //  cout <<"lab thrust phi: " << labThrust.axis.phi() <<endl;

    //    cout <<" Thrust cos theta in lab system: " << cos(labThrust.axis.theta())<<endl;
    //    cout <<" Thrust cos theta in CMS system: " << cos(t.axis.theta())<<endl;
    //    cout <<"thrust px: " << t.axis.x() << " py: "<< t.axis.y()<<" pz: " << t.axis.z()<<endl;

    //

    if(cos(labThrust.axis.theta()) > cuts::maxLabThrustCosTheta || cos(labThrust.axis.theta())<cuts::minLabThrustCosTheta)
      {
	//	cout <<" cut on lab axis " <<endl;
	exitEvent();
	return;
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

    if(rand() % 100 <50)
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
    if(kinematics::thrustMag<cuts::minThrust || abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ|| visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
      {
	bool foundReason=false;
	if(visEnergyOnFile < cuts::minVisEnergy)
	  {
	    //	  cout <<"---------------------------"<<endl<<"cut on vis energy: " << cuts::minVisEnergy <<" found only: " <<visEnergyOnFile <<" GeV" <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }
	if(abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ)
	  {
	    //	  cout <<"---------------------------"<<endl<<"cut on thrust z projection in CM: " << cuts::maxThrustZ <<" found  " <<kinematics::thrustDirCM.z()/kinematics::thrustDirCM.mag() <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }
	if(kinematics::thrustMag<cuts::minThrust)
	  {
	    //	  cout <<"---------------------------"<<endl<<"cut magnitude: " << cuts::minThrust <<" found  " <<kinematics::thrustMag <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }
	if( iChTrks < cuts::minNTracks)
	  {
	    //	    cout <<"-----------------------"<<endl <<" cut on min tracks, need; "<< cuts::minNTracks <<" have: " << iChTrks <<endl;
	    //	  cout <<"-------------------------------"<<endl;
	  foundReason=true;
	  }

	if(!foundReason)
	  {
	  //	  cout <<"exiting event due to some cut" <<endl;
	  }
	exitEvent();
	return;

      }
    //    cout <<"passed event cut " <<endl;
    kinematics::R=1.0;
    //        double R=0.55;
    JetDefinition jet_def(ee_genkt_algorithm,kinematics::R,-1);
    //  JetDefinition jet_def(cambridge_algorithm,R);
    // JetDefinition jet_def(ee_kt_algorithm);
    ClusterSequence cs(fjParticles,jet_def);
    vector<PseudoJet> jets=sorted_by_E(cs.inclusive_jets());
    // cout <<"Clustered with " << jet_def.description() <<endl;
    // cout << " pt y phi modp " <<endl;
    // cout <<"we found " << jets.size()<<" jets " <<endl;
    int numHighEJets=0;

    //    cout <<"----------------------------------"<<endl;
    //    cout <<setw(10)<<" jet # "<<setw(10)<<" Px" <<setw(10)<< "Py"<<setw(10)<<" Pz "<<setw(10)<<"E"<<setw(10)<<" # constituents"<<endl;
    //    cout <<"----------------------------------"<<endl;

    numJets->Fill(jets.size());
    for(unsigned int i=0;i<jets.size();i++)
      {
	jetEnergy->Fill(jets[i].modp());
	//	cout << setw(10)<< i <<  setw(10)<<jets[i].px()<< setw(10)<<jets[i].py()<< setw(10)<<jets[i].pz()<< setw(10)<<jets[i].e()<< setw(10)<< jets[i].constituents().size()<<endl;
	//     cout << "jet " <<i <<": "<<jets[i].perp()<<" " << jets[i].rap() << " " <<jets[i].phi()<<"   " << jets[i].modp()<<endl;//jets[i].theta()<<endl;
	if(jets[i].modp()>2.75)
	  numHighEJets++;
	vector<PseudoJet> constituents=jets[i].constituents();

      }

    //need dijet event...
    if(numHighEJets!=2)
      {
	exitEvent();
	return;
      }

    //switch jets, so that they are not energy ordered...
    int firstJetI=0;
    int secondJetI=1;
    if(kinematics::thrustZReverted)
      {
	//who knows how this is implemented in pseodojet...
	//     std::swap(*jets.begin(),*((jets.begin()++)));
	firstJetI=1;
	secondJetI=0;
      }

    vector<PseudoJet> constituents1=jets[firstJetI].constituents();
    vector<PseudoJet> constituents2=jets[secondJetI].constituents();

    kinematics::jetE1=jets[firstJetI].modp();
    kinematics::jetE2=jets[secondJetI].modp();

    kinematics::jetNumParts1=jets[firstJetI].constituents().size();
    kinematics::jetNumParts2=jets[secondJetI].constituents().size();

    kinematics::jetFluffiness1=AuxFunc::computeFluffiness(jets[firstJetI]);
    kinematics::jetFluffiness2=AuxFunc::computeFluffiness(jets[secondJetI]);


    kinematics::jet1=Hep3Vector(jets[firstJetI].px(),jets[firstJetI].py(),jets[firstJetI].pz());
    //to have the same convention as with thrust we have to flip the direction of the second jet...
    kinematics::jet2=Hep3Vector((-1)*jets[secondJetI].px(),(-1)*jets[secondJetI].py(),(-1)*jets[secondJetI].pz());

    ///-------
    kinematics::dijet[0]=kinematics::jet1;
    kinematics::dijet[1]=kinematics::jet2;

    float diff1=kinematics::dijet[0].angle(kinematics::thrustDirCM);
    float diff2=kinematics::dijet[1].angle(kinematics::thrustDirCM);
    //fill with the one close to the thrust axis

    // cout <<"diff1: " << diff1 << " diff2: " << diff2 <<endl;
    if(diff1<2)
      {
	jetThrustDiff->Fill(diff1);
	//   cout <<"diff: " << diff1 <<endl;
      }
    else
      {
	jetThrustDiff->Fill(diff2);
	//   cout <<"diff: " << diff2 <<endl;
      }
    jetJetDiff->Fill(kinematics::dijet[0].angle(kinematics::dijet[1]));
    // cout <<"angle between jets: " <<kinematics::dijet[0].angle(kinematics::dijet[1])<<endl;

    for(int iJet=0;iJet<2;iJet++)
      {
	vector<PseudoJet> constituents=jets[iJet].constituents();
	numPartInJet->Fill(constituents.size());

      }
    for(int i=0;i<allParticlesBoosted.size();i++)
      {
	float ltheta=AuxFunc::getTheta(kinematics::thrustDirCM,allParticlesBoosted[i]);
	float m_z=2*allPB_E[i]/kinematics::Q;
	m_histos.hEFlowFromThrust->Fill(ltheta,m_z);
      }
    if(evtNr==DEBUG_EVENT)
      cout <<"evt good" <<endl;
    //  cout <<"thrust cms: " << kinematics::thrustDirCM.theta() <<endl;

    //  cout <<"thrust theta, phi: " << kinematics::thrustDirCM.theta() <<" / " << kinematics::thrustDirCM.phi()<<endl;
    //  cout <<"thrust cms after possible flip: " << kinematics::thrustDirCM.theta() <<endl;
    kinematics::thetaEThrust=kinematics::thrustDirCM.angle(kinematics::firstElectronCM.vect());
    //  cout <<"theta - e beam angle: " << kinematics::thetaEThrust <<" thrust theta: "  <<kinematics::thrustDirCM.theta()<<endl;
    //    cout <<"pi0"<<endl;

    int sB=v_allParticles.size();

    //////////

    setParticleProperties();
    findHadronPairs();

#ifdef SAVE_HISTOS
    saveHistos(allParticlesBoosted, allParticlesNonBoosted);
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"about to histos" <<endl;

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"done saving histos" <<endl;
#endif
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"before save tree" <<endl;
    //  cout <<"vor tree" <<endl;



#ifdef XCHECK
    int qCounter=-1;
    (*pXCheck) <<" runNr: "<<runNr <<" eventNr: "<< evtNr;
    (*pXCheck) <<" thrust Mag: " <<kinematics::thrustMag;
    (*pXCheck) <<" thrust dir (CM) theta: "<<     kinematics::thrustDirCM.theta();
    (*pXCheck) <<" phi "<<     kinematics::thrustDirCM.phi()<<endl;

    (*pXCheck) <<" thrust dir (Lab) theta: "<<     kinematics::thrustDirLab.theta();
    (*pXCheck) <<" phi "<<     kinematics::thrustDirLab.phi()<<endl;
    
    for(vector<HadronPair*>::iterator it=v_hadronPairs.begin();it!=v_hadronPairs.end();it++)
    {
      qCounter++;
      HadronPair* hp=(*it);
      //      if(evtNr==1)
	{
	  (*pXCheck) <<" z1: "<< hp->z1 << " z2: "<< hp->z2 << " kT: "<< hp->kT <<" had type1 "<< hp->hadPType1 <<" second had type: "<< hp->hadPType2  <<" first charge: "<< hp->hadCharge1 << " second charge: "<< hp->hadCharge2 <<endl;
	}
      (*pXCheck)<<setprecision(4);
    }
    (*pXCheck) <<endl;
#endif

    saveTree();

    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"dones saving tree" <<endl;

    cleanUp();
    //      cout <<"aft cleaning"<<endl;
    if(evtNr==DEBUG_EVENT || evtNr==DEBUG_EVENT2)
      cout <<"cleaning"<<endl;
   
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
    pTreeSaver->fillWPairData(v_hadronPairs,m_evtInfo);
  }
  void ptSpect::saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted)
  {

    for(int i=0;i<v_vertexR.size();i++)
      {
	m_histos.hVertexR->Fill(v_vertexR[i]);
	m_histos.hVertexZ->Fill(v_vertexZ[i]);
      }
    for(int i=0;i<v_pi0GammaE.size();i++)
      {
	m_histos.hPi0GammaE->Fill(v_pi0GammaE[i]);
      }
    for(int i=0;i<v_gammaE.size();i++)
      {
	m_histos.hGammaE->Fill(v_gammaE[i]);
      }
    for(int i=0;i<v_asyms.size();i++)
      {
	m_histos.hGammaAsym->Fill(v_asyms[i]);
      }
    //for comparison with the thrust from the manager
    Evtcls_hadron_info_Manager& hadronInfo_mgr = Evtcls_hadron_info_Manager::get_manager();
    m_histos.hThrustOnFile->Fill(hadronInfo_mgr.begin()->Thrust());



    for(int i=0;i<v_allParticles.size();i++)
      {
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_allParticles[i]->userInfo());
	//      cout <<"sveing phi: " << pinf.phi <<", theta: " << pinf.theta<<endl;
	if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	  continue;

	m_histos.hz->Fill(pinf.z);
	m_histos.hTheta->Fill(pinf.cmsTheta);
	m_histos.hCosTheta->Fill(cos(pinf.cmsTheta));
	m_histos.hThrustProj->Fill(pinf.thrustProj);



	if(v_allParticles[i]->charge()>0)
	  {

	    m_histos.hThetaPos->Fill(pinf.cmsTheta);
	    m_histos.hZPos->Fill(pinf.z);
	  }
	else
	  {
	    if(v_allParticles[i]->charge()==0)
	      {
		m_histos.hPi0Mass->Fill(dynamic_cast<ParticleInfoMass&>(pinf).mass);

		m_histos.hThetaNeut->Fill(pinf.cmsTheta);
		m_histos.hZNeut->Fill(pinf.z);
	      }
	    else
	      {

		m_histos.hThetaNeg->Fill(pinf.cmsTheta);
		m_histos.hZNeg->Fill(pinf.z);
	      }
	  }
	m_tup->dumpData();
      }



    m_histos.hThrust->Fill(kinematics::thrustMag);
    m_histos.hThrustX->Fill(kinematics::thrustDirCM.x());
    m_histos.hThrustY->Fill(kinematics::thrustDirCM.y());
    m_histos.hThrustZ->Fill(kinematics::thrustDirCM.z());
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
		v_firstHemi.push_back(*it);
	      }

	  }
	else//particle is in the other hemisphere
	  {
	    if((*it)->pType().charge()!=0)
	      {
		//rather pointers...
		v_secondHemi.push_back(*it);
	      }

	  }

      }
  }

  void ptSpect::findHadronPairs()
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
	    hp->firstHadron=*it;
	    hp->secondHadron=*it2;
	    hp->hadCharge=AnaDef::PN;
	    hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());

	    //	  cout <<"R1: " << hp->phiR<<endl;
	    //	  hp->computeThrustTheta(kinematics::thrustDirCM);
	    hp->compute();
	    v_hadronPairs.push_back(hp);
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
    //      allParticlesBoosted.clear();


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
    for(int i=0;i<chargedKCandidates.size();i++){
      delete chargedKCandidates[i];
    }
    chargedKCandidates.clear();
    //    cout <<" done charged K"<< chargedKCandidates.size()<<endl;


    //    cout <<"pi0s: "<< pi0Candidates.size()<<endl;
    for(int i=0;i<pi0Candidates.size();i++){
      delete pi0Candidates[i];
    }
    pi0Candidates.clear();
    //    cout <<" done pi0"<<endl;

    for(int i=0;i<leptonCandidates.size();i++){
      delete leptonCandidates[i];
    }
    leptonCandidates.clear();
    //    cout <<" done leptons"<<endl;

    //    cout <<"Ks candidates: " << KsCandidates.size()<<endl;
    for(int i=0;i<KsCandidates.size();i++){
      delete KsCandidates[i];
    }
    KsCandidates.clear();


    //    cout <<" done Ks"<<endl;
    for(int i=0;i<chargedPiCandidates.size();i++){
      delete chargedPiCandidates[i];
    }
    chargedPiCandidates.clear();
    //    cout <<" done charged Pi"<<endl;
    for(int i=0;i<otherChargedTracks.size();i++){
      delete otherChargedTracks[i];
    }
    otherChargedTracks.clear();
    //    cout <<" charged tracks"<<endl;



    //      cout <<"cleaning up.."<<endl;
    for(vector<HadronPair*>::iterator it=v_hadronPairs.begin();it!=v_hadronPairs.end();it++)
      {
	delete *it;
      }


    //      cout <<"c1"<<endl;
    v_hadronPairs.clear();
    v_firstHemi.clear();
    v_secondHemi.clear();


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
    histoD0Spect->Write();
    histoDStar->Write();
    histoPiSlowMom->Write();

    thetaPhiLab->Write();
    thetaPhiCMS->Write();
    jetThrustDiff->Write();
    jetJetDiff->Write();
    jetEnergy->Write();
    numJets->Write();
    numPartInJet->Write();
    partEnergyInJet->Write();

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
		    histoD0Spect->Fill(d0candidateMass);
		    if(fabs(d0candidateMass-D0mass)<0.1)
		      {
			kinematics::D0Tag=1;
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
				    histoDStar->Fill(dStarcandidateMass);
				    if(fabs(dStarcandidateMass-DStarMass)<0.05)
				      {
					kinematics::DStarTag=1;
					histoPiSlowMom->Fill(allPB[k].mag());
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



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
