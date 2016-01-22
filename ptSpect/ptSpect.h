#ifndef DIHADANA_H
#define DIHADANA_H

#include MDST_H
#include EVTCLS_H
#include "TFile.h"
#include "TH1F.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "particle/Particle.h"

//#include "ptSpect/HadronQuadruple.h"
//#include "ptSpect/HadronPair.h"
#include "ptSpect/EventInfo.h"
#include "ptSpect/DebugHistos.h"
#include "ptSpect/AuxFunc.h"
#include "ptSpect/TreeSaver.h"
#include "ptSpect/mvaVars.h"
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;



using namespace std;
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
class HadronPair;
class HadronQuadruple;
  typedef std::vector<Gen_hepevt*> genhep_vec;
// Module class
class ptSpect : public Module 
{
public:  
  // constructor
  ptSpect( void );

  // destructor
  ~ptSpect( void );

  // initialize
  void init ( int * );
  void printD(bool star=false);
  bool recDStarMC();
  bool recD0MC();
  int goodHadronB() const;
  genhep_vec *getDaughters(const Gen_hepevt &mother);
  // begin_run function
  void begin_run ( BelleEvent*, int* );
  void findDStar(vector<Hep3Vector>& allPB, vector<int>& allPB_Class, vector<int>& allPB_Charge);
  void disp_stat ( const char* ){}
  void saveHistos( vector<Hep3Vector>& v_allParticlesBoosted, vector<Hep3Vector>& v_allParticlesNonBoosted);
  void saveTree();
  void setParticleProperties();
  bool findJetMatch(Hep3Vector& vec, vector<PseudoJet>* const1, vector<PseudoJet>* const2);
  bool recursivePrint(const Gen_hepevt gen_it, string s);
  void findHadronPairs();

  void cleanUp();
  void exitEvent();
  // histogram initialize
  void hist_def ( void );

  // event function
  void event ( BelleEvent*, int* );

  // end_run function
  void end_run ( BelleEvent*, int* );

  //  void other ( int*, BelleEvent*, int* ){}

  // terminate
  void term ( void );
  int test;
  int zNums[4];
  double smpl_;
  char rFileName[200];

  TH2D* thetaPhiLab;
  TH2D* thetaPhiCMS;

  TH1D* jetThrustDiff;
  TH1D* jetEnergy;

  TH1D* numJets;
  TH1D* jetJetDiff;

  TH1D* partEnergyInJet;
  TH1D* numPartInJet;


static int getBin(vector<float>& b1, float value)
{
  int coo1=-1;
  for(int i=0;i<b1.size();i++)
    {
      if(value<=b1[i])
	{
	  coo1=i;
	  break;
	}
    }
  return coo1;
}

  static float getPhi(const Hep3Vector& axis, const Hep3Vector& input)
    {
      return AuxFunc::getPhi(axis,input);
    }
  static float getPhi(const Hep3Vector& axis, const Particle& input)
    {
      return AuxFunc::getPhi(axis,input);
    }
  static float getTheta(const Hep3Vector& axis, const Particle& input)
    {
      return AuxFunc::getTheta(axis,input);
    }
  //reference particle types
  Ptype cPiPlus;
  Ptype cPiNeg;
  Ptype cPiZero;
  Ptype cKPlus;
  Ptype cKNeg;
 protected: 
    vector<Particle*> chargedPiCandidates;
    vector<Particle*> chargedKCandidates;
    vector<Particle*> pi0Candidates;
    vector<Particle*> KsCandidates;
    vector<Particle*> leptonCandidates;
    vector<Particle*> D0Candidates;
    vector<Particle*> chargedDCandidates;
    vector<Particle*> DStarCandidates;

    void reconstructD0();
    void reconstructChargedD();
    void reconstructDStar();
    unsigned doKmFit(Particle &p, double& conflevel, int debug, double mass=0);
    unsigned doKmVtxFit(Particle &p, double& conflevel, int debug);
    unsigned doKmVtxFit2(Particle &p, double& conflevel, int debug, double mass=0);
private:
  //compute distance between decay vertices of quark and antiquark
  float getDecayDist();
  float getDecayDistK();
  float getDecayDistD();
  float getDecayDistPN();
  Hep3Vector getVertex(bool firstHemi);

  TH1D* histoD0Spect;
  TH1D* histoDStar;
  TH1D* histoPiSlowMom;

  TH1D* histoRecDStarSpectToDPi;
  TH1D* histoRecDStarSpectToD0Pi;


  float getTheta(Particle* p1,Particle* p2)
  {
	    double E1,E2;
	    E1=p1->e();
	    E2=p2->e();
	    Hep3Vector mH1=p1->p().vect();
	    Hep3Vector mH2=p2->p().vect();
	    Hep3Vector P_h=p1->p().vect()+p2->p().vect();
	    mH1.rotateZ(-P_h.phi());
	    mH1.rotateY(-P_h.theta());
	    float beta=(mH1.z()+mH2.z())/(E1+E2);
	    float gamma=1/sqrt(1-beta*beta);
	    mH1.setZ(-gamma*beta*E1+gamma*mH1.z());
	    float decayTheta=mH1.theta();
	    return decayTheta;
  }

  float getZ(Particle* p1,Particle* p2)
  {
	    double E1,E2;
	    E1=p1->e();
	    E2=p2->e();
	    float m_z=2*(E1+E2)/kinematics::Q;
	    return m_z;
  }

  //inits the tree and gives the treesaver the mva variable
  void initMvaTree();
  //tree to keep the vars needed for mva to discern charm/uds
  TTree* mvaTree;
  mvaVars m_mvaVars;

  int numKPlusFH;
  int numKMinusFH;
  int numPiPlusFH;
  int numPiMinusFH;

  int numKPlusSH;
  int numKMinusSH;
  int numPiPlusSH;
  int numPiMinusSH;


  int numKPiFH;
  int numPiKFH;
  int numKPiSH;
  int numPiKSH;

  int numKPQ;
  int numPKQ;

  int numPi0FH;
  int numPi0SH;
  bool validRun;

  bool enoughSVDHits(Mdst_charged_Manager::iterator);
  void getDrDz(Mdst_charged_Manager::iterator, int, double&, double&, double&, double&, double&);
  BelleTuple* T_dist;
  BelleTuple* m_tup;
  BelleTuple* m_tupG;
  BelleTuple* m_tupThrust;
  BelleTuple* m_tupThrustPa;
  BelleTuple* m_tupPi0;
  BelleTuple* m_tupEvData;
  vector<HadronPair*> v_hadronPairs;

  vector<Particle*> v_firstHemi;
  vector<Particle*> v_secondHemi;
  vector<Particle*> v_allParticles;

  vector<double> v_vertexR;
  vector<double> v_vertexZ;

  //for histogramms of gamma energy coming from pi0 and not from pi0
  vector<float> v_pi0GammaE;
  vector<float> v_gammaE;
  vector<float> v_asyms;

  vector<HadronQuadruple*> v_hadronQuadruplesPN;
  vector<HadronQuadruple*> v_hadronQuadruplesPNeut;
  vector<HadronQuadruple*> v_hadronQuadruplesNeutN;
  vector<HadronQuadruple*> v_hadronQuadruplesNeutNeut;
  vector<HadronQuadruple*> v_hadronQuadruplesAll;
  //event level info
  EventInfo m_evtInfo;
  TreeSaver* pTreeSaver;
  vector<float> dataF; //float data saved in the tree
  vector<int> dataI; //int data    "

  TFile* m_file;
  DebugHistos m_histos;

  //get phi after the given vector is the z direction

};


extern "C" Module_descr *mdcl_ptSpect()
{
  ptSpect* module = new ptSpect;
  Module_descr* dscr = new Module_descr("ptSpect", module);
  dscr->define_param("smpl", "test parameter", &module->smpl_);
  //hopefully the int is the maximum lenght of the string...
  dscr->define_param("rfname","root file name","S",100,&module->rFileName);
  //registers parameters of ipprofile...
  IpProfile::define_global(dscr);
  return dscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
