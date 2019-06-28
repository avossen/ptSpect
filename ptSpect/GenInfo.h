#ifndef GENINFO_H
#define GENINFO_H

#include "belle.h"
#include <mdst/mdst.h>
#include MDST_H
#include EVTCLS_H
#include MDST_OBS_H
#include HEPEVT_H
#include TRK_H
#include <mdst/Muid_mdst.h>
#include "tuple/BelleTupleManager.h"
#include "ptSpect/AnaConsts.h"
#include "belle.h"
#include "TTree.h"
#include <string>
#include <iomanip>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "ptSpect/HadronPair.h"
#include "particle/Particle.h"
#include "particle/Ptype.h"
#include "ptSpect/Savable.h"
#include "ptSpect/ParticleInfoMass.h"
#include "ptSpect/DebugHistos.h"
#include BELLETDF_H

#include "fastjet/ClusterSequence.hh"
#include <iostream>

using namespace fastjet;
using namespace std;

//taken from Jimmies code...
#define ID_GAMMA 22
#define ID_PI 211
#define ID_PI0 111
#define ID_K 321
#define  ID_P 2212


#if defined(BELLE_NAMESPACE)
namespace Belle {

#endif
  /*
    class that encapsulates generator info and methods to fill it
  */
  Hep3Vector& gretSelf(Hep3Vector& vec);

  class GenInfo
  {

  public:
    GenInfo():cPiPlus("PI+"), cPiNeg("PI-"),cPiZero("PI0"),cKPlus("K+"),cKNeg("K-"), cElectron("e-"), cPositron("e+"), cMuPlus("MU+"), cMuNeg("MU-")
      {
	//	mGZBins.push_back(0.2);
	mGZBins.push_back(0.3);
	mGZBins.push_back(0.4);
	mGZBins.push_back(0.55);
	mGZBins.push_back(0.75);
	mGZBins.push_back(1.1);
      }
      int getIdxFromGeantId(int geantId)
      {
	if(fabs(geantId)==lc_pPlus)
	  {
	    return 4;
	  }
	if(fabs(geantId)==lc_kPlus)
	  {
	    return 3;
	  }
	//pion
	return 2;
      }
  


      int getBin(vector<float>& b1, float value)
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

      vector<float> mGZBins;
    Hep3Vector cmThrust;

    float cmThrustMag;
    Hep3Vector labThrust;
    float thrustPR;
    float vpEnergy;
    float vpPx;
    float vpPy;
    float vpPz;
    int numQuarks;
    float isrPhotonEnergy;
    float quarkAngle;
    float thetaEThrust;


    //the adresses from which the tree should read its data
    static Savable tData;
    int num_f_data;
    int num_fA_data;
    int num_iA_data;
    int num_i_data;

    void setDebugHistos(DebugHistos* m_histos)
      {
	this->m_histos=m_histos;
      }

    void finalize()
    {
      tData.pDataTree->Write();
    }

    void fillInf()
    {
      //      cout <<"compute thrust " <<endl;
      computeGenThrust();
      thrustPR=-log(tan(cmThrust.theta()/2));
    }

    void doAll(vector<Particle*>& ap, bool eventCut)
    {
      //      cout <<"do all.." <<endl;
      //hopefully after thust computation
      //  if(kinematics::thrustMag<cuts::minThrust || abs(kinematics::thrustDirCM.z())/kinematics::thrustDirCM.mag()>cuts::maxThrustZ|| visEnergyOnFile<cuts::minVisEnergy || iChTrks < cuts::minNTracks)
      //      if(cmThrustMag<cuts::minThrust|| abs(cmThrust.z()/cmThrust.mag()> cuts::maxThrustZ))
      //	 return;
      vector<Particle*> v_allParticles;
      vector<Particle*> v_firstHemi;
      vector<Particle*> v_secondHemi;
      vector<HadronPair*> v_hadronPairs;

      getHadronLists(v_allParticles);
      //      cout <<" w/o acceptance " << endl;
      //      cout <<"we have " << v_allParticles.size() <<" particles " <<endl;
      setParticleProperties(v_allParticles, cmThrust,v_firstHemi, v_secondHemi);
      findHadronPairs(v_allParticles, v_hadronPairs, ap, eventCut);   
      //
      fillWPairData(v_hadronPairs);
      //      cout <<"push back event data.. " <<endl; 
      int numF=num_fA_data;
      int numI=num_iA_data;

      //insert here

      tData.dataF.push_back(cmThrust.theta());
      tData.dataF.push_back(cmThrust.phi());
      tData.dataF.push_back(labThrust.theta());
      tData.dataF.push_back(cmThrustMag);
      tData.dataF.push_back(thetaEThrust);
      tData.dataI.push_back(kinematics::runNr);
      tData.dataI.push_back(kinematics::evtNr);

      //      cout <<"save data.. " <<endl;
      if(v_hadronPairs.size()>0)
	{
	  //	  cout <<"saving " << v_hadronPairs.size() <<" gen hadron pairs " <<endl;
	  saveData(tData.dataF,tData.dataI,2*(numI+numF));
	}
      //      cout <<"now clean up. " <<endl;
      //get ready for next call where dataF has to be clean (this is static)
      tData.dataF.clear();
      tData.dataI.clear();
      //and the event data has to be filled afterwards:
      cleanUp(v_allParticles,v_firstHemi,v_secondHemi,v_hadronPairs);
    }
    //reference particle types
    Ptype cPiPlus;
    Ptype cPiNeg;
    Ptype cPiZero;
    Ptype cKPlus;
    Ptype cKNeg;
    Ptype cElectron;
    Ptype cPositron;
    Ptype cMuPlus;
    Ptype cMuNeg;

    void initializeTree()
    {
      if(tData.initialized)
	return;
      num_f_data=0;
      num_fA_data=0;
      num_iA_data=0;
      num_i_data=0;
      //the first time the class is initialized, construct the tree
      tData.pDataTree=new TTree("GenTree","Generated Tree");
      addArrayF("z1_mcWoA");
      addArrayF("z2_mcWoA");

      addArrayF("labTheta1_mcWoA");
      addArrayF("labTheta2_mcWoA");
      
      addArrayF("labPhi1_mcWoA");
      addArrayF("labPhi2_mcWoA");

      addArrayF("cmsTheta1_mcWoA");
      addArrayF("cmsTheta2_mcWoA");

      addArrayF("cmsPhi1_mcWoA");
      addArrayF("cmsPhi2_mcWoA");

      addArrayF("thrustProj1_mcWoA");
      addArrayF("thrustProj2_mcWoA");

      addArrayF("kT_mcWoA");
      //hadron quad level
      //qT
      addArrayF("HadDiffTheta_mcWoA");
      addArrayF("HadDiffPhi_mcWoA");
      addArrayF("qT_mcWoA");

      //hadron pair level seems to be fine

      addArrayI("chargeType_mcWoA");
      addArrayI("particleType_mcWoA");
      addArrayI("chargeType1_mcWoA");
      addArrayI("particleType1_mcWoA");
      addArrayI("chargeType2_mcWoA");
      addArrayI("particleType2_mcWoA");
    
      //missing on event level:
      //thrustTheta
      //thrustThetaLab
      //thetaEThrust

      //fields after arrays to make the filling work....
      addFieldF("thrustTheta_mcWoA");
      addFieldF("thrustPhi_mcWoA");
      addFieldF("thrustThetaLab_mcWoA");
      addFieldF("thrustMag_mcWoA");
      addFieldF("thetaEThrust_mcWoA");

      addFieldI("runNr");
      addFieldI("evtNr");

      tData.initialized=true;
    };
  private:
    void addFieldF(char* fieldname)
    {
      num_f_data++;
      //construct the memory location from which the tree should read the new data field
      float* memLoc=new float;
      tData.treeData.push_back(memLoc);
      tData.pDataTree->Branch(fieldname, memLoc, (fieldname+string("/F")).c_str());
      tData.fieldNamesF.push_back(fieldname);
    };


    void addFieldI(char* fieldname)
    {
      num_i_data++;
      int* memLoc=new int;
      tData.treeData.push_back(memLoc);
      tData.pDataTree->Branch(fieldname,memLoc,(fieldname+string("/I")).c_str());
      tData.fieldNamesI.push_back(fieldname);
    };

    void addArrayI(char* fieldname)
    {
      num_iA_data++;
      //standard lenth, shouldn't be more than that
      int* memLoc=new int[1200];
      int* memLocCounter=new int;
      tData.treeData.push_back(memLocCounter);
      tData.treeData.push_back(memLoc);
      string counterName=string(fieldname)+string("Counter");
      tData.pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
      tData.pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/I")).c_str());
    };

    void addArrayF(char* fieldname)
    {
      num_fA_data++;
      float* memLoc=new float[1200];
      int* memLocCounter=new int;
      tData.treeData.push_back(memLocCounter);
      tData.treeData.push_back(memLoc);
      string counterName=string(fieldname)+string("Counter");
      tData.pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
      tData.pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
    };

    
  void fillWPairData(vector<HadronPair*>& v_pair, bool mcPart=false)
    {
      int numI=-1;
      int numF=-1;
      int counter=0;//offset from one quadruple to the next
      for(vector<HadronPair*>::iterator it=v_pair.begin();it!=v_pair.end();it++)
      {
	HadronPair* pair=*it;
	if(pair==0)
	  {
	    for(int i=0;i<num_fA_data;i++)
	      {
		tData.dataF.push_back(-1);
	      }
	    for(int i=0;i<num_iA_data;i++)
	      {
		tData.dataI.push_back(-1);
	      }
	  }
	else
	  {
	    tData.dataF.push_back(pair->z1);//z1
	    tData.dataF.push_back(pair->z2); //z2

	    float labTheta1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).labTheta;
	    float labTheta2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).labTheta;

	    float labPhi1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).labPhi;
	    float labPhi2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).labPhi;

	    float cmsTheta1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).cmsTheta;
	    float cmsTheta2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).cmsTheta;

	    float cmsPhi1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).cmsPhi;
	    float cmsPhi2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).cmsPhi;
	
	
	    tData.dataF.push_back(labTheta1);
	    tData.dataF.push_back(labTheta2);
	    
	    tData.dataF.push_back(labPhi1);
	    tData.dataF.push_back(labPhi2);

	
	    tData.dataF.push_back(cmsTheta1);
	    tData.dataF.push_back(cmsTheta2);
	    
	    tData.dataF.push_back(cmsPhi1);
	    tData.dataF.push_back(cmsPhi2);

	    tData.dataF.push_back(pair->firstHadron->p3().dot(cmThrust)/(pair->firstHadron->p3().mag()*cmThrust.mag()));
	    tData.dataF.push_back(pair->secondHadron->p3().dot(cmThrust)/(pair->secondHadron->p3().mag()*cmThrust.mag()));
	    tData.dataF.push_back(pair->kT);
	    tData.dataF.push_back(pair->diffTheta);
	    tData.dataF.push_back(pair->diffPhi);

	    tData.dataF.push_back(pair->qT);

	    tData.dataI.push_back(pair->hadCharge);      
	    tData.dataI.push_back(pair->hadPType);//particle type

	    tData.dataI.push_back(pair->hadCharge1);      
	    tData.dataI.push_back(pair->hadPType1);
	    tData.dataI.push_back(pair->hadCharge2);      
	    tData.dataI.push_back(pair->hadPType2);

	  }
	//	cout <<"4"<<endl;      
	numF=tData.dataF.size();  //always the same, 
	numI=tData.dataI.size(); 
	//how many array values do we have? At this point we haven't added event data yet...
	for(int j=0;j<tData.dataF.size();j++)
	  {
	    //	    cout <<"looking at field nr " << j <<endl;
	    if(counter >=1200)
	      {
		cout <<"index q to large" <<endl;
		continue;
	      }
	    //tree data has all fields
	    if(2*j+1 > (tData.treeData.size()-num_f_data))
	    {
	      cout <<"index td to large" << tData.treeData.size() << ", " << 2*j+1<<" num scalar fields: "<< num_f_data <<endl;
	      continue;
	    }
	    ((float*)tData.treeData[2*j+1])[counter]=tData.dataF[j];
	  }
	//tData.dataI has size 4
	//	cout <<" off to integeres.. " <<endl;
	for(int j=0;j<tData.dataI.size();j++)
	  {
	    if(2*(j+numF)+1 > tData.treeData.size())
	      {
		cout <<"index2 td to large" << tData.treeData.size() << ", " << 2*(j+numF)+1 <<endl;
		continue;
	      }

	    ((int*)tData.treeData[2*(j+numF)+1])[counter]=tData.dataI[j];
	  }
	//has to be cleared to save next quad!
       	tData.dataF.clear();
	tData.dataI.clear();
	counter++;
      }
      //	cout <<"5" <<endl;
      counter=v_pair.size();
      //save counter info
      for(int j=0;j<numF;j++)
	{
	  *(int*)tData.treeData[2*j]=counter;
	}
      for(int j=0;j<numI;j++)
	{
	  *(int*)tData.treeData[2*(j+numF)]=counter;
	}
    
      if(numF < 0) //no events, all quad vector empty
	return;

    }

  //saves event data
    void saveData(vector<float>& dataF, vector<int>& dataI, int offset)
    {
      int dataSize=dataF.size()+dataI.size();
      //      cout <<"data size: "<< dataSize <<endl;
      if(dataSize !=tData.treeData.size()-offset)
	{
	  cout << "data size does not match number of branches " <<dataSize <<"+ " << offset << " = " <<tData.treeData.size() <<endl <<flush;
	  exit(0);
	}
      //      cout <<"puting " << dataF.size() <<" floats " << endl;
      for(int i=0;i<dataF.size();i++)
	{
	  *(float*)tData.treeData[i+offset]=dataF[i];
	}
      //      cout <<" putting " << dataI.size() << " integers " << endl;
      for(int i=0;i<dataI.size();i++)
	{
	  *(int*)tData.treeData[i+dataF.size()+offset]=dataI[i];
	}
      //          cout <<"fill..." <<endl;
      tData.pDataTree->Fill();
      //          cout <<"done filling " <<endl;
    };
    //unfortunately just a copy of the ptSpect functions
    void getHadronLists(vector<Particle*>& v_allParticles)
    {
      int motherGenId=0;
      Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
      //      cout <<"iterating" <<endl;
      for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
	{
	  int geantID=abs(gen_it->idhep());//plus is ok, since it is the abs value
	  //for now, take all stable charged particles
	  if(gen_it->isthep()!=1) //isthep==1: stable, ==2, unstable
	    {
	      continue;
	    }
	  //in principle also take leptons... but fine for now. Assume that e id and mu_id is pretty good
	  //pi0 should not be stable...
	  //	  if(!(geantID==lc_pi0 || geantID==lc_piPlus || geantID==lc_kPlus || geantID==lc_pPlus))
	  if(!(geantID==lc_piPlus || geantID==lc_kPlus || geantID==lc_pPlus))
	    continue;

	  //	   cout <<"valid " <<endl;
	  Particle* np=new Particle(*gen_it);
	  //	  cout<<" adding lund: "<< np->lund() << endl;
	  HepLorentzVector boostedVec(gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E());
	  float labTheta=boostedVec.theta();
	  float labPhi=boostedVec.phi();
	  //	  cout <<"got tehta: " << m_theta << " is there a diff with vec3? " << boostedVec.vect().theta()<<endl;
	  boostedVec.boost(kinematics::CMBoost);
	  np->userInfo(ParticleInfoMass()); //gets deleted in destructor of Particle
	  ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(np->userInfo());
	  float m_z=2*boostedVec.t()/kinematics::Q;
	  pinf.motherGenId=motherGenId;
	  //I don't think these fields are used, let's just use [0] and disentangle later with particleId
	  //	  pinf.z[getIdxFromGeantId(geantID)]=m_z;
	  pinf.z[0]=m_z;
	  pinf.z[getIdxFromGeantId(geantID)]=m_z;
	  pinf.boostedMoms[getIdxFromGeantId(geantID)]=boostedVec.vect();
	  pinf.boostedLorentzVec[getIdxFromGeantId(geantID)]=boostedVec;
	  //	  cout <<"geant id is : "<< geantID <<" index: "<< getIdxFromGeantId(geantID) <<endl;
	  //	  cout <<"adding particle with z: "<< m_z <<endl;
	  //need some cutoff so we are not swamped with pairs. This should be the same as in the analysis

	  //  float minCosTheta=-0.511; //cuts for Martin's PID
	  //  float maxCosTheta=0.842; //cuts for Martin's PID

	  //no fiducial cuts based on lab theta of particles for smearing correction
	  if(m_z<0.05)// || cos(labTheta)<-0.511 || cos(labTheta)>0.842)
	    {
	      delete np;
	      continue;
	    }
	  //theta should be the one in the lab system as before...
	  pinf.labTheta=labTheta;
	  pinf.labPhi=labPhi;
	  pinf.cmsTheta=boostedVec.theta();
	  pinf.cmsPhi=boostedVec.phi();

	  Hep3Vector axis=cmThrust;
	  pinf.thrustProj=axis.dot(boostedVec.vect())/(axis.mag()*boostedVec.vect().mag());

	  np->momentum().momentum(boostedVec);
	  v_allParticles.push_back(np);
	}
    }

    void setParticleProperties(vector<Particle*>& v_allParticles, Hep3Vector& mThrustDir,vector<Particle*>& v_firstHemi, vector<Particle*>& v_secondHemi)
    {
      for(vector<Particle*>::const_iterator it=v_allParticles.begin();it!=v_allParticles.end();it++)
	{
	  //	  Hep3Vector axis=jet1;
	  Hep3Vector axis=mThrustDir;
	  //	  if(jet1.dot((*it)->p().vect())<0)
	    //	    axis=jet2;
	  ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	  float phi=AuxFunc::getPhi(axis,**it);

	  //apparently in getHadronlist, the cm theta is saved here:
	  float cmTheta=pinf.cmsTheta;
	  float theta=AuxFunc::getTheta(axis,**it);
	  pinf.phi=phi;
	  //	  pinf.theta=theta;

	  pinf.thrustProj=axis.dot((*it)->p().vect())/(axis.mag()*(*it)->p().vect().mag());	  

	  if(fabs(pinf.thrustProj)<cuts::minThrustProj)
	    {
	      //no cuts...? doch...
	      //	      continue;
	    }
	  if(pinf.thrustProj>0)//one hemisphere
	    {
		v_firstHemi.push_back(*it);

	    }
	  else//particle is in the other hemisphere
	    {
		v_secondHemi.push_back(*it);
	    }
	}
    }

    void findHadronPairs(vector<Particle*>& v_allParticles,vector<HadronPair*>& v_hadronPairs, vector<Particle*>& ap, bool eventCut)
    {
      for(vector<Particle*>::const_iterator it=v_allParticles.begin();it!=v_allParticles.end();it++)
	{
	  for(vector<Particle*>::const_iterator it2=(it+1);it2!=v_allParticles.end();it2++)
	    {
	      ///check if that was already in the rec
	      ////
	      bool foundFirst=false;
	      bool foundSecond=false;
	      //	      cout <<"checking vs " << ap.size()<<" particles " <<endl;
	      for(vector<Particle*>::const_iterator itP=ap.begin();itP!=ap.end();itP++)
		{
		  Gen_hepevt gph;      
		  gph=get_hepevt((*itP)->mdstCharged());
		  if((*it)->genHepevt()==gph)
		    foundFirst=true;
		  if((*it2)->genHepevt()==gph)
		    foundSecond=true;
		}
	      if(foundFirst&& foundSecond)
		{
		  //		  cout <<"found hadron pair already saved " << endl;
		}

	      if(!eventCut && foundFirst&&foundSecond)
		{
		  //only save the ones that were dropped due to acceptance
		  continue;
		}
	      ///
	      //cout <<" is this a new gen hadron pair? " << endl;
	      //back-to-back pair
	      //also have to implement fiducial cuts-->do that in the particle selection when we have the lab momentum



	      if((*it)->p().vect().dot((*it2)->p().vect())<0)
		{
		  //  cout <<"yes it is.." <<endl;
		  //now we have to check if that is already in our reconstruction pairs or should be part of the acceptance correction, otherwise it is already paired with
		  //the accepted ones...
		  HadronPair* hp=new HadronPair();
		  hp->firstHadron=*it;
		  hp->secondHadron=*it2;

		  //done in HadronPair::compute now...
		  //	      nhp->hadCharge=AnaDef::PN;
		  //	      hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
		  //first hemi
		  hp->compute();
		  //		  cout <<" creating pair with z1: "<< hp->z1 <<" z2: "<< hp->z2 <<endl;
		  v_hadronPairs.push_back(hp);


		}
	    }
	}
    }

    void findHadronPairsThrust(vector<Particle*>& v_firstHemi, vector<Particle*>& v_secondHemi, vector<HadronPair*>& v_hadronPairs)
    {

      //changed all this to just this condition:
      //		      if(pinf.boostedMoms[i].dot(pinf2.boostedMoms[j])<0)
      //      cout <<"find hadron pairs with " << v_firstHemiPos.size() << " in pos and " << v_firstHemiNeg.size() <<endl;
      for(vector<Particle*>::const_iterator it=v_firstHemi.begin();it!=v_firstHemi.end();it++)
	{

	  ParticleInfo& pinf=dynamic_cast<ParticleInfo&>((*it)->userInfo());
	  for(vector<Particle*>::const_iterator it2=v_secondHemi.begin();it2!=v_secondHemi.end();it2++)
	    {
	      ParticleInfo& pinf2=dynamic_cast<ParticleInfo&>((*it2)->userInfo());
	      if(pinf.z[0]+pinf2.z[0] < cuts::min2H_Z)
		{
		  //continue;
		}
	      //this doesn't work for non pion pairs
	      if(pinf.z[0]+pinf2.z[0] < 0.001)
		{
		  //	  continue;
		}
	      //now unknowns...
	      HadronPair* hp=new HadronPair();
	      hp->firstHadron=*it;
	      hp->secondHadron=*it2;

	      //done in HadronPair::compute now...
	      //	      hp->hadCharge=AnaDef::PN;
	      //	      hp->hadPType=AuxFunc::getPType((*it)->pType(),(*it2)->pType());
	      //first hemi
	      hp->compute();
	      v_hadronPairs.push_back(hp);
	    }
	}

    }



    void cleanUp(vector<Particle*>& v_allParticles,vector<Particle*>& v_firstHemi, vector<Particle*>& v_secondHemi, vector<HadronPair*>& v_hadronPairs)
    {

      for(vector<HadronPair*>::iterator it=v_hadronPairs.begin();it!=v_hadronPairs.end();it++)
	{
	  delete *it;
	}


      v_hadronPairs.clear();

      //      cout <<"c2"<<endl;

      v_firstHemi.clear();
      v_secondHemi.clear();


      for(int i=0;i<v_allParticles.size();i++)
	{
	  delete v_allParticles[i];
	}
      v_allParticles.clear();
    }


    //mirrors thrust computation for real particles
    void computeGenThrust()
    {
      numQuarks=0;
      isrPhotonEnergy=0.0;
      vector<Hep3Vector> allParticlesBoosted;  
      vector<PseudoJet> fjParticles;
      //in principle easier to just boot the thrust vector I guess...
      vector<Hep3Vector> allParticlesNonBoosted; 
      vector<float> allPB_E;
      Gen_hepevt_Manager& gen_hep_Mgr=Gen_hepevt_Manager::get_manager();
      //      	cout <<"-------------------------new evt--------------------------"<<endl;
      //	cout <<setw(9)<<" ID "<<setw(9)<<"ISTHEP"<<setw(9)<<" LUND_ID  "<<setw(9)<<"Mo first"<<setw(9)<<"    MoLast"<<setw(9)<<"   daFirst"<<setw(9)<<" daLast "<<setw(9)<<" PX  "<<setw(9)<<"PY"<<setw(9)<<"  PZ "<<setw(9)<<"  E "<<setw(9)<<"  M  "<<setw(9)<<" VX "<<setw(9)<<"  VY "<<setw(9)<<"  VZ"<<setw(9)<<"   T " <<endl;
      HepLorentzVector* lv_q1=0;
      HepLorentzVector* lv_q2=0;
      float m_px=0;
      float m_py=0;
      float m_pz=0;
      float m_e=0;

      float m_px2=0;
      float m_py2=0;
      float m_pz2=0;
      float m_e2=0;

      float m_pxVp=0;
      float m_pyVp=0;
      float m_pzVp=0;
      float m_eVp=0;


      HepLorentzVector boostVP(gen_hep_Mgr.begin()->PX(),gen_hep_Mgr.begin()->PY(),gen_hep_Mgr.begin()->PZ(),gen_hep_Mgr.begin()->E());
      boostVP.boost(kinematics::CMBoost);
      vpPx=boostVP.px();
      vpPy=boostVP.py();
      vpPz=boostVP.pz();
      vpEnergy=boostVP.e();

      float pxSum=0.0;
      float pySum=0.0;
      float pzSum=0.0;
      float eSum=0.0;

      //      cout <<endl;
      isrPhotonEnergy=0.0;
      for(Gen_hepevt_Manager::iterator gen_it=gen_hep_Mgr.begin();gen_it!=gen_hep_Mgr.end();gen_it++)
	{

	  //	 	  cout << setw(6)<< gen_it->get_ID() << " |  "  << setw(6)<<gen_it->isthep() << " | " << setw(6)<<gen_it->idhep()  <<" | " <<  setw(6)<< gen_it->moFirst() <<" | " <<setw(6)<< gen_it->moLast() << " | " <<  setw(6)<<gen_it->daFirst() <<" | " <<setw(6)<< gen_it->daLast() <<" | " <<setw(6)<< gen_it->PX() << " | " << setw(6)<< gen_it->PY() <<" | " << setw(6)<< gen_it->PZ() <<" | " <<setw(6)<< gen_it->E() << " | " << setw(6)<< gen_it->M() <<" | " << setw(6)<< gen_it->VX() << " | " << setw(6)<< gen_it->VY() << " | " << setw(6)<< gen_it->VZ() << " | " << setw(6)<< gen_it->T() <<endl;

	  //look for ISR photons. In hepevt these are photons that come from the virtual photon
	  if(gen_it->idhep()==22 && gen_it->mother() && gen_it->mother().idhep()==10022)
	    {
	      isrPhotonEnergy+=gen_it->E();
	    }
	  //	  cout <<"-->isrPhotonEnergy now: " << isrPhotonEnergy <<endl;

	  /*	    if(gen_it->idhep()>94 && gen_it->idhep() <100)
	    cout <<"found special code : " << gen_it->idhep() <<endl;
	    if(gen_it->get_ID()==1)
	    {
	    m_pxVp=gen_it->PX();
	    m_pyVp=gen_it->PY();
	    m_pzVp=gen_it->PZ();
	    m_eVp=gen_it->E();
	    }
	    ///if(gen_it->idhep()!=0 && gen_it->moFirst()==1 && (gen_it->idhep() < 20 || (gen_it->idhep()>20 && gen_it->idhep() <40)  ))//either lepton or gauge boson that came from gamma*
	    //	    	/*       if(gen_it->idhep()!=0 && gen_it->moFirst()==1 && (gen_it->idhep() > 110 || (gen_it->idhep()>1 && gen_it->idhep() <40)  ))//either lepton or gauge boson that came from gamma* */
	  /* 	    	    if(gen_it->idhep()!=0 && gen_it->moFirst()==1) */
	  /* 	      { */
	  /* 		if(gen_it->daLast()==-1 )//either lepton or gauge boson that came from gamma* */
	  /* 		  { */
	  /* 		    m_px+=gen_it->PX(); */
	  /* 		    m_py+=gen_it->PY(); */
	  /* 		    m_pz+=gen_it->PZ(); */
	  /* 		    m_e+=gen_it->E(); */
	  /* 		  } */
	  /* 		else */
	  /* 		  { */
	  /* 		    m_px2+=gen_it->PX(); */
	  /* 		    m_py2+=gen_it->PY(); */
	  /* 		    m_pz2+=gen_it->PZ(); */
	  /* 		    m_e2+=gen_it->E(); */
	  /* 		  } */
	  /* 		  }*\/ */



	   	       if(abs(gen_it->idhep()) <10 && gen_it->idhep()!=0) 
	   		{ 
	   		  if(numQuarks==0) 
			    {
			      if(lv_q1!=0)
				delete lv_q1;
	   		    lv_q1=new HepLorentzVector(gen_it->PX(),gen_it->PY(),gen_it->PZ(), gen_it->E()); 
			    }
	   		  if(numQuarks==1) 
	   		    { 
			      if(lv_q2!=0)
				delete lv_q2;

	   		      lv_q2=new HepLorentzVector(gen_it->PX(),gen_it->PY(),gen_it->PZ(), gen_it->E()); 
	   		      lv_q1->boost(kinematics::CMBoost); 
	   		      lv_q2->boost(kinematics::CMBoost); 
	   		      } 
	   		  numQuarks++; 
		  
	   		} 
	
	  /* 	    /\*	    	    	    if(abs(gen_it->idhep()) <10 && gen_it->idhep()!=0) */
	  /* 	      { */
	  /* 		//		cout <<gen_it->idhep()<<endl; */

	  /* 		//		cout <<"daFirst: " << gen_it->daFirst()<< ", daSecond: " << gen_it->daLast() <<endl; */
	  /* 		//		cout <<"moFirst: " << gen_it->moFirst()<< ", moSecond: " << gen_it->moLast() <<endl; */
	  /* 		//		cout <<"mother: " << gen_it->mother() <<endl; */
	  /* 		if(gen_it->mother()) */
	  /* 		  { */
	  /* 		    //10022 is virtual photon */
	  /* 		    cout <<"motherID " << gen_it->mother().idhep() << " mother addr: " << gen_it->get_ID()<<endl; */
	  /* 		    cout <<"MotherE: " << gen_it->mother().E() <<", px : "<< gen_it->mother().PX() << " , py: " << gen_it->mother().PY() << " pz: " << gen_it->mother().PZ() <<endl; */
	  /*	    cout <<"mother addr: " << gen_it->mother() <<" dafirst : " << gen_it->mother().daFirst() << " last: " << gen_it->mother().daLast() <<endl;
	    if(gen_it->mother().mother()!=0)
	    cout <<"mother has a mother ;-) " <<endl;

	    cout <<"dafirst: " << gen_it->daFirst() <<" last: " << gen_it->daLast()<<endl;
	    }
	    cout <<"E: " << gen_it->E() <<", vx : "<< gen_it->VX() << " , vy: " << gen_it->VY() << " vz: " << gen_it->VZ() << " T: " <<gen_it->T() <<endl;
	    cout <<"px : "<< gen_it->PX() << " , py: " << gen_it->PY() << " pz: " << gen_it->PZ() <<endl;
	    cout <<"quarkCount: " <<quarkCount<< endl;
	    }*/
	  int gId=abs(gen_it->idhep());
	  //	    if(gId==ID_GAMMA || gId==ID_PI || gId==ID_K)
	  //	    if(gen_it->mother()==1 && gen_it->daLast()==-1)
	  //	    if(gen_it->isthep()==2 || gen_it->mother()==1 && gen_it->idhep()<30)
	  //	    if(gen_it->isthep()==1)
	  //	    if(gen_it->mother()==1 && (gen_it->idhep()<10 || (gen_it->idhep()>21)))// ||gen_it->daLast()==0))//gives around 124 mrad, with the ==0 it gives 167
	  //	    if(gen_it->mother()==1 && (gen_it->idhep()<10 || (gen_it->idhep()>21)) && gen_it->daLast()==-1)// ||gen_it->daLast()==0))//gives around 124 mrad, with the ==0 it gives 167



	    //	    if( gen_it->daLast()==-1)
	    
	      HepLorentzVector boostedVec(gen_it->PX(),gen_it->PY(),gen_it->PZ(),gen_it->E());


	      //replace with the question if final state particle...
	      if(gen_it->isthep()==1) //isthep==1: stable, ==2, unstable
		{
		  allParticlesNonBoosted.push_back(boostedVec.vect());
		  boostedVec.boost(kinematics::CMBoost);
	      //	      if(gen_it->idhep()>6 &&(gId==ID_GAMMA || gId==ID_PI || gId==ID_K))
		  fjParticles.push_back(PseudoJet(boostedVec.px(),boostedVec.py(),boostedVec.py(),boostedVec.e()));
		  
		  pxSum+=boostedVec.px();
		  pySum+=boostedVec.py();
		  pzSum+=boostedVec.pz();
		  eSum+=boostedVec.e();

		  //all stable
		  //		  if(gen_it->mother()==1 && (gen_it->idhep()<10 || (gen_it->idhep()>20&& gen_it->idhep() < 26)) && gen_it->daLast()==-1)// ||gen_it->daLast()==0))//gives around 124 mrad, with the ==0 it gives 167
		    {
		      allParticlesBoosted.push_back(boostedVec.vect());
		      allPB_E.push_back(gen_it->E());
		    }
		  
		}


	}

      //      cout <<" sum Px: "<< pxSum << " sum Py: " << pySum <<" sumPz: "<< pzSum << " sumPe: "<< eSum<<endl;
      //	cout << lv_q1->px() <<" " <<lv_q1->py() << " " << lv_q1->pz() <<" " << lv_q1->e() <<endl;
      //	cout << lv_q2->px() <<" " <<lv_q2->py() << " " << lv_q2->pz() <<" " << lv_q2->e() <<endl;
      //	cout <<"p Sum: " << m_px << " " << m_py << " " << m_pz << " " << m_e <<endl;
      //	cout <<"pVp Sum: " << m_pxVp << " " << m_pyVp << " " << m_pzVp << " " << m_eVp <<endl;
      //	    cout <<"num Quarks: " <<numQuarks<<endl;

      quarkAngle= lv_q1->angle(*lv_q2);

      HepLorentzVector boostP(m_px,m_py,m_pz,m_e);
      HepLorentzVector boostP2(m_px2,m_py2,m_pz2,m_e2);
      HepLorentzVector boostPVp(m_pxVp,m_pyVp,m_pzVp,m_eVp);

      boostP.boost(kinematics::CMBoost);
      boostP2.boost(kinematics::CMBoost);
      boostPVp.boost(kinematics::CMBoost);

      //	cout <<"p Sum Boosted: " << boostP.px() << " " << boostP.py() << " " << boostP.pz() << " " << boostP.e()<< " " << abs(boostP.pz())+boostP.e()<<endl;
      //	cout <<"p Sum Boosted2: " << boostP2.px() << " " << boostP2.py() << " " << boostP2.pz() << " " << boostP2.e()<< " " << abs(boostP2.pz())+boostP2.e()<<endl;
      //	cout <<"p SumVp Boosted: " << boostPVp.px() << " " << boostPVp.py() << " " << boostPVp.pz() << " " << boostPVp.e()<< " " << abs(boostPVp.pz())+boostPVp.e()<<endl;
      delete lv_q1;
      delete lv_q2;

    
      ////////////////jet algos....



      ////////////
      cmThrust=thrust(allParticlesBoosted.begin(),allParticlesBoosted.end(),gretSelf);
      Thrust t=thrustall(allParticlesBoosted.begin(),allParticlesBoosted.end(),gretSelf);
      cmThrustMag=t.thru;

      labThrust=thrust(allParticlesNonBoosted.begin(),allParticlesNonBoosted.end(),gretSelf);

      ///well, randomizing the direciton is not a good idea, lets set it in the same dir as the original thrust...
      //      if(rand() % 100 <50)
      if(kinematics::thrustDirCM.angle(cmThrust)>pi/2)
	{
	  cmThrust.setZ((-1)*cmThrust.z());
	  cmThrust.setY((-1)*cmThrust.y());
	  cmThrust.setX((-1)*cmThrust.x());
	}
      thetaEThrust=cmThrust.angle(kinematics::firstElectronCM.vect());
      for(int i=0;i<allParticlesBoosted.size();i++)
	{
	  float m_z=2*allPB_E[i]/kinematics::Q;
	  float theta=AuxFunc::getTheta(cmThrust,allParticlesBoosted[i]);
	  //	  m_histos->hEFlowMC->Fill(theta,m_z);
	}
    }
    DebugHistos* m_histos;
  };



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
