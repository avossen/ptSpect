 #ifndef TREESAVER_H
#define TREESAVER_H
#include "ptSpect/mc.h"  //one central place to put the define mc
#include <mdst/mdst.h>
#include MDST_H
#include EVTCLS_H
#include MDST_OBS_H
#include HEPEVT_H
#include TRK_H

#include "belle.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "tuple/BelleTupleManager.h"
#include "ptSpect/EventInfo.h"
#include "ptSpect/AnaConsts.h"
#include "ptSpect/AuxFunc.h"
#include "ptSpect/mc.h"
#include "ptSpect/GenInfo.h"
#include "ptSpect/ParticleInfoMass.h"
#include <math.h>
//#include "AnaDefs.h"
#include HEPEVT_H
#include BELLETDF_H
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
/*class that facilitates the 
saving to root trees
*/
// class HadronQuadruple;

//#define NUM_F_FIELDS 18 //without the data that is only saved for no mc
//#define NUM_F_FIELDS 22 //without the data that is only saved for no mc
//#define NUM_F_FIELDS 22 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
//#define NUM_F_FIELDS 25 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
////-->with the jet stuff, subtract 8
#define NUM_F_FIELDS 16 //without the data that is only saved for no mc. Once we add qT_mc this will go up by one
//6-->10 with the genIds..
#define NUM_I_FIELDS 6



class TreeSaver
{
public:
  TreeSaver()
  {
    initialize();
    tupF=0;
    tupI=0;
  };
  /*from new mdst.cc*/
    const Gen_hepevt &m_get_hepevt(const Mdst_pi0 &pi0) 
      {
	const int PDG_PI0=111;
	const Gen_hepevt &hepevt1 = gen_level(get_hepevt(pi0.gamma(0)));
	const Gen_hepevt &hepevt2 = gen_level(get_hepevt(pi0.gamma(1)));
	if( hepevt1 && hepevt1.mother() && hepevt2 && hepevt2.mother() )
	  {
	    Gen_hepevt &mother1 = hepevt1.mother();
	    if( mother1.get_ID() == hepevt2.mother().get_ID()&& mother1.idhep() == PDG_PI0 ) 
	      {
		return mother1;
	      }
	  }
	return Gen_hepevt_Manager::get_manager().get_NULL();
    }

    static bool sgn(float f)
      {
	return f>=0;
      }


  //gives the pairs that shuld be written to the tree, called from ptSpect
  void fillWPairData(vector<HadronPair*>& vecHadronPairs,EventInfo& evtInfo)
    {
#ifdef MC

      //      cout <<" gi.fillInf" <<endl;
      //do this now in ptSpect.cc for every event, if we save it here or not
      //      gi.fillInf();
      //      cout <<" gi done " <<endl;
      //      if(sgn(gi.cmThrust.z())!=sgn(kinematics::thrustDirCM.z()))


      //below should not be necessary anymore, since we test in GenInfo if the generated thrust is more than pi/2 away
      //      if(kinematics::thrustZReverted)
      //	{
      //	  gi.cmThrust.setZ((-1)*gi.cmThrust.z()); //so that it is the same as in the computed thrust
      //	  gi.cmThrust.setY((-1)*gi.cmThrust.y()); 
      //	  gi.cmThrust.setX((-1)*gi.cmThrust.x()); 
      //	}
#endif
      int numI=-1;
      int numF=-1;
      int qoffset=0;//offset from one quadruple to the next
      //for the number of hadron quad vectors, e.g. 4

	  if(vecHadronPairs.empty())	  //if one is empty all are
	    {
	      return;
	    }
	  //set this higher?
	  if(vecHadronPairs.size() > 300)  // not so much space in array, looks faulty anyways
	    {
	      cout <<"too many entries " <<endl;
	      return;
	    }
	  //	  cout <<"vecHadronPairsSize " << vecHadronPairs.size()  <<endl;
	  for(int i=0;i<vecHadronPairs.size();i++)
	    {
	      fillWSinglePairData(vecHadronPairs[i]);	 
#ifdef MC
	      //this then also calls fillWQuadrupleData, either with 0 (no correspondence) or mcpart =true
	      getQGenInfo(vecHadronPairs[i]);
#endif
	      //	      cout <<"after getQ " <<endl;
	      //	      cout <<"dataf: " << numF <<" i: " << numI <<endl;
	      numF=dataF.size();  //always the same, 
	      numI=dataI.size(); 
	      //	      	      cout <<"dataf: " << numF <<" i: " << numI <<endl;
	      //dataf 44, then 42 for all, i always 12
	      for(int j=0;j<dataF.size();j++)
		{
		  //		  cout <<"j: "<< j <<endl;
		  if(i+qoffset >=1200)
		    {
		      cout <<"index q to large" <<endl;
		      continue;
		    }
		  if(2*j+1 > treeData.size())
		    {
		      cout <<"index td to large" <<endl;
		      continue;
		    }
		  ((float*)treeData[2*j+1])[i+qoffset]=dataF[j];
		}
	      //	      cout <<"/-1" <<endl;
	      //dataI has size 4
	      //	      cout <<"dataI size: "<< dataI.size()<<endl;
	      for(int j=0;j<dataI.size();j++)
		{
		  //		  cout<<"j2: "<< j <<", treeData size: "<< treeData.size()<<endl;
		  if(2*(j+numF)+1 > treeData.size())
		    {
		      cout <<"index td to large" <<endl;
		      continue;
		    }
		  //		  cout <<"wanting to access treedata " << 2*(j+numF)+1 <<", i+qoffset " << i+qoffset <<endl;
		  ((int*)treeData[2*(j+numF)+1])[i+qoffset]=dataI[j];

		}
	      dataF.clear();
	      dataI.clear();
	    }
	  //	  cout <<"-2 " <<endl;
	  qoffset+=vecHadronPairs.size();
	  //save counter info
	  for(int j=0;j<numF;j++)
	    {
	      *(int*)treeData[2*j]=qoffset;
	    }
	  for(int j=0;j<numI;j++)
	    {
	      *(int*)treeData[2*(j+numF)]=qoffset;
	    }
	  //	  cout <<"2.." <<endl;
      if(numF < 0) //no events, all quad vector empty
	return;
      fillWEvtData(evtInfo);
      //so numI, numF 
      saveData(dataF,dataI,2*(numI+numF));

      dataF.clear();
      dataI.clear();
    };


  //gets the corresponding generator info
  void getQGenInfo(HadronPair* pair)
    {

      vector<Particle*> v_p;
      vector<Gen_hepevt*> v_g;

      Gen_hepevt gph1;
      Gen_hepevt gph2;


      gph1=get_hepevt(pair->firstHadron->mdstCharged());
      gph2=get_hepevt(pair->secondHadron->mdstCharged());


      //id 911 is gamma from background addmixture (I think..)

      v_g.push_back(&gph1);
      v_g.push_back(&gph2);

      for(int i=0;i<v_g.size();i++)
	{
	  //no corresponding mc track or gamma from background admixturem_
	  if(*v_g[i]==0 || (v_g[i])->idhep()==911)
	    {
	      fillWSinglePairData(0);
	      return;
	    }
	}

      /////////////////////////////////
      Gen_hepevt pGph1=gph1;
      Gen_hepevt pGph2=gph2;
      //////////////////////////////////////////////////
      bool validType=true;
      for(int i=0;i<v_g.size();i++)
      {
	//	cout <<"making new particle... "<< endl;
	Particle* np=new Particle(*v_g[i]);

	//	cout <<" done... " << np <<endl;
	v_p.push_back(np);
	HepLorentzVector boostedVec(v_g[i]->PX(),v_g[i]->PY(),v_g[i]->PZ(),v_g[i]->E());
	//this is the theta before boost!
	float labTheta=boostedVec.theta();
	float labPhi=boostedVec.phi();

	boostedVec.boost(kinematics::CMBoost);
	v_p[i]->userInfo(*(new ParticleInfo())); //gets deleted in destructor of Particle
	ParticleInfo& pinf=dynamic_cast<ParticleInfo&>(v_p[i]->userInfo());
	pinf.motherGenId=v_g[i]->mother().idhep();

	float m_z=2*boostedVec.t()/kinematics::Q;
	//	cout <<"z is : " << m_z <<endl;
	int geantID=abs(v_g[i]->idhep());
	//	cout <<"looking at index " << i << " idhep: "<< geantID <<endl;

	//this is mc, so the z is independent of the pid (we know the pid). But because the compute function
	//in HadronPair puts the value of z[0] (pion) as the default, we use it here.
	pinf.z[0]=m_z;

	pinf.z[gi.getIdxFromGeantId(geantID)]=m_z;
	pinf.labTheta=labTheta;
	pinf.labPhi=labPhi;

	pinf.labTheta=labTheta;
	pinf.labPhi=labPhi;

	pinf.cmsTheta=boostedVec.theta();
	pinf.cmsPhi=boostedVec.phi();


	Hep3Vector axis=gi.jet1;
	if(boostedVec.vect().dot(gi.jet1)<0)
	  axis=gi.jet2;
	pinf.thrustProj=axis.dot(boostedVec.vect())/(axis.mag()*boostedVec.vect().mag());

	//turns out that the ones where z is -1 are all the electron/muon 
	//compinations
	if(!(geantID==lc_pi0 || geantID==lc_piPlus || geantID==lc_kPlus || geantID ==lc_pPlus))
	  {
	    //	    cout <<"geantID not pion or kaon:: "<< geantID <<endl;
	    validType=false;
	  }
	//	cout <<" valid Type: "<< validType <<endl;
	v_p[i]->momentum().momentum(boostedVec);

      }
      HadronPair hp;
      ///the below (charge, type ) is anyways overridden by the 'compute' functions...
      if(validType)
	{
	  //	  cout <<"tree saver getting type and charge : "<< endl;
	  hp.hadCharge=AuxFunc::getCharge(v_p[0]->pType(),v_p[1]->pType());
	  hp.hadPType=AuxFunc::getPType(v_p[0]->pType(),v_p[1]->pType());
	  //	  cout <<"got charge: "<< hp.hadCharge <<" type: "<< hp.hadPType<<endl;
	}
      else
	{
	  hp.hadCharge=AnaDef::NA;
	  hp.hadPType=AnaDef::UNKNOWN;

	}

      hp.firstHadron=v_p[0];
      hp.secondHadron=v_p[1];
      if(hp.firstHadron->p().vect()==hp.secondHadron->p().vect())
	{
	  //	  cout <<" hadrons are the same..?" <<endl;
	  validType=false;
	}

      float labTheta1=dynamic_cast<ParticleInfo&>(v_p[0]->userInfo()).labTheta;
      float labTheta2=dynamic_cast<ParticleInfo&>(v_p[1]->userInfo()).labTheta;
      float labPhi1=dynamic_cast<ParticleInfo&>(v_p[0]->userInfo()).labPhi;
      float labPhi2=dynamic_cast<ParticleInfo&>(v_p[1]->userInfo()).labPhi;

      float cmsTheta1=dynamic_cast<ParticleInfo&>(v_p[0]->userInfo()).cmsTheta;
      float cmsTheta2=dynamic_cast<ParticleInfo&>(v_p[1]->userInfo()).cmsTheta;
      float cmsPhi1=dynamic_cast<ParticleInfo&>(v_p[0]->userInfo()).cmsPhi;
      float cmsPhi2=dynamic_cast<ParticleInfo&>(v_p[1]->userInfo()).cmsPhi;

      //for mc should be done with mc thrust dir...
      /*      hp1.computeR(kinematics::thrustDirCM);
      hp1.computeThrustTheta(kinematics::thrustDirCM);
      */
      //      cout <<"valid type: "<< validType <<endl;
      //      cout <<"ts 1 "<<endl;

      bool normalHemi=true;
      if(hp.firstHadron->p().vect().dot(gi.jet1)<0)
	normalHemi=false;

      if(validType)
	{
	  if(normalHemi)
	    hp.compute();
	  else
	    hp.compute();
	}


      if(!validType)
	{
	  hp.z1=-1;
	  hp.z2=-1;
	  hp.z=-1;
	}

      fillWSinglePairData(&hp,true);

      for(int i=0;i<v_p.size();i++)
	{
	  delete v_p[i];
	}
    };

  //get all h-pairs in mc, w/o det acceptance



  //data specific to the event, so no arrays
  void fillWEvtData(EventInfo& evtInfo)
    {
      //      dataF.push_back(-log(tan(kinematics::thrustDirCM.theta()/2)));
      dataF.push_back(kinematics::thrustMag);
      dataF.push_back(kinematics::E_miss);
      dataF.push_back(kinematics::thrustDirCM.theta());
      dataF.push_back(kinematics::thrustDirCM.phi());
      dataF.push_back(kinematics::thrustDirLab.theta());
      dataF.push_back(kinematics::thetaEThrust);
      dataF.push_back(kinematics::jetFluffiness1);
      dataF.push_back(kinematics::jetFluffiness2);

      dataF.push_back(kinematics::jetE1);
      dataF.push_back(kinematics::jetE2);

      dataF.push_back(kinematics::jet1.phi());
      dataF.push_back(kinematics::jet2.phi());

      dataF.push_back(kinematics::jet1.theta());
      dataF.push_back(kinematics::jet2.theta());

#ifdef MC
      float angleToRecThrust=gi.jet1.angle(kinematics::jet1);


      Hep3Vector tmpThrust=gi.cmThrust;      
      if(angleToRecThrust>pi/2)
	{
	  angleToRecThrust=pi-angleToRecThrust;
	  //save to asume that we didn't do the flip
	  tmpThrust.setZ((-1)*tmpThrust.z());
	  tmpThrust.setY((-1)*tmpThrust.y());
	  tmpThrust.setX((-1)*tmpThrust.x());

	}
      float thetaToRecThrust=tmpThrust.theta()-kinematics::thrustDirCM.theta();
      float phiToRecThrust=tmpThrust.phi()-kinematics::thrustDirCM.phi();

      dataF.push_back(tmpThrust.mag());
      dataF.push_back(0.0);
      dataF.push_back(angleToRecThrust);
      dataF.push_back(thetaToRecThrust);
      dataF.push_back(phiToRecThrust);

      //------
      dataF.push_back(gi.fluffiness1);
      dataF.push_back(gi.fluffiness2);

      dataF.push_back(gi.jetE1);
      dataF.push_back(gi.jetE2);

      dataF.push_back(gi.jetPhi1);
      dataF.push_back(gi.jetPhi2);

      dataF.push_back(gi.jetTheta1);
      dataF.push_back(gi.jetTheta2);
      ///-----


      dataF.push_back(gi.vpEnergy);
      dataF.push_back(gi.vpPx);
      dataF.push_back(gi.vpPy);
      dataF.push_back(gi.vpPz);
      dataF.push_back(gi.quarkAngle);
      dataF.push_back(gi.thetaEThrust);
      dataF.push_back(gi.isrPhotonEnergy);
      dataI.push_back(gi.numQuarks);
      dataI.push_back(kinematics::DDecayMC);
      dataI.push_back(kinematics::DStarDecayMC);
      //      cout <<"isr photon energy: " << gi.isrPhotonEnergy <<endl;

#endif
      dataI.push_back(kinematics::runNr);
      dataI.push_back(kinematics::evtNr);
      dataI.push_back(kinematics::jetNumParts1);
      dataI.push_back(kinematics::jetNumParts2);
      dataI.push_back(kinematics::D0Tag);
      dataI.push_back(kinematics::DStarTag);
      dataI.push_back(kinematics::DDecay);
      //           if(kinematics::DStarDecay==2){
	     //      cout <<"tree saver saw dStar decay!!" <<endl<<endl;}
      dataI.push_back(kinematics::DStarDecay);

    };

  void fillWSinglePairData(HadronPair* pair, bool mcPart=false)
    {
//this will only be zero for the mc part, that means we NUM_F_FIELDS should not include fields that are not there in mc (e.g. pion mass)
//--->only include fields that are in the mc part...
      if(pair==0)
	{
	  for(int i=0;i<NUM_F_FIELDS;i++)
	    {
	      dataF.push_back(-1);
	    }
	  for(int i=0;i<NUM_I_FIELDS;i++)
	    {
	      dataI.push_back(-1);
	    }
	}
      else
	{
	  dataF.push_back(pair->z1);//z1
	  dataF.push_back(pair->z2); //z2

	  if(!mcPart)
	    {
	      dataF.push_back(pair->p_PiPi);
	      dataF.push_back(pair->p_PiK);	  
	      dataF.push_back(pair->p_PiP);


	      dataF.push_back(pair->p_KPi);
	      dataF.push_back(pair->p_KK);	  
	      dataF.push_back(pair->p_KP);
	      
	      
	      dataF.push_back(pair->p_PPi);
	      dataF.push_back(pair->p_PK);	  
	      dataF.push_back(pair->p_PP);
	      
	      dataF.push_back(pair->kT_PiPi);
	      dataF.push_back(pair->kT_PiK);	  
	      dataF.push_back(pair->kT_PiP);


	      dataF.push_back(pair->kT_KPi);
	      dataF.push_back(pair->kT_KK);	  
	      dataF.push_back(pair->kT_KP);
	      
	      
	      dataF.push_back(pair->kT_PPi);
	      dataF.push_back(pair->kT_PK);	  
	      dataF.push_back(pair->kT_PP);


	      dataF.push_back(pair->qT_PiPi);
	      dataF.push_back(pair->qT_PiK);	  
	      dataF.push_back(pair->qT_PiP);


	      dataF.push_back(pair->qT_KPi);
	      dataF.push_back(pair->qT_KK);	  
	      dataF.push_back(pair->qT_KP);
	      
	      
	      dataF.push_back(pair->qT_PPi);
	      dataF.push_back(pair->qT_PK);	  
	      dataF.push_back(pair->qT_PP);


	      
	      
	      dataF.push_back(pair->z1_PiPi);
	      dataF.push_back(pair->z1_PiK);
	      dataF.push_back(pair->z1_PiP);


	      dataF.push_back(pair->z1_KPi);
	      dataF.push_back(pair->z1_KK);
	      dataF.push_back(pair->z1_KP);

	      dataF.push_back(pair->z1_PPi);
	      dataF.push_back(pair->z1_PK);
	      dataF.push_back(pair->z1_PP);

	      dataF.push_back(pair->z2_PiPi);
	      dataF.push_back(pair->z2_PiK);
	      dataF.push_back(pair->z2_PiP);

	      dataF.push_back(pair->z2_KPi);
	      dataF.push_back(pair->z2_KK);
	      dataF.push_back(pair->z2_KP);

	      dataF.push_back(pair->z2_PPi);
	      dataF.push_back(pair->z2_PK);
	      dataF.push_back(pair->z2_PP);
	    }

	  float labTheta1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).labTheta;
	  float labTheta2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).labTheta;

	  float labPhi1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).labPhi;
	  float labPhi2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).labPhi;

	  float cmsTheta1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).cmsTheta;
	  float cmsTheta2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).cmsTheta;

	  float cmsPhi1=dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).cmsPhi;
	  float cmsPhi2=dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).cmsPhi;
	

	  if(labTheta1==0 || labTheta2==0 )
	    {
	      //	      cout <<"theta1: " << theta1 << " , theta2:  " << theta2 << " , theta3: " << theta3 << " theta4: " << theta4 <<endl;
	      //4 times theta, 4 times phi
	      for(int i=0;i<8;i++)
		{
		  dataF.push_back(0);
		}
	    }
	  else
	    {

	      dataF.push_back(labTheta1);
	      dataF.push_back(labTheta2);

	      dataF.push_back(labPhi1);
	      dataF.push_back(labPhi2);

	      dataF.push_back(cmsTheta1);
	      dataF.push_back(cmsTheta2);

	      dataF.push_back(cmsPhi1);
	      dataF.push_back(cmsPhi2);

	    }
	
	  dataF.push_back(dynamic_cast<ParticleInfo&>(pair->firstHadron->userInfo()).thrustProj);
	  dataF.push_back(dynamic_cast<ParticleInfo&>(pair->secondHadron->userInfo()).thrustProj);

	  dataF.push_back(pair->kT);
	  dataF.push_back(pair->diffTheta);
	  dataF.push_back(pair->diffPhi);
			  
			 	  
	  //haven't added qT_mc yet, so don't include if mcPart
	  //	  if(!mcPart)// added now
	    dataF.push_back(pair->qT);


	  //	  cout <<"charge: " << vecHQ->hadCharge <<endl;


	  dataI.push_back(pair->hadCharge);      
	  dataI.push_back(pair->hadPType);//particle type

	  dataI.push_back(pair->hadCharge1);      
	  dataI.push_back(pair->hadPType1);

	  dataI.push_back(pair->hadCharge2);      
	  dataI.push_back(pair->hadPType2);

	}
			  
    };
 //std: float datatype
  void addFieldF(char* fieldname)
  {
    //construct the memory location from which the tree should read the new data field
    float* memLoc=new float;
    treeData.push_back(memLoc);
    pDataTree->Branch(fieldname, memLoc, (fieldname+string("/F")).c_str());
    fieldNamesF.push_back(fieldname);
  };


  void addFieldI(char* fieldname)
  {
    int* memLoc=new int;
    treeData.push_back(memLoc);
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("/I")).c_str());
    fieldNamesI.push_back(fieldname);
  };

  void addArrayI(char* fieldname)
  {
    //standard lenth, shouldn't be more than that
    int* memLoc=new int[1200];
    int* memLocCounter=new int;
    treeData.push_back(memLocCounter);
    treeData.push_back(memLoc);
    string counterName=string(fieldname)+string("Counter");
    pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/I")).c_str());
  };

  void addArrayF(char* fieldname)
  {
    float* memLoc=new float[1200];
    int* memLocCounter=new int;
    treeData.push_back(memLocCounter);
    treeData.push_back(memLoc);
    string counterName=string(fieldname)+string("Counter");
    pDataTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
    pDataTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
  };

  void addArrayPi0F(char* fieldname)
    {
      float* memLoc=new float[1200];
      int* memLocCounter=new int;

      float* memLoc2=new float[1200];
      int* memLocCounter2=new int;
      treeDataFalsePi0.push_back(memLocCounter);
      treeDataFalsePi0.push_back(memLoc);

      treeDataRealPi0.push_back(memLocCounter2);
      treeDataRealPi0.push_back(memLoc2);
      string counterName=string(fieldname)+string("Counter");

      pRealPi0Tree->Branch(counterName.c_str(),memLocCounter2,(counterName+string("/I")).c_str());
      pRealPi0Tree->Branch(fieldname,memLoc2,(fieldname+string("[")+counterName+string("]/F")).c_str());

      pFalsePi0Tree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
      pFalsePi0Tree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());
    }

  void addArrayPi0AsymmetryF(char* fieldname)
    {
      float* memLoc=new float[1200];
      int* memLocCounter=new int;
      float* memLoc2=new float[1200];
      int* memLocCounter2=new int;

      treeDataRealPi0Asymmetry.push_back(memLocCounter);
      treeDataRealPi0Asymmetry.push_back(memLoc);

      treeDataFalsePi0Asymmetry.push_back(memLocCounter2);
      treeDataFalsePi0Asymmetry.push_back(memLoc2);

      string counterName=string(fieldname)+string("Counter");
      pRealPi0AsymmetryTree->Branch(counterName.c_str(),memLocCounter,(counterName+string("/I")).c_str());
      pRealPi0AsymmetryTree->Branch(fieldname,memLoc,(fieldname+string("[")+counterName+string("]/F")).c_str());

      pFalsePi0AsymmetryTree->Branch(counterName.c_str(),memLocCounter2,(counterName+string("/I")).c_str());
      pFalsePi0AsymmetryTree->Branch(fieldname,memLoc2,(fieldname+string("[")+counterName+string("]/F")).c_str());

    }


  void saveGenInfo()
  {
#ifdef MC
    gi.fillInf();
    gi.doAll();
#endif
  }

  //my guess as to what is going on here:
  //  dataF, dataI is the event only data (these arrays where deleted before handed to the event-fill function), offset is where the single fields start
  void saveData(vector<float>& dataF, vector<int>& dataI, int offset)
  {
    int dataSize=dataF.size()+dataI.size();
#ifdef MC
    if(dataSize !=treeData.size()-offset)
#else
    if(dataSize !=treeData.size()-offset)
#endif
      {
	/////saving in paw format...  outdated since all the arrays where added...
	/*	if(tupI!=0)
	  {
	    for(int i=0;i<fieldNamesI.size();i++)
	      {
		tupI->column(fieldNamesI[i].c_str(),dataI[i]);
	      }
	    tupI->dumpData();
	  }
	if(tupF!=0)
	  {
	    for(int i=0;i<fieldNamesF.size();i++)
	      {
		tupF->column(fieldNamesF[i].c_str(),dataF[i]);
	      }
	    tupF->dumpData();
	    }*/ //does not work anymore with arrays  :-(
	      /////end save in paw
	cout << "data size does not match number of branches " <<dataSize <<"+ " << offset << " = " <<treeData.size() <<endl <<flush;
	//10+108=122
	exit(0);
      }
    //	cout << "data size does  match number of branches " <<dataSize <<"+ " << offset << " = " <<treeData.size() <<endl <<flush;
    //10+112
    //this works because dataF was cleared before event specific things there saved, 
    for(int i=0;i<dataF.size();i++)
      {
	*(float*)treeData[i+offset]=dataF[i];
      }
    for(int i=0;i<dataI.size();i++)
      {
	*(int*)treeData[i+dataF.size()+offset]=dataI[i];
      }
    //    cout <<"before mc" <<endl;
#ifdef MC
    //    cout <<"doAll" <<endl;
    //maybe this should be done independently from the question if we save an event within the acceptance or not
    //now extra function save gen info, so it is done independently
    //        gi.doAll();
#else

#endif 
    /*    cout <<"pt counter: " <<  *(int*)treeData[50] <<endl;
    for(int i=0 ;i< *(int*)treeData[50] ;i++)
      {
	cout <<"pt: " <<  ((float*)treeData[51])[i] <<endl;
	}*/

    //        cout <<"filling " << endl;



          pDataTree->Fill();

  };

  void setDebugHistos(DebugHistos* d_histos)
  {
    m_histos=d_histos;
    gi.setDebugHistos(d_histos);
  };



private:

  template<class T> void fillTreeEntry(int& counter,T entry)
  {
    counter++;
    *(T*)treeData[counter]=entry;
  }
  template<class T> void fillTreeArrEntry(int& counter,vector<T>& vec)
  {
    counter++;
    *(int*)treeData[counter]=vec.size();
    counter++;
    //cp data from vec in tree array
    for(int i=0;i<vec.size();i++)
      {
	((T*)treeData[counter])[i]=vec[i];
      }

  }

  void initialize()
  {
    if(initialized)
      return;
    //the first time the class is initialized, construct the tree
    pDataTree=new TTree("DataTree","My Transversity Data Tree");
    pRealPi0Tree=new TTree("RealPi0Tree","Real Pi0 Data Tree");
    pFalsePi0Tree=new TTree("FalsePi0Tree","False Pi0 Data Tree");
    pRealPi0AsymmetryTree=new TTree("RealPi0TreeAsymmetry","Real Pi0 Asymmetry Data Tree");
    pFalsePi0AsymmetryTree=new TTree("FalsePi0TreeAsymmetry","False Pi0 Asymmetry Data Tree");
    initialized=true;
#ifdef MC
    gi.initializeTree();
#endif
  };


  static GenInfo gi;
  static TTree* pDataTree;
  static TTree* pRealPi0Tree;
  static TTree* pRealPi0AsymmetryTree;
  static TTree* pFalsePi0Tree;
  static TTree* pFalsePi0AsymmetryTree;
  static bool initialized;
  static vector<float> dataF;
  static vector<int> dataI;

  static vector<float> realPi0_gammaE;
  static vector<float> falsePi0_gammaE;

  static vector<float> realPi0_e9oe25;;
  static vector<float> falsePi0_e9oe25;

  static vector<float> realPi0_mass;
  static vector<float> falsePi0_mass;

  static vector<float> realPi0_gammaAsymmetry;
  static vector<float> falsePi0_gammaAsymmetry;

  //the adresses from which the tree should read its data
  static vector<void*> treeData;
  static vector<void*> treeDataRealPi0;
  static vector<void*> treeDataRealPi0Asymmetry;
  static vector<void*> treeDataFalsePi0;
  static vector<void*> treeDataFalsePi0Asymmetry;
  static vector<string> fieldNamesI;
  static vector<string> fieldNamesF;
  static BelleTuple* tupF;
  static BelleTuple* tupI;
  static DebugHistos* m_histos;
};
#if defined(BELLE_NAMESPACE)
}
#endif
#endif
