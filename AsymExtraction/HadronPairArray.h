#define DEBUG_EVENT 60233

#ifndef HADRON_PAIR_ARRAY_H
#define HADRON_PAIR_ARRAY_H
#include "ReaderBase.h"
#include "TMath.h"
//#define USE_QT
//300 is also the number in the Treesaver, so cannot be larger...
#define Max_ArrSize 300

struct HadronPairArray:public ReaderBase
{
  //cuts:

  float kTCut;
  const float zCut;
  const float zUpperCut;
  const float secondZCut;
  const float hadronTagFiducialCut;
  const bool asymmetryFlag;

  int numPairs;

  float z1[Max_ArrSize];
  float z2[Max_ArrSize];


  float p_PiPi[Max_ArrSize];
  float p_PiK[Max_ArrSize];

  float p_PiP[Max_ArrSize];

  float p_KPi[Max_ArrSize];
  float p_KK[Max_ArrSize];
  float p_KP[Max_ArrSize];

  float p_PPi[Max_ArrSize];
  float p_PK[Max_ArrSize];
  float p_PP[Max_ArrSize];

  float p_PiPi2[Max_ArrSize];
  float p_PiK2[Max_ArrSize];

  float p_PiP2[Max_ArrSize];

  float p_KPi2[Max_ArrSize];
  float p_KK2[Max_ArrSize];
  float p_KP2[Max_ArrSize];

  float p_PPi2[Max_ArrSize];
  float p_PK2[Max_ArrSize];
  float p_PP2[Max_ArrSize];


  float ep_PiPi[Max_ArrSize];
  float ep_PiK[Max_ArrSize];

  float ep_PiP[Max_ArrSize];

  float ep_KPi[Max_ArrSize];
  float ep_KK[Max_ArrSize];
  float ep_KP[Max_ArrSize];

  float ep_PPi[Max_ArrSize];
  float ep_PK[Max_ArrSize];
  float ep_PP[Max_ArrSize];




  float kT_PiPi[Max_ArrSize];
  float kT_PiK[Max_ArrSize];


  float kT_PiP[Max_ArrSize];
  float kT_KPi[Max_ArrSize];

  float kT_KK[Max_ArrSize];
  float kT_KP[Max_ArrSize];
  float kT_PPi[Max_ArrSize];
  float kT_PK[Max_ArrSize];
  float kT_PP[Max_ArrSize];


  float qT_PiPi[Max_ArrSize];
  float qT_PiK[Max_ArrSize];
  float qT_PiP[Max_ArrSize];

  float qT_KPi[Max_ArrSize];
  float qT_KK[Max_ArrSize];
  float qT_KP[Max_ArrSize];
  float qT_PPi[Max_ArrSize];
  float qT_PK[Max_ArrSize];
  float qT_PP[Max_ArrSize];




  float z1_PiPi[Max_ArrSize];
  float z1_PiK[Max_ArrSize];
  float z1_PiP[Max_ArrSize];


  float z1_KPi[Max_ArrSize];
  float z1_KK[Max_ArrSize];
  float z1_KP[Max_ArrSize];

  float z1_PPi[Max_ArrSize];
  float z1_PK[Max_ArrSize];
  float z1_PP[Max_ArrSize];

  float z2_PiPi[Max_ArrSize];
  float z2_PiK[Max_ArrSize];
  float z2_PiP[Max_ArrSize];


  float z2_KPi[Max_ArrSize];
  float z2_KK[Max_ArrSize];
  float z2_KP[Max_ArrSize];


  float z2_PPi[Max_ArrSize];
  float z2_PK[Max_ArrSize];
  float z2_PP[Max_ArrSize];



  float labTheta1[Max_ArrSize];
  float labTheta2[Max_ArrSize];

  float cmsTheta1[Max_ArrSize];
  float cmsTheta2[Max_ArrSize];

  float labPhi1[Max_ArrSize];
  float labPhi2[Max_ArrSize];

  float cmsPhi1[Max_ArrSize];
  float cmsPhi2[Max_ArrSize];

  //sum over all...

  float thrustProj1[Max_ArrSize];
  float thrustProj2[Max_ArrSize];

  TH2D* kinematicSmearingMatrix;


  float kT[Max_ArrSize];

  float qT[Max_ArrSize];

  float hadDiffTheta[Max_ArrSize];
  float hadDiffPhi[Max_ArrSize];


  int chargeType[Max_ArrSize];
  int particleType[Max_ArrSize];


  int chargeType1[Max_ArrSize];
  int particleType1[Max_ArrSize];

  int chargeType2[Max_ArrSize];
  int particleType2[Max_ArrSize];


  int cut[Max_ArrSize];

  // pid dependent z cut is then in MultiPlotter because it has to be applied in each event for each weight differently
  //considering the thrust axis resolution of 0.16+-0.09 rad, a max opening cut of 0.99 is even too large...
  //the pid dependent cuts (z, kT) should be done for the correct hypothesis

 HadronPairArray(TChain* chain, int MCFlag=mcFlagNone):ReaderBase(MCFlag), zCut(0.0),zUpperCut(1.1), secondZCut(0.0), hadronTagFiducialCut(70.2), asymmetryFlag(false),kTCut(100.0)//kTCut(5.31145668)
  {

#ifdef USE_QT
    //qT can have up to sqrt(s), and it can be even bigger if we assume the wrong particles...
    kTCut=35.52;
#endif


    //no chain implies standalone. Cannot branch on the same field twice, this would override
    //    cout <<" do we have a chain? " << chain<<endl;
    if(chain)
      {
	//cout <<"yes " <<endl;
	myChain=chain;
	branchPointersI.push_back(&numPairs);
	branchPointers.push_back(z1);
	branchPointers.push_back(z2);
	//if we just test that it is not mcFlagWoA, the mc pair branches on the same fields, so the data field only gets 0..
	if(mMCFlag==mcFlagNone)
	  {
	    branchPointers.push_back(z1_PiPi);
	    branchPointers.push_back(z1_PiK);
	    branchPointers.push_back(z1_PiP);

	    branchPointers.push_back(z1_KPi);
	    branchPointers.push_back(z1_KK);
	    branchPointers.push_back(z1_KP);

	    branchPointers.push_back(z1_PPi);
	    branchPointers.push_back(z1_PK);
	    branchPointers.push_back(z1_PP);

	    branchPointers.push_back(z2_PiPi);
	    branchPointers.push_back(z2_PiK);
	    branchPointers.push_back(z2_PiP);

	    branchPointers.push_back(z2_KPi);
	    branchPointers.push_back(z2_KK);
	    branchPointers.push_back(z2_KP);

	    branchPointers.push_back(z2_PPi);
	    branchPointers.push_back(z2_PK);
	    branchPointers.push_back(z2_PP);
	    
	    branchPointers.push_back(p_PiPi);
	    branchPointers.push_back(p_PiK);
	    branchPointers.push_back(p_PiP);

	    branchPointers.push_back(p_KPi);
	    branchPointers.push_back(p_KK);
	    branchPointers.push_back(p_KP);

	    branchPointers.push_back(p_PPi);
	    branchPointers.push_back(p_PK);
	    branchPointers.push_back(p_PP);

	    branchPointers.push_back(p_PiPi2);
	    branchPointers.push_back(p_PiK2);
	    branchPointers.push_back(p_PiP2);

	    branchPointers.push_back(p_KPi2);
	    branchPointers.push_back(p_KK2);
	    branchPointers.push_back(p_KP2);

	    branchPointers.push_back(p_PPi2);
	    branchPointers.push_back(p_PK2);
	    branchPointers.push_back(p_PP2);

	    branchPointers.push_back(ep_PiPi);
	    branchPointers.push_back(ep_PiK);
	    branchPointers.push_back(ep_PiP);

	    branchPointers.push_back(ep_KPi);
	    branchPointers.push_back(ep_KK);
	    branchPointers.push_back(ep_KP);

	    branchPointers.push_back(ep_PPi);
	    branchPointers.push_back(ep_PK);
	    branchPointers.push_back(ep_PP);



	    //keep the same address, so we don't have to change
	    //the kt references lateron
	    branchPointers.push_back(kT_PiPi);
	    branchPointers.push_back(kT_PiK);
	    branchPointers.push_back(kT_PiP);
	    
	    branchPointers.push_back(kT_KPi);
	    branchPointers.push_back(kT_KK);
	    branchPointers.push_back(kT_KP);

	    branchPointers.push_back(kT_PPi);
	    branchPointers.push_back(kT_PK);
	    branchPointers.push_back(kT_PP);



	    //different solution: replace kT with qT ifdef USE_QT
//
//	    branchPointers.push_back(qT_PiPi);
//	    branchPointers.push_back(qT_PiK);
//	    branchPointers.push_back(qT_PiP);
//	    
//	    branchPointers.push_back(qT_KPi);
//	    branchPointers.push_back(qT_KK);
//	    branchPointers.push_back(qT_KP);
//
//	    branchPointers.push_back(qT_PPi);
//	    branchPointers.push_back(qT_PK);
//	    branchPointers.push_back(qT_PP);
//
//

	  }


	branchPointers.push_back(labTheta1);
	branchPointers.push_back(labTheta2);

	branchPointers.push_back(labPhi1);
	branchPointers.push_back(labPhi2);
	
	branchPointers.push_back(cmsTheta1);
	branchPointers.push_back(cmsTheta2);

	branchPointers.push_back(cmsPhi1);
	branchPointers.push_back(cmsPhi2);

	branchPointers.push_back(thrustProj1);
	branchPointers.push_back(thrustProj2);

	//if use qt instead of kt, switch them so that the following kt code works on qt instead




	//only for mc...
	//#ifdef USE_QT
	//	branchPointers.push_back(qT);
	//#else
	branchPointers.push_back(kT);
	//#endif

	branchPointers.push_back(hadDiffTheta);
	branchPointers.push_back(hadDiffPhi);

	branchPointersI.push_back(chargeType1);
	branchPointersI.push_back(particleType1);


	branchPointersI.push_back(chargeType2);
	branchPointersI.push_back(particleType2);



	branchPointersI.push_back(chargeType);
	branchPointersI.push_back(particleType);

	if(mMCFlag!=mcFlagNone)
	  {
	  }
	
	branchNamesI.push_back("z1"+addendum+"Counter");
	branchNames.push_back("z1"+addendum);
	branchNames.push_back("z2"+addendum);
	//if we just test that it is not mcFlagWoA, the mc pair branches on the same fields, so the data field only gets 0..
	if(mMCFlag==mcFlagNone)
	  {
	    cout <<"mcflagnone " <<endl;
	    branchNames.push_back("z1_PiPi");
	    branchNames.push_back("z1_PiK");
	    branchNames.push_back("z1_PiP");

	    branchNames.push_back("z1_KPi");
	    branchNames.push_back("z1_KK");
	    branchNames.push_back("z1_KP");

	    branchNames.push_back("z1_PPi");
	    branchNames.push_back("z1_PK");
	    branchNames.push_back("z1_PP");
	    
	    branchNames.push_back("z2_PiPi");
	    branchNames.push_back("z2_PiK");
	    branchNames.push_back("z2_PiP");
	    
	    branchNames.push_back("z2_KPi");
	    branchNames.push_back("z2_KK");
	    branchNames.push_back("z2_KP");
	    
	    branchNames.push_back("z2_PPi");
	    branchNames.push_back("z2_PK");
	    branchNames.push_back("z2_PP");
	    
	    branchNames.push_back("p_PiPi");
	    branchNames.push_back("p_PiK");
	    branchNames.push_back("p_PiP");

	    branchNames.push_back("p_KPi");
	    branchNames.push_back("p_KK");
	    branchNames.push_back("p_KP");
	    
	    branchNames.push_back("p_PPi");
	    branchNames.push_back("p_PK");
	    branchNames.push_back("p_PP");

	    branchNames.push_back("p_PiPi_2");
	    branchNames.push_back("p_PiK_2");
	    branchNames.push_back("p_PiP_2");

	    branchNames.push_back("p_KPi_2");
	    branchNames.push_back("p_KK_2");
	    branchNames.push_back("p_KP_2");
	    
	    branchNames.push_back("p_PPi_2");
	    branchNames.push_back("p_PK_2");
	    branchNames.push_back("p_PP_2");

	    branchNames.push_back("ep_PiPi");
	    branchNames.push_back("ep_PiK");
	    branchNames.push_back("ep_PiP");

	    branchNames.push_back("ep_KPi");
	    branchNames.push_back("ep_KK");
	    branchNames.push_back("ep_KP");
	    
	    branchNames.push_back("ep_PPi");
	    branchNames.push_back("ep_PK");
	    branchNames.push_back("ep_PP");

#ifdef USE_QT
	    branchNames.push_back("qT_PiPi");
	    branchNames.push_back("qT_PiK");
	    branchNames.push_back("qT_PiP");
	    
	    branchNames.push_back("qT_KPi");
	    branchNames.push_back("qT_KK");
	    branchNames.push_back("qT_KP");
	    
	    branchNames.push_back("qT_PPi");
	    branchNames.push_back("qT_PK");
	    branchNames.push_back("qT_PP");

#else	    	    
	    branchNames.push_back("kT_PiPi");
	    branchNames.push_back("kT_PiK");
	    branchNames.push_back("kT_PiP");
	    
	    branchNames.push_back("kT_KPi");
	    branchNames.push_back("kT_KK");
	    branchNames.push_back("kT_KP");
	    
	    branchNames.push_back("kT_PPi");
	    branchNames.push_back("kT_PK");
	    branchNames.push_back("kT_PP");
#endif
	  }

	branchNames.push_back("labTheta1"+addendum);
	branchNames.push_back("labTheta2"+addendum);

	branchNames.push_back("labPhi2"+addendum);
	branchNames.push_back("labPhi2"+addendum);

	branchNames.push_back("cmsTheta1"+addendum);
	branchNames.push_back("cmsTheta2"+addendum);
	branchNames.push_back("cmsPhi2"+addendum);
	branchNames.push_back("cmsPhi2"+addendum);

	branchNames.push_back("thrustProj1"+addendum);
	branchNames.push_back("thrustProj2"+addendum);

	//qt name exists as _mc and _mcWoa (unlike kT which is called just kt in the woa tree)
#ifdef USE_QT
	  branchNames.push_back("qT"+addendum);
#else
	    branchNames.push_back("kT"+addendum);
#endif
	//only for non mc...and woa (so no addendum in the other cases..)
	///-->should be changed now, since we added qT_mc
	//	if(mMCFlag==mcFlagWoA)

	  //#ifndef USE_QT
	  //	  branchNames.push_back("kT"+addendum);
	  //#endif

	  //	else
	  //	  branchNames.push_back("qT");


	branchNames.push_back("HadDiffTheta"+addendum);
	branchNames.push_back("HadDiffPhi"+addendum);
	//pi0 mass doesn't exist as mc
	//		branchNames.push_back("pi0Mass"+sAdd+"1"+addendum);
		//	branchNames.push_back("pi0Mass"+sAdd+"2"+addendum);
	
	branchNamesI.push_back("chargeType1"+addendum);
	branchNamesI.push_back("particleType1"+addendum);

	branchNamesI.push_back("chargeType2"+addendum);
	branchNamesI.push_back("particleType2"+addendum);

	branchNamesI.push_back("chargeType"+addendum);
	branchNamesI.push_back("particleType"+addendum);

	//unfortunately the naming convention was not kept...
	if(mMCFlag==mcFlagWoA)
	  {
	    //	    branchNamesI.push_back("motherGenID"+sAdd+"1"+addendum);
	    //	    branchNamesI.push_back("motherGenID"+sAdd+"2"+addendum);
	  }
	if(mMCFlag==mcFlagMC)
	  {

	  }
	doAllBranching();
      }

    // kinematicSmearingMatrix=new TH2D("kinematicSmearingMatrix","kinematicSmearingMatrix",);

  }

  void afterFill(int evtNr=0, bool print=false)
  {
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"we see " <<numPairs<<endl;
      }
    for(int i=0;i<numPairs;i++)
      {
	cut[i]=0;
	if(isnan(ep_PiPi[i]))
	  {
	    cout <<"eppipi nan" <<endl;
	    cut[i]=1;
	  }

    if(isnan(ep_PiK[i]))
      {
      cout <<"eppik nan" <<endl;
	    cut[i]=1;
	  }


    if(isnan(ep_KPi[i]))
      {
      cout <<"epKpi nan" <<endl;
	    cut[i]=1;
	  }

    if(isnan(ep_KK[i]))
      {
      cout <<"epKK nan" <<endl;
	    cut[i]=1;
	  }


    if(isnan(ep_PiP[i]))
      {
      cout <<"epPiP nan" <<endl;
	    cut[i]=1;
	  }


    if(isnan(ep_PPi[i]))
      {
      cout <<"epPPi nan" <<endl;
	    cut[i]=1;
	  }


    if(isnan(ep_KP[i]))
      {

      cout <<"epKP nan" <<endl;
	    cut[i]=1;
	  }

    if(isnan(ep_PK[i]))
      {
      cout <<"epPK nan" <<endl;
	    cut[i]=1;
	  }



	//	cout <<"looking at p_pipi: "<< p_PiPi[i] << " and " << p_PiPi2[i]<<" e: "<< ep_PiPi[i]<<endl;

    //for the benefit of the xcheck, debug
    if(z1[i]< z2[i])
      {
	float tmp=z2[i];
	z2[i]=z1[i];
	z1[i]=tmp;
      }

    if(z1_PiPi[i]< z2_PiPi[i])
      {
	float tmp=z2_PiPi[i];
	z2_PiPi[i]=z1_PiPi[i];
	z1_PiPi[i]=tmp;
      }

    //    cout <<" kt PiPi: "<< kT_PiPi[i] <<endl;
    //    cout <<"chargeType " << i << ":  " << chargeType[i] << " particle type ; "<< particleType[i] <<endl;
	//	cout <<"hp after fill!" <<endl;


	if(fabs(cmsTheta2[i]-TMath::Pi()/2)>hadronTagFiducialCut)
	  {
	if(print)
	  {
	    cout <<"hadron tag fiducial.." <<endl;
	  }
	    cut[i]=1;
	if(p_PiPi[i]>0 && chargeType[i]==0)
	  {
	    cout <<"cut on fid" <<endl;

	  }
	if(evtNr==DEBUG_EVENT)
	  {
	    cout <<"hadron fid " <<endl;
	  }
	
	cout <<"cut hadon fid cut " <<endl;
	  }
	if(kT_PiPi[i]>kTCut || kT_PiK[i]>kTCut || kT_PiP[i]>kTCut || kT_KPi[i]>kTCut || kT_KK[i]>kTCut || kT_KP[i]>kTCut || kT_PPi[i]>kTCut || kT_PK[i]>kTCut || kT_PPi[i]>kTCut)
	  {
#ifdef USE_QT
	    cout <<"qt defined " <<endl;
#endif
	    //	     cout <<"kt cut, : "<< kT_PiPi[i] <<" pik: "<< kT_PiK[i] <<" kpi: "<< kT_KPi[i] << " KK: " << kT_KK[i] <<" PiP: " << kT_PiP[i] << " KP: "<< kT_KP[i] <<" PPi: "<< kT_PPi[i];
	    //	    cout <<" PK: " << kT_PK[i] << " PP: " << kT_PP[i] <<endl;
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"kt cut " <<endl;
      }


	if(p_PiPi[i]>0 && chargeType[i]==0)
	  {
	    //	  cout <<"cut on kt" <<endl;
	  }
	if(print)
	  {
	    cout <<"cut due to kt cut .." <<endl;
	  }
	    cut[i]=1;
	  }

	////charlotte cuts..
	//	if(cos(labTheta1[i])<-0.511 || cos(labTheta1[i])>=0.842)
	//	       cut[i]=1;

	//	if(cos(labTheta2[i])<-0.511 || cos(labTheta2[i])>=0.842)
	//	       cut[i]=1;

	//	if(cut[i]==0)
	//	  cout <<"hadron pair survived cut, lab1: "<< cos(labTheta1[i]) <<" 2: "<< cos(labTheta2[i])<<endl;
	///

	bool asymFlag=false;
	if(cmsTheta1[i]>TMath::Pi()/2)
	  {
	    if(cmsTheta2[i] < TMath::Pi()- cmsTheta1[i] )
	      asymFlag=true;
	  }
	else
	  {
	    if(cmsTheta2[i] >   TMath::Pi() - cmsTheta1[i] )
	      asymFlag=true;
	  }

	if(asymFlag && asymmetryFlag){


	if(p_PiPi[i]>0 && chargeType[i]==0)
	  {
	    //	  cout <<"cut on asym" <<endl;
	  }
	if(print)
	  {
	    cout <<"cut due to asymmetry cut.." <<endl;
	  }
	  cut[i]=1;

	  //	  cout <<"cut asym " <<endl;
	}

       if(particleType1[i]!=0 || particleType2[i]!=0)
	  {
	    //	    cut[i]=1;
	  }

	if(z2[i]<secondZCut)
	  {
	    //	    cout <<"second z cut " <<endl;
	    if(evtNr==DEBUG_EVENT)
	      {
		cout <<"second z cut d " <<endl;
	      }
	    if(p_PiPi[i]>0 && chargeType[i]==0)
	      {
		//	
	      }
	if(print)
	  {
	    cout <<"second z cut on fid "<< z2[i] <<endl;
	    cout <<"cut due to second z cut.." <<endl;
	  }
	    cut[i]=1;
	  }

	if(p_PiPi[i]>0.5 && z1[i] < 0.1 && z2[i] < 0.1 && chargeType[i]==0)
	  {
	    //	    cout <<"kT pipi: "<< kT_PiPi[i] << " kt: "<< kT[i] <<" z1: " << z1[i] << " z2: " << z2[i] << " z1pipi: " << z1_PiPi[i] <<" z2pipi: " << z2_PiPi[i] <<endl;
	  }

	//////////-----------test: restrict to hadron from uds
	if(mMCFlag==mcFlagMC)
	  {
	    //	    cout <<"mc mother Gen1: " << motherGenId1[i] <<" mother gen 2 : " << motherGenId2[i]<<endl;
	      //	    if(abs(motherGenId1[i]) >3 || abs(motherGenId2[i]) > 3)
	    //

	    //	    cout <<"testing mothergen " <<endl;
	    //	    if(!(motherGenId1[i]==10022 && motherGenId2[i] == 10022))
	      {
		//		cout <<"mother gen id not uds : " << motherGenId1[i] <<" and " << motherGenId2[i] <<endl;
		//		cut[i]=1;
		//		cout <<"done " <<endl;
	      }
	      //	    else
	      {
		//		cout <<"hadron comes from uds..." <<endl;
	      }
	  }


	//in the woa we do not split up into the different combinations and weights, since we know the truth, so do this here with weights 1 and 0
	if(mMCFlag==mcFlagWoA || mMCFlag==mcFlagMC)
	  {
	    //qt <--> kt flip shold have been done already

	    p_PiPi[i]=0.0;
	    p_PiK[i]=0.0;
	    p_PiP[i]=0.0;
	    p_KPi[i]=0.0;
	    p_KK[i]=0.0;
	    p_KP[i]=0.0;
	    p_PPi[i]=0.0;
	    p_PK[i]=0.0;
	    p_PP[i]=0.0;
	    p_PiPi2[i]=0.0;
	    p_PiK2[i]=0.0;
	    p_PiP2[i]=0.0;
	    p_KPi2[i]=0.0;
	    p_KK2[i]=0.0;
	    p_KP2[i]=0.0;
	    p_PPi2[i]=0.0;
	    p_PK2[i]=0.0;
	    p_PP2[i]=0.0;

	    
	    //	    cout <<"looking at particle type: "<< particleType[i] <<endl;
	    if(particleType[i]<0)
	      {
		cout <<"pid < 0" <<endl;

		cut[i]=1;
		continue;
	      }
	    switch(particleType[i])
	      {
	      case PiPi:
		//		cout <<"got pion pair" <<endl;
		p_PiPi[i]=1.0;
		p_PiPi2[i]=1.0;
		kT_PiPi[i]=kT[i];
		z1_PiPi[i]=z1[i];
		z2_PiPi[i]=z2[i];
		break;
	      case PiK:
		p_PiK[i]=1.0;
		p_PiK2[i]=1.0;
		kT_PiK[i]=kT[i];
		z1_PiK[i]=z1[i];
		z2_PiK[i]=z2[i];

		break;
	      case PiP:
		p_PiP[i]=1.0;
		p_PiP2[i]=1.0;
		kT_PiP[i]=kT[i];
		z1_PiP[i]=z1[i];
		z2_PiP[i]=z2[i];

		break;

	      case KPi:
		p_KPi[i]=1.0;
		p_KPi2[i]=1.0;
		kT_KPi[i]=kT[i];
		z1_KPi[i]=z1[i];
		z2_KPi[i]=z2[i];

		break;
	      case KK:
		p_KK[i]=1.0;
		p_KK2[i]=1.0;
		kT_KK[i]=kT[i];
		z1_KK[i]=z1[i];
		z2_KK[i]=z2[i];

		break;
	      case KP:
		p_KP[i]=1.0;
		p_KP2[i]=1.0;
		kT_KP[i]=kT[i];
		z1_KP[i]=z1[i];
		z2_KP[i]=z2[i];

		break;
	      case PPi:
		p_PPi[i]=1.0;
		p_PPi2[i]=1.0;
		kT_PPi[i]=kT[i];
		z1_PPi[i]=z1[i];
		z2_PPi[i]=z2[i];

		break;
	      case PK:
		p_PK[i]=1.0;
		p_PK2[i]=1.0;
		kT_PK[i]=kT[i];
		z1_PK[i]=z1[i];
		z2_PK[i]=z2[i];

		break;
	      case PP:
		p_PP[i]=1.0;
		p_PP2[i]=1.0;
		kT_PP[i]=kT[i];
		z1_PP[i]=z1[i];
		z2_PP[i]=z2[i];

		break;
	      case UNKNOWN:
		//		cout <<"unknown pid " <<endl;
		cut[i]=1;
		continue;
		break;

	      default:
		cout <<"wrong pid in after fill : " << particleType[i] <<" particle type " << particleType[i]<<"charge type : "<< chargeType[i]<<endl;
		cout <<"p_PiPi: " << p_PiPi[i] << " kt_pip : "<< kT_PiPi[i] <<endl;
	      }


	  }
	//	cout <<"good id" <<endl;
	///------------------
	if(z1[i]<zCut || z2[i]<zCut)
	  {
	    //	    cout <<"fail z cut of " << zCut << " : " << z1[i] <<" z2: "<< z2[i] <<endl;
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"fail z cut " <<endl;
      }
	if(p_PiPi[i]>0 && chargeType[i]==0)
	  {
	    //	  cout <<"fail z on fid "<< z1[i] <<" and : " << z2[i]  <<endl;
	  }
	if(print)
	  {
	    cout <<"cut due to z cut (-1).." <<endl;
	  }

	    cut[i]=1;
	  }
	
	if(z1[i]<=0 || z1[i] >1.1|| z2[i]<=0 || z2[i] >1.1)
	  {
	    //	    if(particleType[i]==0 && chargeType[i]==0)
	    //	      cout <<" cut due to wrong z: " << z[i] <<endl;
	    //	    cout <<"wrong z " <<endl;
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"wrong z  " <<endl;
      }

	if(p_PiPi[i]>0 && chargeType[i]==0)
	  {
	    //	  cout <<"wrong z on fid" <<endl;
	  }
	if(print)
	  {
	    cout <<"cut due to z cut (1).." <<endl;
	  }
	    cut[i]=1;
	  }

	if(z1[i] >zUpperCut|| z2[i] >zUpperCut)
	  {
	    //	    cout <<" upper z cut " << endl;
    if(evtNr==DEBUG_EVENT)
      {
	cout <<"upper z " <<endl;
      }
	if(p_PiPi[i]>0 && chargeType[i]==0)
	  {
	    //	  cout <<"upper z on fid" <<endl;
	  }
	if(print)
	  {
	    cout <<"cut due to z upper cut.." <<endl;
	  }
	    cut[i]=1;
	  }





	//	if(fabs(thrustProj1[i])>maxOpeningCut || fabs(thrustProj2[i])>maxOpeningCut)
	//	  cut[i]=1;
	if(isnan(z1[i])|| isnan(z2[i]))
	  {

	   cut[i]=1;
	     cout <<"nan!" <<endl;
	  }

	if(mMCFlag==mcFlagWoA && cut[i]==0)
	  {
	    //	    cout <<"accepted woa event with kt: "<< kT[i] << " z1: " << z1[i] << " z2: "<< z2[i] <<endl;
	  }

	if(cut[i]==1)
	  {
	    if(p_PiPi[i]>0 && chargeType[i]==0)
	      {
		//	      cout <<"cutted on fid" <<endl;
	      }
	    
	    //	    cout <<"Hadron pair was cut .. " <<endl;
	  }	

	if(DEBUG_EVENT==evtNr)
	  {
	    //	    cout <<" pPi: " << p_PiPi[i] << p_PiPi[i] << " chargetype: "<< chargeType[i] <<endl;
	  }

	//pairChargeLikesign==0
		if(p_PiPi[i]>0 && chargeType[i]==0 && !cut[i] && print)
	//	if(print)
	  {

	    //
	    //	    	    cout << "looking at pair pipi with prob: " << p_PiPi[i] <<endl;
	    //	    cout <<" z1: " <<z1[i] <<" z2: "<< z2[i] <<" kT: "<< kT[i] <<endl;
	    //	    	    cout<<"Event " << evtNr;
		    if(cut[i])
		      {
		      cout <<" (cut==1) ";
		      }
		    //		    cout <<" charge type : "<< chargeType[i];
		    //	      	    cout <<std::fixed;
		    //	    	    	    cout.precision(3);
		    //		    	       cout  << " kT " <<kT[i];
		    //		    	       cout <<" z1: "<< z1[i] <<" z2: " << z2[i];
		    //	    	    cout <<endl;
		    //			cout <<"prop pipi: " << p_PiPi[i] << " piK: "<< p_PiK[i] <<" piP: "<< p_PiP[i] << " p_KPi: "<< p_KPi[i] <<" KK: " << p_KK[i] << " KP: ";
			//	cout << p_KP[i]<<endl<< " PPi: " << p_PPi[i] << " PK " << p_PK[i] <<" PP: " << p_PP[i]<<endl;
	  }
	if(p_PiPi[i]>0.5 && z1[i] < 0.1 && z2[i] < 0.1 && chargeType[i]==0)
	  {
	    //	    cout <<"cut? " << cut[i] <<endl;
	  }
      }
    //    cout <<"---------"<<endl;




  }


  //clone this guy...
  HadronPairArray& operator =(HadronPairArray rhs);
  //assign +- randomly, i.e. change direction of R randomly to check for false asymmetries
  void setElement(int pairCounter,HadronPairArray& hp,int index);
  void print();
  protected:
  void setSingleElement(int pairCounter,HadronPairArray& hp,int index);


};


#endif
