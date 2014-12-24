#ifndef HADRON_PAIR_ARRAY_H
#define HADRON_PAIR_ARRAY_H
#include "ReaderBase.h"
#include "TMath.h"

//300 is also the number in the Treesaver, so cannot be larger...
#define Max_ArrSize 300

struct HadronPairArray:public ReaderBase
{
  //cuts:

  const float zCut;
  const float zUpperCut;
  const float secondZCut;

  int numPairs;

  float z1[Max_ArrSize];
  float z2[Max_ArrSize];


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


  //considering the thrust axis resolution of 0.16+-0.09 rad, a max opening cut of 0.99 is even too large...
  HadronPairArray(TChain* chain, int MCFlag=mcFlagNone):ReaderBase(MCFlag), zCut(0.1),zUpperCut(1.4), secondZCut(0.25)
  {
    //no chain implies standalone. Cannot branch on the same field twice, this would override
    //    cout <<" do we have a chain? " << chain<<endl;
    if(chain)
      {
	//cout <<"yes " <<endl;
	myChain=chain;
	branchPointersI.push_back(&numPairs);
	branchPointers.push_back(z1);
	branchPointers.push_back(z2);

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

	branchPointers.push_back(kT);
	//only for mc...
	branchPointers.push_back(qT);

	branchPointers.push_back(hadDiffTheta);
	branchPointers.push_back(hadDiffPhi);


	branchPointersI.push_back(chargeType2);
	branchPointersI.push_back(particleType2);

	branchPointersI.push_back(chargeType1);
	branchPointersI.push_back(particleType1);

	branchPointersI.push_back(chargeType2);
	branchPointersI.push_back(particleType2);

	if(mMCFlag!=mcFlagNone)
	  {
	  }
	
	branchNamesI.push_back("z1"+addendum+"Counter");
	branchNames.push_back("z1"+addendum);
	branchNames.push_back("z2"+addendum);

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

	if(mMCFlag!=mcFlagWoA)
	  branchNames.push_back("kT"+addendum);
	else//woa has not mcWoA ending for kt
	  branchNames.push_back("kT");
	//only for non mc...and woa (so no addendum in the other cases..)

	if(mMCFlag==mcFlagWoA)
	  branchNames.push_back("qT"+addendum);
	else
	  branchNames.push_back("qT");




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
  }

  void afterFill()
  {
    for(int i=0;i<numPairs;i++)
      {
	//	cout <<"hp after fill!" <<endl;
	cut[i]=0;

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

	//	if(asymFlag)
		  //		  cut[i]=1;


	if(z2[i]<secondZCut)
	  cut[i]=1;

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
	//it seems that the WoA has the evtGen codes. So parents don't seem to be quarks but the virutal photon...
	//nope, not true...
	if(mMCFlag==mcFlagWoA)
	  {
	    //	    if(!(motherGenId1[i]==10022 && motherGenId2[i] == 10022))
	      {
		//		cout <<"mother gen id not uds : " << motherGenId1[i] <<" and " << motherGenId2[i] <<endl;
		//				cut[i]=1;
	      }
	      //	    else
	      {
		//		cout <<"hadron comes from uds..." <<endl;
	      }
	  }

	///------------------
	if(z1[i]<zCut || z2[i]<zCut)
	  {
	    cut[i]=1;
	  }

	if(z1[i]<=0 || z1[i] >1.1|| z2[i]<=0 || z2[i] >1.1)
	  {
	    //	    if(particleType[i]==0 && chargeType[i]==0)
	    //	      cout <<" cut due to wrong z: " << z[i] <<endl;
	    cut[i]=1;
	  }

	if(z1[i] >zUpperCut|| z2[i] >zUpperCut)
	  {

	    cut[i]=1;
	  }





	//	if(fabs(thrustProj1[i])>maxOpeningCut || fabs(thrustProj2[i])>maxOpeningCut)
	//	  cut[i]=1;
	if(isnan(z1[i])|| isnan(z2[i]))
	  {
	   cut[i]=1;
	   //	   cout <<"nan!" <<endl;
	  }


      }
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
