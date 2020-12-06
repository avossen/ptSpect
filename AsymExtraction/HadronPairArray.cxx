#include "HadronPairArray.h"

HadronPairArray& HadronPairArray::operator=(HadronPairArray rhs)
  {
    numPairs=rhs.numPairs;
    mMCFlag=rhs.mMCFlag;
    addendum=rhs.addendum;
    for(int i=0;i<rhs.numPairs;i++)
      {
	//use this function so that all and single elements change the same...
	setSingleElement(i,rhs,i);
	/*(	cut[i]=rhs.cut[i];
	z[i]=rhs.z[i];
	mass[i]=rhs.mass[i];
	phiR[i]=rhs.phiR[i];
	phiZero[i]=rhs.phiZero[i];
	thrustProj1[i]=rhs.thrustProj1[i];
	thrustProj2[i]=rhs.thrustProj2[i];
	theta1[i]=rhs.theta1[i];
	theta2[i]=rhs.theta2[i];
	chargeType[i]=rhs.chargeType[i];
	particleType[i]=rhs.particleType[i];
	pi0Sig[i]=rhs.pi0Sig[i];
	pi0Bg[i]=rhs.pi0Bg[i];
	pi0mass1[i]=rhs.pi0mass1[i];
	pi0mass2[i]=rhs.pi0mass2[i];*/
      }
    return *this;
  };



//the the element 'pairCounter' of this pair to the element 'index' of the argument hadpair. This is used for event mixing...

void HadronPairArray::setSingleElement(int pairCounter,HadronPairArray& hp,int index)
{
  z1[pairCounter]=hp.z1[index];
  z2[pairCounter]=hp.z1[index];

  thrustProj1[pairCounter]=hp.thrustProj1[index];
  thrustProj2[pairCounter]=hp.thrustProj2[index];

  labTheta1[pairCounter]=hp.labTheta1[index];
  labTheta2[pairCounter]=hp.labTheta2[index];

  chargeType[pairCounter]=hp.chargeType[index];
  chargeType1[pairCounter]=hp.chargeType1[index];
  chargeType2[pairCounter]=hp.chargeType2[index];
  particleType[pairCounter]=hp.particleType[index];

  particleType1[pairCounter]=hp.particleType1[index];
  particleType2[pairCounter]=hp.particleType2[index];

  cut[pairCounter]=hp.cut[index];

};


void HadronPairArray::setElement(int pairCounter,HadronPairArray& hp,int index)
{
    mMCFlag=hp.mMCFlag;
    addendum=hp.addendum;
    setSingleElement(pairCounter,hp,index);
};

void HadronPairArray::print()
{

  cout <<"Hadron Index | z |  zRatio | mass | phiR |phiZero | thrustProj1 | thrustProj2 | theta1 | theta2 | chargeType | particleType | cut |  hadOpening " <<endl;
  for(int i=0;i<numPairs;i++)
    {
      cout <<i<<" | " << z1[i] << " | " << z2[i] << " | " <<" | " << thrustProj1[i] << " | " <<thrustProj2[i] << " | " <<labTheta1[i] << " | " <<labTheta2[i] << " | " <<chargeType[i] << " | " <<particleType[i] << " | " <<cut[i]  << " | " <<endl;
    }



}

