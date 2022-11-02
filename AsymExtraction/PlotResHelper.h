#ifndef PLOT_RES_HELPER_H_
#define  PLOT_RES_HELPER_H_
#include "PlotResults.h"

//reimplementation of the PlotResults functions so we don't break the serialized objects


PlotResults&  add(PlotResults& lhs, const PlotResults& rhs, float prefactor=1.0)
{
      double intRhsCount=0;
      double intCount=0;
      for(int i=0;i<maxKtBins;i++)
	{
	  cout <<" combining ktbin : " << i <<endl;
	  intRhsCount+=rhs.kTValues[i];
	  intCount+= lhs.kTValues[i];
	  
	  if(lhs.kTValues[i]==0)
	    lhs.kTMeans[i]=rhs.kTMeans[i];
	  if(lhs.kTValues[i]!=0 && rhs.kTValues[i]!=0)
	    {
	      
	      lhs.kTMeans[i]=1.0/(lhs.kTValues[i]+rhs.kTValues[i])*(lhs.kTMeans[i]*lhs.kTValues[i] + rhs.kTMeans[i]*rhs.kTValues[i]);
	    }
	  //just assign the rhs side. Assume that we are not really adding files 
	  //	  if(weakDecayFraction[i]!=0 && rhs.weakDecayFraction[i]!=0)
	  if((lhs.kTValues[i]+rhs.kTValues[i])>0)
	    {
	      //we haven't added to the kT values yet...
	      //total weak decays
	      double origDecayFraction=lhs.weakDecayFraction[i];
	      if(isnan(lhs.weakDecayFraction[i]))
		origDecayFraction=0.0;
	      double rhsWeakDecays=rhs.weakDecayFraction[i];
	      if(origDecayFraction > 1.0 || rhsWeakDecays > 1.0)
		{cout <<"fraction greater one" <<endl;}
	      if(::isnan(rhsWeakDecays))
		{
		  rhsWeakDecays=0.0;
		}
	      //	      total weak decays
	      lhs.weakDecayFraction[i]=origDecayFraction*lhs.kTValues[i]+rhsWeakDecays*rhs.kTValues[i];
	      lhs.weakDecayFraction[i]/=(lhs.kTValues[i]+rhs.kTValues[i]);
	      cout <<"combining fractions: " << origDecayFraction <<" and " << rhsWeakDecays << " to " << lhs.weakDecayFraction[i]<<endl;
	      cout <<"ktVals: "<< lhs.kTValues[i] <<" rhs ktValues: "<< rhs.kTValues[i] <<endl;
	    }
	  lhs.kTValues[i]+=rhs.kTValues[i];
	  lhs.kTValues1[i]+=rhs.kTValues1[i];
	  lhs.kTValues2[i]+=rhs.kTValues2[i];

	  lhs.kTSysUncertainties[i]=sqrt(lhs.kTSysUncertainties[i]*lhs.kTSysUncertainties[i]+rhs.kTSysUncertainties[i]*rhs.kTSysUncertainties[i]);
	  lhs.kTUncertainties[i]=sqrt(lhs.kTUncertainties[i]*lhs.kTUncertainties[i]+rhs.kTUncertainties[i]*rhs.kTUncertainties[i]);
	  lhs.kTUncertainties1[i]=sqrt(lhs.kTUncertainties1[i]*lhs.kTUncertainties1[i]+rhs.kTUncertainties1[i]*rhs.kTUncertainties1[i]);
	  lhs.kTUncertainties2[i]=sqrt(lhs.kTUncertainties2[i]*lhs.kTUncertainties2[i]+rhs.kTUncertainties2[i]*rhs.kTUncertainties2[i]);

	}
      
      if(intCount==0)
	{
	  lhs.meanKinBin1=rhs.meanKinBin1;
	  lhs.meanKinBin2=rhs.meanKinBin2;
	}
      if(intCount!=0 && intRhsCount!=0)
	{
	  lhs.meanKinBin1=1.0/(intCount+intRhsCount)*(lhs.meanKinBin1*intCount+rhs.meanKinBin1*intRhsCount);
	  lhs.meanKinBin2=1.0/(intCount+intRhsCount)*(lhs.meanKinBin2*intCount+rhs.meanKinBin2*intRhsCount);
	}


      lhs.fileType=rhs.fileType;
  lhs.isUds=rhs.isUds;
  lhs.isCharm=rhs.isCharm;
  lhs.isMC=rhs.isMC;
  lhs.calcType=rhs.calcType;
  lhs.binningType=rhs.binningType;
  lhs.chargeBin=rhs.chargeBin;
  lhs.firstKinBin=rhs.firstKinBin;
  lhs.secondKinBin=rhs.secondKinBin;
  lhs.resultIndex=rhs.resultIndex;
  lhs.pidBin=rhs.pidBin;

  return lhs;


}



#endif
