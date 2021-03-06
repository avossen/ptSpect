#ifndef PLOT_RESULTS_H
#define PLOT_RESULTS_H

#include <iostream>
#include <cmath>

#include "TObject.h"
using namespace std;

const int maxKtBins=30;
//#include "TNamed.h"

//structure to save results of fits...
class PlotResults   : public TObject
{
 public:
  enum fileType{uds,charm,eeuu,eess,eecc,tautau,fileTypeEnd};
  float meanKinBin1;
  float meanKinBin2;

  float kTUncertainties[maxKtBins];
  float kTUncertainties1[maxKtBins];
  float kTUncertainties2[maxKtBins];
  float kTSysUncertainties[maxKtBins];
  //in case we have asymmetric uncertainties after unfolding
  float kTSysUncertaintiesLower[maxKtBins];
  float kTValues[maxKtBins];
  //using alternative pid methods
  float kTValues1[maxKtBins];
  float kTValues2[maxKtBins];
  float kTMeans[maxKtBins];

  float weakDecayFraction[maxKtBins];
  

  int numKtValues;


  int exp;
  bool on_res;

  //for mc
  bool isUds;
  bool isCharm;
  bool isMC;
  int fileType;
  

  //1D, 2d, DR
  int calcType;

  int binningType, pidBin,chargeBin, firstKinBin, secondKinBin;
  int resultIndex;
  //needed to create vtable needed for base class...
  PlotResults(): numKtValues(0), mPrint(false)
    {};
  virtual ~PlotResults(){};

  PlotResults& operator +=(const PlotResults& rhs);
  void print();
  void doPrint(bool print=true);

  //no assignment operator or copy constructor since we don't have pointers (would be different with pointers)
 protected:
  bool mPrint;

private:   
  ClassDef(PlotResults,1);
};
inline void PlotResults::doPrint(bool print)
{
  if(print)
    mPrint=print;
}



//PlotResults& PlotResults::operator +=(const PlotResults& rhs);
inline PlotResults& PlotResults::operator +=(const PlotResults& rhs)
{
      double intRhsCount=0;
      double intCount=0;
      for(int i=0;i<maxKtBins;i++)
	{
	  cout <<" combining ktbin : " << i <<endl;
	  intRhsCount+=rhs.kTValues[i];
	  intCount+= kTValues[i];
	  
	  if(kTValues[i]==0)
	    kTMeans[i]=rhs.kTMeans[i];
	  if(kTValues[i]!=0 && rhs.kTValues[i]!=0)
	    {
	      
	      kTMeans[i]=1.0/(kTValues[i]+rhs.kTValues[i])*(kTMeans[i]*kTValues[i] + rhs.kTMeans[i]*rhs.kTValues[i]);
	    }
	  //just assign the rhs side. Assume that we are not really adding files 
	  //	  if(weakDecayFraction[i]!=0 && rhs.weakDecayFraction[i]!=0)
	  if((kTValues[i]+rhs.kTValues[i])>0)
	    {
	      //we haven't added to the kT values yet...
	      //total weak decays
	      double origDecayFraction=weakDecayFraction[i];
	      if(isnan(weakDecayFraction[i]))
		origDecayFraction=0.0;
	      double rhsWeakDecays=rhs.weakDecayFraction[i];
	      if(origDecayFraction > 1.0 || rhsWeakDecays > 1.0)
		{cout <<"fraction greater one" <<endl;}
	      if(::isnan(rhsWeakDecays))
		{
		  rhsWeakDecays=0.0;
		}
	      //	      total weak decays
	      weakDecayFraction[i]=origDecayFraction*kTValues[i]+rhsWeakDecays*rhs.kTValues[i];
	      weakDecayFraction[i]/=(kTValues[i]+rhs.kTValues[i]);
	      cout <<"combining fractions: " << origDecayFraction <<" and " << rhsWeakDecays << " to " << weakDecayFraction[i]<<endl;
	      cout <<"ktVals: "<< kTValues[i] <<" rhs ktValues: "<< rhs.kTValues[i] <<endl;
	    }
	  kTValues[i]+=rhs.kTValues[i];
	  kTValues1[i]+=rhs.kTValues1[i];
	  kTValues2[i]+=rhs.kTValues2[i];

	  kTSysUncertainties[i]=sqrt(kTSysUncertainties[i]*kTSysUncertainties[i]+rhs.kTSysUncertainties[i]*rhs.kTSysUncertainties[i]);
	  kTUncertainties[i]=sqrt(kTUncertainties[i]*kTUncertainties[i]+rhs.kTUncertainties[i]*rhs.kTUncertainties[i]);
	  kTUncertainties1[i]=sqrt(kTUncertainties1[i]*kTUncertainties1[i]+rhs.kTUncertainties1[i]*rhs.kTUncertainties1[i]);
	  kTUncertainties2[i]=sqrt(kTUncertainties2[i]*kTUncertainties2[i]+rhs.kTUncertainties2[i]*rhs.kTUncertainties2[i]);

	}
      
      if(intCount==0)
	{
	  meanKinBin1=rhs.meanKinBin1;
	  meanKinBin2=rhs.meanKinBin2;
	}
      if(intCount!=0 && intRhsCount!=0)
	{
	  meanKinBin1=1.0/(intCount+intRhsCount)*(meanKinBin1*intCount+rhs.meanKinBin1*intRhsCount);
	  meanKinBin2=1.0/(intCount+intRhsCount)*(meanKinBin2*intCount+rhs.meanKinBin2*intRhsCount);
	}


      fileType=rhs.fileType;
  isUds=rhs.isUds;
  isCharm=rhs.isCharm;
  isMC=rhs.isMC;
  calcType=rhs.calcType;
  binningType=rhs.binningType;
  chargeBin=rhs.chargeBin;
  firstKinBin=rhs.firstKinBin;
  secondKinBin=rhs.secondKinBin;
  resultIndex=rhs.resultIndex;
  pidBin=rhs.pidBin;

  return *this;
      
};


//to add to PlotResults. This means we have to incrementaly compute weighted mean etc
inline PlotResults operator+(PlotResults lhs, const PlotResults& rhs)
{
  std::cout <<"not implemented yet!!" <<endl;
  return lhs;
};
inline PlotResults operator-(PlotResults lhs, const PlotResults& rhs)
  {

    return lhs;
  }

#endif
