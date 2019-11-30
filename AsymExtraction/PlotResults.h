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


  int numKtValues;


  int exp;
  bool on_res;

  //for mc
  bool isUds;
  bool isCharm;
  bool isMC;

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
	  intRhsCount+=rhs.kTValues[i];
	  intCount+= kTValues[i];
	  
	  if(kTValues[i]==0)
	    kTMeans[i]=rhs.kTMeans[i];
	  if(kTValues[i]!=0 && rhs.kTValues[i]!=0)
	    {
	      
	      kTMeans[i]=1.0/(kTValues[i]+rhs.kTValues[i])*(kTMeans[i]*kTValues[i] + rhs.kTMeans[i]*rhs.kTValues[i]);
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
