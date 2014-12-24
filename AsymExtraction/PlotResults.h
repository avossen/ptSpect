
#ifndef PLOT_RESULTS_H
#define PLOT_RESULTS_H

#include <cmath>
#include <iostream>

#include "TObject.h"
using namespace std;


//#include "TNamed.h"

//structure to save results of fits...
class PlotResults   : public TObject
{
 public:
  float meanKinBin1;
  float meanKinBin2;

  float kTValues[30];
  float kTMeans[30];

  int numKtValues;


  int exp;
  bool on_res;

  //for mc
  bool isUds;
  bool isCharm;
  bool isMC;

  //1D, 2d, DR
  int calcType;

  int binningType, chargeBin, firstKinBin, secondKinBin;
  int resultIndex;
  //needed to create vtable needed for base class...
  PlotResults(): mPrint(false), numKtValues(0)
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

inline PlotResults& PlotResults::operator +=(const PlotResults& rhs)
{

  return *this;

}


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
