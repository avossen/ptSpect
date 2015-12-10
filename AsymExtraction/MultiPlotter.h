#ifndef MULTI_PLOTTER_H
#define MULTI_PLOTTER_H
#include "NamedExp.h"
#include "MEvent.h"
#include "TwoHadAsymsCommons.h"
#include "TMath.h"
#include "PlotResults.h"
#include "TFile.h"
#include "HadronPairArray.h" 
#include "stdlib.h"
#include <cstdlib>


enum binningType{binType_labTheta_z, binType_ThrustLabTheta_z, binType_z_z,binType_zOnly, binType_labThetaOnly, binType_qTOnly, binType_ThrustOnly,binType_end};
enum pairType{pairChargeInt, pairPN,  pairUnknown,pairTypeEnd};
enum plotType{plotType_2D ,plotType_end};

class MultiPlotter: public ReaderBase, NamedExp//for the normalize angle
{
 public:
  MultiPlotter(const char* filenameBase,string nameAdd, int exNr, bool onRes, bool uds, bool charm,bool mc):NamedExp(filenameBase,nameAdd,exNr,onRes,uds,charm,mc)
    {
      zeroBin=0;
      rFile.mkdir("plotHistos");

      rFile.cd();

      //l      loadBinning(binningKt, binningZ);
      loadBinnings();
      numKtBins=binningKt.size();
      maxKinBins=binningZ.size();
      	  cout <<"loading " << binningZ.size() <<" z bins " << endl;
      if(binningKt.size()>binningZ.size())
	{

	  maxKinBins=numKtBins;
	}
      if(binningLabTheta.size()>maxKinBins)
	{

	  maxKinBins=binningLabTheta.size();
	}
      cout << " loading " <<binningLabTheta.size() << " labT bin " << endl;
      if(binningThrustLabTheta.size()>maxKinBins)
	maxKinBins=binningLabTheta.size();
      if(binningQt.size()>maxKinBins)
	maxKinBins=binningQt.size();
      if(binningCmsThrustTheta.size()>maxKinBins)
	maxKinBins=binningCmsThrustTheta.size();
      if(binningThrust.size()>maxKinBins)
	maxKinBins=binningThrust.size();
      setBinningMap();


      plotResults=new PlotResults[numKinematicBinning*NumCharges*maxKinBins*maxKinBins];
      plotResVect.push_back(plotResults);

      for(unsigned int i=0;i<numKinematicBinning*NumCharges*maxKinBins*maxKinBins;i++)
	{
	  for(vector<PlotResults*>::iterator it=plotResVect.begin();it!=plotResVect.end();it++)
	    {
	      (*it)[i].meanKinBin1=0.0;
	      (*it)[i].meanKinBin2=0.0;
	      (*it)[i].numKtValues=numKtBins;
	      for(unsigned int j=0;j<numKtBins;j++)
		{
		  (*it)[i].kTValues[j]=0.0;
		  (*it)[i].kTMeans[j]=0.0;
		}
	    }
	}


      cout <<"allocating " << numKinematicBinning <<" * " << NumCharge << " * " << maxKinBins <<" * " << maxKinBins <<" * " <<numKtBins <<" * " <<numKtBins<<endl;
      counts=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numKtBins);
      meanValues_kin1=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins);
      meanValues_kin2=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins);
      meanValues_kT=allocateArray<double>(numKinematicBinning,NumCharges,maxKinBins,maxKinBins,numKtBins);

    };

    //
    //    void setFitReuslt();

    void addHadPairArray(HadronPairArray* hq, MEvent& event);
    void setBinningMap();
    void doPlots();
    void savePlots(plotType);

    void setName(string s);
    string getName();
    PlotResults* plotResults;

    vector<PlotResults*> plotResVect;

    //maps the type of binning to two bin pointers
    vector< pair<int*,int*> > binningMap;
    vector< pair<float*,float*> > meanMap;
    vector< pair<int, int> > maxKinMap;
    inline int getResIdx(int binningType,int chargeType, int firstKinBin, int secondKinBin)
      {
	return binningType*NumCharges*maxKinBins*maxKinBins+chargeType*maxKinBins*maxKinBins+firstKinBin*maxKinBins+secondKinBin;
      }
 protected:
    //to reorder arrays used for fitting such that the x values are ascending and do not wrap around
    //should work because y is moved as well
    void reorder(float* mX, float* mY, float* mYErr, int numBins);
    static const int numKinematicBinning;
    static const int NumCharges;

    //to differentiate between different fit results...

    TH1D* hChi2OverNdf;
    TH1D* hOneDVsTwoDA1;
    TH1D* hOneDVsTwoDA2;
    TH1D* hOneDVsTwoDA3;


    //give negative values for first or second bin if they should not be part of the name
    string getBinName(int binningType,int chargeType, int firstBin, int secondBin);
    string getXAxisName(int binningType);
    void loadBinnings();
 public:
    vector<float> binningKt;
    vector<float> binningZ;
    vector<float> binningLabTheta;
    vector<float> binningThrustLabTheta;
    vector<float> binningQt;
    vector<float> binningCmsThrustTheta;
    vector<float> binningThrust;

 protected:
    int zbin1;
    int zbin2;
    int kTBin;
    //always zero
    int zeroBin;
    int labThetaBin1;
    int labThetaBin2;

    int thrustLabThetaBin;
    int cmsThrustThetaBin;
    int cmsThrustPhiBin;
    int multBin;
    int qTBin;
    int thrustBin;


    float z1;
    float z2;
    float kT;

    float labTheta1;
    float labTheta2;
    float thrustLabTheta;


    float cmsThrustTheta;
    float qT;
    float thrust;



    unsigned int numKtBins;
    unsigned int maxKinBins;
    double***** counts;

    double**** meanValues_kin1;
    double**** meanValues_kin2;
    double***** meanValues_kT;
};

inline string MultiPlotter::getName()
{
  return nameAddition;
};

inline void MultiPlotter::setName(string s)
{
  nameAddition=s;
};

#endif
