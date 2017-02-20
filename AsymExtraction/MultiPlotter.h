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
#include "TSVDUnfold.h"

enum binningType{binType_labTheta_z, binType_ThrustLabTheta_z, binType_z_z,binType_zOnly, binType_labThetaOnly, binType_qTOnly, binType_ThrustOnly,binType_end};
//enum binningType{ binType_z_z,binType_zOnly,binType_end};
enum pairType{pairChargeLikesign, pairPN,  pairUnknown,pairTypeEnd};
enum plotType{plotType_2D ,plotType_end};
//enum pidType{PiPi, PiK, PiP, KPi, KK, KP, PPi, PK, PP, pidTypeEnd};

class MultiPlotter: public ReaderBase, NamedExp//for the normalize angle
{
 public:
  MultiPlotter(const char* filenameBase,string nameAdd, int exNr, bool onRes, bool uds, bool charm,bool mc):NamedExp(filenameBase,nameAdd,exNr,onRes,uds,charm,mc)
    {

      maxSmearing=new int*[NumPIDs];
      for(int i=0;i<NumPIDs;i++)
	{
	  maxSmearing[i]=new int[2];
	}
			       
      zeroBin=0;
      rFile.mkdir("plotHistos");

      rFile.cd();
      char buffer[500];
      //l      loadBinning(binningKt, binningZ);
      loadBinnings();
      numKtBins=binningKt.size();
      int maxZBins=binningZ[0].size();
      for(int i=1;i<6;i++)
	{
	  if(binningZ[i].size()>maxZBins)
	    maxZBins=binningZ[i].size();
	}
      maxKinBins=maxZBins;

      kinematicSmearingMatrix=new TH2D***[2];
      xini=new TH1D***[2];
      bini=new TH1D***[2];

      for(int b=0;b<2;b++)
	{

	  kinematicSmearingMatrix[b]=new TH2D**[NumPIDs];
	  xini[b]=new TH1D**[NumPIDs];
	  bini[b]=new TH1D**[NumPIDs];
	  for(int p=0;p<NumPIDs;p++)
	    {
	      pair<int,int> zIdx=pidBin2ZBinningIdx(p);
	      maxSmearing[p][0]=numKtBins*binningZ[zIdx.first].size()*binningZ[zIdx.second].size();
	      maxSmearing[p][1]=numKtBins*binningZ[zIdx.first].size();

	      kinematicSmearingMatrix[b][p]=new TH2D*[NumCharges];
	      xini[b][p]=new TH1D*[NumCharges];
	      bini[b][p]=new TH1D*[NumCharges];
	      for(int c=0;c<NumCharges;c++)
		{
		  sprintf(buffer,"kinematicSmearingMatrix_binning%d_pidBin%d_chargeBin%d",b,p,c);
		  kinematicSmearingMatrix[b][p][c]=new TH2D(buffer,buffer,maxSmearing[p][b],0,maxSmearing[p][b],maxSmearing[p][b],0,maxSmearing[p][b]);
		  sprintf(buffer,"xini_binning%d_pidBin%d_chargeBin%d",b,p,c);
		  xini[b][p][c]=new TH1D(buffer,buffer,maxSmearing[p][b],0,maxSmearing[p][b]);
		  sprintf(buffer,"bini_binning%d_pidBin%d_chargeBin%d",b,p,c);
		  bini[b][p][c]=new TH1D(buffer,buffer,maxSmearing[p][b],0,maxSmearing[p][b]);
		}
	    }
	}
      cout <<"loading " << maxZBins <<" z bins " << endl;
      if(binningKt.size()>maxZBins)
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


      plotResults=new PlotResults[numKinematicBinning*NumPIDs*NumCharges*maxKinBins*maxKinBins];
      plotResVect.push_back(plotResults);

      for(unsigned int i=0;i<numKinematicBinning*NumCharges*NumPIDs*maxKinBins*maxKinBins;i++)
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
      counts=allocateArray<double>(numKinematicBinning,NumPIDs,NumCharges,maxKinBins,maxKinBins,numKtBins);
      meanValues_kin1=allocateArray<double>(numKinematicBinning,NumPIDs,NumCharges,maxKinBins,maxKinBins);
      meanValues_kin2=allocateArray<double>(numKinematicBinning,NumPIDs,NumCharges,maxKinBins,maxKinBins);
      meanValues_kT=allocateArray<double>(numKinematicBinning,NumPIDs,NumCharges,maxKinBins,maxKinBins,numKtBins);

    };

    //
    //    void setFitReuslt();

    void addHadPairArray(HadronPairArray* hq, MEvent& event);
    void addSmearingEntry(HadronPairArray* hq1, HadronPairArray* hq2,bool accSmearing=false);
    void setBinningMap();
    void doPlots();
    void savePlots(plotType);

    //get the histogram that is to be unfolded. This is a 1D histogram binned in z1, (z2) and kT
    TH1D* getHistogram(int binning, int chargeBin, int pidBin);
    //d is a return value
    TH1D* unfold(TH2D* smearingMatrix, TH1D* MC_input,TH1D* MC_out, TH1D* data, TH1D** d);
    TH1D** convertUnfold2Plots(TH1D* input,int binning,  int chargeBin, int pidBin, const char* nameAdd);
    void setName(string s);
    string getName();
    PlotResults* plotResults;

    vector<PlotResults*> plotResVect;

    //maps the type of binning to two bin pointers
    vector< pair<int*,int*> > binningMap;
    vector< pair<float*,float*> > meanMap;
    vector< pair<int, int> >* maxKinMap;
    inline int getResIdx(int binningType,int pidType,int chargeType, int firstKinBin, int secondKinBin)
      {
	return binningType*NumPIDs*NumCharges*maxKinBins*maxKinBins+pidType*NumCharges*maxKinBins*maxKinBins+chargeType*maxKinBins*maxKinBins+firstKinBin*maxKinBins+secondKinBin;
      }

    void saveSmearingMatrix();
    static const int NumCharges;
    static const int NumPIDs;
    pair<int,int> pidBin2ZBinningIdx(int pidBin);
 protected:
    //to reorder arrays used for fitting such that the x values are ascending and do not wrap around
    //should work because y is moved as well
    void reorder(float* mX, float* mY, float* mYErr, int numBins);
    static const int numKinematicBinning;

    //to differentiate between different fit results...

    TH1D* hChi2OverNdf;
    TH1D* hOneDVsTwoDA1;
    TH1D* hOneDVsTwoDA2;
    TH1D* hOneDVsTwoDA3;

    TH2D**** kinematicSmearingMatrix;
    TH1D**** xini;
    TH1D**** bini;


    //give negative values for first or second bin if they should not be part of the name
    string getBinName(int binningType,int pidType,int chargeType, int firstBin, int secondBin);
    string getXAxisName(int binningType);
    void loadBinnings();
 public:
    vector<float> binningKt;
    //make this dependent on particle id
    vector<float>* binningZ;

    vector<float> binningLabTheta;
    vector<float> binningThrustLabTheta;
    vector<float> binningQt;
    vector<float> binningCmsThrustTheta;
    vector<float> binningThrust;

 protected:
    int zbin1;
    int zbin2;
    int kTBin;
    int** maxSmearing;
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
    double****** counts;

    double***** meanValues_kin1;
    double***** meanValues_kin2;
    double****** meanValues_kT;
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
