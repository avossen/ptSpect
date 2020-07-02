#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
#include "MEvent.h"
#include "StyleSetter.h"
//#include "HadronQuadArray.h"
#include "MultiPlotter.h"
#include "TwoHadAsymsCommons.h"
#include "PlotResults.h"
#include "CombPlots.h"
//#define MAX_EVENTS 100

using namespace std;

int main(int argc, char** argv)
{
  if(argc<5)
    {
      cout <<"usage: drawNonQQContrib eeuu eess eecc tautau"<<endl;
      exit(0);
    }
  bool m_useQt=false;
  //should have hadronPairArray where this is defined included
#ifdef USE_QT
  m_useQt=true;
#endif
  cout <<"hallo?"<<endl;
  vector<string> flavor;

  //  flavor.push_back("_uds");
  //  flavor.push_back("_charm");
  flavor.push_back("_all");

  cout <<"argvs..." <<endl;
  
  char* eeuuPath=argv[1];
  char* eessPath=argv[2];
  char* eeccPath=argv[3];
  char* tautauPath=argv[4];
  char* paths[4];
  cout <<"..." <<endl;
  for(int i=0;i<4;i++)
    {
      paths[i]=argv[i+1];
    }
  
  srand(time(NULL));
  enum fileType{uds,charm,eeuu,eess,eecc,tautau,data,fileTypeEnd};

  //reduced asymmetries for charm, uds, all

  setStyleOpts();

  //  vector< pair<string,TChain*> > vFitterNames;
  vector<string> vPlotterNames;


  vector<MultiPlotter*> vPlotters;
  TChain* chUU=0;
  TChain* chSS=0;
  TChain* chCC=0;
    //tautau
  TChain* chTT=0;
  int counter=-1;

  TChain* myChains[4];
  myChains[0]=chUU;
  myChains[1]=chSS;
  myChains[2]=chCC;
  myChains[3]=chTT;

  cout <<"after chains " <<endl;
  string names[4]={"eeUU","eeSS","eeCC","tau-tau"};  
  string plotterName=string("Normal");
  MultiPlotter* plotters[4];
  
  char dataMcNameAdd[100];
  sprintf(dataMcNameAdd,"");
  cout <<"before the loop" <<endl;
  for(int i=0;i<4;i++)
    {
      myChains[i]=new TChain("PlotTree");
      myChains[i]->Add((string(paths[i])+"/"+plotterName+"_*.root").c_str());
      cout <<" adding: "<< (string(paths[i])+"/"+plotterName+"_*.root").c_str() <<endl;
      string fullName=plotterName+string(dataMcNameAdd);
      plotters[i]=new MultiPlotter(m_useQt,const_cast<char*>("."),(fullName+names[i]).c_str(),string(""),0,false,false,false,false,fileTypeEnd);
      plotters[i]->setName(fullName);
      Int_t nevents=myChains[i]->GetEntries();
      PlotResults* plotResults=0;
      myChains[i]->SetBranchAddress("PlotBranch",&plotResults);
      cout <<" running over " << nevents << endl;
      for(long j=0;j<nevents;j++)
	{
	  //	  cout <<"getting entry " << i <<endl;
	  //	  float locW[3]={0.0,0.0,0.0};
	  myChains[i]->GetEntry(j);
	  plotters[i]->plotResults[plotResults->resultIndex]+=(*plotResults);
	}
      cout <<"done with entries.." <<endl;
    }

  //should have now filled all plotters, do the histos

  int binningType=binType_z_z;
  int pidBin=0; //pipi
  int chargeBin=likesign; //same charge

  int maxFirstBin=plotters[0]->maxKinMap[pidBin][binningType].first;
  int maxSecondBin=plotters[0]->maxKinMap[pidBin][binningType].second;
  int resIdx=plotters[0]->getResIdx(binningType,pidBin,chargeBin,0,0);
  int numKtBins=plotters[0]->plotResults[resIdx].numKtValues;
  TGraph** graphs[4];
  for(int i=0;i<4;i++)
    {
      graphs[i]=new TGraph*[maxFirstBin];
    }
  cout <<"maxfirst: "<< maxFirstBin <<" second: " << maxSecondBin <<" resIdx: "<< resIdx <<" numkt: " << numKtBins <<endl;
  char buffer[400];
  double* X=new double[numKtBins];
  double* Y=new double[numKtBins];

  int colors[]={kRed,kBlue,kGreen,kBlack};
  int markerStyles[]={20,21,22,23};

  //max ktVals for each z bin to scale the graphs
  double maxVals[100];
  double maxValsX[100];
  for(int i=0;i<100;i++)
    {
      maxVals[i]=0;
      maxValsX[i]=0;
    }
  for(int i=0;i<4;i++)
    {
      for(int zbin=0;zbin<maxFirstBin;zbin++)
	{
	  cout <<"looking at zbin : "<< zbin <<endl;
	  //diagonal bins
	  resIdx=plotters[0]->getResIdx(binningType,pidBin,chargeBin,zbin,zbin);
	  cout <<"residx2 : "<< resIdx <<endl;
	  int missPoints=0;
	  for(int j=0;j<numKtBins;j++)
	    {
	      //    cout <<"filling x,y num " << j <<" with mean: " << plotters[i]->plotResults[resIdx].kTMeans[j] <<" value: "<< plotters[i]->plotResults[resIdx].kTValues[j]  <<endl;
	      double kTMean=plotters[i]->plotResults[resIdx].kTMeans[j];
	      double kTVal=plotters[i]->plotResults[resIdx].kTValues[j];
	      if(kTVal>maxVals[zbin])
		maxVals[zbin]=kTVal;

	      if(kTMean==0 || isnan(kTMean)|| kTVal==0||isnan(kTVal))
		{
		  missPoints++;
		  continue;
		}
	      if(kTMean>maxVals[zbin])
		maxValsX[zbin]=kTMean;
	      
	      X[j]=plotters[i]->plotResults[resIdx].kTMeans[j];
	      Y[j]=plotters[i]->plotResults[resIdx].kTValues[j];
	    }
	  cout <<"making graphs " << i <<endl;
	  sprintf(buffer,"graph_%s",names[i].c_str());
	  cout <<"name : " << buffer <<endl;
	  graphs[i][zbin]=new TGraph(numKtBins-missPoints,X,Y);
	  graphs[i][zbin]->SetName(buffer);
	  graphs[i][zbin]->SetMarkerColor(colors[i]);
	  //	  graphs[i][zbin]->SetMarkerSize(2);
	  graphs[i][zbin]->SetMarkerStyle(markerStyles[i]);
	}
    }
  TCanvas c("nonQQbarContributions","nonQQBarContributions",1200,800);
  c.Divide(3,4);
  int minCCount=maxFirstBin;
  if(minCCount>12)
    minCCount=12;
  for(int i=0;i<minCCount;i++)
    {
      TVirtualPad* p=c.cd(i+1);
      graphs[0][i]->GetYaxis()->SetRangeUser(0,maxVals[i]*1.5);
      graphs[0][i]->GetXaxis()->SetRangeUser(0,maxValsX[i]*1.2);
      graphs[0][i]->Draw("AP");      
      for(int j=1;j<4;j++)
	{
	  graphs[j][i]->Draw("SAME P");      
	}
    }
  auto legend=new TLegend(0,0,1.0,1,0);
  for(int i=0;i<4;i++)
    {
      legend->AddEntry(graphs[i][0]->GetName(),names[i].c_str(),"lep");
    }
  c.cd(12);
  legend->Draw();
  c.SaveAs("nonQQContribs.png");
  c.SaveAs("nonQQContribs.pdf");
}





