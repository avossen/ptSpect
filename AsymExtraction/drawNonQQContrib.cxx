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

  vector<string> flavor;

  //  flavor.push_back("_uds");
  //  flavor.push_back("_charm");
  flavor.push_back("_all");


  
  char* eeuuPath=argv[1];
  char* eessPath=argv[2];
  char* eeccPath=argv[3];
  char* tautauPath=argv[4];
  char*[4] paths;
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

  
  string names[4]={"UU","SS","CC","TT"};  
  string plotterName=string("Normal");
  MultiPlotter* plotters[4];
  for(int i=0;i<4;i++)
    {
      myChains[i]=new TChain("PlotTree");
      myChains[i]->Add((string(paths[i])+"/"+plotterName+"_*.root").c_str());
      string fullName=plotterName+string(dataMcNameAdd);
      plotters[i]=new MultiPlotter(m_useQt,const_cast<char*>("."),(fullName+names[i]).c_str(),string(""),0,false,false,false,false,fileTypeEnd);
      plotters[i]->setName(fullName);
      Int_t nevents=myChains[i]->GetEntries();
      PlotResults* plotResults=0;
      myChains[i]->SetBranchAddress("PlotBranch",&plotResults);
	  
      for(long i=0;i<nevents;i++)
	{
	  //	  float locW[3]={0.0,0.0,0.0};
	  myChains[i]->GetEntry(i);
	  plotters[i]->plotResults[plotResults->resultIndex]+=(*plotResults);
	}
    }
}





