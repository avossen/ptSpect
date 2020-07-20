#include <stdio.h>
#include <stdlib.h>
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
#include "HadronPairArray.h"
#include "CombPlots.h"
//#define MAX_EVENTS 100

using namespace std;

int main(int argc, char** argv)
{
  char buffer[400];
  enum fileType{udsyes, udsno, charmyes, charmno,fileTypeEnd};  
  bool m_useQt=false;
  //should have hadronPairArray where this is defined included
#ifdef USE_QT
  m_useQt=true;
#endif

  if(argc<5)
    {
      cout <<"usage: drawRadSysContrib yesUDS noUDS yesCharm noCharm"<<endl;
      exit(0);
    }

  vector<string> flavor;

  //  flavor.push_back("_uds");
  //  flavor.push_back("_charm");
  flavor.push_back("_all");

  cout <<"argvs..." <<endl;

  char* yesISRUDS=argv[1];
  char* noISRUDS=argv[2];
  char* yesISRCharm=argv[3];
  char* noISRCharm=argv[4];


  ifstream listYesISRUDS(yesISRUDS);
  ifstream listNoISRUDS(noISRUDS);
  ifstream listYesISRCharm(yesISRCharm);
  ifstream listNoISRCharm(noISRCharm);


  vector<string> vYesISRCharmFilenames;
  vector<string> vNoISRCharmFilenames;
  vector<string> vYesISRUDSFilenames;
  vector<string> vNoISRUDSFilenames;

  vector<string>* fileLists[4];
  fileLists[udsyes]=&vYesISRUDSFilenames;
  fileLists[udsno]=&vNoISRUDSFilenames;
  fileLists[charmyes]=&vYesISRCharmFilenames;
  fileLists[charmno]=&vNoISRCharmFilenames;
  
  string line;
  while(getline(listYesISRUDS,line))
    {
      cout <<" adding " << line <<" to yes " << endl;
      vYesISRUDSFilenames.push_back(line);
    }
  listYesISRUDS.close();
  while(getline(listNoISRUDS,line))
    {
      cout <<" adding " << line <<" to no " << endl;
      vNoISRUDSFilenames.push_back(line);
    }
  listNoISRUDS.close();
    while(getline(listYesISRCharm,line))
    {
      cout <<" adding " << line <<" to yes " << endl;
      vYesISRCharmFilenames.push_back(line);
    }
  listYesISRCharm.close();
  while(getline(listNoISRCharm,line))
    {
      cout <<" adding " << line <<" to no " << endl;
      vNoISRCharmFilenames.push_back(line);
    }
  listNoISRCharm.close();

  int numFiles=vYesISRUDSFilenames.size();
  srand(time(NULL));
  const char** paths[4];
  for(int i=udsyes;i<fileTypeEnd;i++)
    {
      paths[i]=new const char*[numFiles];
      for(int j=0;j<numFiles;j++)
	{
	  paths[i][j]=(*fileLists[i])[j].c_str();
	}
    }

  srand(time(NULL));


  //reduced asymmetries for charm, uds, all

  setStyleOpts();

  //  vector< pair<string,TChain*> > vFitterNames;
  vector<string> vPlotterNames;


  vector<MultiPlotter*> vPlotters;
  int counter=-1;

  TChain** myChains[4];
  cout <<"after chains " <<endl;
  string names[]={"uds","charm"};  
  string plotterName=string("NormalWoA");
  MultiPlotter** plotters[4];
  for(int i=0;i<4;i++)
    {
      myChains[i]=new TChain*[numFiles];
      plotters[i]=new MultiPlotter*[numFiles];
    }
  
  
  char dataMcNameAdd[100];
  sprintf(dataMcNameAdd,"");
  cout <<"before the loop" <<endl;
  for(int i=udsyes;i<fileTypeEnd;i++)
    {
      for(int j=0;j<numFiles;j++)
	{
	  cout <<"accessing chain " << i <<  " file : "<< j <<endl;
	  cout <<" adding: "<< (string(paths[i][j])+"/"+plotterName+"_*.root").c_str() <<endl;
	  myChains[i][j]=new TChain("PlotTree");
	  myChains[i][j]->Add((string(paths[i][j])+"/"+plotterName+"_*.root").c_str());



	  sprintf(buffer,"%s_%s_%i_%i",plotterName.c_str(),dataMcNameAdd,i,j);
	  string fullName=string(buffer);
	  plotters[i][j]=new MultiPlotter(m_useQt,const_cast<char*>("."),(fullName+names[i]).c_str(),string(""),0,false,false,false,false,fileTypeEnd);
	  plotters[i][j]->setName(fullName);
	  Int_t nevents=myChains[i][j]->GetEntries();
	  PlotResults* plotResults=0;
	  myChains[i][j]->SetBranchAddress("PlotBranch",&plotResults);
	  cout <<" running over " << nevents << endl;
	  for(long k=0;k<nevents;k++)
	    {
	  //	  cout <<"getting entry " << i <<endl;
	  //	  float locW[3]={0.0,0.0,0.0};
	      myChains[i][j]->GetEntry(k);
	      plotters[i][j]->plotResults[plotResults->resultIndex]+=(*plotResults);
	    }
	  cout <<"done with entries.." <<endl;
	}
    }

  //should have now filled all plotters, do the histos

  int binningType=binType_z_z;
    int pidBin=0; //pipi
  //    int pidBin=2; //pipi
  int chargeBin=likesign; //same charge

  int maxFirstBin=plotters[0][0]->maxKinMap[pidBin][binningType].first;
  int maxSecondBin=plotters[0][0]->maxKinMap[pidBin][binningType].second;
  int resIdx=plotters[0][0]->getResIdx(binningType,pidBin,chargeBin,0,0);
  int numKtBins=plotters[0][0]->plotResults[resIdx].numKtValues;
  TGraph*** graphs[4];
  TGraph*** graphsWeak[4];
  for(int i=0;i<4;i++)
    {
      graphsWeak[i]=new TGraph**[numFiles];
      graphs[i]=new TGraph**[numFiles];
      for(int j=0;j<numFiles;j++)
	{
	      graphs[i][j]=new TGraph*[maxFirstBin];
	      graphsWeak[i][j]=new TGraph*[maxFirstBin];
	}
    }
    
  cout <<"maxfirst: "<< maxFirstBin <<" second: " << maxSecondBin <<" resIdx: "<< resIdx <<" numkt: " << numKtBins <<endl;

  double* X=new double[numKtBins];
  double* Y=new double[numKtBins];
  double* YWeak=new double[numKtBins];

  int colors[]={kRed,kBlue,kGreen,kBlack,kCyan};
  int markerStyles[]={20,21,22,23,24};

  //max ktVals for each z bin to scale the graphs
  double minVals[100];
  double maxVals[100];
  double maxValsX[100];
  for(int i=0;i<100;i++)
    {
      minVals[i]=1000;
      maxVals[i]=-1000;
      maxValsX[i]=0;
    }
  //0 uds, 1 charm
  cout <<"running over " << numKtBins <<" kt bins " <<endl;
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<numFiles;j++)
	{
	  for(int zbin=0;zbin<maxFirstBin;zbin++)
	    {
	      cout <<"looking at uds/charm: "<< i << " file: "<< j <<", zbin : "<< zbin <<endl;
	      //diagonal bins
	      resIdx=plotters[0][0]->getResIdx(binningType,pidBin,chargeBin,zbin,zbin);
	      cout <<"residx2 : "<< resIdx <<endl;
	      int missPoints=0;

	      for(int k=0;k<numKtBins;k++)
		{
		  //    cout <<"filling x,y num " << j <<" with mean: " << plotters[i]->plotResults[resIdx].kTMeans[j] <<" value: "<< plotters[i]->plotResults[resIdx].kTValues[j]  <<endl;
		  //look at the no radiation (i+1) 
		  double kTMean=plotters[2*i+1][j]->plotResults[resIdx].kTMeans[k];
		  double kTVal=plotters[2*i+1][j]->plotResults[resIdx].kTValues[k];
		  double kTValWRad=plotters[2*i][j]->plotResults[resIdx].kTValues[k];
		  cout <<"ktMean: "<< kTMean <<" ktVal: "<< kTVal <<" ktValWRad: " << kTValWRad <<endl;
		  if(kTMean==0 || isnan(kTMean)|| kTVal==0||isnan(kTVal) )
		    {
		      missPoints++;
		      continue;
		    }

		  //		  double ratio=(kTVal-kTValWRad)/kTVal;
		  //Ralf shows it this way, so do the same for comparison
		  double ratio=kTVal/kTValWRad;
		  float weakRatio=plotters[2*i+1][j]->plotResults[resIdx].weakDecayFraction[k];
		  cout <<"zbin: "<< zbin <<" ktBin: " << k <<" mean: "<< kTMean << "ratio: "<< ratio <<" weakratio: " << weakRatio <<endl;
		  if(ratio>maxVals[zbin])
		    maxVals[zbin]=ratio;
		  if(ratio<minVals[zbin])
		    {
		      minVals[zbin]=ratio;
		      cout <<" setting new minVal for zbin " << zbin <<" ratio: "<< ratio <<" ktbin: "<< k <<" ktmean: "<< kTMean<<endl;
		    }

		  if(kTMean>maxVals[zbin])
		    maxValsX[zbin]=kTMean;
		  
		  X[k]=kTMean;
		  Y[k]=ratio;
		  YWeak[k]=weakRatio;
		}
	      cout <<"making graphs " << i <<endl;
	      sprintf(buffer,"graph_%s_file%i",names[i].c_str(),j);
	      cout <<"name : " << buffer <<endl;
	      graphs[i][j][zbin]=new TGraph(numKtBins-missPoints,X,Y);
	      graphs[i][j][zbin]->SetName(buffer);
	      graphs[i][j][zbin]->SetMarkerColor(colors[j%5]);
	      //	  graphs[i][zbin]->SetMarkerSize(2);
	      graphs[i][j][zbin]->SetMarkerStyle(markerStyles[j%5]);
	      
	      sprintf(buffer,"graphWeak_%s_file%i",names[i].c_str(),j);
	      graphsWeak[i][j][zbin]=new TGraph(numKtBins-missPoints,X,YWeak);
	      graphsWeak[i][j][zbin]->SetName(buffer);
	      graphsWeak[i][j][zbin]->SetMarkerColor(colors[j%5]);
	      //	  graphs[i][zbin]->SetMarkerSize(2);
	      graphsWeak[i][j][zbin]->SetMarkerStyle(markerStyles[j%5]);

	    }
	}
    }

    
  TCanvas c("radSys","radSys",1200,800);
  c.Divide(3,4);
  int minCCount=maxFirstBin;
  if(minCCount>12)
    minCCount=12;

  for(int uc=0;uc<2;uc++)
    {
      for(int i=0;i<minCCount;i++)
	{
	  TVirtualPad* p=c.cd(i+1);
      //      p->SetLogy();

	  
	  graphs[uc][0][i]->GetYaxis()->SetRangeUser(minVals[i]-0.5*fabs(minVals[i]),maxVals[i]+0.5*fabs(maxVals[i]));
	  graphs[uc][0][i]->GetXaxis()->SetRangeUser(0,maxValsX[i]*1.2);
	  graphs[uc][0][i]->Draw("AP");      
	  for(int j=1;j<numFiles;j++)
	    {
	      graphs[uc][j][i]->Draw("SAME P");      
	    }
	}
      auto legend=new TLegend(0,0,1.0,1,0);
      for(int i=0;i<numFiles;i++)
	{
	  legend->AddEntry(graphs[uc][i][0]->GetName(),names[uc].c_str(),"lep");
	}
      c.cd(12);
      legend->Draw();
      if(uc==0)
	sprintf(buffer,"radCorr_uds.png");
      else
	sprintf(buffer,"radCorr_charm.png");
      c.SaveAs(buffer);
    }

  ///weak stuff:

  for(int uc=0;uc<2;uc++)
    {
      for(int i=0;i<minCCount;i++)
	{
	  TVirtualPad* p=c.cd(i+1);
      //      p->SetLogy();
	  graphsWeak[0][0][i]->GetYaxis()->SetRangeUser(0,1.0);
	  graphsWeak[0][0][i]->GetXaxis()->SetRangeUser(0,maxValsX[i]*1.2);
	  graphsWeak[0][0][i]->Draw("AP");      
	  for(int j=1;j<numFiles;j++)
	    {
	      graphsWeak[0][j][i]->Draw("SAME P");      
	    }
	}
      auto legend=new TLegend(0,0,1.0,1,0);
      for(int i=0;i<numFiles;i++)
	{
	  legend->AddEntry(graphsWeak[uc][i][0]->GetName(),names[uc].c_str(),"lep");
	}
      c.cd(12);
      legend->Draw();
      if(uc==0)
	sprintf(buffer,"weakDecay_uds.png");
      else
	sprintf(buffer,"weakDecay_charm.png");
      c.SaveAs(buffer);
    }


  
}





