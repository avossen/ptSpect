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

//#define DEBUG

#ifdef DEBUG
float uds_charm_ratio=1.0;
#else
float uds_charm_ratio=1.634170926;
#endif

//float uds_charm_ratio=1;
using namespace std;

int main(int argc, char** argv)
{

  cout <<"unlikesign: "<< unlikesign <<" likesign: "<<likesign <<endl;
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
  //added the combined as the last index
  string names[]={"uds","charm","nd","nd","combined"};  
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
      cout <<"file type:"<< i <<endl;
      for(int j=0;j<numFiles;j++)
	{
	  cout <<"file number: "<< j <<endl;
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
	      cout <<"combining result index: "<< plotResults->resultIndex <<endl;
	      cout << "binning: "<< plotters[i][j]->getBinningTypeFromIdx(plotResults->resultIndex)<<" "<<plotResults->binningType;
	      cout << " pid: "<< plotters[i][j]->getPidTypeFromIdx(plotResults->resultIndex)<<" " << plotResults->pidBin;
	      cout << " charge: "<< plotters[i][j]->getChargeTypeFromIdx(plotResults->resultIndex)<<" " << plotResults->chargeBin;
	      cout << " first kin: "<< plotters[i][j]->getFirstKinBinFromIdx(plotResults->resultIndex) << " " << plotResults->firstKinBin;
	      cout << " second kin: "<< plotters[i][j]->getSecondKinBinFromIdx(plotResults->resultIndex) << " " << plotResults->secondKinBin;
	      plotters[i][j]->plotResults[plotResults->resultIndex]+=(*plotResults);
	    }
	  cout <<"done with entries.." <<endl;
	  if(i==udsyes || i==udsno)
	    {
	      cout <<"weighting uds " << endl;
	      plotters[i][j]->weight(uds_charm_ratio);
	    }
	  else
	    {
	      cout <<"charm, not weighted " << endl;
	    }
	}
    }
  //should have now filled all plotters, 
  ///calculate ratios 
  int maxCombBin=10000;
  float maxRatio[maxCombBin];
  float minRatio[maxCombBin];

  //only 0 and 1 valid for the getHistogram
  //however, it looks like only z_z (so the binning index 0 is saved in by MultiPlotter)
  //not thetat this is a speciality of 'getHistogram', in the more general case (e.g. used for resIdx, the index
  //for z_z is 6 but in the getHistogram function that gets translated
  
  for(int binningType=0; binningType<1;binningType++)
    {
      for(int pidBin=0;pidBin<3;pidBin++)
	{
	  for(int chargeBin=0;chargeBin<2;chargeBin++)
	    {
	      for(int i=0;i<maxCombBin;i++)
		{
		  maxRatio[i]=-1000.0;
		  minRatio[i]=1000.0;
		}
	      //maxX should be ==maxY
	      int maxX=0;

	      for(int iF=0;iF<numFiles;iF++)
		{

		  ///
		  sprintf(buffer,"bt_%d_pB_%d_cp_%d_if_%d",binningType,pidBin,chargeBin,iF);
		  
		  //histo dimensions should already be correct for that pid
		  TH1D* histoY_UDS=plotters[udsyes][iF]->getHistogram(binningType,chargeBin,pidBin,buffer);
		  TH1D* histoN_UDS=plotters[udsno][iF]->getHistogram(binningType,chargeBin,pidBin,buffer);
		  TH1D* histoY_Charm=plotters[charmyes][iF]->getHistogram(binningType,chargeBin,pidBin,buffer);
		  TH1D* histoN_Charm=plotters[charmno][iF]->getHistogram(binningType,chargeBin,pidBin,buffer);
		  
		  maxX=histoY_UDS->GetNbinsX();
		  if(maxX!=histoN_UDS->GetNbinsX())
		    {
		      cout <<" y and n not the same dimension!" <<endl;
		      cout <<"no: "<< histoN_UDS->GetNbinsX() <<" yes: "<< maxX <<endl;
		      exit(0);
		    }

		  for(int iX=0;iX<histoY_UDS->GetNbinsX();iX++)
		    {
		      cout <<" we have " << histoY_UDS->GetNbinsX() <<" entries " <<endl;
		      if(iX==0)
			{
			  cout <<"uds w rad: "<< histoY_UDS->GetBinContent(iX+1) <<" charm w rad : " << histoY_Charm->GetBinContent(iX+1) <<endl;
			  			  cout <<"uds wo rad: "<< histoN_UDS->GetBinContent(iX+1) <<" charm wo rad : " << histoN_Charm->GetBinContent(iX+1) <<endl;
			}
		      float wISR=histoY_UDS->GetBinContent(iX+1)+histoY_Charm->GetBinContent(iX+1);
		      #ifndef DEBUG
		      float woISR=histoN_UDS->GetBinContent(iX+1)+histoN_Charm->GetBinContent(iX+1);
		      #else
		      	float	      woISR=histoN_Charm->GetBinContent(iX+1);
#endif
		      //		      float woISR=histoN->GetBinContent(iX+1);
		      //    cout <<"w ISR: " << wISR <<" wo: " << woISR <<endl;
		      //at least at high  wISR should be much less

	      
		      if(woISR>0)
			{
			  #ifndef DEBUG
			  float ratio=woISR/wISR;
			  if(iX==0)
			    cout <<"ratio1: "<< ratio <<endl;

			  #else
			  		  float ratio=woISR;
			  #endif
			  if(ratio>maxRatio[iX])
			    {
			      maxRatio[iX]=ratio;
			    }
			  if(ratio<minRatio[iX])
			    {
			      minRatio[iX]=ratio;
			    }
			
			}
		    }
		}

	      ofstream f;
	      char buffer[300];
	      sprintf(buffer,"ISRCorr_binning_%d_charge_%d_pid_%d.txt",binningType,chargeBin,pidBin);
	      f.open(buffer);
	      for(int iX=0;iX<maxX;iX++)
		{
		  f << maxRatio[iX] << " " << minRatio[iX] << " ";
		}

	    }

	}
    }


  ///  
  //do the histos
  int binningType=binType_z_z;
  int pidBin=0; //pipi
  //    int pidBin=2; //pipi
  int chargeBin=likesign; //same charge
  int maxFirstBin=plotters[0][0]->maxKinMap[pidBin][binningType].first;
  int maxSecondBin=plotters[0][0]->maxKinMap[pidBin][binningType].second;
  int resIdx=plotters[0][0]->getResIdx(binningType,pidBin,chargeBin,0,0);
  int numKtBins=plotters[0][0]->plotResults[resIdx].numKtValues;
  //not sure why the dimension is 4 but let's use i==4 (so changed max to 5) for the combined (0,1 corresponds to uc=1, might use 2,3 for mixed,...)
  //changed
  TGraph*** graphs[5];
  TGraph*** graphsWeak[5];
  for(int i=0;i<5;i++)
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
  //there is not really a point for running till 5. We only use 0,1 and 4 signifiying uds,charm and all.
  //it is really not used to index the plotters (udsyes, udsno etc... I guess 2,3 could be with radiation...
  ofstream fweak;
  sprintf(buffer,"WeakFraction_binning_%d_charge_%d_pid_%d.txt",binningType,chargeBin,pidBin);
  fweak.open(buffer);

    for(int i=0;i<5;i++)
    {
      if(i==2 || i==3)
	{
	  //not implemented yet
	  continue;
	}
      //this would be the different tunes
      for(int j=0;j<numFiles;j++)
	{
	  for(int zbin=0;zbin<maxFirstBin;zbin++)
	    {
	      int missPoints=0;
	      for(int zbin2=0;zbin2<maxSecondBin;zbin2++)
	      //	      int zbin2=zbin;
		{
		  cout <<"looking at uds/charm: "<< i << " file: "<< j <<", zbin : "<< zbin <<endl;
		  //diagonal bins
		  resIdx=plotters[0][0]->getResIdx(binningType,pidBin,chargeBin,zbin,zbin2);
		  cout <<"residx2 : "<< resIdx <<endl;


		  for(int k=0;k<numKtBins;k++)
		    {
		      //    cout <<"filling x,y num " << j <<" with mean: " << plotters[i]->plotResults[resIdx].kTMeans[j] <<" value: "<< plotters[i]->plotResults[resIdx].kTValues[j]  <<endl;
		      //look at the no radiation (i+1)
		      double kTMean=0;
		      double kTVal=0;
		      double kTValWRad=0;
		      double weakRatio=0;
		      
		      float kTValU=0;
		      float kTValC=0;
		      float weakRatioU=0;
		      float weakRatioC=0;
		      if(k==0 && i==4)
			{
			  cout <<"first bin uds w isr " << plotters[udsyes][j]->plotResults[resIdx].kTValues[k];
			  cout <<"charm  w isr " << plotters[charmyes][j]->plotResults[resIdx].kTValues[k];
			  cout <<"first bin uds wo isr " << plotters[udsno][j]->plotResults[resIdx].kTValues[k];
			  cout <<"charm  wo isr " << plotters[charmno][j]->plotResults[resIdx].kTValues[k];
			  
			}
		      //i:0 udsyes, 1: udsno ...4: doesn't exist in plotters, I guess introduced here as 'combined'
		      switch(i)
			{
			  //2*i+1 is 3 for case 1, so charmno, for i=0 it would be udsno,
			  //so interpreting i as uds/charm it gives uds/charm no
			case 0:
			case 1:
			  kTMean=plotters[2*i+1][j]->plotResults[resIdx].kTMeans[k];
			  kTVal=plotters[2*i+1][j]->plotResults[resIdx].kTValues[k];
			  kTValWRad=plotters[2*i][j]->plotResults[resIdx].kTValues[k];
			  weakRatio=plotters[2*i+1][j]->plotResults[resIdx].weakDecayFraction[k];
		      break;
			case 4:
			  //just take one of the uds/charm for the means--> [1], [3] corresponds to udsno, charmo
			  kTMean=plotters[udsno][j]->plotResults[resIdx].kTMeans[k];
			  kTValU=plotters[udsno][j]->plotResults[resIdx].kTValues[k];
			  kTValC=plotters[charmno][j]->plotResults[resIdx].kTValues[k];
			  //kTvalues should already be correctly weighted
			  kTVal=kTValU+kTValC;
			  //[0] and [2] are udsyes, charmyes
			  kTValWRad=plotters[udsyes][j]->plotResults[resIdx].kTValues[k]+plotters[charmyes][j]->plotResults[resIdx].kTValues[k];
			  weakRatioU=plotters[udsno][j]->plotResults[resIdx].weakDecayFraction[k];
			  weakRatioC=plotters[charmno][j]->plotResults[resIdx].weakDecayFraction[k];
			  //I guess we could have just added the plotters...
			  weakRatio=(weakRatioU*kTValU+weakRatioC*kTValC)/(kTValU+kTValC);
			  //write weakratio to file
			  fweak <<" " << weakRatio;
			  break;
			default:
			  cout <<"i is wrong " <<endl;
			}
		  
		      cout <<"ktMean: "<< kTMean <<" ktVal: "<< kTVal <<" ktValWRad: " << kTValWRad <<endl;
		      if(kTMean==0 || isnan(kTMean)|| kTVal==0||isnan(kTVal) || kTValWRad==0)
			{
			  cout <<"invalid point " <<endl;
			  
			  //misspoints is only important for the graph where we only do diagonal points..
			  if(zbin==zbin2)
			    {
			      missPoints++;
			    }
			  continue;
		      
			}
		
		      //		  double ratio=(kTVal-kTValWRad)/kTVal;
		      //Ralf shows it this way, so do the same for comparison
		      double ratio=kTVal/kTValWRad;
		      if(k==0 && i==4)
			cout <<"ratio2: "<< ratio <<endl;



		      //for plots, only look at diagonal
		      if(zbin==zbin2)
			{
			  cout <<"zbin: "<< zbin <<" ktBin: " << k <<" mean: "<< kTMean << "ratio: "<< ratio <<" weakratio: " << weakRatio <<endl;			  
			  if(ratio>maxVals[zbin])
			    {
			      cout <<"setting maxVal to " << ratio <<endl;
			      maxVals[zbin]=ratio;
			    }
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
		  
		    }

		}
	      cout <<"making graphs " << i <<endl;
	      sprintf(buffer,"graph_%s_file%i_%i",names[i].c_str(),j,i);
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


  cout <<"minCCount: "<< minCCount <<endl;
  for(int uc=0;uc<5;uc++)
    {
      if(uc==2 || uc==3)
	{
	  //not implemented yet
	  continue;
	}

      for(int i=0;i<minCCount;i++)
	{
	  if(i==6)
	    {
	      cout <<"minVals: "<< minVals[i] <<" max vals : "<< maxVals[i] <<endl;
	    }
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
      if(uc==1)
	sprintf(buffer,"radCorr_charm.png");
      if(uc==4)
	sprintf(buffer,"radCorr_all.png");
      c.SaveAs(buffer);
    }

  ///weak stuff:

  for(int uc=0;uc<5;uc++)
    {
      if(uc==2 || uc==3)
	{
	  //not implemented yet
	  continue;
	}

      for(int i=0;i<minCCount;i++)
	{
	  TVirtualPad* p=c.cd(i+1);
      //      p->SetLogy();
	  graphsWeak[uc][0][i]->GetYaxis()->SetRangeUser(0,1.0);
	  graphsWeak[uc][0][i]->GetXaxis()->SetRangeUser(0,maxValsX[i]*1.2);
	  graphsWeak[uc][0][i]->Draw("AP");      
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
      if(uc==1)
	sprintf(buffer,"weakDecay_charm.png");
      if(uc==4)
	sprintf(buffer,"weakDecay_all.png");

      c.SaveAs(buffer);
    }


  
}





