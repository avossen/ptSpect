
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
#include <stdio.h>
#include <stdlib.h>
//#define MAX_EVENTS 100

using namespace std;

int main(int argc, char** argv)
{
  //calculated from the ratio of woa, no fix mdst mc exp 55
  float uds_charm_ratio=1.634170926;
  
  bool m_useQt=false;
  //should have hadronPairArray where this is defined included
#ifdef USE_QT
  m_useQt=true;
#endif
  vector<string> flavor;

  flavor.push_back("_all");

  set<int> zOnlyResIdx;

  if(argc!=3)
    {
      cout <<"2 arguments needed!!|" <<endl;
      exit(0);
    }

  char* yesISR=argv[1];
  char* noISR=argv[2];

  ifstream listYesISR(yesISR);
  ifstream listNoISR(noISR);

  vector<string> vYesISRFilenames;
  vector<string> vNoISRFilenames;

  string line;
  while(getline(listYesISR,line))
    {
      cout <<" adding " << line <<" to yes " << endl;
      vYesISRFilenames.push_back(line);
    }
  listYesISR.close();
  while(getline(listNoISR,line))
    {
      cout <<" adding " << line <<" to no " << endl;
      vNoISRFilenames.push_back(line);
    }
  cout <<" closing list " <<endl;
  listNoISR.close();
  srand(time(NULL));



  cout <<"got " << vYesISRFilenames.size() <<" and " << vNoISRFilenames.size() <<" filenames " <<endl;
  //assuming that we have to combine uds and charm
  int numPlotters=vYesISRFilenames.size()/2;

  enum flavor{flavUds,flavCharm,flavAll,flavEnd};
  enum onOffRes{res_on,res_off,resEnd};
  //reduced asymmetries for charm, uds, all

  setStyleOpts();

  //  vector< pair<string,TChain*> > vFitterNames;
  vector<string> vPlotterNames;
  vPlotterNames.push_back("Normal");
  string woPlotterName("NormalWoA");

  vector<MultiPlotter*> vPlottersYes;
  vector<MultiPlotter*> vPlottersNo;
  TChain* chAllYes=0;
  TChain* chAllNo=0;
  int counter=-1;

  //  float a[3];
  //  float ea[3];
  cout <<" going through flavors " <<endl;
  for(vector<string>::iterator itFlav=flavor.begin();itFlav!=flavor.end();itFlav++)
    {
      cout <<"flavor: "<< *itFlav<<endl;
      for(int iP=0;iP<numPlotters;iP++)
	{
	  cout << iP << " of " << numPlotters << endl;
      //woA plotter
	  chAllYes=new TChain("PlotTree");
	  chAllNo=new TChain("PlotTree");
	  //assuming that this is for uds and charm
	  chAllYes->Add(vYesISRFilenames[2*iP].c_str());
	  chAllYes->Add(vYesISRFilenames[2*iP+1].c_str());
	  chAllNo->Add(vNoISRFilenames[2*iP].c_str());
	  chAllNo->Add(vNoISRFilenames[2*iP+1].c_str());
  cout <<"adding : "<< vYesISRFilenames[2*iP] <<" and " << vNoISRFilenames[2*iP]<<endl;
  cout <<"also adding : "<< vYesISRFilenames[2*iP+1] <<" and " << vNoISRFilenames[2*iP+1]<<endl;
	  Int_t nevents=chAllYes->GetEntries();
	  if(nevents!=chAllNo->GetEntries())
	    {
	      cout <<"not same number of entries " << endl;
	      exit(0);
	    }
	
	  cout <<"Plotter Name: " << woPlotterName<<endl;
	  char strIp[100];
	  sprintf(strIp,"%d",iP);
	  string fullNameYes=(woPlotterName)+(*itFlav)+"Yes_plotterNum_"+string(strIp);
	  string fullNameNo=(woPlotterName)+(*itFlav)+"No_plotterNum_"+string(strIp);
	  //I guess this doesn't need an output path
	  MultiPlotter* pWoAPlotterYes=new MultiPlotter(m_useQt,const_cast<char*>("."),fullNameYes.c_str(),string(""),0,false,false,false,false);	  
	  MultiPlotter* pWoAPlotterNo=new MultiPlotter(m_useQt,const_cast<char*>("."),fullNameNo.c_str(),string(""),0,false,false,false,false);	  


	  pWoAPlotterYes->setName(fullNameYes);
	  pWoAPlotterNo->setName(fullNameNo);
	  //has to be 0!!
	  PlotResults* plotResultsYes=0;
	  PlotResults* plotResultsNo=0;

	  chAllYes->SetBranchAddress("PlotBranch",&plotResultsYes);
	  chAllNo->SetBranchAddress("PlotBranch",&plotResultsNo);

	  for(long i=0;i<nevents;i++)
	    {
	      //	  float locW[3]={0.0,0.0,0.0};
	      chAllYes->GetEntry(i);
	      chAllNo->GetEntry(i);

	      //for now only continuum:
	      //	      cout <<"onres? "<< endl;
	      //'yes' and 'no' should be synchronized
	      if(plotResultsYes->on_res)
		{
		  cout <<" result is on resonance " <<endl;
		  continue;
		}

	      //check if the  result we are reading right now is compatible with the 
	      //flavor we are looking for
	      //anyways only doing "all" here...
	      if((*itFlav)==string("_charm"))
		{
		  if(!plotResultsYes->isCharm)
		    continue;
		}
	      if((*itFlav)==string("_uds"))
		{
		  if(plotResultsYes->isCharm)
		    continue;
		}

//	      if(plotResultsYes->isUds)
//		{
//		  cout <<"looking at uds, weighting..." <<endl;
//		  plotResultsYes->weight(uds_charm_ratio);
//		}
//	      if(plotResultsNo->isUds)
//		{
//		  	  cout <<"looking at uds, weighting..." <<endl;
//			  plotResultsNo->weight(uds_charm_ratio);
//		}
//	      if(!plotResultsYes->isUds || !plotResultsNo->isUds)
//		{
//		  cout <<" this is charm, no weighting " <<endl;
//		}
	      //the weighting is a function of the plotter. So if we need it, we need to keep the plotter for
	      //uds and charm separate
	      
	      pWoAPlotterYes->plotResults[plotResultsYes->resultIndex]+=(*plotResultsYes);
	      pWoAPlotterNo->plotResults[plotResultsNo->resultIndex]+=(*plotResultsNo);
	    }
	  vPlottersYes.push_back(pWoAPlotterYes);
 	  vPlottersNo.push_back(pWoAPlotterNo);
	}
    }
  
  int maxCombBin=10000;
  float maxRatio[maxCombBin];
  float minRatio[maxCombBin];

  //only 0 and 1 valid for the getHistogram
  //however, it looks like only z_z (so the binning index 0 is saved in by MultiPlotter)
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

	      for(int iP=0;iP<numPlotters;iP++)
		{
		  //histo dimensions should already be correct for that pid
		  TH1D* histoY=vPlottersYes[iP]->getHistogram(binningType,chargeBin,pidBin);
		  TH1D* histoN=vPlottersNo[iP]->getHistogram(binningType,chargeBin,pidBin);
		  
		  maxX=histoY->GetNbinsX();
		  if(maxX!=histoN->GetNbinsX())
		    {
		      cout <<" y and n not the same dimension!" <<endl;
		      cout <<"no: "<< histoN->GetNbinsX() <<" yes: "<< maxX <<endl;
		      exit(0);
		    }

		  for(int iX=0;iX<histoY->GetNbinsX();iX++)
		    {
		      cout <<" we have " << histoY->GetNbinsX() <<" entries " <<endl;
		      float wISR=histoY->GetBinContent(iX+1);
		      float woISR=histoN->GetBinContent(iX+1);
		      //    cout <<"w ISR: " << wISR <<" wo: " << woISR <<endl;
		      //at least at high  wISR should be much less
			      
		      if(woISR>0)
			{
			  float ratio=wISR/woISR;
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

}
