
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

//#define MAX_EVENTS 100


using namespace std;

int main(int argc, char** argv)
{
  set<int> zOnlyResIdx;
  set<int> badOnRes;
  set<int> badCont;
//   badOnRes.insert(49);
//  badOnRes.insert(35);
//
//  badCont.insert(7);
//  badCont.insert(15);
//  badCont.insert(17);
//    badCont.insert(23);
//    badCont.insert(25);
// badCont.insert(27);
// badCont.insert(33);
// badCont.insert(55);

//    badCont.insert(49);
//    badCont.insert(67);
//    badCont.insert(63);
//    badCont.insert(69);
//   //
//    badCont.insert(11);
//   //
//    badCont.insert(13);
//    badCont.insert(35);
//    badCont.insert(47);
//    badCont.insert(43);
   badCont.insert(71);




//    badOnRes.insert(71);
//    badOnRes.insert(55);

  //
    //        int minExp=32;
	  int minExp=0;

  char* rootPath=argv[1];
  srand(time(NULL));
  cout <<"Root path is: " << rootPath <<endl;
  string sRootPath(rootPath);
  int expNumber=-1;

  enum flavor{flavUds,flavCharm,flavAll,flavEnd};
  enum onOffRes{res_on,res_off,resEnd};
  //reduced asymmetries for charm, uds, all


  
  setStyleOpts();

  //  vector< pair<string,TChain*> > vFitterNames;
  vector<string> vPlotterNames;
  vPlotterNames.push_back("Normal");
  vPlotterNames.push_back("NormalWoA");



  vector<MultiPlotter*> vPlotters;
  TChain* chAll=0;
  int counter=-1;

  float a[3];
  float ea[3];

  for(vector<string>::iterator it=vPlotterNames.begin();it!=vPlotterNames.end();it++)
    {
      counter++;
      chAll=new TChain("PlotTree");
      chAll->Add((string(rootPath)+"/"+(*it)+"_*.root").c_str());
            cout <<"adding : "<< (string(rootPath)+"/"+(*it)+"_*.root").c_str() <<endl;
      Int_t nevents=chAll->GetEntries();
      cout <<"Plotter Name: " << *it<<endl;
      MultiPlotter* pPlotter=new MultiPlotter(it->c_str(),string(""),0,false,false,false,false);
      pPlotter->setName(*it);
      vPlotters.push_back(pPlotter);
      //has to be 0!!
      PlotResults* plotResults=0;
      chAll->SetBranchAddress("PlotBranch",&plotResults);
      for(int binningType=binType_labTheta_z; binningType<binType_end;binningType++)
	{
	  for(int chargeBin=0;chargeBin<1;chargeBin++)
	    {
	      for(int firstBin=0;firstBin<pPlotter->maxKinMap[binningType].first;firstBin++)
		{
		  for(int secondBin=0;secondBin<pPlotter->maxKinMap[binningType].second;secondBin++)
		    {
		      int resIdx=pPlotter->getResIdx(binningType,chargeBin,firstBin,secondBin);
		      if(binningType==binType_zOnly && chargeBin==pairChargeInt)
			{
			  //			  cout <<"resIdx is " << resIdx<<" firstBin: "<<firstBin<<" second: "<< secondBin<<endl;
			  zOnlyResIdx.insert(resIdx);
			}
		      //		      pPlotter->plotResults[resIdx];
		      //		      pPlotter->plotResults1D[resIdx];
		    }
		}
	    }
	}


      //enough space for 100 exp * on/off resonance
      float** mX=allocateArray<float>(3,200);
      float** mY=allocateArray<float>(3,200);
      float** mXErr=allocateArray<float>(3,200);
      float** mYErr=allocateArray<float>(3,200);
      float** sumWeights=allocateArray<float>(3,200);

      float** mXOnRes=allocateArray<float>(3,200);
      float** mYOnRes=allocateArray<float>(3,200);

      float** mXErrOnRes=allocateArray<float>(3,200);
      float** mYErrOnRes=allocateArray<float>(3,200);
      float** sumWeightsOnRes=allocateArray<float>(3,200);
      for(int i=0;i<200;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      mX[j][i]=i/2;
	      if(i%2)
		mX[j][i]+=0.5;
	      mXErr[j][i]=0.0;
	    }
	}

      for(long i=0;i<nevents;i++)
	{
	  float locW[3]={0.0,0.0,0.0};
	  chAll->GetEntry(i);

	  if(plotResults->exp<=minExp ||(plotResults->on_res && (badOnRes.find(plotResults->exp)!=badOnRes.end())) ||(!plotResults->on_res && (badCont.find(plotResults->exp)!=badCont.end())))
	    continue;


	  //	  if(plotResults->isCharm)
	  //	    continue;
	  //	  cout <<"result index: "<< plotResults->resultIndex <<endl;
	  cout <<"looking at binning type : "<< plotResults->binningType <<endl;
	  if(binType_labTheta_z == plotResults->binningType)
	    {
	      cout <<"kt bins before : ";
	      for(int i=0;i<10;i++)
		{
		  cout << pPlotter->plotResults[plotResults->resultIndex].kTValues[i] <<" mean: " << pPlotter->plotResults[plotResults->resultIndex].kTMeans[i];
		  cout <<"  rhs mean: " << plotResults->kTMeans[i] <<" ";
		}
	    }
	  //	  cout <<endl;
	  pPlotter->plotResults[plotResults->resultIndex]+=(*plotResults);
	  if(binType_labTheta_z == plotResults->binningType)
	    {
	      //	      cout <<"and then...  ";
	      for(int i=0;i<10;i++)
		{
		  cout << pPlotter->plotResults[plotResults->resultIndex].kTValues[i] <<" mean: " << pPlotter->plotResults[plotResults->resultIndex].kTMeans[i];
	    }
	      cout <<endl;
	    }
	}

      pPlotter->savePlots(plotType_2D);

      for(int i=0;i<200;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      if(sumWeights[j][i]>0)
		{
		  mY[j][i]/=sumWeights[j][i];
		  mYErr[j][i]=1/sqrt(sumWeights[j][i]);
		}
	    }
	}

      //recondition:
      int counter[3]={0,0,0};
      for(int i=0;i<200;i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      if(mYErr[j][i]>0)
		{
		  mYErr[j][counter[j]]=mYErr[j][i];
		  mY[j][counter[j]]=mY[j][i];
		  mX[j][counter[j]]=mX[j][i];
		  //error is zero anyways
		  counter[j]++;
		}
	    }
	}


      //only want to plot for the real results
      if((*it)==string("NormalAccWeighted"))
	{
	  	  cout<<"saving results vs exp" <<endl;
	  TFile tmpFile("resVsExpAccWeighted.root","recreate");
	  TGraphErrors tgA1(counter[0],mX[0],mY[0],mXErr[0],mYErr[0]);
	  TH1D thA1("hIff","hIff",100,-10,10);
	  TH1D thA2("hHand","hHand",100,-10,10);
	  TH1D thA3("hG1T","hG1T",100,-10,10);
	  for(int j=0;j<3;j++)
	    {
	  for(int i=0;i<counter[j];i++)
	    {
	      if(mY[j][i]!=0 && mYErr[j][i]!=0)
		{
		  switch(j)
		    {
		    case 0:
		      thA1.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 1:
		      thA2.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    case 2:
		      thA3.Fill(mY[j][i]/mYErr[j][i]);
		      break;
		    default:
		      break;
		    }
		}
	    }
	    }
	  thA1.Write();
	  thA2.Write();
	  thA3.Write();
	  tgA1.SetName("IffVsExp");
	  TGraphErrors tgA2(counter[1],mX[1],mY[1],mXErr[1],mYErr[1]);
	  tgA2.SetName("HandVsExp");
	  TGraphErrors tgA3(counter[2],mX[2],mY[2],mXErr[2],mYErr[2]);
	  tgA3.SetName("G1TVsExp");
	  tgA1.Write();
	  tgA2.Write();
	  tgA3.Write();


	  tmpFile.Write();
	  tmpFile.Close();
	  //	  cout <<"done " <<endl;
	}

      //only want to plot for the real results
      if((*it)==string("NormalWeighted"))
	{
	  	  cout<<"saving results vs exp" <<endl;
	  TFile tmpFile("resVsExpWeighted.root","recreate");




	  tmpFile.Write();
	  tmpFile.Close();
	  //	  cout <<"done " <<endl;
	}
            //only want to plot for the real results
      if((*it)==string("Normal"))
	{
	  cout<<"saving results vs exp" <<endl;
	  TFile tmpFile("resVsExp.root","recreate");
	  //	  cout <<"done " <<endl;
	}

      //      pPlotter->savePlot(binType_m_m,quadPN);
      //      pPlotter->savePlot(binType_z_z,quadPN);

      delete chAll;
    }

  
};
