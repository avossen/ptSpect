
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

  bool m_useQt=false;
  //should have hadronPairArray where this is defined included
#ifdef USE_QT
  m_useQt=true;
#endif

  vector<string> flavor;

  //  flavor.push_back("_uds");
  //  flavor.push_back("_charm");
  flavor.push_back("_all");

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
  //   badCont.insert(71);

  //    badOnRes.insert(71);
  //    badOnRes.insert(55);



  map<int,int> expCounts;


  //
  int minExp=12;
  //int minExp2=0;
  //  int maxExp=2;
  //  int maxExp2=200;



  int maxExp=35;
  int minExp2=55;
  int maxExp2=69;
  // int minExp=0;

  char* rootPath=argv[1];
  char dataMcNameAdd[100];
  sprintf(dataMcNameAdd,"");
  if(argc>2)
    sprintf(dataMcNameAdd,"%s",argv[2]);

  srand(time(NULL));
  cout <<"Root path is: " << rootPath <<endl;
  string sRootPath(rootPath);
  //  int expNumber=-1;

  enum flavor{flavUds,flavCharm,flavAll,flavEnd};
  enum onOffRes{res_on,res_off,resEnd};
  //reduced asymmetries for charm, uds, all


  
  setStyleOpts();

  //  vector< pair<string,TChain*> > vFitterNames;
  vector<string> vPlotterNames;
  vPlotterNames.push_back("Normal");
  //    vPlotterNames.push_back("NormalWoA");

  vector<MultiPlotter*> vPlotters;
  TChain* chAll=0;
  int counter=-1;

  //  float a[3];
  //  float ea[3];


  for(vector<string>::iterator itFlav=flavor.begin();itFlav!=flavor.end();itFlav++)
    {

      for(vector<string>::iterator it=vPlotterNames.begin();it!=vPlotterNames.end();it++)
	{
	  counter++;
	  chAll=new TChain("PlotTree");
	  chAll->Add((string(rootPath)+"/"+(*it)+"_*.root").c_str());
	  cout <<"adding : "<< (string(rootPath)+"/"+(*it)+"_*.root").c_str() <<endl;
	  Int_t nevents=chAll->GetEntries();
	  cout <<"Plotter Name: " << *it<<endl;
	  string fullName=(*it)+string(dataMcNameAdd)+(*itFlav);
	  //I guess this doesn't need an output path
	  MultiPlotter* pPlotter=new MultiPlotter(m_useQt,const_cast<char*>(""),fullName.c_str(),string(""),0,false,false,false,false);

	  cout <<" setting plotter name to:" << fullName <<endl;
	  pPlotter->setName(fullName);
	  vPlotters.push_back(pPlotter);
	  //has to be 0!!
	  PlotResults* plotResults=0;
	  chAll->SetBranchAddress("PlotBranch",&plotResults);
	  for(int binningType=binType_labTheta_z; binningType<binType_end;binningType++)
	    {
	      for(int pidBin=0;pidBin<3;pidBin++)
		{
		  for(int chargeBin=0;chargeBin<1;chargeBin++)
		    {
		      for(int firstBin=0;firstBin<pPlotter->maxKinMap[pidBin][binningType].first;firstBin++)
			{
			  for(int secondBin=0;secondBin<pPlotter->maxKinMap[pidBin][binningType].second;secondBin++)
			    {
			      int resIdx=pPlotter->getResIdx(binningType,pidBin,chargeBin,firstBin,secondBin);
			      if(binningType==binType_zOnly && chargeBin==unlikesign)
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
	    }

	  //enough space for 100 exp * on/off resonance
	  float** mX=allocateArray<float>(3,200);
	  float** mY=allocateArray<float>(3,200);
	  float** mXErr=allocateArray<float>(3,200);
	  float** mYErr=allocateArray<float>(3,200);
	  float** sumWeights=allocateArray<float>(3,200);

	  //      float** mXOnRes=allocateArray<float>(3,200);
	  //      float** mYOnRes=allocateArray<float>(3,200);

	  //      float** mXErrOnRes=allocateArray<float>(3,200);
	  //      float** mYErrOnRes=allocateArray<float>(3,200);
	  //      float** sumWeightsOnRes=allocateArray<float>(3,200);
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
	      //	  float locW[3]={0.0,0.0,0.0};
	      chAll->GetEntry(i);

	      if(plotResults->exp<minExp ||(plotResults->on_res && (badOnRes.find(plotResults->exp)!=badOnRes.end())) ||(!plotResults->on_res && (badCont.find(plotResults->exp)!=badCont.end())))
		continue;
	      if(plotResults->exp>maxExp && plotResults->exp < minExp2)
		continue;
	      if(plotResults->exp>maxExp2)
		continue;

	      //entry doesn't exist yet
	      if(expCounts.find(plotResults->exp)==expCounts.end())
		{
		  expCounts[plotResults->exp]=0;
		}


	      //for now only continuum:
	      cout <<"onres? "<< endl;
	      if(plotResults->on_res)
		{
		  cout <<" result is on resonance " <<endl;
		  continue;
		}
	      if((*itFlav)==string("_charm"))
		{
		  if(!plotResults->isCharm)
		    continue;
		}
	      if((*itFlav)==string("_uds"))
		{
		  if(plotResults->isCharm)
		    continue;
		}

	      //	  cout <<"result index: "<< plotResults->resultIndex <<endl;
	      cout << " looking at exp: "<< plotResults->exp <<" itFlav: " << *itFlav <<" isCharm?  " << plotResults->isCharm<<endl;

	      if(binType_zOnly == plotResults->binningType && (*itFlav)==string("_all"))
		{
		  for(int k=0;k<30;k++)
		    {
		      expCounts[plotResults->exp]+=(int)plotResults->kTValues[k];
		    }

		}
	      if(binType_labTheta_z == plotResults->binningType)
		{
		  cout <<"kt bins before : ";
		  for(int i=0;i<10;i++)
		    {
		      cout << pPlotter->plotResults[plotResults->resultIndex].kTValues[i] <<" uncertainties: "<< pPlotter->plotResults[plotResults->resultIndex].kTUncertainties[i]<<" mean: " << pPlotter->plotResults[plotResults->resultIndex].kTMeans[i];
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
		      cout << pPlotter->plotResults[plotResults->resultIndex].kTValues[i] <<" uncertainties "<< pPlotter->plotResults[plotResults->resultIndex].kTUncertainties[i] <<" mean: " << pPlotter->plotResults[plotResults->resultIndex].kTMeans[i];
		    }
		  cout <<endl;
		}
	    }

	  pPlotter->savePlots(plotType_2D);
	  ///get smearing matrix, xini, bini
	  TDirectory* dir=gDirectory;
	  TFile* smearingFile=new TFile("smearing.root");
	  //      for(int c=0;c<MultiPlotter::NumCharges;c++)
	  TCanvas cnvs;
	  //z1_z2 binning (b==0) and onlyZ
	  //      for(int b=0;b<2;b++)
	  for(int b=0;b<1;b++)
	    {

	      for(int c=0;c<2;c++)
		{
		  //just do the pik
		  //	  for(int p=0;p<9;p++)
		  for(int p=1;p<2;p++)
		    //	  for(int p=0;p<MultiPlotter::NumPIDs;p++)
		    {
		      char buffer[500];
		      sprintf(buffer,"kinematicSmearingMatrix_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      TH2D* smearingMatrix=(TH2D*)smearingFile->Get(buffer);
		      if(smearingMatrix->Integral()<10)
			{
			  cout <<" c: "<< c <<" p: "<< p << " integral: " << smearingMatrix->Integral()<<endl;
			  continue;
			}
		      smearingMatrix->Draw("colz");
		      sprintf(buffer,"debug_smM_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      sprintf(buffer,"xini_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      cout <<"trying to load " << buffer <<endl;

		      TH1D* xini=(TH1D*)smearingFile->Get(buffer);
		      xini->Draw();
		      sprintf(buffer,"debug_xini_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      sprintf(buffer,"bini_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      cout <<"trying to load " << buffer <<endl;
		      TH1D* bini=(TH1D*)smearingFile->Get(buffer);
		      bini->Draw();
		      sprintf(buffer,"debug_bini_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      //get combined z/kT histogram for this charge, pid bin
		      TH1D* combinedHisto=pPlotter->getHistogram(b,c,p);
		      combinedHisto->Draw();
		      sprintf(buffer,"debug_combinedH_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      TH1D** d=new (TH1D*);

		      //no unfolding for now
		      //	      TH1D* output=(TH1D*)combinedHisto->Clone("sth");
		      TH1D* output=pPlotter->unfold(smearingMatrix,xini,bini,combinedHisto,d);
		      //for closure test, bini is output....
		      //	      for(int t=0;t<combinedHisto->GetNbinsX();t++)
		      //		{
		      //		  combinedHisto->SetBinContent(t+1,bini->GetBinContent(t+1)+1);
		      //		}
		      //	      TH1D* output=pPlotter->unfold(smearingMatrix,xini,bini,bini,d);
		      output->Draw();
		      sprintf(buffer,"debug_unfoldedH_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);

		      //	      (*d)->Draw();
		      //	      sprintf(buffer,"debug_D_pid%d_charge_%d.png",p,c);
		      //	      cnvs.SaveAs(buffer);
		      //

		      TH1D** sepKtZHistos_mcInput=pPlotter->convertUnfold2Plots(xini,b,c,p,"mcInput");
		      TH1D** sepKtZHistos_mcOutput=pPlotter->convertUnfold2Plots(bini,b,c,p,"mcOut");
		      TH1D** sepKtZHistos=pPlotter->convertUnfold2Plots(output,b,c,p,"dataUnfold");
		      TH1D** sepKtZHistosDataInput=pPlotter->convertUnfold2Plots(combinedHisto,b,c,p,"dataInput");
		      //    TH1D*** sepAllKtZHistos=pPlotter->convertAllUnfold2Plots(output,b,c,p,"allDataInput");
		      TH1D*** sepAllKtZHistos=pPlotter->convertAllUnfold2Plots(combinedHisto,b,c,p,"allDataInput");

		      TCanvas cnvs2;
		      //need one canvas for the legend
		      TLegend leg(0.0,0,1.0,1.0);
		      pair<int,int> zIdx=pPlotter->pidBin2ZBinningIdx(p);
		      int maxZ1=pPlotter->binningZ[zIdx.first].size();
		      int maxZ2=pPlotter->binningZ[zIdx.second].size();
		      int minMaxZ=maxZ1;
		      saveToTxt(b,c,p,maxZ1,maxZ2,pPlotter->getNumKtBins(),sepAllKtZHistos);
		      if(maxZ2<maxZ1)
			minMaxZ=maxZ2;
		      if(minMaxZ>5)
			cnvs2.Divide(3,3);
		      else
			cnvs2.Divide(2,3);
		      cout <<" got maxZ1: "<< maxZ1 << " maxZ2: "<< maxZ2 <<" minMaxZ: "<< minMaxZ <<endl;


		      for(int iZ=0;iZ<minMaxZ;iZ++)
			{
			  cout <<"iZ: "<< iZ <<endl;
			  TVirtualPad* pad=cnvs2.cd(iZ+1);
			  sepKtZHistos_mcInput[iZ]->SetMarkerColor(kRed);
			  sepKtZHistos_mcInput[iZ]->SetMarkerStyle(21);
			  sepKtZHistos_mcOutput[iZ]->SetMarkerColor(kGreen);
			  sepKtZHistos_mcOutput[iZ]->SetMarkerStyle(12);
			  sepKtZHistosDataInput[iZ]->SetMarkerStyle(34);
			  sepKtZHistosDataInput[iZ]->SetMarkerColor(kBlack);
			  sepKtZHistos[iZ]->SetMarkerColor(kBlue);
			  sepKtZHistos[iZ]->SetMarkerStyle(23);
			  //		  		  sepKtZHistos[iZ]->Draw("A P  E1");
			  //		  		  sepKtZHistos_mcInput[iZ]->Draw("SAME P E1");
			  sepKtZHistos[iZ]->Draw("SAME P  E1");
			  sepKtZHistos_mcInput[iZ]->Draw("A P E1");
			  sepKtZHistos_mcOutput[iZ]->Draw("SAME P E1");		  

			  sepKtZHistosDataInput[iZ]->Draw("SAME P  E1");
			  pad->Update();
			  Double_t x1,x2,y1,y2;
			  pad->GetRangeAxis(x1,y1,x2,y2);
			  pad->DrawFrame(x1,y1,x2,y2*1.3);
			  sepKtZHistos_mcInput[iZ]->Draw("SAME P E1");
			  sepKtZHistos_mcOutput[iZ]->Draw("SAME P E1");		  
			  sepKtZHistos[iZ]->Draw("SAME P  E1");
			  sepKtZHistosDataInput[iZ]->Draw("SAME P  E1");
			  if(iZ==0)
			    {
			      leg.AddEntry(sepKtZHistos_mcInput[iZ],"MC Input","lep");
			      leg.AddEntry(sepKtZHistos_mcOutput[iZ],"MC Output","lep");
			      leg.AddEntry(sepKtZHistosDataInput[iZ],"Data Input","lep");
			      leg.AddEntry(sepKtZHistos[iZ],"MC Unfolded","lep");
			    }
			}
		      cnvs2.cd(minMaxZ+1);
		      leg.Draw();
		      sprintf(buffer,"unfoldedResult_binning_%d_pid_%d_charge_%d.png",b,p,c);
		      cnvs2.SaveAs(buffer);

		      //	      sprintf("");

		    }
		}
	    }
	  dir->cd();
	  smearingFile->Close();
	  ///save the returned plots... put them on the same plot etc..

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

	  for(map<int,int>::iterator itEx=expCounts.begin();itEx!=expCounts.end();itEx++)
	    {
	      cout <<"exp " << itEx->first <<" : " << itEx->second <<endl;
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
    }
  
};
