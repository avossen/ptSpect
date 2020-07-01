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
  if(argc<2)
    {
      cout <<"usage: CombPlots rootPath smearingFileName rawSmearingFileName"<<endl;
      exit(0);
    }
  //  bool closureTest=true;
  bool closureTest=false;
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
  int minExp=0;
  //int minExp2=0;
  //  int maxExp=2;
  //  int maxExp2=200;

  int maxExp=135;
  int minExp2=0;
  int maxExp2=169;
  // int minExp=0;

  char* rootPath=argv[1];
  char smearingFileName[400];
  char smearingFileRawName[400];
  if(argc>2)
    {
      sprintf(smearingFileName,"%s",argv[2]);
    }
  else{
    sprintf(smearingFileName,"%s","smearing.root");
      }
  if(argc>3)
    {
      sprintf(smearingFileRawName,argv[3]);
    }
  else
    {
      sprintf(smearingFileRawName,"%s","smearing.root");
    }

  char dataMcNameAdd[100];
  sprintf(dataMcNameAdd,"");
  if(argc>4)
    sprintf(dataMcNameAdd,"%s",argv[4]);

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
  string woPlotterName("NormalWoA");

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
	  MultiPlotter* pPlotter=new MultiPlotter(m_useQt,const_cast<char*>("."),fullName.c_str(),string(""),0,false,false,false,false,fileTypeEnd);

	  cout <<" setting plotter name to:" << fullName <<endl;
	  pPlotter->setName(fullName);
	  vPlotters.push_back(pPlotter);
	  //has to be 0!!
	  PlotResults* plotResults=0;

	  chAll->SetBranchAddress("PlotBranch",&plotResults);
	  for(int binningType=binType_zOnly; binningType<binType_end;binningType++)
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
	  //this is for the run stability plots
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
	

	  ///////////up to hear semi-useful preparation-------


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

	      //check if the  result we are reading right now is compatible with the 
	      //flavor we are looking for
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
	      if(binType_z_z == plotResults->binningType)
		{
		  cout <<"kt bins before : ";
		  for(int i=0;i<10;i++)
		    {
		      		      cout << pPlotter->plotResults[plotResults->resultIndex].kTValues[i] <<" uncertainties: "<< pPlotter->plotResults[plotResults->resultIndex].kTUncertainties[i]<<" mean: " << pPlotter->plotResults[plotResults->resultIndex].kTMeans[i];
		      cout <<"  rhs mean: " << plotResults->kTMeans[i] <<" ";
		    }
		}

	      //
	      //    above is hardly ever used
	      ///
	      //	  cout <<endl;
	      //I guess this is where the magic happens
	      ////
	      ///------>>>MAGIC<<<<----------
	      ///
	      pPlotter->plotResults[plotResults->resultIndex]+=(*plotResults);
	      ////
	      ///
	      ///--------->>>> END MAGIC <<<<<------------
	      //
	      if(binType_z_z == plotResults->binningType)
		{
		  //	      cout <<"and then...  ";
		  for(int i=0;i<10;i++)
		    {
		      cout << pPlotter->plotResults[plotResults->resultIndex].kTValues[i] <<" uncertainties "<< pPlotter->plotResults[plotResults->resultIndex].kTUncertainties[i] <<" mean: " << pPlotter->plotResults[plotResults->resultIndex].kTMeans[i];
		    }
		  cout <<endl;
		}
	    }

	  //////////////now we have all the results in our plotter, time to do something with it


	
	  cout <<"save plot " <<endl;
	  pPlotter->savePlots(plotType_2D);
	  cout <<"print debug!" <<endl;
	  pPlotter->printDebug(plotType_2D);
	  cout <<"done " <<endl;
	  ///get smearing matrix, xini, bini
	  TDirectory* dir=gDirectory;

	  //this could be created by just adding all smering files
	  TFile* myOutputFile=new TFile("unfoldingOut.root","RECREATE");
	  TFile* smearingFile=new TFile(smearingFileName);
	  TFile* smearingFileRaw=new TFile(smearingFileRawName);
	  cout <<"loaded smearingFile: "<< smearingFileName <<" and raw file: "<< smearingFileRawName <<endl;
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
		  //		  for(int p=1;p<2;p++)
		  for(int p=0;p<9;p++)
		    //	  for(int p=0;p<MultiPlotter::NumPIDs;p++)
		    {
		      char buffer[500];
		      char buffer2[500];
		      sprintf(buffer,"kinematicSmearingMatrix_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      TH2D* smearingMatrix=(TH2D*)smearingFile->Get(buffer);

		      //check for max and average
		      double max=-1;
		      double avg=0;
		      
		      int smDimX=smearingMatrix->GetNbinsX();
		      int smDimY=smearingMatrix->GetNbinsY();

		      for(int ix=0;ix<smearingMatrix->GetNbinsX();ix++)
			{
			  for(int iy=0;iy<smearingMatrix->GetNbinsY();iy++)
			    {
			      double cont=smearingMatrix->GetBinContent(ix+1,iy+1);
			      if(cont>max)
				max=cont;
			      avg+=cont;
			      if(ix+1==iy+1)
				{
				  //			      smearingMatrix->SetBinContent(ix+1,iy+1,1.0);
				}
			      else
				{
				  //			      smearingMatrix->SetBinContent(ix+1,iy+1,0.0);
				}
			      
			    }
			}

		      cout <<"sm max: "<< max<<endl;
		      for(int ix=0;ix<smearingMatrix->GetNbinsX();ix++)
			{
			  for(int iy=0;iy<smearingMatrix->GetNbinsY();iy++)
			    {
			      double v=smearingMatrix->GetBinContent(ix+1,iy+1);
			      //			      smearingMatrix->SetBinContent(ix+1,iy+1,v/max);
			    }
			}




		      cout <<"smearing matrix average: " << avg/((double)smDimX*smDimY)<<" max: "<< max <<" avg from int: "<< smearingMatrix->Integral()/((double)smDimX*smDimY)<<endl;
		      if(smearingMatrix->Integral()<10)
			{
			  cout <<" c: "<< c <<" p: "<< p << " integral: " << smearingMatrix->Integral()<<endl;
			  continue;
			}

		      ////////nothing has happened so far other than us loading the smearing matrix


		      cnvs.SetLogz(true);
		      smearingMatrix->Draw("colz");
		      sprintf(buffer,"debug_smM_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      cnvs.SetLogz(false);
		      sprintf(buffer,"xini_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      cout <<"trying to load " << buffer <<endl;



		      //ormalize?


		      cnvs.SetLogy();
		      TH1D* xini=(TH1D*)smearingFileRaw->Get(buffer);
		      xini->Draw();
		      sprintf(buffer,"debug_xini_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      sprintf(buffer,"debug_xini_binning%d_pid%d_charge_%d.pdf",b,p,c);
		      cnvs.SaveAs(buffer);

		      sprintf(buffer,"bini_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      cout <<"trying to load " << buffer <<endl;
		      TH1D* bini=(TH1D*)smearingFile->Get(buffer);
		      bini->Draw();
		      sprintf(buffer,"debug_bini_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      sprintf(buffer,"debug_bini_binning%d_pid%d_charge_%d.pdf",b,p,c);
		      cnvs.SaveAs(buffer);


		      sprintf(buffer,"backgroundCounts_binning%d_pidBin%d_chargeBin%d",b,p,c);
		      TH1D* bgCounts=(TH1D*)smearingFile->Get(buffer);
		      sprintf(buffer,"debug_backgroundCounts_binning%d_pidBin%d_chargeBin%d.png",b,p,c);
		      bgCounts->Draw();
		      sprintf(buffer,"debug_backgroundCounts_binning%d_pidBin%d_chargeBin%d.pdf",b,p,c);
		      bgCounts->Draw();
		      cnvs.SaveAs(buffer);
		      //subtract the backgrund
		      //		      cout <<"subtracting background" <<endl;--->now in xini
		      ////		      bini->Add(bgCounts,-1);

		      cout <<"getting combined histo for b: "<< b <<" c: "<< c <<" p: "<< p <<endl;
		      //get combined z/kT histogram for this charge, pid bin
		      TH1D* combinedHisto=pPlotter->getHistogram(b,c,p,"");
		      TH1D* combinedHistoSys=pPlotter->getHistogram(b,c,p,"sys",1);
		      TH1D* combinedHistoUpperSys=pPlotter->getHistogram(b,c,p,"sysUp");
		      TH1D* combinedHistoLowerSys=pPlotter->getHistogram(b,c,p,"sysDown");


		      combinedHisto->Draw();
		      sprintf(buffer,"debug_combinedH_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      sprintf(buffer,"debug_combinedH_binning%d_pid%d_charge_%d.pdf",b,p,c);
		      cnvs.SaveAs(buffer);
		      TH1D** d=new (TH1D*);
		      TH2D** statCov=new TH2D*;
		      TH2D** sysCov=new TH2D*;
		      TH2D** mcStatCov=new TH2D*;

		      max=-1;
		      avg=0;
		      for(int ix=0;ix<bini->GetNbinsX();ix++)
			{
			  double cont=bini->GetBinContent(ix+1);
			  if(cont>max)
			    max=cont;
			  avg+=cont;
			}
		      cout <<"bini max: "<< max<<" avg: " << avg/((double)bini->GetNbinsX()) <<" count: "<< avg <<" integral: " << bini->Integral()<<endl;

		      max=-1;
		      avg=0;
		      for(int ix=0;ix<xini->GetNbinsX();ix++)
			{
			  double cont=xini->GetBinContent(ix+1);
			  if(cont>max)
			    max=cont;
			  avg+=cont;
			}
		      cout <<"xini max: "<< max<<" avg: " << avg/((double)xini->GetNbinsX()) <<" count: "<< avg <<" integral: " << xini->Integral()<<endl;

		      //no unfolding for now
		      //	      TH1D* output=(TH1D*)combinedHisto->Clone("sth");
		      TH1D* output;
		      TH1D* outputHighSys;
		      TH1D* outputLowSys;
		      cout <<"writing pidCorrected" <<endl;

		      sprintf(buffer,"xini_binning%d_pid%d_charge_%d.txt",b,p,c);
		      MultiPlotter::printMatrix(xini,buffer);
		      sprintf(buffer,"bini_binning%d_pid%d_charge_%d.txt",b,p,c);
		      MultiPlotter::printMatrix(bini,buffer);
		      sprintf(buffer,"pidCorrected_binning%d_pid%d_charge_%d.txt",b,p,c);
		      MultiPlotter::printMatrix(combinedHisto,buffer);
		      sprintf(buffer,"pidCorrectedUncert_binning%d_pid%d_charge_%d.txt",b,p,c);
		      MultiPlotter::printMatrix(combinedHisto,buffer,true);
		      sprintf(buffer,"pidSys_binning%d_pid%d_charge_%d.txt",b,p,c);
		      MultiPlotter::printMatrix(combinedHistoSys,buffer);


		      sprintf(buffer,"_binning%d_pid%d_charge_%d.png",b,p,c);
		      if(!closureTest)
			{
			  cout <<"unfolding w/o closure test .." <<endl;
			  output=pPlotter->unfold(smearingMatrix,xini,bini,combinedHisto,combinedHistoSys,d,statCov,mcStatCov,sysCov,buffer);
			  cout <<"done .. got  "<< output->GetNbinsX() <<" dim matrix " <<endl;
			  cout <<"printing unfolding output " <<endl;
			  sprintf(buffer,"unfolded_binning%d_pid%d_charge_%d.txt",b,p,c);
			  MultiPlotter::printMatrix(output,buffer);
			  sprintf(buffer,"statCov_binning%d_pid%d_charge_%d.txt",b,p,c);
			  MultiPlotter::printMatrix(*statCov,buffer);
			  sprintf(buffer,"mcStatCov_binning%d_pid%d_charge_%d.txt",b,p,c);
			  MultiPlotter::printMatrix(*mcStatCov,buffer);
			  sprintf(buffer,"sysCov_binning%d_pid%d_charge_%d.txt",b,p,c);
			  MultiPlotter::printMatrix(*sysCov,buffer);
			  //			  outputHighSys=pPlotter->unfold(smearingMatrix,xini,bini,combinedHistoUpperSys,combinedHistoSys,d,statCov,mcStatCov,sysCov,buffer);
			  //			  outputLowSys=pPlotter->unfold(smearingMatrix,xini,bini,combinedHistoLowerSys,combinedHistoSys,d,statCov,mcStatCov,sysCov,buffer);
		      //for closure test, bini is output....
			}
		      else
			{
			  for(int t=0;t<combinedHisto->GetNbinsX();t++)
			    {
			      combinedHisto->SetBinContent(t+1,bini->GetBinContent(t+1));
			    }
			  cout <<"unfolding .." <<endl;
			   output=pPlotter->unfold(smearingMatrix,xini,bini,combinedHisto,combinedHistoSys,d,statCov,mcStatCov,sysCov,buffer);
			   cout <<"done .. got  "<< output->GetNbinsX() <<" dim matrix " <<endl;
			}
		      ///--->save debugs
		      ///-->just for tmp
		       //		      output=combinedHisto;
		      output->Draw();
		      sprintf(buffer,"debug_unfoldedH_binning%d_pid%d_charge_%d.png",b,p,c);
		      cnvs.SaveAs(buffer);
		      sprintf(buffer,"debug_unfoldedH_binning%d_pid%d_charge_%d.pdf",b,p,c);
		      cnvs.SaveAs(buffer);
		      ////---->now we should save the unfolded result into the plotResults again and save the whole plotter
		      pPlotter->setHistogram(b,c,p,combinedHisto,combinedHistoUpperSys,combinedHistoLowerSys);

		      /////----->
		      //	      (*d)->Draw();
		      //	      sprinntf(buffer,"debug_D_pid%d_charge_%d.png",p,c);
		      //	      cnvs.SaveAs(buffer);
		      //

		      //this is only one array with the diagonal.....

		      TH1D** sepKtZHistos_mcInput=pPlotter->convertUnfold2Plots(xini,b,c,p,"mcInput");
		      cout <<" got sepkt 1" <<endl;
		      TH1D** sepKtZHistos_mcOutput=pPlotter->convertUnfold2Plots(bini,b,c,p,"mcOut");
		      cout <<" got sepkt 2" <<endl;
		      TH1D** sepKtZHistos=pPlotter->convertUnfold2Plots(output,b,c,p,"dataUnfold");
		      cout <<" got sepkt 3" <<endl;
		      TH1D** sepKtZHistosDataInput=pPlotter->convertUnfold2Plots(combinedHisto,b,c,p,"dataInput");
		      //    TH1D*** sepAllKtZHistos=pPlotter->convertAllUnfold2Plots(output,b,c,p,"allDataInput");

		      cout <<" got sepkt " <<endl;
		      //these are all z bins...
		      myOutputFile->cd();
		      TH1D*** sepAllKtZHistos=pPlotter->convertAllUnfold2Plots(combinedHisto,b,c,p,"dataInput");
		      TH1D*** sepAllKtZHistosUnfolded=pPlotter->convertAllUnfold2Plots(output,b,c,p,"unfolded");
		      cout <<" converted" <<endl;
		      TCanvas cnvs2;
		      //need one canvas for the legend
		      TLegend leg(0.0,0,1.0,1.0);
		      pair<int,int> zIdx=pPlotter->pidBin2ZBinningIdx(p);



		      int maxZ1=pPlotter->binningZ[zIdx.first].size();
		      int maxZ2=pPlotter->binningZ[zIdx.second].size();
		      cout <<" directory " <<endl;

		      //zs can be different due to PID (kaon, pion, proton have different numbers of z bins);
		      TDirectory* tempDir=gDirectory;

		      for(int iz1=0;iz1<maxZ1;iz1++)
			{
			  for(int iz2=0;iz2<maxZ2;iz2++)
			    {
			      sepAllKtZHistos[iz1][iz2]->Write();
			      sepAllKtZHistosUnfolded[iz1][iz2]->Write();
			    }
			}
		      tempDir->cd();
		      cout<<" before out to txt" <<endl;

		      int minMaxZ=maxZ1;
		      saveToTxt(b,c,p,maxZ1,maxZ2,pPlotter->getNumKtBins(),sepAllKtZHistos,"input");
		      saveToTxt(b,c,p,maxZ1,maxZ2,pPlotter->getNumKtBins(),sepAllKtZHistosUnfolded,"output");
		      cout <<"saved to txt" <<endl;
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
		      sprintf(buffer,"unfoldedResult_binning_%d_pid_%d_charge_%d.pdf",b,p,c);
		      cnvs2.SaveAs(buffer);
		    }
		}
	    }
	  dir->cd();
	  smearingFile->Close();
	  myOutputFile->Write();
	  myOutputFile->Close();
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
