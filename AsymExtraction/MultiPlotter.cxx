#include "MultiPlotter.h"
#include "PlotResults.h"


//compute means and fill plotResults
void MultiPlotter::doPlots()
{

  for(int bt=binType_labTheta_z; bt<binType_end;bt++)
    {
      //      cout <<"looking at bin type " << bt <<endl;
      //      for(int chargeBin=0;chargeBin<NumCharges;chargeBin++)
      for(int chargeBin=0;chargeBin<2;chargeBin++)
	{
	  for(int pidBin=0;pidBin<NumPIDs;pidBin++)
	    {
	      for(int firstBin=0;firstBin<maxKinMap[pidBin][bt].first;firstBin++)
		{
		  for(int secondBin=0;secondBin<maxKinMap[pidBin][bt].second;secondBin++)
		    {
		      double locCount=0;
		      int resIdx=getResIdx(bt,pidBin,chargeBin,firstBin,secondBin);
		      for(unsigned int ktBin=0;ktBin<numKtBins;ktBin++)
			{
			  locCount+=counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  plotResults[resIdx].kTValues[ktBin]=counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  plotResults[resIdx].kTMeans[ktBin]=meanValues_kT[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]/counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  if(bt==binType_labTheta_z && chargeBin==0)
			    {
			      cout <<"kTbin: " <<ktBin << " kt count: "<< plotResults[resIdx].kTValues[ktBin] << " mean kt: "<< plotResults[resIdx].kTMeans[ktBin]<<endl;
			    }
			}
		      plotResults[resIdx].meanKinBin1=meanValues_kin1[bt][pidBin][chargeBin][firstBin][secondBin]/locCount;
		      plotResults[resIdx].meanKinBin1=meanValues_kin2[bt][pidBin][chargeBin][firstBin][secondBin]/locCount;
		      plotResults[resIdx].firstKinBin=firstBin;
		      plotResults[resIdx].secondKinBin=secondBin;
		      if(bt==binType_labTheta_z && chargeBin==0)
			{
			  cout <<"resIdx: " << resIdx <<" mean 1 : " << plotResults[resIdx].meanKinBin1;
			  cout <<",  mean2: " << plotResults[resIdx].meanKinBin2 <<",  firstBin: " << firstBin <<" second: " << secondBin;
			  
			}

		      plotResults[resIdx].pidBin=pidBin;
		      plotResults[resIdx].chargeBin=chargeBin;
		      plotResults[resIdx].binningType=bt;
		      plotResults[resIdx].exp=m_expNr;
		      plotResults[resIdx].on_res=m_onRes;
		      plotResults[resIdx].isUds=m_uds;
		      plotResults[resIdx].isCharm=m_charm;
		      plotResults[resIdx].isMC=m_mc;
		      plotResults[resIdx].resultIndex=resIdx;		  
		    }
		  
		}
	    }
	  
	}
    }
}

void MultiPlotter::loadBinnings()
{

  binningZ=new vector<float>[6];

  //pions
  binningZ[0].push_back(0.1);
  binningZ[0].push_back(0.15);
  binningZ[0].push_back(0.2);
  binningZ[0].push_back(0.3);
  ////  binningZ.push_back(0.7);
  binningZ[0].push_back(0.4); //// 
  //  binningZ.push_back(0.5);
  //  binningZ.push_back(0.6);
  binningZ[0].push_back(2.6);

  //the rest
  for(int i=1;i<6;i++)
    {
      binningZ[i].push_back(0.2);
      binningZ[i].push_back(0.3);
      binningZ[i].push_back(0.4);
      ////  binningZ.push_back(0.7);
      binningZ[i].push_back(0.5); //// 
  //  binningZ.push_back(0.5);
  //  binningZ.push_back(0.6);
      binningZ[i].push_back(2.6);


    }


  binningKt.push_back(0.15);
  binningKt.push_back(0.3);
  binningKt.push_back(0.45);
   binningKt.push_back(0.6);
  binningKt.push_back(0.75);
    binningKt.push_back(0.9);
  binningKt.push_back(1.05);
  //  binningKt.push_back(1.5);
  binningKt.push_back(2.0);
  binningKt.push_back(10000);

  //    binningLabTheta.push_back(0.9);
  //  binningLabTheta.push_back(1.1);

  ////  binningLabTheta.push_back(0.5);
  ////  binningLabTheta.push_back(0.9);
  ////  binningLabTheta.push_back(1.4);
  ////  binningLabTheta.push_back(1.8);
  
  ////  binningLabTheta.push_back(2.5);
  binningLabTheta.push_back(8.0);

  binningThrustLabTheta.push_back(0.5);
  binningThrustLabTheta.push_back(0.9);
  binningThrustLabTheta.push_back(1.4);
  binningThrustLabTheta.push_back(1.8);

  //  binningThrustLabTheta.push_back(2.5);
  binningThrustLabTheta.push_back(8.0);




  binningQt.push_back(0.5);
  binningQt.push_back(1);
  binningQt.push_back(2);
  binningQt.push_back(3);
  binningQt.push_back(4);
  binningQt.push_back(5);
  binningQt.push_back(100);

  //is from 0 to pi
  binningCmsThrustTheta.push_back(1.3);
  binningCmsThrustTheta.push_back(1.7);
  binningCmsThrustTheta.push_back(5.0);

  binningThrust.push_back(0.85);
  binningThrust.push_back(0.9);
  binningThrust.push_back(1.85);



};


void MultiPlotter::saveSmearingMatrix()
{
  rFile.cd();
  //numCharges is 3, but only the first two should be populated (likesign, unlikesign)
  for(int b=0;b<2;b++)
    {
  for(int c=0;c<NumCharges;c++)
    {
      for(int p=0;p<NumPIDs;p++)
	{
	  kinematicSmearingMatrix[b][p][c]->Write();
      
	  xini[b][p][c]->Write();
	  bini[b][p][c]->Write();
	}
    }
    }
}


//convert the convuoluted histogram we get from unfolding into 'regular' plots
//also put in the binwidth factors that were not used for the unfolding
//for the z_z binning, just put out the parallel bins for now...
TH1D** MultiPlotter::convertUnfold2Plots(TH1D* input, int binning,  int chargeBin, int pidBin, const char* nameAdd)
{
  char buffer[300];

  pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);
  int maxZBin=(binningZ[zIdx.first].size()>binningZ[zIdx.second].size()) ? binningZ[zIdx.first].size() : binningZ[zIdx.second].size();
  int minMaxZBin=(binningZ[zIdx.first].size()<binningZ[zIdx.second].size()) ? binningZ[zIdx.first].size() : binningZ[zIdx.second].size();

  //for z1 binning, just take number of z1 bins. For z1/z2 where we do the diagonal, take smaller of the two
  int locMaxZBin=binningZ[zIdx.first].size();
  if(binning==0)
    {
      locMaxZBin=minMaxZBin;
    }

  TH1D** ret=new TH1D*[locMaxZBin];
  //this should be the smaller of the two since we are doing the diagonal bins
  for(int zBin=0;zBin<locMaxZBin;zBin++)
    {
      sprintf(buffer,"un_convert_binning_%d_cBin_%d_pBin_%d_zBin_%d_%s",binning,chargeBin,pidBin,zBin,nameAdd);
      ret[zBin]=new TH1D(buffer,buffer,numKtBins,0,numKtBins);
    }

  Double_t value;
  //      int recBin=z1Bin1*numKtBins+kTBin1;
  for(int zBin=0;zBin<locMaxZBin;zBin++)
    {
      for(int kTBin=0;kTBin<binningKt.size();kTBin++)
	{
	  float binWidthFactor=1.0;
	  if(0==kTBin)
	    {
	      binWidthFactor=binningKt[0];
	    }
	  else
	    {
	      binWidthFactor=binningKt[kTBin]-binningKt[kTBin-1];
	    }
	  //for the last bin, it doesn't make sense to divide by 1000 or so...
	  binWidthFactor > 1.0 ?  (binWidthFactor=1.0) : true ;
	  binWidthFactor<=0 ?   (binWidthFactor=1.0) : true;
	  //the first one is z2, so the array size is multiplied with the max z1 bns
	  ///see : 
	  int combBin=zBin*binningZ[zIdx.first].size()*numKtBins + zBin*numKtBins + kTBin;
	  if(binning==1)//onlyZ
	    combBin=zBin*numKtBins+kTBin;
	  //	  cout<<endl <<"combBin : " << combBin<<endl;
	  if(combBin>=input->GetNbinsX())
	    {
	      cout <<"convert, zBin " << combBin <<" greater than " << input->GetNbinsX()<<endl;
	    }
	  //	  cout <<"convert unfold.. binWidthFactor: " << binWidthFactor<<endl;
	  value=input->GetBinContent(combBin+1);
	  //	  cout <<"value first " << value <<endl;
	  value=input->GetBinContent(combBin+1)/binWidthFactor;
	  //	  cout <<"now: "<< value<<endl;
	  ret[zBin]->SetBinContent(kTBin+1,value);
	}
    }
  return ret;
}

//mc_input is what is called xini in the tsvdunfold docu, MC_out is what is called bini
TH1D* MultiPlotter::unfold(TH2D* smearingMatrix, TH1D* MC_input,TH1D* MC_out, TH1D* data, TH1D** d)
{
  //  TH1D* xini=(TH1D*)file2.Get("xini");
  //  TH1D* bini=(TH1D*)file2.Get("bini");
  
  //  TH1D* d=(TH1D*)file.Get("bini");
  //  TH1D* k=(TH1D*)bini->Clone();
  
  //  for(Int_t i=0;i<1000;i++)
  //    {
  //      d->Fill(i%105);
  //
  //    }
  char buffer[300];
  sprintf(buffer,"statcof_%s",data->GetName());
  TH2D* statcovMatrix=new TH2D(buffer,buffer,data->GetNbinsX(),0,data->GetNbinsX(),data->GetNbinsX(),0,data->GetNbinsX());
  cout <<" input has " << MC_input->GetNbinsX();
  cout <<" bins and output " << MC_out->GetNbinsX() <<" data: " << data->GetNbinsX() <<" smearing matrix " << smearingMatrix->GetNbinsX() <<" x " << smearingMatrix->GetNbinsY() <<endl;
  for(int i =0;i<MC_input->GetNbinsX();i++)
    {
      cout <<MC_input->GetBinContent(i+1) <<" ";
    }
  cout <<endl;
  for(int i =0;i<MC_out->GetNbinsX();i++)
    {
      cout <<MC_out->GetBinContent(i+1) <<" ";

    }
  cout <<endl<<endl;
  for(int i =0;i<data->GetNbinsX();i++)
    {
      cout <<data->GetBinContent(i+1) <<" error: " << data->GetBinError(i+1)<<" ";
  for(int j =0;j<data->GetNbinsX();j++)
    {
      statcovMatrix->SetBinContent(i+1,j+1,data->GetBinError(i+1)*data->GetBinError(j+1));
    }
 
    }
  cout <<endl<<endl;
  for(int i =0;i<smearingMatrix->GetNbinsX();i++)
    {

      for(int j =0;j<smearingMatrix->GetNbinsY();j++)
	{
	  cout <<smearingMatrix->GetBinContent(i+1,j+1) <<" " ;	  
	}
      cout <<endl;
    }
  cout <<endl;

  int rank=MC_input->GetNbinsX();
  TSVDUnfold* f=new TSVDUnfold(data,statcovMatrix,MC_out,MC_input,smearingMatrix);
  f->SetNormalize(false);

  TH1D* ret=f->Unfold(rank);
  TH2D* uadetcov=f->GetAdetCovMatrix(100);
  TH2D* utaucov= f->GetXtau();
  utaucov->Add(uadetcov);
  sprintf(buffer,"unfold_from_%s",data->GetName());
  cout <<" unfolding " << buffer <<endl;
  TH1D* ret2=(TH1D*)ret->Clone(buffer);
  for(int i=0;i<ret2->GetNbinsX();i++)
    {
      ret2->SetBinError(i,TMath::Sqrt(utaucov->GetBinContent(i,i)));
    }
  (*d)=f->GetD();
  for(int i=0;i<ret->GetNbinsX();i++)
    {
      //      cout <<"data: "<< data->GetBinContent(i+1) <<" unfolded: " << ret->GetBinContent(i+1) <<" mc out: "<< MC_out->GetBinContent(i+1) <<" mc in: "<< MC_input->GetBinContent(i+1)<<endl;

    }



  sprintf(buffer,"debug_D_from%s.png",MC_input->GetName());
  cout <<"d: bins: "<< (*d)->GetNbinsX()<<endl;
  TCanvas c;
  (*d)->Draw();
  c.SetLogy();
  c.SaveAs(buffer);
  return ret2;
}

//void MultiPlotter::getIntAsymmetry(float a[3], float ea[3],int binningType,int chargeType, bool save1D)

//for now only for the zOnly binning--> changed to z1,z2 and zOnly (binning argument)
//don't do the bin width normalization since we also don't do it for the xini, bini
TH1D* MultiPlotter::getHistogram(int binning, int chargeType, int pidType)
{
  PlotResults* m_plotResults=plotResults;
  PlotResults* loc_plotResults=0;
  char buffer[200];
  char buffer1[200];
  float mX[50];
  float mY[50];
  float mXErr[50];
  float mYErr[50];

  int binningType=binType_z_z;
  if(1==binning)
    binningType=binType_zOnly;

  string binName=getBinName(binningType,pidType,chargeType,-1,-1);
  sprintf(buffer,"%s",binName.c_str());
  sprintf(buffer1,"histo_%s",buffer);
  TH1D* ret=new TH1D(buffer1,buffer1,maxSmearing[pidType][binning],0,maxSmearing[pidType][binning]);
  //  cout <<" getH: " << maxSmearing <<endl;
    cout <<"saving graph for " << binName <<" buffer; " << buffer<<endl;
  for(int i=0;i<maxKinMap[pidType][binningType].first;i++)
    {
      //      for(int j=0;j<maxKinMap[binningType].second;j++)
      for(int j=0;j<maxKinMap[pidType][binningType].second;j++)
	{
	  double normFactor=1.0;
	  double maxVal=-1.0;
	  double maxValNorm=-1.0;
	  int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
	  for(unsigned int iKtBin=0;iKtBin<numKtBins;iKtBin++)
	    {
	      float binWidthFactor=1.0;
	      if(0==iKtBin)
		{
		  binWidthFactor=binningKt[0];
		}
	      else
		{
		  binWidthFactor=binningKt[iKtBin]-binningKt[iKtBin-1];
		}
	      /////	      binWidthFactor > 1.0 ?  (binWidthFactor=1.0) : true ;
	      ////	      binWidthFactor<=0 ?   (binWidthFactor=1.0) : true;

	      binWidthFactor=1.0;

	      //normalize so that the final points have the same maximum
	      if(maxValNorm< m_plotResults[resIdx].kTValues[iKtBin]/binWidthFactor)
		{
		  maxVal=m_plotResults[resIdx].kTValues[iKtBin];
		  maxValNorm=m_plotResults[resIdx].kTValues[iKtBin]/binWidthFactor;
		}

	    }
	  loc_plotResults=&m_plotResults[resIdx];
	  //	  normFactor=1.0/maxValNorm;
	  normFactor=1.0;

	  for(unsigned int iKtBin=0;iKtBin<numKtBins;iKtBin++)
	    {
	      float binWidthFactor=1.0;
	      if(0==iKtBin)
		{
		  binWidthFactor=binningKt[0];
		}
	      else
		{
		  binWidthFactor=binningKt[iKtBin]-binningKt[iKtBin-1];
		}
	      //for the last bin, it doesn't make sense to divide by 1000 or so...
	      ////	      binWidthFactor > 1.0 ?  (binWidthFactor=1.0) : true ;
	      ////	      binWidthFactor<=0 ?   (binWidthFactor=1.0) : true;

	      binWidthFactor=1.0;
	      //this is a bit confusing: becuase the z_z binning saves as pair z1, z2, but in this case j is the index of z1, we have to take this into account
	      int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
	      //	      cout <<"looking at index:" << resIdx<<endl;
	      mX[iKtBin]=m_plotResults[resIdx].kTMeans[iKtBin];
	      //	      cout <<"mX["<<iKtBin <<"] " << mX[iKtBin]<<endl;
	      if((iKtBin>0)&& mX[iKtBin]<=mX[iKtBin-1]) 
		{
		  cout <<"MultiPlotter::getHistogram wanting to set X["<<iKtBin<<"] to: " << mX[iKtBin] <<" but the one before is: " << mX[iKtBin-1] <<endl;
		  mX[iKtBin]=mX[iKtBin-1]+0.1;
		}
	      //	  cout <<"setting x: " << mX[iKtBin] <<endl;
	      mXErr[iKtBin]=0.0;
	      mY[iKtBin]=m_plotResults[resIdx].kTValues[iKtBin]*normFactor/binWidthFactor;
	      Double_t binContent=m_plotResults[resIdx].kTValues[iKtBin]*normFactor/binWidthFactor;
	      //have to add one due to histos counting from 1
	      if(binning==0)//for the z_z binning, the z1, z2 used for getResIdx are switchted...
		ret->SetBinContent(j*maxKinMap[pidType][binningType].first*numKtBins+i*numKtBins+iKtBin+1,binContent);
	      else
		ret->SetBinContent(j*numKtBins+iKtBin+1,binContent);
	      //	      cout <<" getH, bin number: j: "<< j << " iKtBin  "<<iKtBin <<" bin number: " << j*numKtBins+iKtBin << " content: " << binContent <<endl;
	      mYErr[iKtBin]=sqrt(m_plotResults[resIdx].kTValues[iKtBin])*normFactor/binWidthFactor;
	      //	      cout <<"mY["<<iKtBin <<"] " << mY[iKtBin]<<endl;
	      //	      cout <<"mYErr["<<iKtBin <<"] " << mYErr[iKtBin]<<endl;
	      //	  cout <<"y: " << mY[iKtBin] << ", " << mYErr[iKtBin] <<endl;
	    }
	}
    }
  return ret;
}



void MultiPlotter::savePlots( plotType mPlotType)
{

  PlotResults* m_plotResults=plotResults;
  PlotResults* loc_plotResults=0;

  rFile.cd();
  TTree *tree = new TTree("PlotTree","PlotTree");
  tree->Branch("PlotBranch","PlotResults",&loc_plotResults,32000,99);

  //  int numKinBin1=0;
  //  int numKinBin2=0;
  char buffer[200];
  char buffer1[200];
  float mX[50];
  float mY[50];
  float mXErr[50];
  float mYErr[50];
  for(int binningType=binType_labTheta_z; binningType<binType_end;binningType++)
    {
      for(int pidType=0;pidType<9;pidType++)
	{
	  for(int chargeType=0;chargeType<2;chargeType++)
	    {
	      string binName=getBinName(binningType,pidType,chargeType,-1,-1);
	      sprintf(buffer,"%s",binName.c_str());
	      cout <<"saving graph for " << binName <<" buffer; " << buffer<<endl;
	      for(int i=0;i<maxKinMap[pidType][binningType].first;i++)
		{
	      for(int j=0;j<maxKinMap[pidType][binningType].second;j++)
		{
		  //		  cout <<" bin: " << i << ", " << j << endl;
		  
		  
		  double normFactor=1.0;
		  double maxVal=-1.0;
		  double maxValNorm=-1.0;
		  int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
		  for(unsigned int iKtBin=0;iKtBin<numKtBins;iKtBin++)
		    {
		      float binWidthFactor=1.0;
		      if(0==iKtBin)
			{
			  binWidthFactor=binningKt[0];
			}
		      else
			{
			  binWidthFactor=binningKt[iKtBin]-binningKt[iKtBin-1];
			}
		      binWidthFactor > 1.0 ?  (binWidthFactor=1.0) : true ;
		      binWidthFactor<=0 ?   (binWidthFactor=1.0) : true;
		      //normalize so that the final points have the same maximum
		      if(maxValNorm< m_plotResults[resIdx].kTValues[iKtBin]/binWidthFactor)
			{
			  maxVal=m_plotResults[resIdx].kTValues[iKtBin];
			  maxValNorm=m_plotResults[resIdx].kTValues[iKtBin]/binWidthFactor;
			}

		    }
		  loc_plotResults=&m_plotResults[resIdx];
		  tree->Fill();
		  normFactor=1.0/maxValNorm;
	  

		  for(unsigned int iKtBin=0;iKtBin<numKtBins;iKtBin++)
		    {
		      float binWidthFactor=1.0;
		      if(0==iKtBin)
			{
			  binWidthFactor=binningKt[0];
			}
		      else
			{
			  binWidthFactor=binningKt[iKtBin]-binningKt[iKtBin-1];
			}
		      //for the last bin, it doesn't make sense to divide by 1000 or so...
		      binWidthFactor > 1.0 ?  (binWidthFactor=1.0) : true ;
		      binWidthFactor<=0 ?   (binWidthFactor=1.0) : true;
		      int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
		      //	      cout <<"looking at index:" << resIdx<<endl;
		      mX[iKtBin]=m_plotResults[resIdx].kTMeans[iKtBin];
		      //	      cout <<"mX["<<iKtBin <<"] " << mX[iKtBin]<<endl;
		      if((iKtBin>0)&& mX[iKtBin]<=mX[iKtBin-1]) 
			{
			  cout <<"MultiPlotter::saveGraph, wanting to set X["<<iKtBin<<"] to: " << mX[iKtBin] <<" but the one before is: " << mX[iKtBin-1] <<endl;
			  mX[iKtBin]=mX[iKtBin-1]+0.1;
			}
		      //	  cout <<"setting x: " << mX[iKtBin] <<endl;
		      mXErr[iKtBin]=0.0;
		      mY[iKtBin]=m_plotResults[resIdx].kTValues[iKtBin]*normFactor/binWidthFactor;
		      mYErr[iKtBin]=sqrt(m_plotResults[resIdx].kTValues[iKtBin])*normFactor/binWidthFactor;
		      //	      cout <<"mY["<<iKtBin <<"] " << mY[iKtBin]<<endl;
		      //	      cout <<"mYErr["<<iKtBin <<"] " << mYErr[iKtBin]<<endl;
		      //	  cout <<"y: " << mY[iKtBin] << ", " << mYErr[iKtBin] <<endl;
		    }
		  rFile.cd();
		  TGraphErrors graph(numKtBins,mX,mY,mXErr,mYErr);
		  sprintf(buffer1,"%s_ptSpect_%s_bin%d_%d",nameAddition.c_str(),buffer,i,j);
		  graph.SetName(buffer1);
		  graph.SetTitle(buffer1);
		  graph.GetYaxis()->SetTitle("normalized counts [arb. units]");
		  graph.GetXaxis()->SetTitle("kT [GeV]");
		  cout <<"saved as " << buffer1 <<endl;
		  graph.Write();

		}
	    }
	}
      //make sure this is saved...
    }
    }
  rFile.Write();
}


pair<int,int> MultiPlotter::pidBin2ZBinningIdx(int pidBin)
{
  pair<int,int> ret;
  ret.first=0;
  ret.second=0;
  if(pidBin< 0 || pidBin > PP)
    {
      cout <<"pidbin2zbinninningidx: wrong index" <<endl;
      return ret;
    }
   switch(pidBin)
       {
       case PiPi:
	 //already set to 0,0
	 return ret;
	 break;

       case PiK:
	 ret.second=1;
	 return ret;
	 break;

       case PiP:
	 ret.second=2;
	 return ret;
	 break;


	 break;
       case KPi:
	 ret.first=1;
	 return ret;

	 break;

       case KK:
	 ret.first=1;
	 ret.second=1;
	 return ret;
	 break;

       case KP:
	 ret.first=1;
	 ret.second=2;
	 return ret;
	 break;
       case PPi:
	 ret.first=2;
	 return ret;
	 break;

       case PK:
	 ret.first=2;
	 ret.second=1;
	 return ret;

	 break;

       case PP:
	 ret.first=2;
	 ret.second=2;
	 return ret;
	 break;

       }
   return ret;
 
}
//if accSmearing the number of pairs in the array might be different
//hp1 is the measured, hp2 the mc one, so we only check on hp1 if it is cut
void MultiPlotter::addSmearingEntry(HadronPairArray* hp1, HadronPairArray* hp2, bool accSmearing)
 {
  //needed for the mean computation...
   if(hp1->numPairs!=hp2->numPairs)
     {
       cout <<" smearing  hadron pairs not the same size: "<< hp1->numPairs <<", to " << hp2->numPairs <<endl;
       return;
     }
  for(int i=0;i<hp1->numPairs;i++)
    {
      if(hp1->cut[i] )
	{
	  //	  	  cout <<"hadron pair cut" <<endl;
	  continue;
	}
      else
	{
	  //	  cout <<"hadron pair survived " <<endl;
	}


      int chargeBin=hp2->chargeType[i];

      //      int particleBin1=hp->particleType1[i];
      //      int particleBin2=hp->particleType2[i];
      //      int particleBin=hp->particleType[i];

      //this is from MC, so there is no PID smearing, second hp is MC truth so we take that
      int pidBin=hp2->particleType[i];

      pair<int,int> zIdx2=pidBin2ZBinningIdx(pidBin);
      //      cout <<"pidBin: " << pidBin <<" chargeBin: " << chargeBin <<endl;
      //probably no matching particle
      if(pidBin<0) 
	continue;

     int kTBin1=getBin(binningKt,hp1->kT[i]);
     int kTBin2=getBin(binningKt,hp2->kT[i]);

     //       cout <<" kt data " << hp1->kT[i] << " kit mc: "<< hp2->kT[i] << " data PiPi " << hp1->kT_PiPi[i] << " mc " << hp2->kT_PiPi[i] <<endl;


     int  z2Bin1=getBin(binningZ[zIdx2.first],hp2->z1[i]);
     int z2Bin2=getBin(binningZ[zIdx2.second],hp2->z2[i]);
     int  z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1[i]);
     int z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2[i]);

     switch(pidBin)
       {
       case PiPi:
	 kTBin1=getBin(binningKt,hp1->kT_PiPi[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_Pi[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_Pi[i]);
	 break;

       case PiK:
	 kTBin1=getBin(binningKt,hp1->kT_PiK[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_Pi[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_K[i]);

	 break;

       case PiP:
	 kTBin1=getBin(binningKt,hp1->kT_PiP[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_Pi[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_P[i]);

	 break;
       case KPi:
	 kTBin1=getBin(binningKt,hp1->kT_KPi[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_K[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_Pi[i]);

	 break;

       case KK:
	 kTBin1=getBin(binningKt,hp1->kT_KK[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_K[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_K[i]);

	 break;

       case KP:
	 kTBin1=getBin(binningKt,hp1->kT_KP[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_K[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_P[i]);

	 break;
       case PPi:
	 kTBin1=getBin(binningKt,hp1->kT_PPi[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_P[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_Pi[i]);

	 break;

       case PK:
	 kTBin1=getBin(binningKt,hp1->kT_PK[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_P[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_K[i]);

	 break;

       case PP:
	 kTBin1=getBin(binningKt,hp1->kT_PP[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_P[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_P[i]);

	 break;

       }

      //      cout <<"kt: " << kT <<" bin: "<< kTBin<<endl;

     //     cout <<" first (mc) z: "<< hp1->z1[i] <<endl;
     //          cout <<" first (mc) z2: "<< hp1->z2[i] <<endl;
     
     //          cout <<" second (mc) z: "<< hp2->z1[i] <<endl;
	  //          cout <<" second (mc) z2: "<< hp2->z2[i] <<endl;
      ///let's only use the first z bin...


  pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);

     int numZBins1=binningZ[zIdx.first].size();
     int numZBins2=binningZ[zIdx.second].size();
      int recBin0=z1Bin2*numZBins1*numKtBins+z1Bin1*numKtBins+kTBin1;
      int iniBin0=z2Bin2*numZBins1*numKtBins+z2Bin1*numKtBins+kTBin2;
      int recBin1=z1Bin1*numKtBins+kTBin1;
      int iniBin1=z2Bin1*numKtBins+kTBin2;

      //      cout <<"ini bin: " << iniBin <<" recBin: " << recBin <<" z1Bin1: " << z1Bin1 << " kTBin1: " << kTBin1 << " z2bin1: " << z2Bin1 << " kTBin2: " << kTBin2<<endl;
      //      cout <<"pidBin: "<< pidBin <<" chargeBin " << chargeBin <<" iniBin: "<< iniBin <<endl;
      xini[0][pidBin][chargeBin]->Fill(iniBin0);
      bini[0][pidBin][chargeBin]->Fill(recBin0);
      xini[1][pidBin][chargeBin]->Fill(iniBin1);
      bini[1][pidBin][chargeBin]->Fill(recBin1);

      //true observable on the y axis, reconstrubted on the x axis
      kinematicSmearingMatrix[0][pidBin][chargeBin]->Fill(recBin0,iniBin0);
      kinematicSmearingMatrix[1][pidBin][chargeBin]->Fill(recBin1,iniBin1);
    }
 }



void MultiPlotter::addHadPairArray(HadronPairArray* hp, MEvent& event)
{
  //    cout <<"filling with " << hq->numHadQuads << " quads " << endl;
  //needed for the mean computation...
  this->thrustLabTheta=event.thrustThetaLab;
  this->cmsThrustTheta=event.thrustThetaCMS;

  this->thrust=event.Thrust;
  normalizeAngle(cmsThrustTheta);

  //  cout <<"multiplicity: " << multiplicity <<" mult bin: " <<multBin <<endl;
  thrustLabThetaBin=getBin(binningThrustLabTheta,thrustLabTheta);
  ///  cout <<"thrustLabTheta: "<< thrustLabTheta <<endl;
  //  cout <<"thrustLabtheta bin: " << thrustLabThetaBin <<endl;
  cmsThrustThetaBin=getBin(binningCmsThrustTheta,cmsThrustTheta);
  //  cout <<"got cmsThrust theta: " << cmsThrustTheta <<" bin: " << cmsThrustThetaBin <<" phi: " << cmsThrustPhi <<" bin: " << cmsThrustPhiBin<<endl;

  for(int i=0;i<hp->numPairs;i++)
    {
      if(hp->cut[i] )
	{
	  //	  	  cout <<"hadron pair cut" <<endl;
	  continue;
	}
      else
	{
	  //	  cout <<"hadron pair survived " <<endl;
	}

      //      int chargeBin1=hp->chargeType1[i];
      //      int chargeBin2=hp->chargeType2[i];
      int chargeBin=hp->chargeType[i];
      for(int p =PiPi;p<UNKNOWN;p++)
	{

      int pidBin=(int)p;

      //      int particleBin1=hp->particleType1[i];
      //      int particleBin2=hp->particleType2[i];
      //      int particleBin=hp->particleType[i];

      float weight=0.0;
      //don'te care for now...

      //      cout <<" pid: " << p <<endl;

      switch(p)
	{
	case PiPi:
	  this->z1=hp->z1_Pi[i];
	  this->z2=hp->z2_Pi[i];
	  this->kT=hp->kT_PiPi[i];
	  weight=hp->p_PiPi[i];
	  break;
	case PiK:
	  this->z1=hp->z1_Pi[i];
	  this->z2=hp->z2_K[i];
	  this->kT=hp->kT_PiK[i];
	  weight=hp->p_PiK[i];
	  break;
	case PiP:
	  this->z1=hp->z1_Pi[i];
	  this->z2=hp->z2_P[i];
	  this->kT=hp->kT_PiP[i];
	  weight=hp->p_PiP[i];
	  break;
	case KPi:
	  this->z1=hp->z1_K[i];
	  this->z2=hp->z2_Pi[i];
	  this->kT=hp->kT_KPi[i];
	  weight=hp->p_KPi[i];
	  break;
	case KK:
	  this->z1=hp->z1_K[i];
	  this->z2=hp->z2_K[i];
	  this->kT=hp->kT_KK[i];
	  weight=hp->p_KK[i];
	  break;
	case KP:
	  this->z1=hp->z1_K[i];
	  this->z2=hp->z2_P[i];
	  this->kT=hp->kT_KP[i];
	  weight=hp->p_KP[i];
	  break;

	case PPi:
	  this->z1=hp->z1_K[i];
	  this->z2=hp->z2_Pi[i];
	  this->kT=hp->kT_KPi[i];
	  weight=hp->p_PPi[i];
	  break;
	case PK:
	  this->z1=hp->z1_K[i];
	  this->z2=hp->z2_K[i];
	  this->kT=hp->kT_KK[i];

	  weight=hp->p_PK[i];
	  break;
	case PP:
	  this->z1=hp->z1_K[i];
	  this->z2=hp->z2_P[i];
	  this->kT=hp->kT_KP[i];
	  weight=hp->p_PP[i];
	  break;


	default:
	  cout <<"wrong pid " << endl;
	  this->z1=hp->z1[i];
	  this->z2=hp->z2[i];
	  this->kT=hp->kT[i];
	  weight =0.0;
	  cout <<"done with default " << endl;
	}

      this->qT=hp->qT[i];

      this->labTheta1=hp->labTheta1[i];
      this->labTheta2=hp->labTheta2[i];

      thrustBin=getBin(binningThrust,event.Thrust);

      qTBin=getBin(binningQt,qT);
      kTBin=getBin(binningKt,kT);
      //      cout <<"kt: " << kT <<" bin: "<< kTBin<<endl;

      pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);
      zbin1=getBin(binningZ[zIdx.first],this->z1);
      //      cout <<"zbin1: "<< zbin1 <<endl;
      zbin2=getBin(binningZ[zIdx.second],this->z2);

      labThetaBin1=getBin(binningLabTheta,hp->labTheta1[i]);
      labThetaBin2=getBin(binningLabTheta,hp->labTheta2[i]);
      //      cout <<"getting mass: " << hq->hp1.mass[i] <<endl;
	
      for(int bt=binType_labTheta_z; bt<binType_end;bt++)
	{
	  int firstBin=*(binningMap[bt].first);
	  int secondBin=*(binningMap[bt].second);
	  float firstKin=*(meanMap[bt].first);
	  float secondKin=*(meanMap[bt].second);
	  //	  cout << "firstBin " << firstBin <<" secondBin: "<< secondBin << " firstKin " << firstKin <<" secondKin " << secondKin <<endl;
	  //	  cout <<"chargeBin : " << chargeBin <<" bt: " << bt <<" pidBin: "<< pidBin <<endl;
	  if(bt<0 || chargeBin <0 || firstBin<0 || secondBin < 0)
	    {
	      //this gets called for all the woa events because they don't have thrust phi saved, 
	      //so let's exclude the binning with phi
	      //	      if(bt==binType_ThrustThetaPhi || bt==binType_ThrustPhiTheta)
	      //		continue;
	      //	      cout<<" hadOpen 1: " << hadronOpen1<<" and 2: " << hadronOpen2
	      //	      	      cout <<"bt: " << bt << " chargeBin: " << chargeBin << "firstBin: " << firstBin <<" second " << secondBin <<endl;
	      //		      cout <<"thrust: "<< this->thrust <<endl;
	      //	      cout <<"thrust theta: " <<this->cmsThrustTheta <<" thrust phi: "<< this->cmsThrustPhi<<endl;
	      //	      cout <<" z: " << this->z1 <<" z2: " << this->z2 << endl;
	      continue;
	    }

	  //	  chargeBin=0;

	  if(bt==binType_ThrustLabTheta_z)
	    {
	      //	      cout <<"filling thrusttheta bin: "<< firstBin << " zb in : "<< secondBin <<endl;
	    }
	  //	    cout <<"bt: " << bt <<" chargeBin: " << chargeBin<< " firstBin: " << firstBin << " second: " << secondBin <<" kt: "<< kTBin <<endl;
	  counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=weight;
	  meanValues_kin1[bt][pidBin][chargeBin][firstBin][secondBin]+=(weight*firstKin);
	  meanValues_kin2[bt][pidBin][chargeBin][firstBin][secondBin]+=(weight*secondKin);
	  meanValues_kT[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=(weight*kT);
	}
	}
    }
};

void MultiPlotter::setBinningMap()
{


  //add hoc hack of course with the problem that the max z bin is not pid dependent. So there might be empty bins
  
  maxKinMap=new vector< pair<int, int> >[NumPIDs];


  for(int bt=binType_labTheta_z; bt<binType_end;bt++)
    {
      switch(bt)
	{
	case binType_labTheta_z:
	  binningMap.push_back(pair<int*, int* >(&(this->labThetaBin1), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->labTheta1),&(this->z1)));
	  break;
	case binType_ThrustLabTheta_z:
	  binningMap.push_back(pair<int*, int* >(&(this->thrustLabThetaBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->thrustLabTheta),&(this->z1)));
	  break;
	case binType_z_z:
	  binningMap.push_back(pair<int*, int* >(&(this->zbin1), &(this->zbin2)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->z2)));
	  break;
	case binType_zOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->z1) ));
	  break;


	case binType_labThetaOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->labThetaBin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->labTheta1),&(this->labTheta1) ));
	  break;

	case binType_qTOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->qTBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->qT),&(this->qT) ));
	  break;


	case binType_ThrustOnly:
	  binningMap.push_back(pair<int*,int* > (&(this->zeroBin),&(this->thrustBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->thrust),&(this->thrust) ));
	  break;

	default:
	  cout <<"binning not recognized!!"<<endl;
	  exit(0);
	}

    }



  //we have to do maxKinMap extra, because here the z binning changes for the pids
  for(int pidBin=0;pidBin<NumPIDs;pidBin++)
    {
      pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);
      int maxZ1=zIdx.first;
      int maxZ2=zIdx.second;
      
      
      for(int bt=binType_labTheta_z; bt<binType_end;bt++)
	{
	  switch(bt)
	    {
	    case binType_labTheta_z:
	      maxKinMap[pidBin].push_back(pair<int,int>(binningLabTheta.size(),maxZ1));
	      break;
	    case binType_ThrustLabTheta_z:
	  maxKinMap[pidBin].push_back(pair<int,int>(binningThrustLabTheta.size(),maxZ1));
	  break;
	    case binType_z_z:
	      maxKinMap[pidBin].push_back(pair<int,int>(maxZ1,maxZ2));
	      break;
	    case binType_zOnly:
	      maxKinMap[pidBin].push_back(pair<int,int>(1,maxZ1));
	      break;
	      
	    case binType_labThetaOnly:
	      maxKinMap[pidBin].push_back(pair<int,int>(1,binningLabTheta.size()));
	      break;
	      
	    case binType_qTOnly:
	      maxKinMap[pidBin].push_back(pair<int,int>(1,binningLabTheta.size()));
	      break;
	      
	      
	    case binType_ThrustOnly:
	      maxKinMap[pidBin].push_back(pair<int,int>(1,binningThrust.size()));
	      break;
	      
	    default:
	      cout <<"binning not recognized!!"<<endl;
	      exit(0);
	    }
	  
	}
    }
}




////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


string MultiPlotter::getXAxisName(int binningType)
{
  string ret;
  switch(binningType)
    {

    case binType_z_z:
      ret+="z";
      break;

    case binType_labTheta_z:
      ret+="z";
      break;

    case binType_zOnly:
      ret+="z";
      break;

    case binType_labThetaOnly:
      ret+="lab #theta";
      break;
    case binType_qTOnly:
      ret+="Q_{T} [GeV]";
      break;
    case binType_ThrustOnly:
      ret+="Thrust";
      break;

    default:
      ret+="not_rec_";
      cout <<"wrong binning !!!" <<endl;
    }

  return ret;
}

//give negative values for first or second bin if they should not be part of the name
 string MultiPlotter::getBinName(int binningType,int pidType,int chargeType, int firstBin, int secondBin)
{
  string ret;
  switch(binningType)
    {
    case binType_z_z:
      ret+="z_z_";
      break;

    case binType_labTheta_z:
      ret+="labT_z_";
      break;

    case binType_ThrustLabTheta_z:
      ret+="labThrustTheta_z_";
      break;

    case binType_zOnly:
      ret+="onlyZ_";
      break;

    case binType_labThetaOnly:
      ret+="onlyLabT_";
      break;
    case binType_qTOnly:
      ret+="onlyQt_";
      break;

    case binType_ThrustOnly:
      ret+="onlyThrust";
      break;

    default:
      ret+="not_rec_";
      cout <<"wrong binning !!!" <<endl;
    }


  //enum pidType{PiPi, PiK, PiP, KPi, KK, KP, PPi, PK, PP, pidTypeEnd};
  switch(pidType)
    {

    case PiPi:
      ret+="_PiPi_";
      break;

    case PiK:
      ret+="_PiK_";
      break;

    case PiP:
      ret+="_PiP_";
      break;

    case KPi:
      ret+="_KPi_";
      break;

    case KK:
      ret+="_KK_";
      break;

    case KP:
      ret+="_KP_";
      break;

    case PPi:
      ret+="_PPi_";
      break;

    case PK:
      ret+="_PK_";
      break;

    case PP:
      ret+="_PP_";
      break;

    default:
      ret+="_notPID_";
      cout <<"wrong pid binning " <<endl;
    }

  switch(chargeType)
    {
    case pairChargeLikesign:
      ret+="ChargeLikeSign_";
      break;
    case pairPN:
      ret+="PN_";
      break;
    default:
      break;
      
    }

  char buffer[10];
  sprintf(buffer,"%d",firstBin);

  if(firstBin>=0)
    ret=ret+"bin_"+buffer;
  sprintf(buffer,"%d",secondBin);

  if(secondBin>=0)
    ret=ret+"bin_"+buffer;
  //  ret=ret+"_"+buffer;
  return ret;
}




void MultiPlotter::reorder(float* mX, float* mY, float* mYErr, int numBins)
{
  cout <<"reordering " <<endl;
  float tmpX[100];
  float tmpY[100];
  float tmpEY[100];
  int firstXBin=0;
  float minXVal=10000000;
  for(int i=0;i<numBins;i++)
    {
      tmpX[i]=mX[i];
      tmpY[i]=mY[i];
      tmpEY[i]=mYErr[i];
      if(mX[i]<minXVal)
	{
	  minXVal=mX[i];
	  firstXBin=i;
	}
    }
  for(int i=firstXBin;i<firstXBin+numBins;i++)
    {
      //wrap around
      int counter=i%numBins;

      mX[i-firstXBin]=tmpX[counter];
      mY[i-firstXBin]=tmpY[counter];
      mYErr[i-firstXBin]=tmpEY[counter];
      cout <<"mX["<<i-firstXBin <<"] from tmpX["<<counter<<"] is " <<  mX[i-firstXBin]<<endl;
    }
}


const int MultiPlotter::numKinematicBinning=7;
const int MultiPlotter::NumCharges=3;
const int MultiPlotter::NumPIDs=10;
//there is also an unknown flag which e.g. is used for all the electron/muon combinations
