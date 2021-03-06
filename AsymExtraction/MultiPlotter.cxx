#include "MultiPlotter.h"
#include "PlotResults.h"
string MultiPlotter::getParticlePairName(int p)
{
   switch(p)
       {
       case PiPi:
	 return "PiPi";
       case PiK:
	 return "PiK";
       case PiP:
	 return "PiP";
       case KPi:
	 return "KPi";
       case KK:
	 return "KK";

       case KP:
	 return "KP";
       case PPi:
	 return "PPi";

       case PK:
	 return "PK";

       case PP:
	 return "PP";

       }
   return "unknown";


}

//compute means and fill plotResults
void MultiPlotter::doPlots(bool print)
{
  for(int bt=binType_z_z; bt<binType_end;bt++)
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
			  if(print)
			    {
			      cout <<"looking at bt: " <<  bt <<" pidBin: "<< pidBin<<" firstbin: "<< firstBin;
			      cout <<" secondBin: "<< secondBin <<" kTBin : " << ktBin<<endl;
			    }

			  locCount+=counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  plotResults[resIdx].kTValues[ktBin]=counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  if(print)
			    {
			  cout <<"kt vals: "<< plotResults[resIdx].kTValues[ktBin];
			    }
			  plotResults[resIdx].kTValues1[ktBin]=counts1[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  plotResults[resIdx].kTValues2[ktBin]=counts2[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  plotResults[resIdx].kTUncertainties[ktBin]=sqrt(uncertainties[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]);
			  cout << " sys Uncert: "<< sysUncertainties[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]<<endl;


			  if(sysUncertainties[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]<0)
			    {
			      cout <<"sysUncertainties to small " <<endl;
			    }
			  plotResults[resIdx].kTSysUncertainties[ktBin]=sqrt(sysUncertainties[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]);
			  if(print)
			    {
			      cout <<"using counts: "<< counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			      cout <<", uncertainties: " << uncertainties[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			      cout <<", sqrt: "<< sqrt(uncertainties[bt][pidBin][chargeBin][firstBin][secondBin][ktBin])<<endl;
			    }

			  plotResults[resIdx].kTMeans[ktBin]=meanValues_kT[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]/counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  plotResults[resIdx].weakDecayFraction[ktBin]=weakDecayFraction[bt][pidBin][chargeBin][firstBin][secondBin][ktBin]/counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin];
			  //			  cout <<" weakFraction: " << weakDecayFraction[bt][pidBin][chargeBin][firstBin][secondBin][ktBin] << " counts: "<< counts[bt][pidBin][chargeBin][firstBin][secondBin][ktBin] <<" xxx: " <<  plotResults[resIdx].weakDecayFraction[ktBin] <<endl;
			  if(bt==binType_z_z && chargeBin==0)
			    {
			      if(print)
				{
				  cout <<"resIdx: "<< resIdx<<endl;
				  cout <<"kTbin: " <<ktBin << " kt count: "<< plotResults[resIdx].kTValues[ktBin] << " mean kt: "<< plotResults[resIdx].kTMeans[ktBin]<<endl;
				}
			    }
			}//end kt bins
		      if(locCount>0)
			{
			  plotResults[resIdx].meanKinBin1=meanValues_kin1[bt][pidBin][chargeBin][firstBin][secondBin]/locCount;
			  plotResults[resIdx].meanKinBin2=meanValues_kin2[bt][pidBin][chargeBin][firstBin][secondBin]/locCount;
			}
		      else{
			plotResults[resIdx].meanKinBin1=0;
			plotResults[resIdx].meanKinBin2=0;
		      }


		      plotResults[resIdx].firstKinBin=firstBin;
		      plotResults[resIdx].secondKinBin=secondBin;
		      if(bt==binType_z_z && chargeBin==0)
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
		      plotResults[resIdx].fileType=m_ftype;
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
  //binningZ[0].push_back(0.05);
  binningZ[0].push_back(0.1);
  binningZ[0].push_back(0.125);
  binningZ[0].push_back(0.15);
  binningZ[0].push_back(0.175);
  binningZ[0].push_back(0.2);
  binningZ[0].push_back(0.25);
  binningZ[0].push_back(0.3);
  ////  binningZ.push_back(0.7);
  binningZ[0].push_back(0.4); //// 
  binningZ[0].push_back(0.5);
  binningZ[0].push_back(0.6);
  binningZ[0].push_back(1.11);

  //the rest
  for(int i=1;i<6;i++)
    {
      binningZ[i].push_back(0.2);
      binningZ[i].push_back(0.25);
      binningZ[i].push_back(0.3);
      binningZ[i].push_back(0.4);
      ////  binningZ.push_back(0.7);
      binningZ[i].push_back(0.5); //// 
      binningZ[i].push_back(0.6); //// 
  //  binningZ.push_back(0.5);
  //  binningZ.push_back(0.6);
      binningZ[i].push_back(1.11);
    }



  if(useQt)
    {
      binningKt.push_back(1.0);
      binningKt.push_back(1.5);
      binningKt.push_back(2.0);
      binningKt.push_back(2.5);
      binningKt.push_back(3.0);
      binningKt.push_back(4.0);
      binningKt.push_back(6.0);
      binningKt.push_back(8.0);
      //  binningKt.push_back(1.5);
      binningKt.push_back(50.3);

    }
  else
    {
      binningKt.push_back(0.23);
      binningKt.push_back(0.35);
      binningKt.push_back(0.45);
      binningKt.push_back(0.55);
      binningKt.push_back(0.65);
      binningKt.push_back(0.78);
      binningKt.push_back(0.96);
      binningKt.push_back(2.5);
      binningKt.push_back(5000.3);
    }

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

  binningQt.push_back(1.0);
  binningQt.push_back(1.5);
  binningQt.push_back(2.0);
  binningQt.push_back(2.5);
  binningQt.push_back(3.0);
  binningQt.push_back(4.0);
  binningQt.push_back(6.0);
  binningQt.push_back(8.0);
  binningQt.push_back(100);

  //is from 0 to pi
  binningCmsThrustTheta.push_back(1.3);
  binningCmsThrustTheta.push_back(1.7);
  binningCmsThrustTheta.push_back(5.0);

  binningThrust.push_back(0.55);
  binningThrust.push_back(0.7);
  binningThrust.push_back(0.75);

  binningThrust.push_back(1.85);



};

void MultiPlotter::saveXini()
{
 rFile.cd();
 //numCharges is 3, but only the first two should be populated (likesign, unlikesign)
 for(int b=0;b<2;b++)
   {
     for(int c=0;c<NumCharges;c++)
       {
	 for(int p=0;p<NumPIDs;p++)
	   {
	     xini[b][p][c]->Write();
	}
       }
   }
 rFile.Write();
}
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
	  backgroundCounts[b][p][c]->Write();
	  kinematicSmearingMatrix[b][p][c]->Write();
      
	  xini[b][p][c]->Write();
	  bini[b][p][c]->Write();
	}
    }
    }
  rFile.Write();
}


//convert the convuoluted histogram we get from unfolding into 'regular' plots
//also put in the binwidth factors that were not used for the unfolding
//for the z_z binning, just put out the diagonal bins for now...
//
TH1D** MultiPlotter::convertUnfold2Plots(TH1D* input, int binning,  int chargeBin, int pidBin, const char* nameAdd)
{
  char buffer[300];
  float binWidthFactorZ1=1.0;
  float binWidthFactorZ2=1.0;
  float lowerZ1Cut=0.1;
  float lowerZ2Cut=0.1;
  if(pidBin==PiPi || pidBin== PiK || pidBin== PiP)
    lowerZ1Cut=0.05;
  if(pidBin==PiPi || pidBin== KPi || pidBin== PPi)
    lowerZ2Cut=0.05;

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

  Double_t value=0.0;
  Double_t valueUncertainty=0.0;
  //      int recBin=z1Bin1*numKtBins+kTBin1;
  for(int zBin=0;zBin<locMaxZBin;zBin++)
    {
	  if(1==binning)
	    {
	      //.first is only one bin, so 
	      binWidthFactorZ2=1.0-lowerZ2Cut;
      
	      //for the single z bin, the second counter is the z1 for some reason
	      if(zBin==0)
		binWidthFactorZ1=binningZ[zIdx.first][0]-lowerZ1Cut;
	      else
		binWidthFactorZ1=binningZ[zIdx.first][zBin]-binningZ[zIdx.first][zBin-1];
	    }
	  else
	    {
	      if(zBin==0)
		binWidthFactorZ1=binningZ[zIdx.first][0]-lowerZ1Cut;
	      else
		binWidthFactorZ1=binningZ[zIdx.first][zBin]-binningZ[zIdx.first][zBin-1];
	      if(zBin==0)
		binWidthFactorZ2=binningZ[zIdx.second][0]-lowerZ2Cut;
	      else
		binWidthFactorZ2=binningZ[zIdx.second][zBin]-binningZ[zIdx.second][zBin-1];
	    }

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
	  binWidthFactor*=binWidthFactorZ1*binWidthFactorZ2;
	  binWidthFactor > 100.0 ?  (binWidthFactor=1.0) : true ;
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
	  valueUncertainty=input->GetBinError(combBin+1);
	  //	  cout <<"value first " << value <<endl;
	  value=input->GetBinContent(combBin+1)/binWidthFactor;
	  valueUncertainty=input->GetBinError(combBin+1)/binWidthFactor;
	  //	  cout <<"now: "<< value<<endl;
	  ret[zBin]->SetBinContent(kTBin+1,value);
	  ret[zBin]->SetBinError(kTBin+1,valueUncertainty);
	}
    }
  return ret;
}



//return all z1,z2
TH1D*** MultiPlotter::convertAllUnfold2Plots(TH1D* input, int binning,  int chargeBin, int pidBin, const char* nameAdd)
{
  char buffer[300];

  pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);
  int maxZ1Bin = binningZ[zIdx.first].size();
  int maxZ2Bin = binningZ[zIdx.second].size();
  int maxZBin=(binningZ[zIdx.first].size()>binningZ[zIdx.second].size()) ? binningZ[zIdx.first].size() : binningZ[zIdx.second].size();
  

  float binWidthFactorZ1=1.0;
  float binWidthFactorZ2=1.0;
  float lowerZ1Cut=0.1;
  float lowerZ2Cut=0.1;
  if(pidBin==PiPi || pidBin== PiK || pidBin== PiP)
    lowerZ1Cut=0.05;
  if(pidBin==PiPi || pidBin== KPi || pidBin== PPi)
    lowerZ2Cut=0.05;




  TH1D*** ret=new TH1D**[maxZ1Bin];
  //this should be the smaller of the two since we are doing the diagonal bins
  for(int zBin1=0;zBin1<maxZ1Bin;zBin1++)
    {
      ret[zBin1]=new TH1D*[maxZ2Bin];
      for(int zBin2=0;zBin2<maxZ2Bin;zBin2++)
	{
	  sprintf(buffer,"un_convert_binning_%d_cBin_%d_pBin_%d_zBin1_%d_zBin2_%d_%s",binning,chargeBin,pidBin,zBin1,zBin2,nameAdd);
	  ret[zBin1][zBin2]=new TH1D(buffer,buffer,numKtBins,0,numKtBins);
	}
    }

  
  Double_t value=0.0;
  Double_t valueUncert=0.0;
  //      int recBin=z1Bin1*numKtBins+kTBin1;
  for(int zBin1=0;zBin1<maxZ1Bin;zBin1++)
    {
      for(int zBin2=0;zBin2<maxZ2Bin;zBin2++)
	{
	  if(1==binning)
	    {
	      //.first is only one bin, so 
	      binWidthFactorZ2=1.0-lowerZ2Cut;
      
	      //for the single z bin, the second counter is the z1 for some reason
	      if(zBin2==0)
		binWidthFactorZ1=binningZ[zIdx.first][0]-lowerZ1Cut;
	      else
		binWidthFactorZ1=binningZ[zIdx.first][zBin2]-binningZ[zIdx.first][zBin2-1];
	    }
	  else
	    {
	      if(zBin1==0)
		binWidthFactorZ1=binningZ[zIdx.first][0]-lowerZ1Cut;
	      else
		binWidthFactorZ1=binningZ[zIdx.first][zBin1]-binningZ[zIdx.first][zBin1-1];
	      if(zBin2==0)
		binWidthFactorZ2=binningZ[zIdx.second][0]-lowerZ2Cut;
	      else
		binWidthFactorZ2=binningZ[zIdx.second][zBin2]-binningZ[zIdx.second][zBin2-1];
	    }


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
	  binWidthFactor > 100.0 ?  (binWidthFactor=1.0) : true ;
	  binWidthFactor<=0 ?   (binWidthFactor=1.0) : true;
	  binWidthFactor*=binWidthFactorZ1*binWidthFactorZ2;
	  //the first one is z2, so the array size is multiplied with the max z1 bns
	  ///see : 
	  int combBin=zBin1*binningZ[zIdx.second].size()*numKtBins + zBin2*numKtBins + kTBin;
	  if(binning==1)//onlyZ
	    combBin=zBin1*numKtBins+kTBin;
	  //	  cout<<endl <<"combBin : " << combBin<<endl;
	  if(combBin>=input->GetNbinsX())
	    {
	      cout <<"convert all, zBin " << combBin <<" greater than " << input->GetNbinsX()<<endl;
	    }
	  //	  cout <<"convert unfold.. binWidthFactor: " << binWidthFactor<<endl;
	  value=input->GetBinContent(combBin+1);
	  valueUncert=input->GetBinError(combBin+1);
	  //	  cout <<"value first " << value <<endl;
	  value=input->GetBinContent(combBin+1)/binWidthFactor;
	  valueUncert=input->GetBinError(combBin+1)/binWidthFactor;
	  //	  cout <<"now: "<< value<<endl;
	  ret[zBin1][zBin2]->SetBinContent(kTBin+1,value);
	  ret[zBin1][zBin2]->SetBinError(kTBin+1,valueUncert);
	}
    }
    }
  return ret;
}

//mc_input is what is called xini in the tsvdunfold docu, MC_out is what is called bini
//needs stat uncertainties to construct bcov-->right now in 'statcov' but that takes it from the binerror. At some point this has to be set and also
//propagated from data1 to data --->actually is set as uncertainty on data1, so all we have to do is propagate...
//return three cov matrices: uncert due to limited MC statistics (getatdetcovmatrix), stat uncert (getxtau) and PID: propagate pid uncert (sys)
//with getunfoldcov matrix
TH1D* MultiPlotter::unfold(TH2D* smearingMatrix1, TH1D* MC_input1,TH1D* MC_out1, TH1D* data1, TH1D* sys1, TH1D** d, TH2D** statCov, TH2D** mcStatCov, TH2D** sysCov, const char* name)
{

  cout <<"unfolding with data dim: "  << data1->GetNbinsX() <<endl;

  vector<int> indicesThatIStillHave;
  vector<int> indicesThatIRemove;
  indicesThatIStillHave.push_back(350);
  indicesThatIRemove.push_back(449);
  indicesThatIRemove.push_back(907);
  indicesThatIRemove.push_back(1015);
  
  int countThreshold=0;
  vector<int> lowCountRows;
  for(int ix=0;ix<smearingMatrix1->GetNbinsX();ix++)
    {
      int countX=0;
      int countY=0;
      for(int iy=-1;iy<smearingMatrix1->GetNbinsY()+1;iy++)
	{
	  countX+=smearingMatrix1->GetBinContent(iy+1,ix+1);

	  if(find(indicesThatIStillHave.begin(),indicesThatIStillHave.end(),ix)!=indicesThatIStillHave.end())
		    {
		      if(smearingMatrix1->GetBinContent(iy+1,ix+1)>0)
			cout <<" found entry in "<<ix+1 <<", " << iy+1 << ": " << smearingMatrix1->GetBinContent(ix+1,iy+1) <<endl;
		    }

	}
      //this row has low counts, let's see if the column is small too
      //      	  cout <<"index : " << ix <<" rowSum: " << countX <<endl;
      if(find(indicesThatIStillHave.begin(),indicesThatIStillHave.end(),ix)!=indicesThatIStillHave.end())
	{
	  cout <<"index found : " << ix <<" rowSum: " << countX <<endl;
	}
      if(find(indicesThatIRemove.begin(),indicesThatIRemove.end(),ix)!=indicesThatIRemove.end())
	{
	  cout <<"removed index : " << ix <<" rowSum: " << countX <<endl;
	}
      //just check every column
      //      if(countX<=countThreshold)
      if(true)
	{
	  for(int iy=-1;iy<smearingMatrix1->GetNbinsY()+1;iy++)
	    {
	      //switch x and y to check column
	      countY+=smearingMatrix1->GetBinContent(ix+1,iy+1);
	      	  if(find(indicesThatIStillHave.begin(),indicesThatIStillHave.end(),ix)!=indicesThatIStillHave.end())
		    {
		      if(smearingMatrix1->GetBinContent(ix+1,iy+1)>0)
			cout <<" found entry in "<<ix+1 <<", " << iy+1 << ": " << smearingMatrix1->GetBinContent(ix+1,iy+1) <<endl;
		    }

		  
	    }
	  if(find(indicesThatIRemove.begin(),indicesThatIRemove.end(),ix)!=indicesThatIRemove.end())
	      {
		cout <<"removed index : " << ix <<" columnSum: " << countY <<endl;
	      }

	  if(find(indicesThatIStillHave.begin(),indicesThatIStillHave.end(),ix)!=indicesThatIStillHave.end())
	    {
	      cout <<"index : " << ix <<" otherDim: " << countY <<endl;
	    }

	  if(countY<=countThreshold)
	    {
	      lowCountRows.push_back(ix);
	    }
	}
      
    }
  //check for diagonal elements
  for(int ix=0;ix<smearingMatrix1->GetNbinsX();ix++)
    {
      	if(find(indicesThatIRemove.begin(),indicesThatIRemove.end(),ix)!=indicesThatIRemove.end())
	      {
		cout <<"removed index : " << ix <<" diagonal " << smearingMatrix1->GetBinContent(ix+1,ix+1) <<endl;
	      }

    if(smearingMatrix1->GetBinContent(ix+1,ix+1)<=countThreshold)
      {
	if(find(lowCountRows.begin(),lowCountRows.end(),ix)==lowCountRows.end())
	  {
	    lowCountRows.push_back(ix);
	  }

      }
    }
  
  //check for mcinput ==0
  for(int ix=0;ix<MC_input1->GetNbinsX();ix++)
    {
      	if(find(indicesThatIRemove.begin(),indicesThatIRemove.end(),ix)!=indicesThatIRemove.end())
	      {
		cout <<"removed index : " << ix <<" mc input " << MC_input1->GetBinContent(ix+1) <<endl;
	      }

      if(MC_input1->GetBinContent(ix+1)<=countThreshold)
	{
	  if(find(lowCountRows.begin(),lowCountRows.end(),ix)==lowCountRows.end())
	    {
	      lowCountRows.push_back(ix);
	    }
	  }
    }

  
  cout <<"found " << lowCountRows.size()<<" low count rows " <<endl;
  for(int i=0;i<lowCountRows.size();i++)
    {
      cout <<lowCountRows[i]<<", ";
    }
  cout <<endl;


  //smearingMatrix

  //  TH1D* xini=(TH1D*)file2.Get("xini");
  //  TH1D* bini=(TH1D*)file2.Get("bini");
  
  //  TH1D* d=(TH1D*)file.Get("bini");
  //  TH1D* k=(TH1D*)bini->Clone();
  
  //  for(Int_t i=0;i<1000;i++)
  //    {
  //      d->Fill(i%105);
  //
  //    }

  //  int maxBin=99*10;
  cout <<"bins in the beginning: "<< smearingMatrix1->GetNbinsX();
  int initialDimension=smearingMatrix1->GetNbinsX();
  int maxBin=initialDimension-lowCountRows.size();
  cout <<"maxBin: " << maxBin <<endl;
  cout <<" low countsRows size: "<< lowCountRows.size() <<" maxBins: "<< maxBin<<endl;
  TH2D* smearingMatrix=new TH2D("tmpSm","tmpSm",maxBin,0,maxBin,maxBin,0,maxBin);
  TH1D* MC_input=new TH1D("tmpMCIn","tmpMCIn",maxBin,0,maxBin);
  TH1D* MC_out=new TH1D("tmpMCOut","tmpMCOut",maxBin,0,maxBin);
  TH1D* data=new TH1D("data","data",maxBin,0,maxBin);
  TH1D* sys=new TH1D("tmpSys","tmpSys",maxBin,0,maxBin);


  int i2=0;
  //  for(int i=0;i<smearingMatrix->GetNbinsX();i++)
  for(int i=0;i<smearingMatrix1->GetNbinsX();i++)
    {
      //not one of the low count rows
      if(find(lowCountRows.begin(),lowCountRows.end(),i)==lowCountRows.end())
	{
	  MC_input->SetBinContent(i2+1,MC_input1->GetBinContent(i+1));
	  MC_input->SetBinError(i2+1,MC_input1->GetBinError(i+1));
	  MC_out->SetBinContent(i2+1,MC_out1->GetBinContent(i+1));
	  MC_out->SetBinError(i2+1,MC_out1->GetBinError(i+1));
	  sys->SetBinContent(i2+1,sys1->GetBinContent(i+1));
	  sys->SetBinError(i2+1,sys1->GetBinError(i+1));

	  data->SetBinContent(i2+1,data1->GetBinContent(i+1));
	  data->SetBinError(i2+1,data1->GetBinError(i+1));
	  int j2=0;
	  for(int j=0;j<smearingMatrix1->GetNbinsY();j++)
	    { 
	      if(find(lowCountRows.begin(),lowCountRows.end(),j)==lowCountRows.end())
		{
		  smearingMatrix->SetBinContent(i2+1,j2+1,smearingMatrix1->GetBinContent(i+1,j+1));
		  j2++;
		}
	    }
	  i2++;
	}
    }

  char buffer[300];
  sprintf(buffer,"statcov_%s%s",data->GetName(),name);
  TH2D* statcovMatrix=new TH2D(buffer,buffer,data->GetNbinsX(),0,data->GetNbinsX(),data->GetNbinsX(),0,data->GetNbinsX());
  sprintf(buffer,"syscov_%s%s",data->GetName(),name);
  TH2D* syscovMatrix=new TH2D(buffer,buffer,data->GetNbinsX(),0,data->GetNbinsX(),data->GetNbinsX(),0,data->GetNbinsX());
    cout <<" input has " << MC_input->GetNbinsX();
  cout <<" bins and output " << MC_out->GetNbinsX() <<" data: " << data->GetNbinsX() <<" smearing matrix " << smearingMatrix->GetNbinsX() <<" x " << smearingMatrix->GetNbinsY() <<endl;
  for(int i =0;i<MC_input->GetNbinsX();i++)
    {
      //      cout <<MC_input->GetBinContent(i+1) <<" ";
    }
  cout <<endl;
  for(int i =0;i<MC_out->GetNbinsX();i++)
    {
      //      cout <<MC_out->GetBinContent(i+1) <<" ";

    }
  cout <<endl<<endl;
  for(int i =0;i<data->GetNbinsX();i++)
    {
      cout <<data->GetBinContent(i+1) <<" error: " << data->GetBinError(i+1)<<" ";
      ///--->
      //only diagonal
      cout <<"setting statcov " << i+1 <<", " << i+1<<": " << data->GetBinError(i+1)*data->GetBinError(i+1);
      cout <<" fraction data: " << data->GetBinError(i+1)/data->GetBinContent(i+1)<<endl;
      cout <<"setting syscov " << i+1 <<", " << i+1 << ": " << sys->GetBinError(i+1)*sys->GetBinError(i+1);
      cout <<"should probably be  " << i+1 <<", " << i+1 << ": " << sys->GetBinContent(i+1)*sys->GetBinContent(i+1) <<" (from square of " << sys->GetBinContent(i+1) <<") ";
      //comparing with charlotte, it seems that sys->GetBinContent gives the correct sys error. The error on that is
      //just the squar root
      cout <<" fraction data: " << sys->GetBinError(i+1)/data->GetBinContent(i+1)<<endl;
      cout <<" or better: " << sys->GetBinContent(i+1)/data->GetBinContent(i+1)<<endl;
      statcovMatrix->SetBinContent(i+1,i+1,data->GetBinError(i+1)*data->GetBinError(i+1));
      syscovMatrix->SetBinContent(i+1,i+1,sys->GetBinContent(i+1)*sys->GetBinContent(i+1));
      //      for(int j =0;j<data->GetNbinsX();j++)
	{
	  //--->	  statcovMatrix->SetBinContent(i+1,j+1,data->GetBinError(i+1)*data->GetBinError(j+1));

	}
    }

 


  cout <<endl<<"smearingMatrix: " <<endl;
  double scalingFactor=1.0;
  //loop seems for scaling which we don't do anyways
////  for(int j =0;j<smearingMatrix->GetNbinsY();j++)
////    {
////      MC_out->SetBinContent(j+1,MC_out->GetBinContent(j+1)/scalingFactor);
////      MC_input->SetBinContent(j+1,MC_input->GetBinContent(j+1)/scalingFactor);
////      data->SetBinContent(j+1,data->GetBinContent(j+1)/scalingFactor);
////      statcovMatrix->SetBinContent(j+1,j+1,statcovMatrix->GetBinContent(j+1,j+1)/scalingFactor);
////      for(int i =0;i<smearingMatrix->GetNbinsX();i++)
////	{
////	  //	  cout <<smearingMatrix->GetBinContent(i+1,j+1) <<" " ;	  
////	  smearingMatrix->SetBinContent(i+1,j+1,smearingMatrix->GetBinContent(i+1,j+1)/scalingFactor);
////
////	  if(i==j)
////	    {
////	      //      smearingMatrix->SetBinContent(i+1,j+1,50.0);
////	    }
////	  else
////	    {
////	      //	    smearingMatrix->SetBinContent(i+1,j+1,0.0);
////	    }
////	}
////      //      cout <<endl;
////    }
  cout <<"smMatrixEnd"<<endl;
  TCanvas cT("c1","c1",0,0,2000,2000);
  cT.SetLogz();
  smearingMatrix->Draw("colz");
  cT.SaveAs("smTest123.png");

  cT.SaveAs("smTest123.root");

  cout <<endl;
  cout <<"--->"<<endl;
  int rank=MC_input->GetNbinsX();
  //  rank-=100;
  cout <<"input statcov: " << endl;
  for(int i=0;i<statcovMatrix->GetNbinsX();i++)
    {
      for(int j=0;j<statcovMatrix->GetNbinsY();j++)
	{
	  cout<< statcovMatrix->GetBinContent(i+1,j+1)<<" ";
	}
      cout <<endl;
    }

  TSVDUnfold* f=new TSVDUnfold(data,statcovMatrix,MC_out,MC_input,smearingMatrix);
  f->SetNormalize(false);
  cout <<"input statcov: " << endl;


  for(int i=0;i<statcovMatrix->GetNbinsX();i++)
    {
      for(int j=0;j<statcovMatrix->GetNbinsY();j++)
	{
	  cout << statcovMatrix->GetBinContent(i+1,j+1)<<" ";
	}
      cout <<endl;
    }
  

  
  //  f->SetNormalize(true);
    //  cout<<"normalized, using rank--> "<< rank <<endl;
  TH1D* ret=f->Unfold(rank);
  cout <<"use " <<rank <<" ranks " <<endl;

  cout <<"4"<<endl;
    TH2D* utaucov= f->GetXtau();
      (*statCov)=utaucov;
  //  (*statCov)=f->GetUnfoldCovMatrix(statcovMatrix,1200);
  cout <<"output statcov: " << endl;
  for(int i=0;i<(*statCov)->GetNbinsX();i++)
    {
      for(int j=0;j<(*statCov)->GetNbinsY();j++)
	{
	  cout << "("<<i+1<<","<<j+1<<": " <<(*statCov)->GetBinContent(i+1,j+1)<<" ";
	}
      cout <<endl;
    }
    cout <<"1 testing statcov 2,1: "<< (*statCov)->GetBinContent(2,1)<<endl;
  TH2D* uadetcov=f->GetAdetCovMatrix(10);
  (*mcStatCov)=uadetcov;
  cout <<"statcov dim: "<< (*statCov)->GetNbinsX()<<endl;
  ///  (*sysCov)=f->GetUnfoldCovMatrix(syscovMatrix,10);
  (*sysCov)=f->GetUnfoldCovMatrix(syscovMatrix,10);

  //return separate, -->if we would add here, that would also affect our return value, since we only saved the pointer above...
  ///    utaucov->Add(uadetcov);

  sprintf(buffer,"unfold_from_%s%s",data->GetName(),name);
  cout <<" unfolding " << buffer <<endl;
  TH1D* ret2=(TH1D*)ret->Clone(buffer);
  for(int i=0;i<ret2->GetNbinsX();i++)
    {
      //            ret2->SetBinError(i,TMath::Sqrt(utaucov->GetBinContent(i,i)));
    }

  cout <<"2 testing statcov 2,1: "<< (*statCov)->GetBinContent(2,1)<<endl;
  (*d)=f->GetD();
  for(int i=0;i<ret->GetNbinsX();i++)
    {
      //      cout <<"data: "<< data->GetBinContent(i+1) <<" unfolded: " << ret->GetBinContent(i+1) <<" mc out: "<< MC_out->GetBinContent(i+1) <<" mc in: "<< MC_input->GetBinContent(i+1)<<endl;
    }
  sprintf(buffer,"debug_D_from%s%s.png",MC_input->GetName(),name);
  cout <<"d: bins: "<< (*d)->GetNbinsX()<<endl;
  TCanvas c;
  (*d)->Draw();
  c.SetLogy();
  c.SaveAs(buffer);
  sprintf(buffer,"%s_ret3%s",ret2->GetTitle(),name);
  TH1D* ret3=new TH1D(buffer,buffer,initialDimension,ret2->GetBinLowEdge(1),ret2->GetBinLowEdge(ret2->GetNbinsX())+ret2->GetBinWidth(ret2->GetNbinsX()));
  cout <<"created ret3 with dim: "<< initialDimension<<" low edge: " << ret2->GetBinLowEdge(1) <<" high: " << ret2->GetBinLowEdge(ret2->GetNbinsX())+ret2->GetBinWidth(ret2->GetNbinsX()) <<endl;

  float lowerEdge=ret2->GetBinLowEdge(1);
  float upperEdge=ret2->GetBinLowEdge(ret2->GetNbinsX())+ret2->GetBinWidth(ret2->GetNbinsX());
  cout <<"3 testing statcov 2,1: "<< (*statCov)->GetBinContent(2,1)<<endl;
  sprintf(buffer,"%s_%s_statCov",data->GetName(),name);
  TH2D* statCovFullsize=new TH2D(buffer,buffer,initialDimension,lowerEdge,upperEdge,initialDimension,lowerEdge,upperEdge);
  sprintf(buffer,"%s_%s_mcStatCov",data->GetName(),name);
  TH2D* mcStatCovFullsize=new TH2D(buffer,buffer,initialDimension,lowerEdge,upperEdge,initialDimension,lowerEdge,upperEdge);
  sprintf(buffer,"%s_%s_sysStatCov",data->GetName(),name);
  TH2D* sysCovFullsize=new TH2D(buffer,buffer,initialDimension,lowerEdge,upperEdge,initialDimension,lowerEdge,upperEdge);


  int redCount=0;
  int redCount2=0;
  cout <<" original unfold: "<<endl;
  for(int i=0;i<ret2->GetNbinsX();i++)
    {
      cout <<i<<"-->"<<ret2->GetBinContent(i+1)<<", "<<endl;
    }
  cout <<"4 testing statcov 2,1: "<< (*statCov)->GetBinContent(2,1)<<endl;
  cout <<"filling output statcov" <<endl;
  redCount=0;
  for(int i=0;i<ret3->GetNbinsX();i++)
    {
      //not dropped
      if(find(lowCountRows.begin(),lowCountRows.end(),i)==lowCountRows.end())
	{
	  cout <<i+1 <<"("<<redCount+1<<")-->"<<ret2->GetBinContent(redCount+1)<<endl;
	  ret3->SetBinContent(i+1,ret2->GetBinContent(redCount+1));
	  redCount2=0;
	  for(int j=0;j<ret3->GetNbinsX();j++)
	    {
	      if(find(lowCountRows.begin(),lowCountRows.end(),j)==lowCountRows.end())
		{
		  if((i==1&& j==0) || (i==0 &&j==1) || (i==0&& j==0))
		      cout <<"5 testing statcov 2,1: "<< (*statCov)->GetBinContent(2,1)<<endl;
		    
		  cout <<i+1<<","<< j+1 <<"("<<redCount+1<<","<<redCount2+1<<")-->"<<(*statCov)->GetBinContent(redCount+1,redCount2+1)<<endl;
		  statCovFullsize->SetBinContent(i+1,j+1,(*statCov)->GetBinContent(redCount+1,redCount2+1));
		  sysCovFullsize->SetBinContent(i+1,j+1,(*sysCov)->GetBinContent(redCount+1,redCount2+1));
		  mcStatCovFullsize->SetBinContent(i+1,j+1,(*mcStatCov)->GetBinContent(redCount+1,redCount2+1));
		  redCount2++;
		}
	    }
	  cout <<endl;
	  redCount++;
	}
    }
  (*statCov)=statCovFullsize;
  (*sysCov)=sysCovFullsize;
  (*mcStatCov)=mcStatCovFullsize;

  return ret3;
}



void MultiPlotter::setHistogram(int binning, int chargeType, int pidType, TH1D* histo, TH1D* upperSys, TH1D* lowerSys)
{
  int binningType=binType_z_z;
  if(1==binning)
    binningType=binType_zOnly;

 PlotResults* m_plotResults=plotResults;


  for(int i=0;i<maxKinMap[pidType][binningType].first;i++)
    {
      //      for(int j=0;j<maxKinMap[binningType].second;j++)
      for(int j=0;j<maxKinMap[pidType][binningType].second;j++)
	{
	  int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
	  for(unsigned int iKtBin=0;iKtBin<numKtBins;iKtBin++)
	    {

	  float val=0;
	  float eVal=0;

	  float sysUp=0;
	  float sysDown=0;
	      if(binning==0)//for the z_z binning, the z1, z2 used for getResIdx are switchted...
		{
		 val= histo->GetBinContent(j*maxKinMap[pidType][binningType].first*numKtBins+i*numKtBins+iKtBin+1);
		 eVal= histo->GetBinError(j*maxKinMap[pidType][binningType].first*numKtBins+i*numKtBins+iKtBin+1);

		 sysUp=upperSys->GetBinContent(j*maxKinMap[pidType][binningType].first*numKtBins+i*numKtBins+iKtBin+1);
		 sysDown=lowerSys->GetBinContent(j*maxKinMap[pidType][binningType].first*numKtBins+i*numKtBins+iKtBin+1);


		}
	      else
		{
		  val=histo->GetBinContent(j*numKtBins+iKtBin+1);
		  eVal=histo->GetBinError(j*numKtBins+iKtBin+1);
		  sysUp=upperSys->GetBinContent(j*numKtBins+iKtBin+1);
		  sysDown=lowerSys->GetBinContent(j*numKtBins+iKtBin+1);
		}

	      m_plotResults[resIdx].kTValues[iKtBin]=val;
	      m_plotResults[resIdx].kTUncertainties[iKtBin]=eVal;
	      //let's hope that the systematic uncertainties are symmetricxb

	      m_plotResults[resIdx].kTSysUncertainties[iKtBin]=sysUp;
	      m_plotResults[resIdx].kTSysUncertaintiesLower[iKtBin]=sysDown;
	}

    }

    }
}
//void MultiPlotter::getIntAsymmetry(float a[3], float ea[3],int binningType,int chargeType, bool save1D)

//for now only for the zOnly binning--> changed to z1,z2 and zOnly (binning argument)
//don't do the bin width normalization since we also don't do it for the xini, bini

//pidcorrection:0 just regular ktValues, 1: ktValues1 (data) , 2: ktValues2 (mc)
TH1D* MultiPlotter::getHistogram(int binning, int chargeType, int pidType, const char* nameAdd, int getSys, int pidCorrection)

{
  //  cout <<"getting histo for binning: " << binning <<" charge: "<< chargeType <<" pidType: " << pidType <<endl;
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
  sprintf(buffer,"%s_%s",binName.c_str(),nameAdd);
  sprintf(buffer1,"histo_%s_%s",buffer,nameAdd);
  TH1D* ret=new TH1D(buffer1,buffer1,maxSmearing[pidType][binning],0,maxSmearing[pidType][binning]);
  //  cout <<" getH: " << maxSmearing <<endl;
    cout <<"saving graph for " << binName <<" buffer; " << buffer<<endl;
  for(int i=0;i<maxKinMap[pidType][binningType].first;i++)
    {
      //      for(int j=0;j<maxKinMap[binningType].second;j++)
      for(int j=0;j<maxKinMap[pidType][binningType].second;j++)
	{
	  pair<int,int> zIdx=pidBin2ZBinningIdx(pidType);
	  int maxZ1=binningZ[zIdx.first].size();
	  int maxZ2=binningZ[zIdx.second].size();

	  double normFactor=1.0;
	  double maxVal=-1.0;
	  double maxValNorm=-1.0;
	  int resIdx=getResIdx(binningType,pidType,chargeType,i,j);

	  float binWidthFactorZ1=1.0;
	  float binWidthFactorZ2=1.0;
	  float lowerZ1Cut=0.1;
	  float lowerZ2Cut=0.1;
	  if(pidType==PiPi || pidType== PiK || pidType== PiP)
	    lowerZ1Cut=0.05;
	  if(pidType==PiPi || pidType== KPi || pidType== PPi)
	    lowerZ2Cut=0.05;
	      
	  if(1==binning)
	    {
	      //.first is only one bin, so the other factor is the whole range minus the cutoff
		binWidthFactorZ2=1.0-lowerZ2Cut;


	      //for the single z bin, the second counter is the z1 for some reason
	      if(j==0)
		binWidthFactorZ1=binningZ[zIdx.first][0]-lowerZ1Cut;
	      else
		binWidthFactorZ1=binningZ[zIdx.first][j]-binningZ[zIdx.first][j-1];
	    }
	  else
	    {
	      if(i==0)
		binWidthFactorZ1=binningZ[zIdx.first][0]-lowerZ1Cut;
	      else
		binWidthFactorZ1=binningZ[zIdx.first][i]-binningZ[zIdx.first][i-1];
	      if(j==0)
		binWidthFactorZ2=binningZ[zIdx.second][0]-lowerZ2Cut;
	      else
		binWidthFactorZ2=binningZ[zIdx.second][j]-binningZ[zIdx.second][j-1];

	    }


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
	      /////	      binWidthFactor > 100.0 ?  (binWidthFactor=1.0) : true ;
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
	      binWidthFactor*=(binWidthFactorZ1*binWidthFactorZ2);
	      //and the z binning


	      //for the last bin, it doesn't make sense to divide by 1000 or so...
	      ////	      binWidthFactor > 100.0 ?  (binWidthFactor=1.0) : true ;
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

	      //--->pidco
	      float mKTValues=m_plotResults[resIdx].kTValues[iKtBin];
	      if(pidCorrection==1)
		{
		 mKTValues= m_plotResults[resIdx].kTValues1[iKtBin];
		}
	      if(pidCorrection==2)
		{
		 mKTValues= m_plotResults[resIdx].kTValues2[iKtBin];
		}
	      
	      //	      mY[iKtBin]=m_plotResults[resIdx].kTValues[iKtBin]*normFactor/binWidthFactor;
	      mY[iKtBin]=mKTValues*normFactor/binWidthFactor;
	      //	      Double_t binContent=m_plotResults[resIdx].kTValues[iKtBin]*normFactor/binWidthFactor;
	      	      Double_t binContent=mKTValues*normFactor/binWidthFactor;
	      Double_t binUncertainty=m_plotResults[resIdx].kTUncertainties[iKtBin]*normFactor/binWidthFactor;
	      Double_t binSysUncertainty=m_plotResults[resIdx].kTSysUncertainties[iKtBin]*normFactor/binWidthFactor;

	      //have to add one due to histos counting from 1
	      if(binning==0)//for the z_z binning, the z1, z2 used for getResIdx are switchted...
		{
		  if(!getSys)
		    {
		      ret->SetBinContent(i*maxKinMap[pidType][binningType].second*numKtBins+j*numKtBins+iKtBin+1,binContent);
		      ret->SetBinError(i*maxKinMap[pidType][binningType].second*numKtBins+j*numKtBins+iKtBin+1,binUncertainty);
		    }
		  else
		    {
		      ret->SetBinContent(i*maxKinMap[pidType][binningType].second*numKtBins+j*numKtBins+iKtBin+1,binSysUncertainty);
		    }
		}
	      else
		{
		  if(!getSys)
		    {
		      ret->SetBinContent(j*numKtBins+iKtBin+1,binContent);
		      ret->SetBinError(j*numKtBins+iKtBin+1,binUncertainty);
		    }
		  else{
		    ret->SetBinContent(j*numKtBins+iKtBin+1,binSysUncertainty);
		  }
		}
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


void MultiPlotter::printMatrix(TH1D* histo, const char* filename, bool saveUncert)
{
  cout <<" printing matrix of size " << histo->GetNbinsX() <<endl;
  double count=0;
  int numZBins=0;
  ofstream of;
  of.open(filename);
  for(int i=0;i<histo->GetNbinsX();i++)
    {
      if((i%9)==0)
	{
	  //	  of <<"new z bin "<<numZBins <<endl;
	  numZBins++;
	}

      if(saveUncert)
	{
	of << histo->GetBinError(i+1)<<" ";
	count+=histo->GetBinError(i+1);
	}
      else
	{
	of << histo->GetBinContent(i+1)<<" ";
	count+=histo->GetBinContent(i+1);
	}

    }
  of<<endl;
  of<<" total: " << count<<endl;
  of.close();

}


void MultiPlotter::printMatrix(TH2D* histo, const char* filename, bool saveUncert)
{
  ofstream of;
  of.open(filename);
  for(int i=0;i<histo->GetNbinsX();i++)
    {
      for(int j=0;j<histo->GetNbinsY();j++)
	{
	  if(saveUncert)
	    of << histo->GetBinError(i+1,j+1)<<" ";
	  else
	    of << histo->GetBinContent(i+1,j+1)<<" ";
	}
      //      cout <<endl;
    }
  of<<endl;
  of.close();

}



//just print a textfile with the bin content for the cross check
void MultiPlotter::printDebug(plotType mPlotType, int chargeType)
{
  cout <<"print debug, weightSum "<<weightSum <<endl;
  PlotResults* m_plotResults=plotResults;

  //  int numKinBin1=0;
  //  int numKinBin2=0;
  int total=0;
  char buffer[200];
  char buffer1[200];
  int binningType=binType_z_z;
  int pidType=PiPi;
  //  int chargeType=pairChargeLikesign;
  //    int chargeType=pairPN;
  //  for(int binningType=binType_labTheta_z; binningType<binType_end;binningType++)
    {
      //      for(int pidType=0;pidType<9;pidType++)
	{
	  //  for(int chargeType=0;chargeType<2;chargeType++)
	    {
	      string binName=getBinName(binningType,pidType,chargeType,-1,-1);
	      sprintf(buffer,"%s",binName.c_str());
	      for(int i=0;i<maxKinMap[pidType][binningType].first;i++)
		{
		  for(int j=0;j<maxKinMap[pidType][binningType].second;j++)
		    {
		      //		  cout <<" bin: " << i << ", " << j << endl;
		      double normFactor=1.0;
		      double maxVal=-1.0;
		      double maxValNorm=-1.0;
		      int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
		      //		  normFactor=1.0/maxValNorm;
		      normFactor=1.0;

		      for(unsigned int iKtBin=0;iKtBin<numKtBins;iKtBin++)
			{
			  float binWidthFactor=1.0;
			  
			  //for the last bin, it doesn't make sense to divide by 1000 or so...
			  int resIdx=getResIdx(binningType,pidType,chargeType,i,j);
			  //	      cout <<"looking at index:" << resIdx<<endl;
			  
			  //	  cout <<"setting x: " << mX[iKtBin] <<endl;
			  cout <<" " <<  m_plotResults[resIdx].kTValues[iKtBin];
			  cout <<" " <<  m_plotResults[resIdx].kTValues1[iKtBin];
			  cout <<" " <<  m_plotResults[resIdx].kTValues2[iKtBin];
			  cout <<" " <<  m_plotResults[resIdx].kTUncertainties[iKtBin];
			  cout <<" " <<  m_plotResults[resIdx].kTSysUncertainties[iKtBin];
			  total+=m_plotResults[resIdx].kTValues[iKtBin];

			}
		      //		  graph.GetYaxis()->SetTitle("normalized counts [arb. units]");
		      
		    }
		}
	    }
	  //make sure this is saved...
	}
    }

    cout <<endl<<" total: " << total <<endl;
}


void MultiPlotter::savePlots( plotType mPlotType, bool print)
{
  cout <<"in save plots " << endl;
  PlotResults* m_plotResults=plotResults;
  PlotResults* loc_plotResults=0;
  cout <<" accessing file " << endl;
  rFile.cd();
  cout <<" makeing tree.. " << endl;
  TTree *tree = new TTree("PlotTree","PlotTree");
  cout <<"about to branch" <<tree <<endl;
  tree->Branch("PlotBranch","PlotResults",&loc_plotResults,32000,99);

  cout <<"save plot after branch" <<endl;
  //  int numKinBin1=0;
  //  int numKinBin2=0;
  char buffer[200];
  char buffer1[200];
  float mX[50];
  float mY[50];
  float mXErr[50];
  float mYErr[50];
  for(int binningType=binType_z_z; binningType<binType_end;binningType++)
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
		      binWidthFactor > 100.0 ?  (binWidthFactor=1.0) : true ;
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
		  //		  normFactor=1.0/maxValNorm;
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
		      binWidthFactor > 100.0 ?  (binWidthFactor=1.0) : true ;
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
		  //		  graph.GetYaxis()->SetTitle("normalized counts [arb. units]");

		  graph.GetYaxis()->SetTitle("counts / [GeV]");

		  if(useQt)
		    graph.GetXaxis()->SetTitle("q_{T} [GeV]");
		  else
		    graph.GetXaxis()->SetTitle("k_{T} [GeV]");

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
  if(pidBin< 0 || pidBin > UNKNOWN)
    {
      cout <<"pidbin2zbinninningidx: wrong index" << pidBin <<endl;
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
//called the argument hp2 to be consistent with the mc entry of addSmearingEntry
void MultiPlotter::addXiniEntry(HadronPairArray* hp2)
{
  //  cout <<"adding xini entry for " << hp2->numPairs <<endl;
  for(int i=0;i<hp2->numPairs;i++)
    {
      int chargeBin=hp2->chargeType[i];
      //      int particleBin1=hp->particleType1[i];
      //      int particleBin2=hp->particleType2[i];
      //      int particleBin=hp->particleType[i];

      //this is from MC, so there is no PID smearing, second hp is MC truth so we take that
      int pidBin=hp2->particleType[i];
      //cut on dot product should not be necessary for MC truth
      if(hp2->cut[i] || pidDependentCut(hp2->z1[i],hp2->z2[i],hp2->kT[i],pidBin))
	continue;
      pair<int,int> zIdx2=pidBin2ZBinningIdx(pidBin);
      //            cout <<"pidBin: " << pidBin <<" chargeBin: " << chargeBin <<endl;
      //probably no matching particle

      if(pidBin<0) 
	continue;

     int kTBin2=getBin(binningKt,hp2->kT[i]);
     //       cout <<" kt data " << hp1->kT[i] << " kit mc: "<< hp2->kT[i] << " data PiPi " << hp1->kT_PiPi[i] << " mc " << hp2->kT_PiPi[i] <<endl;
     int  z2Bin1=getBin(binningZ[zIdx2.first],hp2->z1[i]);
     int z2Bin2=getBin(binningZ[zIdx2.second],hp2->z2[i]);
      //      cout <<"kt: " << kT <<" bin: "<< kTBin<<endl;

     //     cout <<" first (mc) z: "<< hp1->z1[i] <<endl;
     //          cout <<" first (mc) z2: "<< hp1->z2[i] <<endl;
     
     //          cout <<" second (mc) z: "<< hp2->z1[i] <<endl;
	  //          cout <<" second (mc) z2: "<< hp2->z2[i] <<endl;
      ///let's only use the first z bin...
     pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);

     int numZBins1=binningZ[zIdx.first].size();
     int numZBins2=binningZ[zIdx.second].size();
     //ini,rec bin 0 is for z_z binning, ini, rec 1 for zOnly
     int iniBin0=z2Bin1*numZBins2*numKtBins+z2Bin2*numKtBins+kTBin2;
     int iniBin1=z2Bin1*numKtBins+kTBin2;

      //      cout <<"ini bin: " << iniBin <<" recBin: " << recBin <<" z1Bin1: " << z1Bin1 << " kTBin1: " << kTBin1 << " z2bin1: " << z2Bin1 << " kTBin2: " << kTBin2<<endl;
      //      cout <<"pidBin: "<< pidBin <<" chargeBin " << chargeBin <<" iniBin: "<< iniBin <<endl;
      xini[0][pidBin][chargeBin]->Fill(iniBin0);
      xini[1][pidBin][chargeBin]->Fill(iniBin1);


            if(pidBin==0 && chargeBin==0)
	      {
		cout <<std::fixed;
		cout.precision(3);
		//		cout <<"Generated : z1: "<< hp2->z1[i] << " z2: "<< hp2->z2[i] << " sys: " << sys<<endl;
	  //	  cout <<"Generated : z1: "<< hp2->z1[i] << " z2: "<< hp2->z2[i] << " kT: "<< hp2->kT[i] << " bin: ";
	  //	  cout <<z2Bin1<<"*"<<numZBins2<<"*" << numKtBins<<"+"<<z2Bin2<<"*"<<numKtBins<<"+"<<kTBin2<<"="<< iniBin0<<endl;
	  //	  cout << "lab theta hadron1: " << hp2->labTheta1[i] << " hadron2: "<< hp2->labTheta2[i] << " lab phi 1: "<< hp2->labPhi1[i] << ", phi2: "<< hp2->labPhi2[i] <<endl;
//	  //	  cout << "cms theta hadron1: " << hp2->cmsTheta1[i] << " hadron2: "<< hp2->cmsTheta2[i] << " cms phi 1: "<< hp2->cmsPhi1[i] << " "<< hp2->cmsPhi2[i] <<endl;
//

	      }
	      

//      if(pidBin==0 && chargeBin==0 && (iniBin0==7 || iniBin0==1084 || iniBin0==1083|| iniBin0==1080))
//	{
//	  cout <<std::fixed;
//	  cout.precision(3);
//	  cout <<"Generated : z1: "<< hp2->z1[i] << " z2: "<< hp2->z2[i] << " kT: "<< hp2->kT[i] << " bin: ";
//	  cout <<z2Bin1<<"*"<<numZBins2<<"*" << numKtBins<<"+"<<z2Bin2<<"*"<<numKtBins<<"+"<<kTBin2<<"="<< iniBin0<<endl;
//	  cout << "lab theta hadron1: " << hp2->labTheta1[i] << " hadron2: "<< hp2->labTheta2[i] << " lab phi 1: "<< hp2->labPhi1[i] << ", phi2: "<< hp2->labPhi2[i] <<endl;
//	  //	  cout << "cms theta hadron1: " << hp2->cmsTheta1[i] << " hadron2: "<< hp2->cmsTheta2[i] << " cms phi 1: "<< hp2->cmsPhi1[i] << " "<< hp2->cmsPhi2[i] <<endl;
//
//
//	}


    }
}
//if accSmearing the number of pairs in the array might be different
//hp1 is the measured, hp2 the mc one, so we only check on hp1 if it is cut
void MultiPlotter::addSmearingEntry(HadronPairArray* hp1, HadronPairArray* hp2, bool accSmearing)
 {
   bool mcCut=false;
   bool accCut=false;
  //needed for the mean computation...
   if(evtNr==DEBUG_EVENT)
     {
       cout <<"looking at " << hp1->numPairs <<" hadron pairs" <<endl;
     }
   if(hp1->numPairs!=hp2->numPairs)
     {
       cout <<" smearing  hadron pairs not the same size: "<< hp1->numPairs <<", to " << hp2->numPairs <<endl;
       return;
     }
   else
     {
       //       cout <<"same size " <<endl;
     }
   for(int i=0;i<hp1->numPairs;i++)
    {
      //reset for each pair
      mcCut=false;
      accCut=false;

      ////-----> need a check here if hadron pair 2 (the mc) was cut or has -1 entries (no match found)
      ///that would mean that the entry is 'background'
      int pidBin=hp2->particleType[i];
      int chargeBin=hp2->chargeType[i];
      //      int particleBin1=hp->particleType1[i];
      //      int particleBin2=hp->particleType2[i];
      //      int particleBin=hp->particleType[i];


      //      cout <<"adding smearing entry for pid bin " << pidBin << " chargetype: "<< chargeBin<<endl;
      pair<int,int> zIdx2=pidBin2ZBinningIdx(pidBin);
      //      cout <<"pidBin: " << pidBin <<" chargeBin: " << chargeBin <<endl;
      //probably no matching particle
      if(evtNr==DEBUG_EVENT)
	{
	  cout << " pid bin: " << pidBin <<" charge bin: " << chargeBin  <<endl;
	}

      if(pidBin<0) 
	continue;




     int kTBin1=getBin(binningKt,hp1->kT[i]);
     int kTBin2=getBin(binningKt,hp2->kT[i]);

     //       cout <<" kt data " << hp1->kT[i] << " kit mc: "<< hp2->kT[i] << " data PiPi " << hp1->kT_PiPi[i] << " mc " << hp2->kT_PiPi[i] <<endl;
     int  z2Bin1=getBin(binningZ[zIdx2.first],hp2->z1[i]);
     int z2Bin2=getBin(binningZ[zIdx2.second],hp2->z2[i]);
     int  z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1[i]);
     int z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2[i]);

     float pidZ1=hp2->z1[i];
     float pidZ2=hp2->z2[i];
     float pidKt=hp2->kT[i];

      //this is from MC, so there is no PID smearing, second hp is MC truth so we take that

      //accepted hadron pair cut, should still be in xini
      if(hp1->cut[i])
	{
	  //	  	  cout <<"hadron pair cut" <<endl;
	  //	  continue;
	  accCut=true;
	}
      else
	{
	  //	  cout <<"hadron pair survived " <<endl;
	}


      //reconstructed pair but no mc pair
      //the hp2 will be cut due to the z cut, so the test for -1 is superfluous
      if(hp2->cut[i] || hp2->z1[i]==-1 || hp2->dp[i]>0)
	{
	  mcCut=true; 
	}


     switch(pidBin)
       {
       case PiPi:
	 if(hp1->z1_PiPi[i] > 1.01 || hp1->z2_PiPi[i]>1.01)
	   {
	     cout <<"indeed a z that is too high: " << hp1->z1_PiPi[i] <<" or " << hp2->z2_PiPi[i] <<endl;
	     continue;
	   }
	 kTBin1=getBin(binningKt,hp1->kT_PiPi[i]);
	 pidZ1=hp1->z1_PiPi[i];
	 pidZ2=hp1->z2_PiPi[i];
	 pidKt=hp1->kT_PiPi[i];
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_PiPi[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_PiPi[i]);

	 break;

       case PiK:
	 if(hp1->z1_PiK[i] > 1.01 || hp1->z2_PiK[i]>1.01)
	   continue;
	 kTBin1=getBin(binningKt,hp1->kT_PiK[i]);
	 pidZ1=hp1->z1_PiK[i];
	 pidZ2=hp1->z2_PiK[i];
	 pidKt=hp1->kT_PiK[i];
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_PiK[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_PiK[i]);

	 break;

       case PiP:
	 if(hp1->z1_PiP[i] > 1.01 || hp1->z2_PiP[i]> 1.01)
	   continue;
	 pidZ1=hp1->z1_PiP[i];
	 pidZ2=hp1->z2_PiP[i];
	 pidKt=hp1->kT_PiP[i];

	 kTBin1=getBin(binningKt,hp1->kT_PiP[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_PiP[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_PiP[i]);

	 break;
       case KPi:
	 if(hp1->z1_KPi[i] > 1.01 || hp1->z2_KPi[i]>1.01)
	   continue;

	 pidZ1=hp1->z1_KPi[i];
	 pidZ2=hp1->z2_KPi[i];
	 pidKt=hp1->kT_KPi[i];

	 kTBin1=getBin(binningKt,hp1->kT_KPi[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_KPi[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_KPi[i]);

	 break;

       case KK:
	 if(hp1->z1_KK[i] > 1.01 || hp1->z2_KK[i]> 1.01)
	   continue;
	 pidZ1=hp1->z1_KK[i];
	 pidZ2=hp1->z2_KK[i];
	 pidKt=hp1->kT_KK[i];

	 kTBin1=getBin(binningKt,hp1->kT_KK[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_KK[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_KK[i]);

	 break;

       case KP:
	 if(hp1->z1_KP[i] > 1.01 || hp1->z2_KP[i] > 1.01)
	   continue;

	 pidZ1=hp1->z1_KP[i];
	 pidZ2=hp1->z2_KP[i];
	 pidKt=hp1->kT_KP[i];

	 kTBin1=getBin(binningKt,hp1->kT_KP[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_KP[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_KP[i]);

	 break;
       case PPi:
	 if(hp1->z1_PPi[i] > 1.01 || hp1->z2_PPi[i]>1.01)
	   continue;

	 pidZ1=hp1->z1_PPi[i];
	 pidZ2=hp1->z2_PPi[i];
	 pidKt=hp1->kT_PPi[i];

	 kTBin1=getBin(binningKt,hp1->kT_PPi[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_PPi[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_PPi[i]);

	 break;

       case PK:
	 if(hp1->z1_PK[i] > 1.01 || hp1->z2_PK[i]> 1.01)
	   continue;

	 pidZ1=hp1->z1_PK[i];
	 pidZ2=hp1->z2_PK[i];
	 pidKt=hp1->kT_PK[i];

	 kTBin1=getBin(binningKt,hp1->kT_PK[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_PK[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_PK[i]);

	 break;

       case PP:
	 if(hp1->z1_PP[i] > 1.01 || hp1->z2_PP[i]>1.01)
	   continue;

	 pidZ1=hp1->z1_PP[i];
	 pidZ2=hp1->z2_PP[i];
	 pidKt=hp1->kT_PP[i];
	 kTBin1=getBin(binningKt,hp1->kT_PP[i]);
	 z1Bin1=getBin(binningZ[zIdx2.first],hp1->z1_PP[i]);
	 z1Bin2=getBin(binningZ[zIdx2.second],hp1->z2_PP[i]);

	 break;

       }
     if(pidDependentCut(pidZ1,pidZ2,pidKt,pidBin)|| cutOnDotProduct(hp1,i,pidBin))
       {
	 accCut=true;
       }
     //for mc, the z1, z2 already use the correct pid, might be a good idea to 
     //check on the dot product, but mc dot product not yet in tree (should be though)
     if(pidDependentCut(hp2->z1[i],hp2->z2[i],hp2->kT[i],pidBin))
       {
	 mcCut=true;
       }


     //     cout <<"kt (data): " << hp1->kT[i] << " pid: " << pidKt << " bin: "<< kTBin1<<endl;
     //     cout <<"kt (mc): " << hp2->kT[i] <<" bin: "<< kTBin2<<endl;


     //     cout <<" first (data) z: "<< hp1->z1[i] <<endl;
     //     cout <<" firpst (data) z2: "<< hp1->z2[i] <<endl;
     
     //     cout <<" second  (mc) z: "<< hp2->z1[i] << " pid: " << pidZ1 <<endl;
     //     cout <<" second  (mc) z2: "<< hp2->z2[i] << " pid: " << pidZ2 <<endl;
      ///let's only use the first z bin...
     pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);

     int numZBins1=binningZ[zIdx.first].size();
     int numZBins2=binningZ[zIdx.second].size();
     //ini,rec bin 0 is for z_z binning, ini, rec 1 for zOnly
     int recBin0=z1Bin1*numZBins2*numKtBins+z1Bin2*numKtBins+kTBin1;
     int iniBin0=z2Bin1*numZBins2*numKtBins+z2Bin2*numKtBins+kTBin2;
     int recBin1=z1Bin1*numKtBins+kTBin1;
     int iniBin1=z2Bin1*numKtBins+kTBin2;


     if(z1Bin1==z1Bin2==z2Bin1==z2Bin2 && chargeBin==0 && pidBin==0)
       {
	 //	 cout <<"looking at diagonal : " <<endl;
	 //	 cout <<"ini bin0: " << iniBin0 <<" recBin0: " << recBin0 <<" z1Bin1: " << z1Bin1 << " z1Bin2: "<< z1Bin2 <<" kTBin1: " << kTBin1 << " z2bin1: " << z2Bin1 <<" z2Bin2: "<< z2Bin2<< " kTBin2: " << kTBin2<<endl;
	 //	 cout <<"kt (data): " << hp1->kT[i] << " pid: " << pidKt << " bin: "<< kTBin1<<endl;
	 //	 cout <<"kt (mc): " << hp2->kT[i] <<" bin: "<< kTBin2<<endl;
	 //	 cout <<"ini bin0: " << iniBin0 <<" recBin0: " << recBin0 <<" z1Bin1: " << z1Bin1 << " z1Bin2: "<< z1Bin2 <<" kTBin1: " << kTBin1 << " z2bin1: " << z2Bin1 <<" z2Bin2: "<< z2Bin2<< " kTBin2: " << kTBin2<<endl;
	 //	 cout <<"ini bin1: " << iniBin1 <<" recBin1: " << recBin1 <<" z1Bin1: " << z1Bin1 << " kTBin1: " << kTBin1 << " z2bin1: " << z2Bin1 << " kTBin2: " << kTBin2<<endl;
	 //	 cout <<"pidBin: "<< pidBin <<" chargeBin " << chargeBin <<endl<<endl;;
	 //          cout <<" first (data) z: "<< hp1->z1[i] <<endl;
	  //          cout <<" first (data) z2: "<< hp1->z2[i] <<endl;
     
	  //          cout <<" second  (mc) z: "<< hp2->z1[i] << " pid: " << pidZ1 <<endl;
	  //          cout <<" second  (mc) z2: "<< hp2->z2[i] << " pid: " << pidZ2 <<endl;
       }
     //if mc is cut and data is not. In principle we could also expect TSVD to deal with background, but its easy to remove and not clear if TSVD handels it correctly
     if(mcCut && !accCut)
       {
	 backgroundCounts[0][pidBin][chargeBin]->Fill(recBin0);
	 backgroundCounts[1][pidBin][chargeBin]->Fill(recBin1);
	 //background events shouldn't contribute to smearing, since we subtract the bg
	 //before we apply the smearing inversion
	 continue;
       }

     //if mc is not cut 
     if(!mcCut)
       {
	 //we fill xini even if the bini is cut. TSVDUnfold should then deal with the acc correction
	 xini[0][pidBin][chargeBin]->Fill(iniBin0);
	 xini[1][pidBin][chargeBin]->Fill(iniBin1);
       }
     //since we already bailed out for mcCut && !accCut before, this is implicitly for !accCut && !mcCut
      if(!accCut)
	{
	  bini[0][pidBin][chargeBin]->Fill(recBin0);
	  bini[1][pidBin][chargeBin]->Fill(recBin1);

	  //true observable on the y axis, reconstrubted on the x axis
	}
      if(!accCut && !mcCut)
	{
	  if(recBin1==0 && iniBin1==0 && pidBin==0)
	    {
	      //	      cout <<"rec 1D both zero, recBin0: " << recBin0 << " iniBin0: "<< iniBin0 <<endl;
	      //	      cout << "lowest Diag entry 2D: "<< 	  kinematicSmearingMatrix[0][pidBin][chargeBin]->GetBinContent(1,1);
	      //	      cout << "lowest Diag entry 1D: "<< 	  kinematicSmearingMatrix[1][pidBin][chargeBin]->GetBinContent(1,1);
	    }
	  //for x-check with charlotte
	  //	  if(pidBin==0 && chargeBin==0 && iniBin0==3 && recBin0==4)
	  if(pidBin==0 && chargeBin==0 )//&& iniBin0==102 && recBin0==3)
	    {
//	     	      cout <<std::fixed;
//	     	      cout.precision(5);
//	     	      cout <<" MC Truth: z1: "<< hp2->z1[i] << " z2: "<< hp2->z2[i] << " kT: "<< hp2->kT[i] << " bin: ";
//	     	      cout <<z2Bin1<<"*"<<numZBins2<<"*" << numKtBins<<"+"<<z2Bin2<<"*"<<numKtBins<<"+"<<kTBin2<<"="<< iniBin0<<endl;
//	     	      cout <<"Reconstructed: z1: " << pidZ1 <<" z2: "<< pidZ2<< " kT: "<< pidKt << " bin: ";
//		      cout <<z1Bin1<<"*"<<numZBins2<<"*"<<numKtBins<<"+"<<z1Bin2<<"*"<<numKtBins<<"+"<<kTBin1<<"="<<recBin0<<endl;
	    }


	  kinematicSmearingMatrix[0][pidBin][chargeBin]->Fill(recBin0,iniBin0);
	  kinematicSmearingMatrix[1][pidBin][chargeBin]->Fill(recBin1,iniBin1);
	  if(recBin1==0 && iniBin1==0)
	    {
	      //	      cout << " after fill lowest Diag entry 2D: "<< 	  kinematicSmearingMatrix[0][pidBin][chargeBin]->GetBinContent(1,1);
	      //	      cout << "lowest Diag entry 1D: "<< 	  kinematicSmearingMatrix[1][pidBin][chargeBin]->GetBinContent(1,1)<<endl;
	    }
	}
    }
 }


bool MultiPlotter::cutOnDotProduct(HadronPairArray* hp, int index, int pidBin)
{
  
  
  switch(pidBin)
    {
    case PiPi:
      if(hp->dp_PiPi[index]<0)
	return false;
      break;
    case PiK:
      if(hp->dp_PiK[index]<0)
	return false;
      break;

    case PiP:
      if(hp->dp_PiP[index]<0)
	return false;
      break;

    case KPi:
      if(hp->dp_KPi[index]<0)
	return false;
      break;

    case KK:
      if(hp->dp_KK[index]<0)
	return false;
      break;

    case KP:
      if(hp->dp_KP[index]<0)
	return false;
      break;

    case PPi:
      if(hp->dp_PPi[index]<0)
	return false;
      break;

    case PK:
      if(hp->dp_PK[index]<0)
	return false;
      break;

    case PP:
      if(hp->dp_PP[index]<0)
	return false;
      break;

    default:
      //      cout <<"wrong pid in pid dependent cuts!!!: "<< pidBin<<endl;
      return true;

    }

  return true;
}




bool MultiPlotter::pidDependentCut(float z1, float z2, float kT, int pidBin )
{
  float z1Cut=0.05;
  float z2Cut=0.05;

  float z1UpperCut=1.1;
  float z2UpperCut=1.1;

  switch(pidBin)
    {
    case PiPi:
      z1Cut=zCutPi;
      z2Cut=zCutPi;
      break;
    case PiK:
      z1Cut=zCutPi;
      z2Cut=zCutPK;
      break;

    case PiP:
      z1Cut=zCutPi;
      z2Cut=zCutPK;
      break;

    case KPi:
      z1Cut=zCutPK;
      z2Cut=zCutPi;
      break;

    case KK:
      z1Cut=zCutPK;
      z2Cut=zCutPK;
      break;

    case KP:
      z1Cut=zCutPK;
      z2Cut=zCutPK;
      break;

    case PPi:
      z1Cut=zCutPK;
      z2Cut=zCutPi;
      break;

    case PK:
      z1Cut=zCutPK;
      z2Cut=zCutPK;
      break;

    case PP:
      z1Cut=zCutPK;
      z2Cut=zCutPK;
      break;

    default:
      //      cout <<"wrong pid in pid dependent cuts!!!: "<< pidBin<<endl;
      return true;

    }

  if(useQt && kT> 12.0)
    return true;
  if(!useQt && kT > 5.31145668)
    return true;
  if(z1 < z1Cut || z2< z2Cut || z1> z1UpperCut || z2> z2UpperCut) 
    {
      return true;
    }
  return false;
}


void MultiPlotter::addHadPairArray(HadronPairArray* hp, MEvent& event,bool print)
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
  if(numEvts>0)
    {
      //            cout <<"hadPairCount: " << hadPairCount <<" num evts: " << numEvts <<endl;
      //            cout <<"pair per event: " << (float)hadPairCount/(float)numEvts <<endl;
    }
  bool alreadyInc=false;

  //for xcheck
  int numWeak=0;
  int numStrong=0;
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

      /////           	  cout <<"event nr: " << event.evtNr <<endl;

      if(hp->z1[i]>0.1 && hp->z2[i]>0.1 && hp->particleType[i]==0 && hp->chargeType[i]==1)
	{
	  if(!alreadyInc)
	    {
	      alreadyInc=true;
	      numEvts++;
	    }
	  hadPairCount++;
	}
      if(print)
	{
	  cout << std::fixed;
	  cout <<setprecision(4);
	  cout <<endl;
	  	  cout <<"event nr: " << event.evtNr <<endl;
	  if(hp->chargeType[i]==0)
	    {
	      	    cout << " same-charged pair";
	    }
	  else
	    {
	      cout <<" opposite-charged pair";
	    }
	}
      for(int p =PiPi;p<UNKNOWN;p++)
	{

	  int pidBin=(int)p;

      //      int particleBin1=hp->particleType1[i];
      //      int particleBin2=hp->particleType2[i];
      //      int particleBin=hp->particleType[i];

	  float weight=0.0;

      //for the two alternative PID
	  float weight1=0.0;
	  float weight2=0.0;
      //to keep track of the uncertainty on the PID
	  float sys=0.0;
	  //don'te care for now...
	  
	  //      cout <<" pid: " << p <<endl;


      switch(p)
	{
	case PiPi:
	  this->z1=hp->z1_PiPi[i];
	  this->z2=hp->z2_PiPi[i];
	  this->kT=hp->kT_PiPi[i];
	  this->dotProduct=hp->dp_PiPi[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p) || dotProduct > 0.0)
	    continue;

	  weight1=hp->p_PiPi1[i];
	  weight2=hp->p_PiPi2[i];
	  //weight=(weight1+weight2)/2;

	  weight=hp->p_PiPi[i];
	  //
	  //
	  ////for xcheck, do unweighted
	  //if(hp->particleType[i]==0)
	  //	    weight=1;
	  //	  else
	  //	    weight=0;
	  
	  //cout <<"looking at weight1: "<< weight1 <<" weight2: " << weight2 <<" weight: "<< weight <<endl;
	  sys=hp->ep_PiPi[i];
	  //	  cout <<"sys is : "<< sys <<endl;
	  break;
	case PiK:
	  this->z1=hp->z1_PiK[i];
	  this->z2=hp->z2_PiK[i];
	  this->kT=hp->kT_PiK[i];
	  this->dotProduct=hp->dp_PiK[i];
	  if(pidDependentCut(z1,z2,kT,p) || dotProduct>0.0)
	    continue;



	  weight1=hp->p_PiK1[i];
	  weight2=hp->p_PiK2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_PiK[i];
	  sys=hp->ep_PiK[i];

	  //	  cout <<"adding weight: " << weight <<endl;
	  break;
	case PiP:
	  this->z1=hp->z1_PiP[i];
	  this->z2=hp->z2_PiP[i];
	  this->kT=hp->kT_PiP[i];
	  this->dotProduct=hp->dp_PiP[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p) || dotProduct > 0.0)
	    continue;

	  weight1=hp->p_PiP1[i];
	  weight2=hp->p_PiP2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_PiP[i];
	  sys=hp->ep_PiP[i];

	  break;
	case KPi:
	  this->z1=hp->z1_KPi[i];
	  this->z2=hp->z2_KPi[i];
	  this->kT=hp->kT_KPi[i];
	  this->dotProduct=hp->dp_KPi[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p)|| dotProduct >0.0)
	     continue;

	  weight1=hp->p_KPi1[i];
	  weight2=hp->p_KPi2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_KPi[i];
	  sys=hp->ep_KPi[i];

	  break;
	case KK:
	  this->z1=hp->z1_KK[i];
	  this->z2=hp->z2_KK[i];
	  this->kT=hp->kT_KK[i];
	  this->dotProduct=hp->dp_KK[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p) || dotProduct > 0.0)
	    continue;
	  weight1=hp->p_KK1[i];
	  weight2=hp->p_KK2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_KK[i];
	  sys=hp->ep_KK[i];

	  break;
	case KP:
	  this->z1=hp->z1_KP[i];
	  this->z2=hp->z2_KP[i];
	  this->kT=hp->kT_KP[i];
	  this->dotProduct=hp->dp_KP[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p) || dotProduct > 0.0)
	    continue;
	  weight1=hp->p_KP1[i];
	  weight2=hp->p_KP2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_KP[i];

	  sys=hp->ep_KP[i];

	  break;

	case PPi:
	  this->z1=hp->z1_PPi[i];
	  this->z2=hp->z2_PPi[i];
	  this->kT=hp->kT_KPi[i];
	  this->dotProduct=hp->dp_PPi[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p) || dotProduct > 0.0)
	    continue;
	  weight1=hp->p_PPi1[i];
	  weight2=hp->p_PPi2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_PPi[i];
	  sys=hp->ep_PPi[i];

	  break;
	case PK:
	  this->z1=hp->z1_PK[i];
	  this->z2=hp->z2_PK[i];
	  this->kT=hp->kT_PK[i];
	  this->dotProduct=hp->dp_PK[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p) || dotProduct >0.0)
	    continue;
	  weight1=hp->p_PK1[i];
	  weight2=hp->p_PK2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_PK[i];
	  sys=hp->ep_PK[i];


	  break;
	case PP:
	  this->z1=hp->z1_PP[i];
	  this->z2=hp->z2_PP[i];
	  this->kT=hp->kT_PP[i];
	  this->dotProduct=hp->dp_PP[i];
	  if(pidDependentCut(this->z1,this->z2,this->kT,p)|| dotProduct > 0.0)
	    continue;
	  weight1=hp->p_PP1[i];
	  weight2=hp->p_PP2[i];
	  //	  weight=(weight1+weight2)/2;
	  weight=hp->p_PP[i];
 	  sys=hp->ep_PP[i];

	  break;


	default:
	  cout <<"wrong pid " << endl;
	  this->z1=hp->z1[i];
	  this->z2=hp->z2[i];
	  this->kT=hp->kT[i];
	  weight =0.0;
	  weight1 =0.0;
	  weight2=0.0;
	  sys=0.0;

	  cout <<"done with default " << endl;
	}
      ////correct for the wrong factor we put in the trees
      sys*=(2*0.25*0.25);

      //doesn't matter, since the qt value is now in the kt field
      //but I guess this qT is still from some old code where one bins in qT (e.g asymmetries)
      this->qT=hp->qT[i];

      this->labTheta1=hp->labTheta1[i];
      this->labTheta2=hp->labTheta2[i];

      thrustBin=getBin(binningThrust,event.Thrust);
      //doesn't matter, since the qt value is now in the kt field
      qTBin=getBin(binningQt,qT);

      kTBin=getBin(binningKt,kT);
      //n      if(z1<0.1 && z2<0.1)
      //	cout <<"z1: "<< z1 <<" z2: "<< z2 <<" kt: " << kT <<" bin: "<< kTBin<<endl;

      pair<int,int> zIdx=pidBin2ZBinningIdx(pidBin);
      zbin1=getBin(binningZ[zIdx.first],this->z1);
      //      cout <<"zbin1: "<< zbin1 <<endl;
      zbin2=getBin(binningZ[zIdx.second],this->z2);

      labThetaBin1=getBin(binningLabTheta,hp->labTheta1[i]);
      labThetaBin2=getBin(binningLabTheta,hp->labTheta2[i]);
      //      cout <<"getting mass: " << hq->hp1.mass[i] <<endl;
	  


    if(print)
      {
	switch(p)
	  {

	  case PiPi:
	    if(hp->flip_PiPi[i])
	      cout <<" z1: " << hp->z2_PiPi[i] <<" z2: " << hp->z1_PiPi[i] << " kT: " << hp->kT_PiPi[i];
	    else
	      cout <<" z1: " << hp->z1_PiPi[i] <<" z2: " << hp->z2_PiPi[i] << " kT: " << hp->kT_PiPi[i]<<endl;;


	    cout <<" (combined/data/mc) probabilities for ";
	    cout << "PiPi: " << hp->p_PiPi[i] <<" / " << hp->p_PiPi1[i] << " / " << hp->p_PiPi2[i];
	    break;
	  case PiK:
	    //check if the assignment is flipped that would produce this
	    if(hp->flip_KPi[i])
	      cout << "PiK: " << hp->p_KPi[i] <<" / " << hp->p_KPi[i] << " / " << hp->p_KPi2[i];
	    else
	      cout << "PiK: " << hp->p_PiK[i] <<" / " << hp->p_PiK[i] << " / " << hp->p_PiK2[i];
	    break;


	  case PiP:
	    if(hp->flip_PPi[i])
	      cout << "PiP: " << hp->p_PPi[i] <<" / " << hp->p_PPi1[i] << " / " << hp->p_PPi2[i];
	    else
	      cout << "PiP: " << hp->p_PiP[i] <<" / " << hp->p_PiP1[i] << " / " << hp->p_PiP2[i];
	    break;

	  case KPi:
	    if(hp->flip_PiK[i])
	      cout << "KPi: " << hp->p_PiK[i] <<" / " << hp->p_PiK1[i] << " / " << hp->p_PiK2[i];
	    else
	      cout << "KPi: " << hp->p_KPi[i] <<" / " << hp->p_KPi1[i] << " / " << hp->p_KPi2[i];
	    break;

	  case KK:
	    cout << "KK: " << hp->p_KK[i] <<" / " << hp->p_KK[i] << " / " << hp->p_KK2[i];
	    break;

	  case KP:
	    //	    cout <<"dp: "<<hp->dp_KP[i] <<" dp PK: " << hp->dp_PK[i]<< " z1 kp: " << hp->z1_KP[i] <<" z2 kp: " << hp->z2_KP[i] << " kt kP: "<< hp->kT_KP[i] <<endl;
	    if(hp->flip_PK[i])
	      cout << "KP (flipped): " << hp->p_PK[i] <<" / " << hp->p_PK1[i] << " / " << hp->p_PK2[i];
	    else
	      cout << "KP: " << hp->p_KP[i] <<" / " << hp->p_KP1[i] << " / " << hp->p_KP2[i];
	    break;

	  case PPi:
	    if(hp->flip_PiP[i])
	      cout << "PPi: " << hp->p_PiP[i] <<" / " << hp->p_PiP1[i] << " / " << hp->p_PiP2[i];
	    else
	      cout << "PPi: " << hp->p_PPi[i] <<" / " << hp->p_PPi1[i] << " / " << hp->p_PPi2[i];
	    break;

	  case PK:
	    cout <<"dp: "<<hp->dp_PK[i] <<" dp PK: " << hp->dp_PK[i]<< " z1 pk: " << hp->z1_PK[i] <<" z2 pk: " << hp->z2_PK[i] << " kt kP: "<< hp->kT_PK[i] <<endl;
	    if(hp->flip_KP[i])
	      cout << "PK (flipped): " << hp->p_KP[i] <<" / " << hp->p_KP[i] << " / " << hp->p_KP2[i];
	    else
	      cout << "PK: " << hp->p_PK[i] <<" / " << hp->p_PK[i] << " / " << hp->p_PK2[i];
	    
	    break;


	  case PP:
	    cout << "PP: " << hp->p_PP[i] <<" / " << hp->p_PP1[i] << " / " << hp->p_PP2[i];
	    break;
	  }

	cout <<endl;
	cout <<" combined sys uncert for ";
	switch(p)
	  {
	  case PiPi:
	    cout <<" PiPi: " << hp->ep_PiPi[i];
	    break;
	  case PiK:
	    if(hp->flip_KPi[i])
	      cout <<" PiK: " << hp->ep_KPi[i];
	    else
	      cout <<" PiK: " << hp->ep_PiK[i];
	    break;

	  case PiP:
	    if(hp->flip_PPi[i])
	      cout <<" PiP: " << hp->ep_PPi[i];
	    else
	      cout <<" PiP: " << hp->ep_PiP[i];
	    break;

	  case KPi:
	    if(hp->flip_PiK[i])
	      cout <<" KPi: " << hp->ep_PiK[i];
	    else
	      cout <<" KPi: " << hp->ep_KPi[i];
	    break;

	  case KK:
	    cout <<" KK: " << hp->ep_KK[i];
	    break;

	  case KP:
	    if(hp->flip_PK[i])
	      cout <<" KP: " << hp->ep_PK[i];
	    else
	      cout <<" KP: " << hp->ep_KP[i];
	    break;

	  case PPi:
	    if(hp->flip_PiP[i])
	      cout <<" PPi: " << hp->ep_PiP[i];
	    else
	      cout <<" PPi: " << hp->ep_PPi[i];
	    break;

	  case PK:
	    if(hp->flip_KP[i])
	      cout <<" PK: " << hp->ep_KP[i];
	    else
	      cout <<" PK: " << hp->ep_PK[i];
	    break;

	  case PP:
	    cout <<" PP: " << hp->ep_PP[i];
	  }

	cout <<endl;
      }


	  if(print)
	    {
	      if(hp->particleType[i]==PiPi && p==PiPi)
		{
		  		  cout <<" found pipi, z1 pipi: " << hp->z1_PiPi[i]<<" z2: "<< hp->z2_PiPi[i];
		  //		  cout <<" kt pipi: " << hp->kT_PiPi[i] <<" weight: " << hp->p_PiPi[i] <<endl;
		  //		  cout <<"hadron pType: " << hp->particleType[i];
		  		  cout <<" z1: " << hp->z1[i] <<", z2: " <<hp->z2[i];
		  //		  cout <<" kt: "<< hp->kT[i]<<endl;
		  

		  }



	    }
	
      for(int bt=binType_z_z; bt<binType_end;bt++)
	{
	  int firstBin=*(binningMap[bt].first);
	  int secondBin=*(binningMap[bt].second);
	  float firstKin=*(meanMap[bt].first);
	  float secondKin=*(meanMap[bt].second);

	  if(hp->particleType[i]==PiPi && p==PiPi && print)
	    {
	  	  cout << "firstBin " << firstBin <<" secondBin: "<< secondBin << " firstKin " << firstKin <<" secondKin " << secondKin <<endl;
	  	  cout <<"chargeBin : " << chargeBin <<" bt: " << bt <<" pidBin: "<< pidBin <<endl;
	    }
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
	  //	  	    cout <<"bt: " << bt <<" chargeBin: " << chargeBin<< " firstBin: " << firstBin << " second: " << secondBin <<" kt: "<< kTBin <<endl;
	  //		    cout <<"weight: "<< weight <<endl;



	  
	  /////for x-check, looking at PiPi weights
	  //	  if(hp->chargeType[i]==0 && p==PiPi && this->z1 <0.1 && this->z2 < 0.1 && this->kT < 0.45)
	  if(hp->chargeType[i]==0 && hp->particleType[i]==PiPi&&  p==PiPi && bt==binType_z_z)
	    {
	      weightSum+=weight;
	      cout << std::fixed;
	      cout <<setprecision(4);

	      

	      //	      cout <<"event nr: " << event.evtNr << " weight: " << weight <<endl;
	      //	      cout <<"event nr: " << event.evtNr << 
	      //	            cout << "accepted z1: " << this->z1 <<" z2: "<< this->z2 << " sys: " << sys <<endl;
	      //  cout << "event nr: "<< event.evtNr-1 <<" z1: " << hp->z1[i] <<" z2: "<< hp->z2[i] << " cms theta1: " << hp->cmsTheta1[i]<<" cms theta2: " <<hp->cmsTheta2[i]<< "  weak: " << hp->isWeakDecay[i]<<endl;

	      if(hp->isWeakDecay[i])
		numWeak++;
	      else
		numStrong++;


	      //	      cout << weight <<endl;
	      //		      	      	      	      	      	      cout << "accepted z1: " << this->z1 <<" z2: "<< this->z2 << " kT: "<< this->kT<< " weight1: "<<weight1 <<" weight2: "<< weight2 << " weight: "<< weight<< endl;
	      //		      	                  cout <<"originally identified as " << getParticlePairName(hp->particleType[i]) <<" theta cms1: "<< hp->cmsTheta1[i] <<" theta cms2: "<< hp->cmsTheta2[i]<<endl;
	      //      	            		      cout <<"theta lab1: "<< hp->labTheta1[i] <<", theta lab2: "<< hp->labTheta2[i] <<endl;
	      //	                  cout <<"particle type1: "<< hp->particleType1[i] << " particle type2: "<< hp->particleType2[i]<<endl;

	    }


	  /////
	  //	  cout <<"increasing counts: "<< counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]<<endl;
	  counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=weight;
	  if(weight==0)
	    {
	      //	      cout <<"weight zero for bt: "<< bt <<" pid: "<< pidBin <<" charge: "<< chargeBin <<" k1: " <<firstBin <<" k2: "<< secondBin <<" Ktbin: "<< kTBin <<endl; 
	    }

	  
	  //	  cout <<" now: " <<  counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]<<endl;
	  counts1[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=weight1;
	  counts2[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=weight2;
	  uncertainties[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=(weight*weight);

	  if(isnan(sys*sys))
	    {
	      cout <<"sys*sys is nan: "<< sys<<endl;
	      sys=0;
	    }

	  //	  sysUncertainties[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=(sys*sys);
	  //the factor is because in the trees I have a factor of 0.5 for each term, but in reality that should be 0.25^2, since there is a factor of 1/4 in front of the square root, so we have to do 1/0.5 * 0.25^2
	  sysUncertainties[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=sys;
	  //	  	  cout <<"adding sys: "<< sys*sys <<" adding to " << sysUncertainties[bt][pidBin][chargeBin][firstBin][secondBin][kTBin] <<endl;

	  meanValues_kin1[bt][pidBin][chargeBin][firstBin][secondBin]+=(weight*firstKin);
	  meanValues_kin2[bt][pidBin][chargeBin][firstBin][secondBin]+=(weight*secondKin);
	  //	  cout <<"filling mean val for kt bin : " << kTBin <<" with kT: " << kT << " weight: " << weight <<endl;
	  meanValues_kT[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=(weight*kT);

	  //need to multiply by weight since otherwise even in the woa case we are adding weight 0 events
	  weakDecayFraction[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]+=weight*hp->isWeakDecay[i];
	  //	  cout <<"weakDecayFracinevent: " << weakDecayFraction[bt][pidBin][chargeBin][firstBin][secondBin]<<" " <<endl;
	  if(weakDecayFraction[bt][pidBin][chargeBin][firstBin][secondBin][kTBin]>counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin])
	    {
	      //	      	      cout <<" weight: " << weight <<endl;
	      //	      	      cout << " weak greater counts " << weakDecayFraction[bt][pidBin][chargeBin][firstBin][secondBin][kTBin] << " is greater than " << counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin] <<endl;
	    }
	   

	  
	  if(chargeBin==pairChargeLikesign && pidBin==PiPi && bt==binType_z_z)
	    {
	      if(print)
		{
		      int resIdx=getResIdx(bt,pidBin,chargeBin,firstBin,secondBin);
		      cout <<"residx: "<< resIdx<<endl;
	      	      cout <<"filling in bin " << firstBin*(maxKinMap[pidBin][bt].second*numKtBins)+secondBin*numKtBins+kTBin<<endl;
		      cout <<" firstBin:" << firstBin <<" secondBin: " << secondBin <<" kTBin: "<< kTBin<<endl;
		      cout <<"counts : " << 	  counts[bt][pidBin][chargeBin][firstBin][secondBin][kTBin];
		      cout <<" counts2 : " << 	  counts2[bt][pidBin][chargeBin][firstBin][secondBin][kTBin];
		      cout <<" mean1 : "<< 	  meanValues_kin1[bt][pidBin][chargeBin][firstBin][secondBin];
		      cout <<" mean2 : "<< 	  meanValues_kin2[bt][pidBin][chargeBin][firstBin][secondBin];
		      
		}
	    }
	}
	}
    }
    if(numWeak>0 || numStrong>0)
      {
  //    cout <<"event nr: " << event.evtNr-1 <<" num weak: "<< numWeak << " num strong: "<< numStrong<<endl;
      }


};

void MultiPlotter::setBinningMap()
{


  //add hoc hack of course with the problem that the max z bin is not pid dependent. So there might be empty bins
  
  maxKinMap=new vector< pair<int, int> >[NumPIDs];


  //have to run over all
  for(int bt=binType_qTOnly; bt<binType_end;bt++)
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
      int maxZ1=binningZ[zIdx.first].size();
      int maxZ2=binningZ[zIdx.second].size();
      
      
      for(int bt=binType_qTOnly; bt<binType_end;bt++)
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
void MultiPlotter::weight(float weightFactor)
{
  for(int bt=binType_z_z; bt<binType_end;bt++)
    {
      for(int chargeBin=0;chargeBin<2;chargeBin++)
	{
	  for(int pidBin=0;pidBin<NumPIDs;pidBin++)
	    {
	      for(int firstBin=0;firstBin<maxKinMap[pidBin][bt].first;firstBin++)
		{
		  for(int secondBin=0;secondBin<maxKinMap[pidBin][bt].second;secondBin++)
		    {
		      int resIdx=getResIdx(bt,pidBin,chargeBin,firstBin,secondBin);
		      for(unsigned int ktBin=0;ktBin<numKtBins;ktBin++)
			{
			  plotResults[resIdx].kTValues[ktBin]*=weightFactor;
			  plotResults[resIdx].kTValues1[ktBin]*=weightFactor;
			  plotResults[resIdx].kTValues2[ktBin]*=weightFactor;
			  //shouldn't have to weight the weakDecayFraction since we are already weighting the kTValue
			  //
			}
		    }
		}
	    }
	}
    }

}

const int MultiPlotter::numKinematicBinning=7;
const int MultiPlotter::NumCharges=3;
const int MultiPlotter::NumPIDs=10;
//there is also an unknown flag which e.g. is used for all the electron/muon combinations
