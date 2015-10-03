#include "MultiPlotter.h"
#include "PlotResults.h"


//compute means and fill plotResults
void MultiPlotter::doPlots()
{

  for(int bt=binType_labTheta_z; bt<binType_end;bt++)
    {
      //      cout <<"looking at bin type " << bt <<endl;
      //      for(int chargeBin=0;chargeBin<NumCharges;chargeBin++)
      for(int chargeBin=0;chargeBin<1;chargeBin++)
	{
	  for(int firstBin=0;firstBin<maxKinMap[bt].first;firstBin++)
	    {
	      for(int secondBin=0;secondBin<maxKinMap[bt].second;secondBin++)
		{
		  double locCount=0;
		  int resIdx=getResIdx(bt,chargeBin,firstBin,secondBin);
		  for(int ktBin=0;ktBin<numKtBins;ktBin++)
		    {
		      locCount+=counts[bt][chargeBin][firstBin][secondBin][ktBin];
		      plotResults[resIdx].kTValues[ktBin]=counts[bt][chargeBin][firstBin][secondBin][ktBin];
		      plotResults[resIdx].kTMeans[ktBin]=meanValues_kT[bt][chargeBin][firstBin][secondBin][ktBin]/counts[bt][chargeBin][firstBin][secondBin][ktBin];
		      if(bt==binType_labTheta_z && chargeBin==0)
			{
			  cout <<"kTbin: " <<ktBin << " kt count: "<< plotResults[resIdx].kTValues[ktBin] << " mean kt: "<< plotResults[resIdx].kTMeans[ktBin]<<endl;
			}
		    }
		  plotResults[resIdx].meanKinBin1=meanValues_kin1[bt][chargeBin][firstBin][secondBin]/locCount;
		  plotResults[resIdx].meanKinBin1=meanValues_kin2[bt][chargeBin][firstBin][secondBin]/locCount;
		  plotResults[resIdx].firstKinBin=firstBin;
		  plotResults[resIdx].secondKinBin=secondBin;
		  if(bt==binType_labTheta_z && chargeBin==0)
		    {
		      cout <<"resIdx: " << resIdx <<" mean 1 : " << plotResults[resIdx].meanKinBin1;
		      cout <<",  mean2: " << plotResults[resIdx].meanKinBin2 <<",  firstBin: " << firstBin <<" second: " << secondBin;

		    }

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


void MultiPlotter::loadBinnings()
{



  /////  binningZ.push_back(0.2);
  binningZ.push_back(0.3);
  binningZ.push_back(0.4);
  binningZ.push_back(0.5);
  ////  binningZ.push_back(0.7);
  binningZ.push_back(0.6); //// 
  binningZ.push_back(1.5);

  binningKt.push_back(0.1);
  binningKt.push_back(0.15);
  binningKt.push_back(0.2);
   binningKt.push_back(0.25);
  binningKt.push_back(0.3);
    binningKt.push_back(0.35);
  binningKt.push_back(0.4);
    binningKt.push_back(0.45);
  binningKt.push_back(0.5);
    binningKt.push_back(0.55);
  binningKt.push_back(0.6);
    binningKt.push_back(0.65);
  binningKt.push_back(0.7);
    binningKt.push_back(0.75);
  binningKt.push_back(0.8);
  binningKt.push_back(0.9);
    binningKt.push_back(1.0);
    binningKt.push_back(1.25);
    binningKt.push_back(1.5);
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





//void MultiPlotter::getIntAsymmetry(float a[3], float ea[3],int binningType,int chargeType, bool save1D)

void MultiPlotter::savePlots( plotType mPlotType)
{

  PlotResults* m_plotResults=plotResults;
  PlotResults* loc_plotResults=0;

  rFile.cd();
  TTree *tree = new TTree("PlotTree","PlotTree");
  tree->Branch("PlotBranch","PlotResults",&loc_plotResults,32000,99);


  int numKinBin1=0;
  int numKinBin2=0;
  char buffer[200];
  char buffer1[200];
  float mX[50];
  float mY[50];
  float mXErr[50];
  float mYErr[50];
  for(int binningType=binType_labTheta_z; binningType<binType_end;binningType++)
    {
      for(int chargeType=0;chargeType<1;chargeType++)
	{

	  string binName=getBinName(binningType,chargeType,-1,-1);
	  sprintf(buffer,"%s",binName.c_str());
	  cout <<"saving graph for " << binName <<" buffer; " << buffer<<endl;
	  for(int i=0;i<maxKinMap[binningType].first;i++)
	    {
	      for(int j=0;j<maxKinMap[binningType].second;j++)
		{
		  //		  cout <<" bin: " << i << ", " << j << endl;

		  
		  double normFactor=1.0;
		  double maxVal=-1.0;
		  double maxValNorm=-1.0;
		  int resIdx=getResIdx(binningType,chargeType,i,j);
		  for(int iKtBin=0;iKtBin<numKtBins;iKtBin++)
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
	  

		  for(int iKtBin=0;iKtBin<numKtBins;iKtBin++)
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
		      int resIdx=getResIdx(binningType,chargeType,i,j);
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
  rFile.Write();
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

      int chargeBin1=hp->chargeType1[i];
      int chargeBin2=hp->chargeType2[i];
      int chargeBin=hp->chargeType[i];

      int particleBin1=hp->particleType1[i];
      int particleBin2=hp->particleType2[i];
      int particleBin=hp->particleType[i];

      //don'te care for now...

      this->z1=hp->z1[i];
      this->z2=hp->z2[i];
      this->qT=hp->qT[i];
      this->kT=hp->kT[i];

      this->labTheta1=hp->labTheta1[i];
      this->labTheta2=hp->labTheta2[i];

      thrustBin=getBin(binningThrust,event.Thrust);

      qTBin=getBin(binningQt,qT);
      kTBin=getBin(binningKt,kT);
      //      cout <<"kt: " << kT <<" bin: "<< kTBin<<endl;

      zbin1=getBin(binningZ,hp->z1[i]);
      //      cout <<"zbin1: "<< zbin1 <<endl;
      zbin2=getBin(binningZ,hp->z2[i]);

      labThetaBin1=getBin(binningLabTheta,hp->labTheta1[i]);
      labThetaBin2=getBin(binningLabTheta,hp->labTheta2[i]);
      //      cout <<"getting mass: " << hq->hp1.mass[i] <<endl;

      for(int bt=binType_labTheta_z; bt<binType_end;bt++)
	{
	  int firstBin=*(binningMap[bt].first);
	  int secondBin=*(binningMap[bt].second);
	  float firstKin=*(meanMap[bt].first);
	  float secondKin=*(meanMap[bt].second);

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

	  double weight=1.0;
	  chargeBin=0;

	  if(bt==binType_ThrustLabTheta_z)
	    {
	      //	      cout <<"filling thrusttheta bin: "<< firstBin << " zb in : "<< secondBin <<endl;
	    }
	  //	  Cout <<"bt: " << bt <<" chargeBin: " << chargeBin<< " firstBin: " << firstBin << " second: " << secondBin <<" kt: "<< kTBin <<endl;
	  counts[bt][chargeBin][firstBin][secondBin][kTBin]+=weight;
	  meanValues_kin1[bt][chargeBin][firstBin][secondBin]+=(weight*firstKin);
	  meanValues_kin2[bt][chargeBin][firstBin][secondBin]+=(weight*secondKin);
	  meanValues_kT[bt][chargeBin][firstBin][secondBin][kTBin]+=(weight*kT);
	}
    }
};

void MultiPlotter::setBinningMap()
{

  for(int bt=binType_labTheta_z; bt<binType_end;bt++)
    {
      switch(bt)
	{
	case binType_labTheta_z:
	  binningMap.push_back(pair<int*, int* >(&(this->labThetaBin1), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->labTheta1),&(this->z1)));
	  maxKinMap.push_back(pair<int,int>(binningLabTheta.size(),binningZ.size()));
	  break;
	case binType_ThrustLabTheta_z:
	  binningMap.push_back(pair<int*, int* >(&(this->thrustLabThetaBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->thrustLabTheta),&(this->z1)));
	  maxKinMap.push_back(pair<int,int>(binningThrustLabTheta.size(),binningZ.size()));
	  break;
	case binType_z_z:
	  binningMap.push_back(pair<int*, int* >(&(this->zbin1), &(this->zbin2)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->z2)));
	  maxKinMap.push_back(pair<int,int>(binningZ.size(),binningZ.size()));
	  break;


	case binType_zOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->zbin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->z1),&(this->z1) ));
	  maxKinMap.push_back(pair<int,int>(1,binningZ.size()));
	  break;


	case binType_labThetaOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->labThetaBin1)));
	  meanMap.push_back(pair<float*, float*>(&(this->labTheta1),&(this->labTheta1) ));
	  maxKinMap.push_back(pair<int,int>(1,binningLabTheta.size()));
	  break;

	case binType_qTOnly:
	  binningMap.push_back(pair<int*, int* > (  &(this->zeroBin), &(this->qTBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->qT),&(this->qT) ));
	  maxKinMap.push_back(pair<int,int>(1,binningLabTheta.size()));
	  break;


	case binType_ThrustOnly:
	  binningMap.push_back(pair<int*,int* > (&(this->zeroBin),&(this->thrustBin)));
	  meanMap.push_back(pair<float*, float*>(&(this->thrust),&(this->thrust) ));
	  maxKinMap.push_back(pair<int,int>(1,binningThrust.size()));
	  break;

	default:
	  cout <<"binning not recognized!!"<<endl;
	  exit(0);
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
string MultiPlotter::getBinName(int binningType,int chargeType, int firstBin, int secondBin)
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

  switch(chargeType)
    {
    case pairChargeInt:
      ret+="ChargeInt_";
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
