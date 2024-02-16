#include "TRandom3.h"
#include "TF1.h"
#include "TFitResult.h"
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
#include "PlotResHelper.h"
//#define MAX_EVENTS 100

using namespace std;

void loadPIDMatrix();
void initPIDMatrix();
#define  pb 17   //(number of momentum bins)
#define  thb 9 //  (number of theta bins)

float ****pidMatrixPositive;
float ****pidMatrixNegative;

float ****pidMatrixPositive2;
float ****pidMatrixNegative2;

float ****pidUncertPositive;
float ****pidUncertNegative;

float ****pidUncertPositive2;
float ****pidUncertNegative2;


//from AnaDef
static int Pos=0;
static int Neg=1;

static int Pion=0;
static int Kaon=1;
static int Proton=2;


//this is the belle convention that is followed in the matrices. However, in the hadron pair, the pion starts at 0, since it is an enum that starts at 0
#define pionIdx 2
#define kaonIdx 3
#define protonIdx 4




float ****pidMatrixPositiveTrial;
float ****pidMatrixNegativeTrial;

float ****pidMatrixPositive2Trial;
float ****pidMatrixNegative2Trial;


float getPProduct(float p1_1, float p1_2, float p2_1, float p2_2)
{
  return (p1_1+p2_1)/2.0 * (p1_2+p2_2)/2.0;
}



static const int numKinematicBinning=2;

    static const int numMomBins=pb+1;
  static const int numThetaBins=thb+1;
  static const int numPIDs=10;
//  static const int numPIDs=3;
  static const int numSinglePIDs=5;
static const int numChargeBins=3;

static int numKtBins=9;
static int numZBins=11; //max

static int numIterations=200;

int main(int argc, char** argv)
{

  //or initialized with   time(NULL));
  TRandom3 randInst(0);
  char buffer[300];
  char* pidFile=argv[1];
  initPIDMatrix();
  loadPIDMatrix();

  double****** counts;
  int***** histoCheck=allocateArray<int>(numPIDs,numChargeBins,numZBins,numZBins,numKtBins);
  TH1D****** histos=allocHistos1DPointers(numPIDs, numChargeBins, numZBins, numZBins, numKtBins);
  //  TH1D****** histos=
  TH2D****** covMatrices=allocHistos(numPIDs,numChargeBins, numZBins,numZBins, numKtBins);
  
    //  pseudocode
  
  //init histogram for each pid, charge, z1, z2 and kt bin
  //covariance matrix for each z1,z2,kt (not across, pid, charge, since the matrix uncertainties are independent for those (probably not true, but that is what we have)
  
  //load tree
  TFile file(pidFile);
  TTree* tree=(TTree*)file.Get("pidSysTree");
  //  cout <<" makeing tree.. " << endl;
  //  TTree *tree = new TTree("PlotTree","PlotTree");
  //  cout <<"about to branch" <<tree <<endl;
  //  tree->Branch("PlotBranch","PlotResults",&loc_plotResults,32000,99);

  int zBin1;
  int zBin2;
  int ptBin;
  int pBin1;
  int pBin2;
  int thetaBin1;
  int thetaBin2;
  int chargeBin;
  int pidBin;
  int idAsBin1;
  int idAsBin2;
  int charge1;
  int charge2;
  int pid1;
  int pid2;
  float ptVal;
  
  tree->SetBranchAddress("zBin1",&zBin1);
  tree->SetBranchAddress("zBin2",&zBin2);
  tree->SetBranchAddress("ptBin",&ptBin);
  tree->SetBranchAddress("ptVal",&ptVal);
  tree->SetBranchAddress("pBin1",&pBin1);
  tree->SetBranchAddress("pBin2",&pBin2);
  tree->SetBranchAddress("thetaBin1",&thetaBin1);
  tree->SetBranchAddress("thetaBin2",&thetaBin2);
  tree->SetBranchAddress("chargeBin",&chargeBin);
  tree->SetBranchAddress("pidBin",&pidBin);
  tree->SetBranchAddress("idAsBin1",&idAsBin1);
  tree->SetBranchAddress("idAsBin2",&idAsBin2);
  tree->SetBranchAddress("charge1",&charge1);
  tree->SetBranchAddress("charge2",&charge2);
  tree->SetBranchAddress("pid1",&pid1);
  tree->SetBranchAddress("pid2",&pid2);
  
  srand(time(NULL));
  int randInt=rand()%200;

  float rndNum=1.0-randInt/100.0;

  //sample random number from gauss
  
  //for each pid bin,
  //for each charge bin
  //for each z1, z2
  //for each kT
  //get weight + uncertainty*gaussValue
  //put in histo
  //add to covariance matrix

  //done done done

  //normalize gausians, correlation matrix
  //get sigma etc..

  Long64_t nevents=tree->GetEntries();
  cout <<" we have " << nevents <<" entries " << endl;
  for(int i=0;i<numIterations;i++)
    {

    for(int iP=0;iP<numMomBins;iP++)
      {
      	for(int iT=0;iT<numThetaBins;iT++)
	  {
	  	    for(int iPID1=0;iPID1<numSinglePIDs;iPID1++)
		      {
		      	    for(int iPID2=0;iPID2<numSinglePIDs;iPID2++)
			      {
				randInt=rand()%200;
				rndNum=1.0-randInt/100.0;
				////				cout <<"reading matrix for iP : " << iP << " iT: " << iT <<" pid1: " << iPID1 <<" pid2: "<< iPID2 <<endl;
				//				pidMatrixPositiveTrial[iP][iT][iPID1][iPID2]=pidMatrixPositive[iP][iT][iPID1][iPID2]+rndNum*pidUncertPositive[iP][iT][iPID1][iPID2];
				////				cout <<" pidMatrixPos before:  "<< pidMatrixPositive[iP][iT][iPID1][iPID2] <<" uncert: " << pidUncertPositive[iP][iT][iPID1][iPID2] <<" randInst: "<< pidMatrixPositiveTrial[iP][iT][iPID1][iPID2] <<endl;
				pidMatrixPositiveTrial[iP][iT][iPID1][iPID2]=randInst.Gaus(pidMatrixPositive[iP][iT][iPID1][iPID2],pidUncertPositive[iP][iT][iPID1][iPID2]);
				///				cout <<" pidMatrixPos: "<< pidMatrixPositive[iP][iT][iPID1][iPID2] <<" uncert: " << pidUncertPositive[iP][iT][iPID1][iPID2] <<" randInst: "<< pidMatrixPositiveTrial[iP][iT][iPID1][iPID2] <<endl;
				/////				cout <<"rndNum " << rndNum <<" pidMatrix before: " << pidMatrixPositive[iP][iT][iPID1][iPID2] <<" uncert: " << pidUncertPositive[iP][iT][iPID1][iPID2] <<" trial matrix: "<< pidMatrixPositiveTrial[iP][iT][iPID1][iPID2] <<endl;
				randInt=rand()%200;
				rndNum=1.0-randInt/100.0;
				///				pidMatrixPositive2Trial[iP][iT][iPID1][iPID2]=pidMatrixPositive[iP][iT][iPID1][iPID2]+rndNum*pidUncertPositive[iP][iT][iPID1][iPID2];
				pidMatrixPositive2Trial[iP][iT][iPID1][iPID2]=randInst.Gaus(pidMatrixPositive2[iP][iT][iPID1][iPID2],pidUncertPositive2[iP][iT][iPID1][iPID2]);

				
				randInt=rand()%200;
				rndNum=1.0-randInt/100.0;
				///				pidMatrixNegativeTrial[iP][iT][iPID1][iPID2]=pidMatrixNegative[iP][iT][iPID1][iPID2]+rndNum*pidUncertNegative[iP][iT][iPID1][iPID2];
				pidMatrixNegativeTrial[iP][iT][iPID1][iPID2]=randInst.Gaus(pidMatrixNegative[iP][iT][iPID1][iPID2],pidUncertNegative[iP][iT][iPID1][iPID2]);
				randInt=rand()%200;
				rndNum=1.0-randInt/100.0;
				////				pidMatrixNegative2Trial[iP][iT][iPID1][iPID2]=pidMatrixNegative[iP][iT][iPID1][iPID2]+rndNum*pidUncertNegative[iP][iT][iPID1][iPID2];
				pidMatrixNegative2Trial[iP][iT][iPID1][iPID2]=randInst.Gaus(pidMatrixNegative2[iP][iT][iPID1][iPID2],pidUncertNegative2[iP][iT][iPID1][iPID2]);
				
			      }
		      }
	  }
      }

    
 


    
      //run over the tree for each iteration
      //allocateArray(all bins);
      counts=allocateArray<double>(numKinematicBinning,numPIDs,numChargeBins,numZBins,numZBins,numKtBins);
    
      //fill with the weights from the tree
      for(Long64_t i=0;i<nevents;i++)
	{
	  if(i> nevents/3)
	    break;
	  tree->GetEntry(i);
	  //positive or negative?
	  cout <<"lookign at pid bin: " << pidBin <<endl;
	  float weight1_1=0;
	  float weight1_2=0;
	  for(int hypo1=pionIdx;hypo1<=protonIdx;hypo1++)
	    {
	      if(charge1==Pos)
		{
		  //there is a shift of 2 between the pid in the hadron pair and in the matrix
		  //the pid of the hadron pair is the 'idAs'
		  weight1_1+=pidMatrixPositiveTrial[pBin1][thetaBin1][hypo1][pid1+2];
		  weight1_2+=pidMatrixPositive2Trial[pBin1][thetaBin1][hypo1][pid1+2];
		}
	      else
		{
		  weight1_1+=pidMatrixNegativeTrial[pBin1][thetaBin1][hypo1][pid1+2];
		  weight1_2+=pidMatrixNegative2Trial[pBin1][thetaBin1][hypo1][pid1+2];
		}
	    }
	  float weight2_1=0;
	  float weight2_2=0;
	  for(int hypo2=pionIdx;hypo2<=protonIdx;hypo2++)
	    {
	      if(charge2==Pos)
		{
		  //there is a shift of 2 between the pid in the hadron pair and in the matrix
		  weight2_1=+pidMatrixPositiveTrial[pBin2][thetaBin2][hypo2][pid2+2];
		  weight2_2=+pidMatrixPositive2Trial[pBin2][thetaBin2][hypo2][pid2+2];
		}
	      else
		{
		  weight2_1=+pidMatrixNegativeTrial[pBin2][thetaBin2][hypo2][pid2+2];
		  weight2_2=+pidMatrixNegative2Trial[pBin2][thetaBin2][hypo2][pid2+2];
		}
	    }
	  counts[0][pidBin][chargeBin][zBin1][zBin2][ptBin]+=getPProduct(weight1_1,weight1_2,weight2_1,weight2_2);
	}
      //insert into histogram
      

    
      for(int pid=0;pid<numPIDs;pid++)
	{
	  for(int chargeBin=0;chargeBin<numChargeBins;chargeBin++)
	    {
	      for(int z1Bin=0;z1Bin<numZBins;z1Bin++)
		{
		  for(int z2Bin=0;z2Bin<numZBins;z2Bin++)
		    {
		      for(int ktBin=0;ktBin<numKtBins;ktBin++)
			{

			  
			  //cout <<"Filling histo with count : "<< counts[0][pid][chargeBin][z1Bin][z2Bin][ktBin] <<endl;
			  //haven't created this particular histogram yet
			  if(histoCheck[pid][chargeBin][z1Bin][z2Bin][ktBin]==0)
			    {
			      cout <<" have note already allocated histo for " << pid <<" charge: " << chargeBin <<" z1Bin: " << z1Bin << " z2Bin: " << z2Bin << " ktBin: " << ktBin <<endl;
			      histoCheck[pid][chargeBin][z1Bin][z2Bin][ktBin]=1;
			      float locCounts=counts[0][pid][chargeBin][z1Bin][z2Bin][ktBin];
			      int upperBin=1.5*locCounts;
			      int lowerBin=0.5*locCounts;
			      cout <<" locCounts: " << locCounts<< "upper bin: " << upperBin << " lowerBin: " << lowerBin <<endl;
			      sprintf(buffer,"histo_pid%d_chargeBin%d_z1Bin%d_z2Bin%d_ktBin%d",pid,chargeBin,z1Bin,z2Bin,ktBin);
			      cout <<"creating histo " << buffer <<endl;
			      histos[pid][chargeBin][z1Bin][z2Bin][ktBin]=new TH1D(buffer,buffer,100,lowerBin, upperBin);
			      cout <<" done..." <<endl;
			    }
			  else
			    {
			      //			      cout <<" should have already allocated histo for " << pid <<" charge: " << chargeBin <<" z1Bin: " << z1Bin << " z2Bin: " << z2Bin << " ktBin: " << ktBin <<endl;
			    }

			    histos[pid][chargeBin][z1Bin][z2Bin][ktBin]->Fill(counts[0][pid][chargeBin][z1Bin][z2Bin][ktBin]);
			  
			  //for each of the entries, but the count into the histogram (get one entry per iteration)
			}
		    }
		}
	    }
	}
    }

  TFile outFile("pidSysOut.root","RECREATE");
    for(int pid=0;pid<numPIDs;pid++)
	{
	  for(int chargeBin=0;chargeBin<numChargeBins;chargeBin++)
	    {
	      for(int z1Bin=0;z1Bin<numZBins;z1Bin++)
		{
		  for(int z2Bin=0;z2Bin<numZBins;z2Bin++)
		    {
		      for(int ktBin=0;ktBin<numKtBins;ktBin++)
			{
			  cout <<" looking at pid: " << pid <<"  chargeBin: " << chargeBin <<" z1Bin; " << z1Bin << " z2Bin: " << z2Bin <<" ktBin: "<< ktBin <<endl;
		  //fit with gaussian

			 histos[pid][chargeBin][z1Bin][z2Bin][ktBin]->Fit("gaus");
			 cout <<" done with fit " << endl;
			  //parameters should be amplitude, mean and sigma
			  TF1 *fitFunc = histos[pid][chargeBin][z1Bin][z2Bin][ktBin]->GetFunction("gaus");

			  cout <<" get fit func " <<endl;
			  double meanCountsFromFit=0;
			  if(fitFunc)
			    {
			      meanCountsFromFit=fitFunc->GetParameter(1);
			    }
			  else
			    {
			      cout <<" fit failed ... " << endl;
			      cout <<" histo has: " << histos[pid][chargeBin][z1Bin][z2Bin][ktBin]->GetEntries() << " entries " << endl;
			      continue;
			    }
			  if(meanCountsFromFit==0)
			    {
			      cout <<" mean counts is 0" <<endl;
			      continue;
			    }
			  double sigma= fitFunc->GetParameter(2);
			  double sigmaErr=fitFunc->GetParError(2);
			  cout <<"extract sigma: " << sigma <<" err: "<< sigmaErr << " peak position: " << meanCountsFromFit <<  " rel sigma: " << sigma/meanCountsFromFit << endl;

			  histos[pid][chargeBin][z1Bin][z2Bin][ktBin]->SetDirectory(&outFile);
			  histos[pid][chargeBin][z1Bin][z2Bin][ktBin]->Write();
			}
		    }
		}
	    }
	}
    outFile.Write();
  
//now get the width etc
  
}

void initPIDMatrix()
{

    pidMatrixPositive=new float***[numMomBins];
    pidMatrixNegative=new float***[numMomBins];
    pidMatrixPositive2=new float***[numMomBins];
    pidMatrixNegative2=new float***[numMomBins];
    pidMatrixPositiveTrial=new float***[numMomBins];
    pidMatrixNegativeTrial=new float***[numMomBins];
    pidMatrixPositive2Trial=new float***[numMomBins];
    pidMatrixNegative2Trial=new float***[numMomBins];

    pidUncertPositive=new float***[numMomBins];
    pidUncertNegative=new float***[numMomBins];
    pidUncertPositive2=new float***[numMomBins];
    pidUncertNegative2=new float***[numMomBins];

    for(int i=0;i<numMomBins;i++)
      {
	pidMatrixPositive[i]=new float**[numThetaBins];
	pidMatrixNegative[i]=new float**[numThetaBins];
	pidMatrixPositive2[i]=new float**[numThetaBins];
	pidMatrixNegative2[i]=new float**[numThetaBins];

	pidMatrixPositiveTrial[i]=new float**[numThetaBins];
	pidMatrixNegativeTrial[i]=new float**[numThetaBins];
	pidMatrixPositive2Trial[i]=new float**[numThetaBins];
	pidMatrixNegative2Trial[i]=new float**[numThetaBins];

	
	pidUncertPositive[i]=new float**[numThetaBins];
	pidUncertNegative[i]=new float**[numThetaBins];
	pidUncertPositive2[i]=new float**[numThetaBins];
	pidUncertNegative2[i]=new float**[numThetaBins];


	for(int j=0;j<numThetaBins;j++)
	  {
	    pidMatrixPositive[i][j]=new float*[numPIDs];
	    pidMatrixNegative[i][j]=new float*[numPIDs];
	    pidMatrixPositive2[i][j]=new float*[numPIDs];
	    pidMatrixNegative2[i][j]=new float*[numPIDs];

	    pidMatrixPositiveTrial[i][j]=new float*[numPIDs];
	    pidMatrixNegativeTrial[i][j]=new float*[numPIDs];
	    pidMatrixPositive2Trial[i][j]=new float*[numPIDs];
	    pidMatrixNegative2Trial[i][j]=new float*[numPIDs];

	    
	    pidUncertPositive[i][j]=new float*[numPIDs];
	    pidUncertNegative[i][j]=new float*[numPIDs];
	    pidUncertPositive2[i][j]=new float*[numPIDs];
	    pidUncertNegative2[i][j]=new float*[numPIDs];

	    for(int k=0;k<numPIDs;k++)
	      {
		pidMatrixPositive[i][j][k]=new float[numPIDs];
		pidMatrixNegative[i][j][k]=new float[numPIDs];
		pidMatrixPositive2[i][j][k]=new float[numPIDs];
		pidMatrixNegative2[i][j][k]=new float[numPIDs];

		pidMatrixPositiveTrial[i][j][k]=new float[numPIDs];
		pidMatrixNegativeTrial[i][j][k]=new float[numPIDs];
		pidMatrixPositive2Trial[i][j][k]=new float[numPIDs];
		pidMatrixNegative2Trial[i][j][k]=new float[numPIDs];


		
		pidUncertPositive[i][j][k]=new float[numPIDs];
		pidUncertNegative[i][j][k]=new float[numPIDs];
		pidUncertPositive2[i][j][k]=new float[numPIDs];
		pidUncertNegative2[i][j][k]=new float[numPIDs];


		for(int l=0;l<numPIDs;l++)
		  {
		    //this will be overwritten by Martin's matrices. So set this to zero so we realize that we loaded incorrect values..
		    if(k==l)
		      {
			pidMatrixPositive[i][j][k][l]=0.0;
			pidMatrixNegative[i][j][k][l]=0.0;
			pidMatrixPositive2[i][j][k][l]=0.0;
			pidMatrixNegative2[i][j][k][l]=0.0;

			pidMatrixPositiveTrial[i][j][k][l]=0.0;
			pidMatrixNegativeTrial[i][j][k][l]=0.0;
			pidMatrixPositive2Trial[i][j][k][l]=0.0;
			pidMatrixNegative2Trial[i][j][k][l]=0.0;


		      }
		  }
	      }
	  }
      }

}
  void  loadPIDMatrix()
  {
    //   TFile* fpid = new TFile("newpid.root","read");

    TFile* fpid = new TFile("~/Documents/myProjects/belle/ptSpect/invertedpidmatrices_setb061810I_inv2_realdataalways.root","read");
    TFile* fpid2 = new TFile("~/Documents/myProjects/belle/ptSpect/invertedpidmatrices_setb061810I_inv1_MConlyatlooseends.root","read");

    char matrix_name[300];
    char uncert_minus_name[300];
    char uncert_plus_name[300];
    if (fpid->IsZombie() || fpid2->IsZombie()) {
      printf("File code.root does not exist.\n");
      return;
    }
  
    // reading all matrices together, as the code crashes if I try to open the file too many times...
    for (Int_t u = 0; u < pb; u++)
      for (Int_t v = 0; v < thb; v++)
	for (Int_t w = 0; w < 2; w++){
	  //w is charge ( 0  negative, 1 positive);
	  sprintf(matrix_name,"invanalyticmatrix_u%d_v%d_w%d",u,v,w); // invanalyticmatrix_u16_v0_w0 does not exists!!! 
	  sprintf(uncert_plus_name,"invuncertrejplmatrix_u%d_v%d_w%d",u,v,w); 
	  sprintf(uncert_minus_name,"invuncertrejmimatrix_u%d_v%d_w%d",u,v,w); 
	  //u16 starts from v3 (invanalyticmatrix_u16_v3_w0)
	  // corresponds to z > 1 (smearing?) additional uncertainties would be 1-2
	  // up to 9 July 2015! int k= u*8*2+v*2+w; using 8 instead of 9!
	  //    int k= u*thb*2+v*2+w;
	  if(u<16 || v>2 ){
	    //       cout <<" u:  " << u << " v: " << v << " w: " << w <<endl;
	    TMatrixD mat = *(TMatrixD*)fpid->Get(matrix_name);
	    TMatrixD mat2 = *(TMatrixD*)fpid2->Get(matrix_name);
	    TMatrixD matUncertPos = *(TMatrixD*)fpid->Get(uncert_plus_name);
	    TMatrixD matUncertNeg = *(TMatrixD*)fpid->Get(uncert_minus_name);

	    TMatrixD matUncertPos2 = *(TMatrixD*)fpid2->Get(uncert_plus_name);
	    TMatrixD matUncertNeg2 = *(TMatrixD*)fpid2->Get(uncert_minus_name);

	    for (Int_t i = 0; i <= 4; i++)
	      for (Int_t j = 0; j <= 4; j++){

		//	   cout <<"mat: " << mat(i,j) <<" mat2: "<< mat2(i,j) <<endl;
		//	   cout <<"pos uncert: " << matUncertPos(i,j) <<" neg: "<< matUncertNeg(i,j) <<endl;
		//	   	   cout <<"pos uncert2: " << matUncertPos2(i,j) <<" neg: "<< matUncertNeg2(i,j) <<endl;
		float symUncert=matUncertPos(i,j);
		float symUncert2=matUncertPos2(i,j);
		if(isnan(matUncertPos(i,j)))
		  {
		    cout <<" matuncert pos nan"<<endl;
		    //to force using the neg uncert
		    symUncert=0;
		  }
		if(isnan(matUncertNeg(i,j)))
		  {
		    cout <<" matuncert neg nan"<<endl;
		  }
		else
		  {
		    if(fabs(matUncertNeg(i,j))>symUncert)
		      symUncert=fabs(matUncertNeg(i,j));
		  }
		if(isnan(matUncertNeg2(i,j)))
		  {
		    cout <<" matuncert neg2 nan"<<endl;
		  }
		else
		  {
		    if(fabs(matUncertNeg2(i,j))>symUncert2)
		      symUncert2=fabs(matUncertNeg2(i,j));
		  }

		//not sure if that makes sense here since we do the averaging later on again
		//		if(fabs(symUncert2)<100.0)
		//		  symUncert=sqrt((symUncert*symUncert+symUncert2*symUncert2))/2;
		if(isnan(symUncert))
		  cout<<"resulting uncertainty still nan!!!!!!!!"<<endl<<endl;


		//now we got the symmetric uncertainties

		//at this point symUncert is the data uncertainty if there is a problem with MC, otherwise the geometric mean
		//

		if(w==1)
		  {
#ifdef noPID
		    if(i==j)
		      {
			pidMatrixPositive[u][v][i][j]=1;
			pidMatrixPositive2[u][v][i][j]=1;
			pidUncertPositive[u][v][i][j]=0.0;
			pidUncertPositive2[u][v][i][j]=0.0;
		      }
		    else
		      {
			pidMatrixPositive[u][v][i][j]=0;
			pidMatrixPositive2[u][v][i][j]=0;
			pidUncertPositive[u][v][i][j]=0.0;
			pidUncertPositive2[u][v][i][j]=0.0;
		      }
#else
		    pidMatrixPositive[u][v][i][j]=mat(i,j);
		    pidMatrixPositive2[u][v][i][j]=mat2(i,j);
		    //			    cout <<"loading " << mat(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" positive " <<endl;
		    //		    cout <<"loading 2" << mat2(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" positive " <<endl;

		    if(fabs(pidMatrixPositive2[u][v][i][j])>100.0 || fabs(symUncert2)>100.0)
		      {
			symUncert2=symUncert;
			pidMatrixPositive2[u][v][i][j]=pidMatrixPositive[u][v][i][j];
			//the uncert asseignment doesnt' matter, since we use the average between the two anyways
			//			symUncert2=symUncert;
		      }

		    if(isnan(pidMatrixPositive2[u][v][i][j]) || isnan(matUncertPos2(i,j)) || isnan(matUncertNeg2(i,j)))
		      {
			symUncert2=symUncert;
		      }
		    //		    symUncert=sqrt(symUncert*symUncert+symUncert2*symUncert2)/2;
		    pidUncertPositive[u][v][i][j]=symUncert;
		    pidUncertPositive2[u][v][i][j]=symUncert2;
		    if(isnan(symUncert))
		      {
			cout <<"poblem with pos symUncert!!!" <<endl;
		      }
		    if(symUncert>100.0)
		      {
			cout <<" uncert > 100 " <<endl;
		      }	   

#endif
		    //	     cout <<"loading positive "<< i <<", " << j << " "  << mat(i,j) <<endl;




		  }
		else
		  {
#ifdef noPID
		    if(i==j)
		      {
			pidMatrixNegative[u][v][i][j]=1;
			pidMatrixNegative2[u][v][i][j]=1;
			pidUncertNegative[u][v][i][j]=0.0;
			pidUncertNegative2[u][v][i][j]=0.0;
		      }
		    else
		      {
			pidMatrixNegative[u][v][i][j]=0;
			pidMatrixNegative2[u][v][i][j]=0;
			pidUncertNegative[u][v][i][j]=0.0;
			pidUncertNegative2[u][v][i][j]=0.0;
		      }
#else
		    //		    cout <<"loading " << mat(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" negative " <<endl;
		    //		    cout <<"symUncert: "<< symUncert <<endl;
		    //		    cout <<"loading 2 " << mat2(i,j) << " for u: "<< u <<" v: " << v << ", i : "<< i << " j: " << j <<" negative " <<endl;
		    pidMatrixNegative[u][v][i][j]=mat(i,j);
		    pidMatrixNegative2[u][v][i][j]=mat2(i,j);
		    if(fabs(pidMatrixNegative2[u][v][i][j])>100.0 || fabs(symUncert2)>100.0) 
		      {
			pidMatrixNegative2[u][v][i][j]=pidMatrixNegative[u][v][i][j];
			//			symUncert2=symUncert;
		      }
		    if(isnan(pidMatrixNegative2[u][v][i][j]) || isnan(matUncertPos2(i,j)) || isnan(matUncertNeg2(i,j)))
		      {
			pidMatrixNegative2[u][v][i][j]=pidMatrixNegative[u][v][i][j];
			symUncert2=symUncert;
		      }

		    //		    symUncert=sqrt(symUncert*symUncert+symUncert2*symUncert2)/2;

		    if(isnan(symUncert))
		      cout <<"problem with neg symUncert!!!" <<endl;
		    pidUncertNegative[u][v][i][j]=symUncert;
		    pidUncertNegative2[u][v][i][j]=symUncert2;

		  
#endif

		    //	     cout <<"loading negative " << i << ", " << j << " "  << mat(i,j) <<endl;




		  }

		//        matrix[k][i][j] = mat(i,j);
		//if(k==72) cout << i << " " << j <<" " <<  mat(i,j) << endl;
	      }
	  }
	}
    fpid->Close();
    cout << std::fixed;
    cout << setprecision(4);
    for(Int_t w=0;w<2;w++)
      {
	for (Int_t u = 0; u < pb; u++)
	  {

	    for (Int_t v = 0; v < thb; v++)
	      {
		cout <<"momentum bin " << u;
		cout <<", cosTheta bin " <<v<<endl;
		for (Int_t i = 0; i <= 4; i++)
		  {
		    for (Int_t j = 0; j <= 4; j++)
		      {
			if(w==0)
			  cout <<"data: " << pidMatrixNegative[u][v][i][j] <<" mc: " << pidMatrixNegative2[u][v][i][j] <<" ";
			else
			  cout <<"data: " <<pidMatrixPositive[u][v][i][j]<< " mc: "<< pidMatrixPositive2[u][v][i][j] <<" ";
			//			pidUncertNegative[u][v][i][j];
		      }
		    cout <<endl;
		  }
		cout <<"symmetrized uncertainties  "<<endl;
		for (Int_t i = 0; i <= 4; i++)
		  {
		    for (Int_t j = 0; j <= 4; j++)
		      {
			if(w==0)
			  cout <<"data: " << pidUncertNegative[u][v][i][j] <<" mc: "<< pidUncertNegative2[u][v][i][j]<< " " ;
			else
			  cout <<"data: " << pidUncertPositive[u][v][i][j] <<" mc: "<< pidUncertPositive2[u][v][i][j]<<" ";
			//			pidUncertNegative[u][v][i][j];
		      }
		    cout <<endl;
		  }

		cout <<endl;
	      }
	  }
      }

  }


