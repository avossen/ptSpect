#include <iostream>
#include <fstream>
#include <cstdlib>
#include "MEvent.h"
#include "StyleSetter.h"
#include "HadronPairArray.h"
#include "MultiPlotter.h"
#include "TwoHadAsymsCommons.h"


//#define MAX_EVENTS 100000


using namespace std;

int main(int argc, char** argv)
{
  cout <<" Computing asymmetries... with argument: " << argv[1] <<endl;
  //  TFile tstFile;
  //  TTree myTree;

  //use this for all mc so we have WoA data

  Bool_t mcData=false;
  //  kMCFlags isMC=mcFlagNone;
  //  kMCFlags isMC=mcFlagWoA;
  char* rootPath=argv[1];
  char* basePath=argv[2];
  char* argMCFlag=argv[3];
  kMCFlags isMC=mcFlagNone;


  if(argc==4 && string(argMCFlag).find("mc")!=string::npos)
    {
      cout <<"using mc info" <<endl;
      isMC=mcAsData;
    }

  srand(time(NULL));
  cout <<"Root path is: " << rootPath <<endl;
  cout <<"base path is: "<< basePath <<endl;
  string sRootPath(rootPath);
  if(sRootPath.find_last_of('/')==sRootPath.length()-1)
    {
      sRootPath.erase(sRootPath.find_last_of('/'));
    }
  size_t found=sRootPath.find_last_of('/');
  string folderName=sRootPath.substr(found+1);
  bool onResonance=false;
  bool isUds=false;
  bool isCharm=false;

  cout <<"folder Name: "<< folderName <<endl;

  if(folderName.find("on_resonance")!=string::npos)
    onResonance=true;
  if(folderName.find("MC")!=string::npos)
    mcData=true;
  if(folderName.find("uds")!=string::npos)
    isUds=true;
  if(folderName.find("charm")!=string::npos)
    isCharm=true;

  char buffer[100];
  if(isCharm)
    sprintf(buffer,"invMass_charm");
  else
    sprintf(buffer,"invMass_uds");


  if(isCharm)
    sprintf(buffer,"z1_charm");
  else
    sprintf(buffer,"z1_uds");
  TH1D z1(buffer,buffer,500,0,1);
  if(isCharm)
    sprintf(buffer,"z2_charm");
  else
    sprintf(buffer,"z2_uds");
  TH1D z2(buffer,buffer,500,0,1);

  //for the data..
  size_t numPos=folderName.find("ex");
  int expNumber=-1;
  if(numPos!=string::npos)
    {
      char tmpNum[3];
      tmpNum[0]=folderName[numPos+2];
      cout <<"first num: " << tmpNum[0] << ", pos: " << numPos <<endl;
      tmpNum[1]=folderName[numPos+3];
      cout <<"sec num: " << tmpNum[1] << ", pos: " << numPos <<endl;
      tmpNum[2]='\n';
      expNumber=atoi(tmpNum);
    }
  cout <<"experiment : " << expNumber <<endl;
  //sometimes leads to crash?? (in setting marker size... strange...)
  //  setStyleOpts();
  TChain* chAll;
  TChain* chWoA=0;

  if(mcFlagWoA==isMC)
    chAll=new TChain("GenTree");
  else
    {
      cout <<"define datatree" <<endl;
      chAll=new TChain("DataTree");

      if(isMC==mcAsData)
	chWoA=new TChain("GenTree");
    }
  
  if(mcData)
    {
      chAll->Add((string(rootPath)+"/*.root").c_str());
      //      chAll->Add((string(rootPath)+"/*.root").c_str());
    }
  else
    {
      //      cout <<" real data.. adding: " << (string(rootPath)+"/*.root")<<endl;
      chAll->Add((string(rootPath)+"/*.root").c_str());
    }
  if(chWoA)
    {
      cout <<"adding woa files.."<<endl;
      //we have it all in one folder now...
      chWoA->Add((string(rootPath)+"/*.root").c_str());
      //      chWoA->Add((string(rootPath)+"/uds/*.root").c_str());
      //      chWoA->Add((string(rootPath)+"/charm/*.root").c_str());
    }
  //  cout <<" had quad  " <<endl;

  //only one class is actually reading the tree
  kMCFlags dataMCFlag=mcFlagNone;
  //    kMCFlags dataMCFlag=mcFlagNone;
  if(isMC==mcFlagMC)
    dataMCFlag=mcFlagMC;


  cout <<"dataMCFlag: "<< dataMCFlag <<endl;
  if(dataMCFlag==mcFlagMC)
    cout <<"mc Flag!! " <<endl;
  HadronPairArray hadPair(chAll,dataMCFlag);
  hadPair.zOrdered=true;

  kMCFlags hadMCFlag=dataMCFlag;
  if(isMC==mcAsData)
    {
      cout <<" mcAsData" <<endl;
      hadMCFlag=mcFlagMC;
    }

  HadronPairArray* hadPairMC;

  if(isMC!=mcFlagNone)
    {
      hadPairMC=new HadronPairArray(chAll,hadMCFlag);
      hadPairMC->followFlip=true;
      hadPairMC->relatedHP=&hadPair;
      hadPairMC->zOrdered=false;
    }

  cout <<"done with had pairs....." <<endl;
  kMCFlags mEventDataMCFlag=dataMCFlag;
  if(isMC==mcAsData)
    {
      mEventDataMCFlag=mcAsData;
    }
  cout <<"mevent mc flag: " <<mEventDataMCFlag << " data flag: "<< dataMCFlag <<endl;
  MEvent myEvent(chAll,mEventDataMCFlag);
  //  MEvent myEventMC(chAll,mcFlagMC);
  cout <<" 1 " << endl;


  cout <<" 2 " << endl;
  MEvent* pMyEventWoA;
  HadronPairArray* pHadPairWoA;
  if(chWoA) 
    {
      pMyEventWoA=new MEvent(chWoA,mcFlagWoA);
      pHadPairWoA=new HadronPairArray(chWoA,mcFlagWoA);
      pHadPairWoA->zOrdered=true;
    }


  //  MEvent myEventWoA(chWoA,mcFlagWoA);

  cout << "done event " << endl;
  cout <<"how many? "<<endl;
  Int_t nevents=chAll->GetEntries();
  cout <<"we have " << nevents <<endl;
  //NUM_PHI_BINS angular bins
  stringstream ss;
  ss <<"_ex"<<expNumber;
  if(onResonance)
    ss<<"_onRes_";
  else
    ss<<"_continuum_";
  if(isUds)
    ss <<"_uds_";
  if(isCharm)
    ss <<"_charm_";
  if(mcData && !isUds && !isCharm)
    ss <<"_mcAll_";
  else 
    ss<<"_data_";

  bool m_useQt=false;
#ifdef USE_QT
  m_useQt=true;
  cout<<"setting m_useQt" <<endl;
#endif

  cout <<"m_useqt: " << m_useQt <<endl;
  MultiPlotter smearingPlotter(m_useQt,const_cast<char*>(basePath),const_cast<char*>("smearingPlotter"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData );
  //this one is essentially to save xini for the raw pythia
  MultiPlotter smearingPlotterRaw(m_useQt,const_cast<char*>(basePath),const_cast<char*>("smearingPlotterRaw"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData );
  MultiPlotter plotter(m_useQt,const_cast<char*>(basePath),const_cast<char*>("Normal"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData);
  //  MultiPlotter plotterMC(const_cast<char*>("NormalMC"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData);
  MultiPlotter plotterWoA(m_useQt,const_cast<char*>(basePath),const_cast<char*>("NormalWoA"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData);
  //  MultiPlotter fitPi0SigMinusMix(const_cast<char*>("fitPi0SigMinusMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);
  //  MultiPlotter fitPi0BgMinusMix(const_cast<char*>("fitPi0BgMinusMix"),ss.str(),expNumber,onResonance,isUds,isCharm,mcData,NUM_PHI_BINS);

#ifdef USE_QT
  smearingPlotter.useQt=true;
  smearingPlotterRaw.useQt=true;
  plotter.useQt=true;
  plotterWoA.useQt=true;
#endif

  plotter.setName("Normal");
  //  plotterMC.setName("NormalMC");
  plotterWoA.setName("NormalWoA");
  if(chWoA)
    {
      Int_t neventsWoA=chWoA->GetEntries();
      cout <<"we have " << neventsWoA <<" WoA events" <<endl;

      for(long i=0;i<neventsWoA;i++)
	{
#ifdef MAX_EVENTS
	  if(i>MAX_EVENTS)
	    break;
#endif
	  if(!(i%10000))
	    cout <<"processing woa event nr " << i << " of " << neventsWoA << "(" << 100*i/(float)neventsWoA<< "% )"<<endl;

	  chWoA->GetEntry(i);
	  pMyEventWoA->afterFill();
	  
	  if(pMyEventWoA->cutEvent)
	    {
	      continue;
	    }
	  //  cout <<"woa after fill " <<endl;
	  pHadPairWoA->afterFill();
	  //  cout <<"done " <<endl;
	  //	  	  cout <<"adding woa had quad to plotter... " <<endl;
	  plotterWoA.addHadPairArray(pHadPairWoA, *pMyEventWoA);
	  //	  cout <<"runNumber: " << pMyEventWoA->runNr <<" eventNumber: "<< pMyEventWoA->evtNr <<endl;
	  smearingPlotterRaw.evtNr=pMyEventWoA->evtNr;
	  smearingPlotterRaw.addXiniEntry(pHadPairWoA);

	}
    }
  if(chWoA)
    plotterWoA.doPlots();

  for(long i=0;i<nevents;i++)
    {
#ifdef MAX_EVENTS
      if(i>MAX_EVENTS)
	break;
#endif
      if(!(i%10000))
	cout <<"processing event nr " << i << " of " << nevents << "(" << 100*i/(float)nevents<< "% )"<<endl;
      chAll->GetEntry(i);

      myEvent.afterFill();
      if(myEvent.cutEvent)
	{
	  continue;
	}

      //      cout <<"normal quad after fill" <<endl;
      //           cout<<" data after fill " <<endl;
      //      cout <<"had pair data " <<endl;


      //      cout <<"reconstructed pair: " <<endl;
      hadPair.afterFill(myEvent.evtNr,false);
      //           cout <<"mc after fill " <<endl;
      if(isMC!=mcFlagNone)
	{
	  //	  	    cout <<"filling mc pair " <<endl;
	  //	  	  cout <<endl<<"generated pair: " << endl;
	  hadPairMC->afterFill(myEvent.evtNr,false);
	  //	    cout <<"done " <<endl;
	}
      //      cout <<endl<<endl;
      for(int i=0;i<hadPair.numPairs;i++)
	{
	  z1.Fill(hadPair.z1[i]);
	  z2.Fill(hadPair.z2[i]);
	}

      plotter.addHadPairArray(&hadPair, myEvent);
      if(isMC!=mcFlagNone)
	{

	  //for x-check with Charlotte
	  //	  cout <<"runNumber: " << myEvent.runNr <<" eventNumber: "<< myEvent.evtNr <<endl;
	  smearingPlotter.evtNr=myEvent.evtNr;
	  smearingPlotter.addSmearingEntry(&hadPair,hadPairMC);
	}
    }
  if(isMC!=mcFlagNone)
    {
      smearingPlotterRaw.saveXini();
      smearingPlotter.saveSmearingMatrix();
    }


  if(isCharm)
    sprintf(buffer,"z1_charm.root");
  else
    sprintf(buffer,"z1_uds.root");
  z1.SaveAs(buffer);
  if(isCharm)
    sprintf(buffer,"z2_charm.root");
  else
    sprintf(buffer,"z2_uds.root");
  z2.SaveAs(buffer);
  if(isCharm)
    sprintf(buffer,"invMass_charm");
  else
    sprintf(buffer,"invMass_uds");

  plotter.doPlots();

  cout <<" printing debug " <<endl;
  plotter.printDebug(plotType_2D);
  cout <<"done " <<endl;
  plotter.savePlots(plotType_2D);
  cout <<"now woa..." <<endl;
  if(chWoA)
    plotterWoA.savePlots(plotType_2D);

  //  plotterWoA.savePlot(0,0,);

  vector<MultiPlotter*> myPlotterArray;
  myPlotterArray.push_back(&plotter);

  if(chWoA)
    {
      myPlotterArray.push_back(&plotterWoA);

    }
     
  for(vector<MultiPlotter*>::iterator it=myPlotterArray.begin();it!=myPlotterArray.end();it++)
    {
      //      (*it)->savePlots();   //was save asymmetries...

      //      (*it)->savePlot(binType_z_z,quadPN);


      //      for(plotType pt=plotType_1D;pt<plotType_end;pt=(plotType)((int)pt+1))
	{
	  //	  (*it)->savePlot(binType_m_m,quadPN,pt);
	  //	  (*it)->savePlot(binType_z_z,quadPN,pt);
	
	}
    }


};

//  LocalWords:  endl
