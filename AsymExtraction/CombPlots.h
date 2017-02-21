#ifndef __COMP_PLOTS_H__
#define __COMP_PLOTS_H__

#include <iostream>
#include <fstream>
#include "TH1D.h"

using namespace std;

void saveToTxt(int b,int c,int p,int maxZ1,int maxZ2,int numKtBins,TH1D*** sepKtZHistos)
{
  if(maxZ1==0 || maxZ2==0|| numKtBins==0)
    {
      cout <<"saveToTxt problem with zero" <<endl;
      return;
    }

  ofstream f;
  char buffer[300];
  sprintf(buffer,"out_binning_%d_charge_%d_pid_%d.txt",b,c,p);
  f.open(buffer);
  f <<" " << p <<" " << c << " " <<b <<" " << maxZ1 <<" " << maxZ2 <<" " << numKtBins<<endl;
  for(int iZ1=0;iZ1<maxZ1;iZ1++)
    {
      for(int iZ2=0;iZ2<maxZ2;iZ2++)
	{
	   for(int iKt=0;iKt<numKtBins;iKt++) 
	     {
	       f<<sepKtZHistos[iZ1][iZ2]->GetBinContent(iKt+1) << " ";
	     }
	   cout <<endl;
	   //error bars not so clear...
	   for(int iKt=0;iKt<numKtBins;iKt++) 
	     {
	       f<<sqrt(sepKtZHistos[iZ1][iZ2]->GetBinContent(iKt+1)) << " ";
	     }
	   cout <<endl;
	}
    }




}

#endif
