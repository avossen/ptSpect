void dumpHistoToTxt(const char* filename)
{
  TFile* file=new TFile(filename);
  cout <<"opening " << filename <<endl;
  TH2D* histo=(TH2D*)file->Get("kinematicSmearingMatrix_binning0_pidBin0_chargeBin0");
  int numXBins=histo->GetNbinsX();
  int numYBins=histo->GetNbinsY();
  //  cout <<std::fixed;
  //  cout.precision(3);
  for(int i=0;i<numXBins;i++)
    {
      for(int j=0;j<numYBins;j++)
	{
	  cout <<" "<<histo->GetBinContent(i+1,j+1);
	}
      cout <<endl;
    }
}
