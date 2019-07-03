void dumpXiniToTxt(const char* filename)
{
  TFile* file=new TFile(filename);
  cout <<"opening " << filename <<endl;
  TH1D* histo=(TH1D*)file->Get("xini_binning0_pidBin0_chargeBin0");
  int numXBins=histo->GetNbinsX();
  //  cout <<std::fixed;
  //  cout.precision(3);
  for(int i=0;i<numXBins;i++)
    {
	  cout <<" "<<histo->GetBinContent(i+1);
    }
  cout <<endl;
}
