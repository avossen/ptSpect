
TGraphErrors* getRatioPlot(TGraphErrors* g1, TGraphErrors* g2)
{
  Double_t* x1=g1->GetX();
  Double_t* y1=g1->GetY();
  Double_t* y2=g2->GetY();

  Double_t* ey1=g1->GetEY();
  Double_t* ey2=g2->GetEY();

  float x[200];
  float y[200];
  float ey[200];
  float ex[200];

  for(int i=0;i<g1->GetN();i++)
    {
      x[i]=x1[i];
      cout <<" n: "<< g1->GetN()<<endl;
      cout <<"i : " << i << " y: "<< y1[i] << " y2: " << y2[i] <<endl;
      ex[i]=0;
      if(y1[i]>0 && y2[i]>0)
	{
	  y[i]=y1[i]/y2[i];
	  cout <<"y new: "<< y[i] <<endl;
	  //	  ey[i]=sqrt(ey1[i]*ey1[i]/(y1[i]*y1[i]) + ey2[i]*ey2[i]/(y2[i]*y2[i]))*y[i];
	  ey[i]=0.0;
	}
      else
	{
	  y[i]=0;
	  ey[i]=0;
	}
    }
  cout <<endl;
  for(int i=0;i<g1->GetN();i++)
    {
      cout <<"x: " << x[i] << ", y: " <<y[i] << ", ex: " << ex[i] << " ey: " << ey[i] <<endl;
    }

  TGraphErrors* gRatio=new TGraphErrors(g1->GetN(),x,y,ex,ey);
  gRatio->SetMarkerSize(2);
  gRatio->SetMarkerStyle(24);
  return gRatio;

}



void doCompSpectra(char* filename, char* filenameWoA)
{
  gROOT->SetStyle("Plain");
  TFile file(filename);
  TFile filew(filenameWoA);
  cout <<"reading file: "<< filename<<endl;


  TGraphErrors* tgBin0_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin0_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin0_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin0_3=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_3");
  TGraphErrors* tgBin0_4=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_4");


  TGraphErrors* tgBin0w_0=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin0w_1=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin0w_2=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin0w_3=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_3");
  TGraphErrors* tgBin0w_4=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_4");




  TGraphErrors* tgBin0_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin0_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin0_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin0_3=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_3");
  TGraphErrors* tgBin0_4=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_4");



  TGraphErrors* tgBin0w_0=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin0w_1=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin0w_2=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin0w_3=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_3");
  TGraphErrors* tgBin0w_4=(TGraphErrors*)filew.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_4");



  TGraphErrors* z_zBin0_0=(TGraphErrors*)file.Get("Normal_ptSpect_z_z_ChargeInt__bin0_0");
  TGraphErrors* z_zBin1_1=(TGraphErrors*)file.Get("Normal_ptSpect_z_z_ChargeInt__bin1_1");
  TGraphErrors* z_zBin2_2=(TGraphErrors*)file.Get("Normal_ptSpect_z_z_ChargeInt__bin2_2");
  TGraphErrors* z_zBin3_3=(TGraphErrors*)file.Get("Normal_ptSpect_z_z_ChargeInt__bin3_3");
  TGraphErrors* z_zBin4_4=(TGraphErrors*)file.Get("Normal_ptSpect_z_z_ChargeInt__bin4_4");


  TGraphErrors* z_zBinw0_0=(TGraphErrors*)filew.Get("Normal_ptSpect_z_z_ChargeInt__bin0_0");
  TGraphErrors* z_zBinw1_1=(TGraphErrors*)filew.Get("Normal_ptSpect_z_z_ChargeInt__bin1_1");
  TGraphErrors* z_zBinw2_2=(TGraphErrors*)filew.Get("Normal_ptSpect_z_z_ChargeInt__bin2_2");
  TGraphErrors* z_zBinw3_3=(TGraphErrors*)filew.Get("Normal_ptSpect_z_z_ChargeInt__bin3_3");
  TGraphErrors* z_zBinw4_4=(TGraphErrors*)filew.Get("Normal_ptSpect_z_z_ChargeInt__bin4_4");




  z_zBin0_0->SetMarkerSize(2);
  z_zBin0_0->SetMarkerStyle(28);
  z_zBin0_0->SetMarkerColor(kBlue);


  z_zBin1_1->SetMarkerSize(2);
  z_zBin1_1->SetMarkerStyle(28);
  z_zBin1_1->SetMarkerColor(kBlue);

  z_zBin2_2->SetMarkerSize(2);
  z_zBin2_2->SetMarkerStyle(28);
  z_zBin2_2->SetMarkerColor(kBlue);

  z_zBin3_3->SetMarkerSize(2);
  z_zBin3_3->SetMarkerStyle(28);
  z_zBin3_3->SetMarkerColor(kBlue);

  z_zBin4_4->SetMarkerSize(2);
  z_zBin4_4->SetMarkerStyle(28);
  z_zBin4_4->SetMarkerColor(kBlue);


  z_zBinw0_0->SetMarkerSize(2);
  z_zBinw0_0->SetMarkerStyle(24);
  z_zBinw0_0->SetMarkerColor(kGreen);


  z_zBinw1_1->SetMarkerSize(2);
  z_zBinw1_1->SetMarkerStyle(24);
  z_zBinw1_1->SetMarkerColor(kGreen);

  z_zBinw2_2->SetMarkerSize(2);
  z_zBinw2_2->SetMarkerStyle(24);
  z_zBinw2_2->SetMarkerColor(kGreen);

  z_zBinw3_3->SetMarkerSize(2);
  z_zBinw3_3->SetMarkerStyle(24);
  z_zBinw3_3->SetMarkerColor(kGreen);

  z_zBinw4_4->SetMarkerSize(2);
  z_zBinw4_4->SetMarkerStyle(24);
  z_zBinw4_4->SetMarkerColor(kGreen);


  TCanvas c("c","c",0,0,1600,1200);
  c.Divide(3,2);
  c.cd(1);
  z_zBin0_0->Draw("AP");
  z_zBinw0_0->Draw("SAME P");
  c.cd(2);
  z_zBin1_1->Draw("AP");
  z_zBinw1_1->Draw("SAME P");
  c.cd(3);
  z_zBin2_2->Draw("AP");
  z_zBinw2_2->Draw("SAME P");
  c.cd(4);
  z_zBin3_3->Draw("AP");
  z_zBinw3_3->Draw("SAME P");
  c.cd(5);
  z_zBin4_4->Draw("AP");
  z_zBinw4_4->Draw("SAME P");

  c.SaveAs("cCompHadSpectra.png");

  }

