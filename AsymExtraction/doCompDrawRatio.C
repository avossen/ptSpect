
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
  
      cout <<" n: "<< g1->GetN()<<endl;
      cout <<"i : " << i << " y: "<< y1[i] << " y2: " << y2[i] <<endl;
      x[i]=x1[i];
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
  return gRatio;

}


void doCompDrawRatio(char* filename, char* filenameWoA)
{
  gROOT->SetStyle("Plain");
  TFile file(filename);
  TFile filew(filenameWoA);
  cout <<"reading file: "<< filename<<endl;

  TGraphErrors* tgBin0_0=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin1_0=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin1_0");
  TGraphErrors* tgBin2_0=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin2_0");
  TGraphErrors* tgBin3_0=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin3_0");
  TGraphErrors* tgBin4_0=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin4_0");


  TGraphErrors* tgBin0w_0=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin1w_0=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin1_0");
  TGraphErrors* tgBin2w_0=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin2_0");
  TGraphErrors* tgBin3w_0=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin3_0");
  TGraphErrors* tgBin4w_0=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin4_0");

  TGraphErrors* tgBin0_1=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin1_1=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin1_1");
  TGraphErrors* tgBin2_1=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin2_1");
  TGraphErrors* tgBin3_1=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin3_1");
  TGraphErrors* tgBin4_1=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin4_1");


  TGraphErrors* tgBin0w_1=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin1w_1=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin1_1");
  TGraphErrors* tgBin2w_1=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin2_1");
  TGraphErrors* tgBin3w_1=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin3_1");
  TGraphErrors* tgBin4w_1=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin4_1");
  

  TGraphErrors* tgBin0_2=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin1_2=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin1_2");
  TGraphErrors* tgBin2_2=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin2_2");
  TGraphErrors* tgBin3_2=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin3_2");
  TGraphErrors* tgBin4_2=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin4_2");


  TGraphErrors* tgBin0w_2=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin1w_2=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin1_2");
  TGraphErrors* tgBin2w_2=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin2_2");
  TGraphErrors* tgBin3w_2=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin3_2");
  TGraphErrors* tgBin4w_2=(TGraphErrors*)filew.Get("NormalWoA_ptSpect_labThrustTheta_z_ChargeInt__bin4_2");


  TGraphErrors* ratio2_0=getRatioPlot(tgBin2_0,tgBin2w_0);
  TGraphErrors* ratio2_1=getRatioPlot(tgBin2_1,tgBin2w_1);
  TGraphErrors* ratio2_2=getRatioPlot(tgBin2_2,tgBin2w_2);

  ratio2_0->SetMarkerSize(2);
  ratio2_1->SetMarkerSize(2);
  ratio2_2->SetMarkerSize(2);

  ratio2_0->SetMarkerStyle(28);
  ratio2_1->SetMarkerStyle(28);
  ratio2_2->SetMarkerStyle(28);

  ratio2_0->SetMarkerColor(kBlue);
  ratio2_1->SetMarkerColor(kBlue);
  ratio2_2->SetMarkerColor(kBlue);


  tgBin3_0->SetMarkerSize(2);
  tgBin3_1->SetMarkerSize(2);
  tgBin3_2->SetMarkerSize(2);

  tgBin3w_0->SetMarkerSize(2);
  tgBin3w_1->SetMarkerSize(2);
  tgBin3w_2->SetMarkerSize(2);


  //  tgBin1->SetMarkerSize(2);
  //  tgBin2->SetMarkerSize(2);
  //    tgBin3->SetMarkerSize(2);
  //  tgBin4->SetMarkerSize(2);
  //
  //  tgBin0w->SetMarkerSize(2);
  //  tgBin1w->SetMarkerSize(2);
  //  tgBin2w->SetMarkerSize(2);
  //   tgBin3w->SetMarkerSize(2);
  //  tgBin4w->SetMarkerSize(2);


  tgBin0_0->SetMarkerStyle(24);
  tgBin1_0->SetMarkerStyle(24);
  tgBin2_0->SetMarkerStyle(24);
  tgBin3_0->SetMarkerStyle(24);
  tgBin4_0->SetMarkerStyle(24);


  tgBin0w_0->SetMarkerStyle(28);
  tgBin1w_0->SetMarkerStyle(28);
  tgBin2w_0->SetMarkerStyle(28);
  tgBin3w_0->SetMarkerStyle(28);
  tgBin4w_0->SetMarkerStyle(28);

  tgBin0_0->SetMarkerColor(kBlack);
  tgBin1_0->SetMarkerColor(kRed);
  tgBin2_0->SetMarkerColor(kBlue);
  tgBin3_0->SetMarkerColor(kGreen);
  tgBin4_0->SetMarkerColor(kMagenta);

  tgBin0w_0->SetMarkerColor(kBlack);
  tgBin1w_0->SetMarkerColor(kRed);
  tgBin2w_0->SetMarkerColor(kBlue);
  tgBin3w_0->SetMarkerColor(kGreen);
  tgBin4w_0->SetMarkerColor(kMagenta);



  ////---

  tgBin0_1->SetMarkerStyle(24);
  tgBin1_1->SetMarkerStyle(24);
  tgBin2_1->SetMarkerStyle(24);
  tgBin3_1->SetMarkerStyle(24);
  tgBin4_1->SetMarkerStyle(24);


  tgBin0w_1->SetMarkerStyle(28);
  tgBin1w_1->SetMarkerStyle(28);
  tgBin2w_1->SetMarkerStyle(28);
  tgBin3w_1->SetMarkerStyle(28);
  tgBin4w_1->SetMarkerStyle(28);

  tgBin0_1->SetMarkerColor(kBlack);
  tgBin1_1->SetMarkerColor(kRed);
  tgBin2_1->SetMarkerColor(kBlue);
  tgBin3_1->SetMarkerColor(kGreen);
  tgBin4_1->SetMarkerColor(kMagenta);

  tgBin0w_1->SetMarkerColor(kBlack);
  tgBin1w_1->SetMarkerColor(kRed);
  tgBin2w_1->SetMarkerColor(kBlue);
  tgBin3w_1->SetMarkerColor(kGreen);
  tgBin4w_1->SetMarkerColor(kMagenta);


  ///---

  tgBin0_2->SetMarkerStyle(24);
  tgBin1_2->SetMarkerStyle(24);
  tgBin2_2->SetMarkerStyle(24);
  tgBin3_2->SetMarkerStyle(24);
  tgBin4_2->SetMarkerStyle(24);


  tgBin0w_2->SetMarkerStyle(28);
  tgBin1w_2->SetMarkerStyle(28);
  tgBin2w_2->SetMarkerStyle(28);
  tgBin3w_2->SetMarkerStyle(28);
  tgBin4w_2->SetMarkerStyle(28);

  tgBin0_2->SetMarkerColor(kBlack);
  tgBin1_2->SetMarkerColor(kRed);
  tgBin2_2->SetMarkerColor(kBlue);
  tgBin3_2->SetMarkerColor(kGreen);
  tgBin4_2->SetMarkerColor(kMagenta);

  tgBin0w_2->SetMarkerColor(kBlack);
  tgBin1w_2->SetMarkerColor(kRed);
  tgBin2w_2->SetMarkerColor(kBlue);
  tgBin3w_2->SetMarkerColor(kGreen);
  tgBin4w_2->SetMarkerColor(kMagenta);




  TCanvas c("c","c",0,0,1400,1000);
  c.Divide(2,2);
  c.cd(1);
  cout <<"?"<<endl;

  ratio2_0->Draw("AP");
  //          tgBin0_0->Draw("SAME P");
  //        tgBin1_0->Draw("SAME P");

  //        tgBin3_0->Draw("SAME P");
  //	    tgBin3_0->Draw("AP");
	    //      tgBin4_0->Draw("SAME P");


    //      tgBin4w_0->Draw("SAME P");

    c.cd(2);

    ratio2_1->Draw("AP");
    //          tgBin1_1->Draw("SAME P");

    //           tgBin0_1->Draw("SAME P");
    //            tgBin3_1->Draw("SAME P");
	    //        tgBin3_1->Draw("AP");
    //        tgBin4_1->Draw("SAME P");
    
    //          tgBin0w_1->Draw("SAME P");
    //          tgBin1w_1->Draw("SAME P");

    //      tgBin3w_1->Draw("SAME P");
    //      tgBin4w_1->Draw("SAME P");
    
  //
  c.cd(3);
  //
  ratio2_2->Draw("AP");
    //        tgBin0_2->Draw("SAME P");
    //        tgBin1_2->Draw("SAME P");
  
    //        tgBin3_2->Draw("SAME P");
	//      tgBin3_2->Draw("AP");
    //        tgBin4_2->Draw("SAME P");
  
    //          tgBin0w_2->Draw("SAME P");
    //          tgBin1w_2->Draw("SAME P");
    //      tgBin2w_2->Draw("SAME P");
      //        tgBin3w_2->Draw("SAME P");
      //        tgBin4w_2->Draw("SAME P");
  //
  c.SaveAs("cCompThrustRatio.png");
  //c.SaveAs("cCompThrust.jpg");
  //  c.SaveAs("cCompThrust.pdf");
  //c.SaveAs("cCompThrust.C");
  cout <<"done.."<<endl;
}

