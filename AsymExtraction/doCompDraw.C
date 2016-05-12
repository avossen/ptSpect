
TGraphErrors* getRatioPlot(TGraphErrors* g1, TGraphErrors* g2)
{
  float* x=g1->GetX();
  float* y1=g1->GetY();
  float* y2=g2->GetY();

  float* ey1=g1->GetEY();
  float* ey2=g2->GetEY();


  float y[200];
  float ey[200];
  float ex[200];

  for(int i=0;i<g1->GetN();i++)
    {
      ex[i]=0;
      if(y1[i]>0 && y2[i]>0)
	{
	  y[i]=y1[i]/y2[i];
	  ey[i]=sqrt(ey1[i]*ey1[i]/(y1[i]*y1[i]) + ey2[i]*ey2[i]/(y2[i]*y2[i]))*y[i];
	}
      else
	{
	  y[i]=0;
	  ey[i]=0;
	}
    }

  TGraphErrors* gRatio=new TGraphErrors(g1->GetN(),x,y,ex,ey);
  return gRatio;

}


void doCompDraw(char* filename, char* filenameWoA)
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

      tgBin2_0->Draw("AP");
  //          tgBin0_0->Draw("SAME P");
  //        tgBin1_0->Draw("SAME P");

  //        tgBin3_0->Draw("SAME P");
  //	    tgBin3_0->Draw("AP");
	    //      tgBin4_0->Draw("SAME P");


    tgBin2w_0->Draw("SAME P");
    //        tgBin0w_0->Draw("SAME P");
    //        tgBin1w_0->Draw("SAME P");

    //      tgBin3w_0->Draw("SAME P");
    //      tgBin4w_0->Draw("SAME P");

    c.cd(2);

    tgBin2_1->Draw("AP");
    //          tgBin1_1->Draw("SAME P");

    //           tgBin0_1->Draw("SAME P");
    //            tgBin3_1->Draw("SAME P");
	    //        tgBin3_1->Draw("AP");
    //        tgBin4_1->Draw("SAME P");
    
    //          tgBin0w_1->Draw("SAME P");
    //          tgBin1w_1->Draw("SAME P");
    tgBin2w_1->Draw("SAME P");
    //      tgBin3w_1->Draw("SAME P");
    //      tgBin4w_1->Draw("SAME P");
    
  //
  c.cd(3);
  //
    tgBin2_2->Draw("AP");
    //        tgBin0_2->Draw("SAME P");
    //        tgBin1_2->Draw("SAME P");
  
    //        tgBin3_2->Draw("SAME P");
	//      tgBin3_2->Draw("AP");
    //        tgBin4_2->Draw("SAME P");
  
    //          tgBin0w_2->Draw("SAME P");
    //          tgBin1w_2->Draw("SAME P");
      tgBin2w_2->Draw("SAME P");
      //        tgBin3w_2->Draw("SAME P");
      //        tgBin4w_2->Draw("SAME P");
  //
  c.SaveAs("cCompThrust.png");
  //c.SaveAs("cCompThrust.jpg");
  //  c.SaveAs("cCompThrust.pdf");
  //c.SaveAs("cCompThrust.C");
  cout <<"done.."<<endl;
}

