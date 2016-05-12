void doDrawHad(char* filename)
{
  gROOT->SetStyle("Plain");
  TFile file(filename);
  cout <<"reading file: "<< filename<<endl;

  TGraphErrors* tgBin0_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin1_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin1_0");
  TGraphErrors* tgBin2_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin2_0");
  TGraphErrors* tgBin3_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin3_0");
  TGraphErrors* tgBin4_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin4_0");

  TGraphErrors* tgBin0_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_1");
  TGraphErrors* tgBin1_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin1_1");
  TGraphErrors* tgBin2_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin2_1");
  TGraphErrors* tgBin3_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin3_1");
  TGraphErrors* tgBin4_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin4_1");

  TGraphErrors* tgBin0_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin1_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin1_2");
  TGraphErrors* tgBin2_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin2_2");
  TGraphErrors* tgBin3_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin3_2");
  TGraphErrors* tgBin4_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin4_2");



  //    tgBin0->SetMarkerSize(2);
  //  tgBin1->SetMarkerSize(2);
  //  tgBin2->SetMarkerSize(2);
  //  tgBin3->SetMarkerSize(2);
  //  tgBin4->SetMarkerSize(2);


  tgBin0_0->SetMarkerStyle(20);
  tgBin1_0->SetMarkerStyle(21);
  tgBin2_0->SetMarkerStyle(22);
  tgBin3_0->SetMarkerStyle(23);
  tgBin4_0->SetMarkerStyle(24);

  tgBin0_1->SetMarkerStyle(20);
  tgBin1_1->SetMarkerStyle(21);
  tgBin2_1->SetMarkerStyle(22);
  tgBin3_1->SetMarkerStyle(23);
  tgBin4_1->SetMarkerStyle(24);

  tgBin0_2->SetMarkerStyle(20);
  tgBin1_2->SetMarkerStyle(21);
  tgBin2_2->SetMarkerStyle(22);
  tgBin3_2->SetMarkerStyle(23);
  tgBin4_2->SetMarkerStyle(24);



  tgBin0_0->SetMarkerColor(kBlack);
  tgBin1_0->SetMarkerColor(kRed);
  tgBin2_0->SetMarkerColor(kBlue);
  tgBin3_0->SetMarkerColor(kGreen);
  tgBin4_0->SetMarkerColor(kMagenta);

  tgBin0_1->SetMarkerColor(kBlack);
  tgBin1_1->SetMarkerColor(kRed);
  tgBin2_1->SetMarkerColor(kBlue);
  tgBin3_1->SetMarkerColor(kGreen);
  tgBin4_1->SetMarkerColor(kMagenta);


  tgBin0_2->SetMarkerColor(kBlack);
  tgBin1_2->SetMarkerColor(kRed);
  tgBin2_2->SetMarkerColor(kBlue);
  tgBin3_2->SetMarkerColor(kGreen);
  tgBin4_2->SetMarkerColor(kMagenta);

  TCanvas c("c","c",0,0,1400,1000);
  c.Divide(2,2);
  c.cd(1);
  tgBin2_0->Draw("AP");
  tgBin1_0->Draw("SAME P");
  tgBin0_0->Draw("SAME P");
  tgBin3_0->Draw("SAME P");
  tgBin4_0->Draw("SAME P");

  c.cd(2);
  tgBin2_1->Draw("AP");
  tgBin1_1->Draw("SAME P");
  tgBin0_1->Draw("SAME P");
  tgBin3_1->Draw("SAME P");
  tgBin4_1->Draw("SAME P");

  c.cd(3);
  tgBin2_1->Draw("AP");
  tgBin1_1->Draw("SAME P");
  tgBin0_1->Draw("SAME P");
  tgBin3_1->Draw("SAME P");
  tgBin4_1->Draw("SAME P");


  c.SaveAs("cHad.png");

  }

