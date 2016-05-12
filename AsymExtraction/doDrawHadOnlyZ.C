void doDrawHadOnlyZ(char* filename)
{
  gROOT->SetStyle("Plain");
  TFile file(filename);
  cout <<"reading file: "<< filename<<endl;

  TGraphErrors* tgBin0_0=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_0");

  TGraphErrors* tgBin0_1=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_1");


  TGraphErrors* tgBin0_2=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_2");
  TGraphErrors* tgBin0_3=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_3");
  TGraphErrors* tgBin0_4=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_4");
  TGraphErrors* tgBin0_5=(TGraphErrors*)file.Get("Normal_ptSpect_labT_z_ChargeInt__bin0_5");



  //    tgBin0->SetMarkerSize(2);
  //  tgBin1->SetMarkerSize(2);
  //  tgBin2->SetMarkerSize(2);
  //  tgBin3->SetMarkerSize(2);
  //  tgBin4->SetMarkerSize(2);


  tgBin0_0->SetMarkerStyle(20);
  tgBin0_1->SetMarkerStyle(21);
  tgBin0_2->SetMarkerStyle(22);
  tgBin0_3->SetMarkerStyle(23);
  tgBin0_4->SetMarkerStyle(24);
  tgBin0_5->SetMarkerStyle(25);




  tgBin0_0->SetMarkerColor(kBlack);
  tgBin0_1->SetMarkerColor(kRed);
  tgBin0_2->SetMarkerColor(kBlue);
  tgBin0_3->SetMarkerColor(kGreen);
  tgBin0_4->SetMarkerColor(kMagenta);
  tgBin0_5->SetMarkerColor(kYellow);


  TCanvas c("c","c",0,0,1400,1000);

  tgBin0_0->Draw("AP");
  tgBin0_1->Draw("SAME P");
  tgBin0_2->Draw("SAME P");
  tgBin0_3->Draw("SAME P");
  tgBin0_4->Draw("SAME P");
  tgBin0_5->Draw("SAME P");



  c.SaveAs("cHadOnlyZ.png");

  }

