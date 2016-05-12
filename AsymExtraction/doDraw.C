void doDraw(char* filename)
{
  gROOT->SetStyle("Plain");
  TFile file(filename);
  cout <<"reading file: "<< filename<<endl;

  TGraphErrors* tgBin0=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin0_0");
  TGraphErrors* tgBin1=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin1_0");
  TGraphErrors* tgBin2=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin2_0");
  TGraphErrors* tgBin3=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin3_0");
  TGraphErrors* tgBin4=(TGraphErrors*)file.Get("Normal_ptSpect_labThrustTheta_z_ChargeInt__bin4_0");

  tgBin0->SetMarkerSize(2);
  tgBin1->SetMarkerSize(2);
  tgBin2->SetMarkerSize(2);
  tgBin3->SetMarkerSize(2);
  tgBin4->SetMarkerSize(2);


  tgBin0->SetMarkerStyle(20);
  tgBin1->SetMarkerStyle(21);
  tgBin2->SetMarkerStyle(22);
  tgBin3->SetMarkerStyle(23);
  tgBin4->SetMarkerStyle(24);

  tgBin0->SetMarkerColor(kBlack);
  tgBin1->SetMarkerColor(kRed);
  tgBin2->SetMarkerColor(kBlue);
  tgBin3->SetMarkerColor(kGreen);
  tgBin4->SetMarkerColor(kYellow);

  TCanvas c;

  tgBin0->Draw("AP");
  tgBin1->Draw("SAME P");
  tgBin2->Draw("SAME P");
  tgBin3->Draw("SAME P");
  tgBin4->Draw("SAME P");


  c.SaveAs("c.png");

  }

