{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Nov 14 21:13:43 2014) by ROOT version5.28/00c
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(9);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetMarkerSize(2);
   gre->SetPoint(0,0.1317298,1.027204);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,0.2775689,1.006796);
   gre->SetPointError(1,0,0);
   gre->SetPoint(2,0.4225202,0.9825228);
   gre->SetPointError(2,0,0);
   gre->SetPoint(3,0.5652312,0.9533765);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,0.7134338,0.9364967);
   gre->SetPointError(4,0,0);
   gre->SetPoint(5,0.8768423,0.9090653);
   gre->SetPointError(5,0,0);
   gre->SetPoint(6,1.109578,0.8705665);
   gre->SetPointError(6,0,0);
   gre->SetPoint(7,1.526904,0.9975388);
   gre->SetPointError(7,0,0);
   gre->SetPoint(8,nan,0);
   gre->SetPointError(8,0,0);
   gre->Draw("ap");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
