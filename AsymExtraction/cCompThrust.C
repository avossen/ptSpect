{
//=========Macro generated from canvas: c/c
//=========  (Sat Jul 12 15:16:39 2014) by ROOT version5.28/00c
   TCanvas *c = new TCanvas("c", "c",0,0,1400,1000);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c_1
   TPad *c_1 = new TPad("c_1", "c_1",0.01,0.51,0.49,0.99);
   c_1->Draw();
   c_1->cd();
   c_1->Range(-0.1164395,-0.151461,1.404769,1.363149);
   c_1->SetFillColor(0);
   c_1->SetBorderMode(0);
   c_1->SetBorderSize(2);
   c_1->SetFrameBorderMode(0);
   c_1->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin0_0");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin0_0");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetMarkerSize(2);
   gre->SetPoint(0,0.1370952,0.3814433);
   gre->SetPointError(0,0,0.06270889);
   gre->SetPoint(1,0.2767749,0.5979381);
   gre->SetPointError(1,0,0.07851312);
   gre->SetPoint(2,0.4259005,0.8350515);
   gre->SetPointError(2,0,0.0927835);
   gre->SetPoint(3,0.5723587,1);
   gre->SetPointError(3,0,0.1015346);
   gre->SetPoint(4,0.7180457,0.6907216);
   gre->SetPointError(4,0,0.08438507);
   gre->SetPoint(5,0.8831728,0.742268);
   gre->SetPointError(5,0,0.08747713);
   gre->SetPoint(6,1.151234,0.7731959);
   gre->SetPointError(6,0,0.08928097);
   gre->SetPoint(7,nan,0);
   gre->SetPointError(7,0,0);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin0_01__1 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin0_01__1","Normal_ptSpect_not_rec_ChargeInt__bin0_0",100,0.03568128,1.252648);
   Normal_ptSpect_not_rec_ChargeInt__bin0_01__1->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_01__1->SetMaximum(1.211688);
   Normal_ptSpect_not_rec_ChargeInt__bin0_01__1->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_01__1->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_01__1->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin0_01__1->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin0_01);
   
   gre->Draw("ap");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin1_0");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin1_0");
   gre->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1294957,0.4);
   gre->SetPointError(0,0,0.05547002);
   gre->SetPoint(1,0.2719184,0.8923077);
   gre->SetPointError(1,0,0.08284869);
   gre->SetPoint(2,0.4184492,0.976923);
   gre->SetPointError(2,0,0.0866879);
   gre->SetPoint(3,0.5692582,1);
   gre->SetPointError(3,0,0.08770581);
   gre->SetPoint(4,0.7187214,0.7384616);
   gre->SetPointError(4,0,0.07536892);
   gre->SetPoint(5,0.898324,0.4846154);
   gre->SetPointError(5,0,0.0610558);
   gre->SetPoint(6,1.203328,0.3692308);
   gre->SetPointError(6,0,0.05329387);
   gre->SetPoint(7,nan,0);
   gre->SetPointError(7,0,0);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin1_02__2 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin1_02__2","Normal_ptSpect_not_rec_ChargeInt__bin1_0",100,0.0221125,1.310711);
   Normal_ptSpect_not_rec_ChargeInt__bin1_02__2->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_02__2->SetMaximum(1.196476);
   Normal_ptSpect_not_rec_ChargeInt__bin1_02__2->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_02__2->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_02__2->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin1_02__2->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin1_02);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin2_0");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin2_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1321255,0.5210701);
   gre->SetPointError(0,0,0.001575784);
   gre->SetPoint(1,0.2785937,0.9067606);
   gre->SetPointError(1,0,0.002078715);
   gre->SetPoint(2,0.4244567,1);
   gre->SetPointError(2,0,0.002182974);
   gre->SetPoint(3,0.5678643,0.8102332);
   gre->SetPointError(3,0,0.00196496);
   gre->SetPoint(4,0.7168022,0.4270731);
   gre->SetPointError(4,0,0.001426592);
   gre->SetPoint(5,0.885142,0.2617526);
   gre->SetPointError(5,0,0.001116848);
   gre->SetPoint(6,1.157475,0.1580104);
   gre->SetPointError(6,0,0.0008677436);
   gre->SetPoint(7,1.531191,0.002177777);
   gre->SetPointError(7,0,0.0001018721);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin2_03__3 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin2_03__3","Normal_ptSpect_not_rec_ChargeInt__bin2_0",100,0,1.671098);
   Normal_ptSpect_not_rec_ChargeInt__bin2_03__3->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_03__3->SetMaximum(1.102194);
   Normal_ptSpect_not_rec_ChargeInt__bin2_03__3->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_03__3->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_03__3->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin2_03__3->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin2_03);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin3_0");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin3_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1322966,0.5148906);
   gre->SetPointError(0,0,0.00074435);
   gre->SetPoint(1,0.2786025,0.8983256);
   gre->SetPointError(1,0,0.000983188);
   gre->SetPoint(2,0.4246067,1);
   gre->SetPointError(2,0,0.001037337);
   gre->SetPoint(3,0.5674884,0.7958033);
   gre->SetPointError(3,0,0.0009253853);
   gre->SetPoint(4,0.7168631,0.4191647);
   gre->SetPointError(4,0,0.0006716021);
   gre->SetPoint(5,0.8852012,0.2600521);
   gre->SetPointError(5,0,0.0005289929);
   gre->SetPoint(6,1.15132,0.1482917);
   gre->SetPointError(6,0,0.0003994645);
   gre->SetPoint(7,1.530534,0.001605492);
   gre->SetPointError(7,0,4.156462e-05);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin3_04__4 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin3_04__4","Normal_ptSpect_not_rec_ChargeInt__bin3_0",100,0,1.670358);
   Normal_ptSpect_not_rec_ChargeInt__bin3_04__4->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_04__4->SetMaximum(1.100985);
   Normal_ptSpect_not_rec_ChargeInt__bin3_04__4->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_04__4->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_04__4->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin3_04__4->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin3_04);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin4_0");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin4_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1321293,0.5232435);
   gre->SetPointError(0,0,0.001123073);
   gre->SetPoint(1,0.2784429,0.9095691);
   gre->SetPointError(1,0,0.001480723);
   gre->SetPoint(2,0.4244992,1);
   gre->SetPointError(2,0,0.001552587);
   gre->SetPoint(3,0.5673434,0.7912363);
   gre->SetPointError(3,0,0.001381049);
   gre->SetPoint(4,0.7165446,0.4116216);
   gre->SetPointError(4,0,0.000996105);
   gre->SetPoint(5,0.8849229,0.2533609);
   gre->SetPointError(5,0,0.0007814943);
   gre->SetPoint(6,1.151118,0.1425755);
   gre->SetPointError(6,0,0.000586244);
   gre->SetPoint(7,1.530489,0.00164639);
   gre->SetPointError(7,0,6.299737e-05);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin4_05__5 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin4_05__5","Normal_ptSpect_not_rec_ChargeInt__bin4_0",100,0,1.670325);
   Normal_ptSpect_not_rec_ChargeInt__bin4_05__5->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_05__5->SetMaximum(1.10155);
   Normal_ptSpect_not_rec_ChargeInt__bin4_05__5->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_05__5->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_05__5->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin4_05__5->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin4_05);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_0");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_0");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1242837,0.25);
   gre->SetPointError(0,0,0.1767767);
   gre->SetPoint(1,0.294622,0.875);
   gre->SetPointError(1,0,0.3307189);
   gre->SetPoint(2,0.4288761,1);
   gre->SetPointError(2,0,0.3535534);
   gre->SetPoint(3,0.5773602,0.625);
   gre->SetPointError(3,0,0.2795085);
   gre->SetPoint(4,0.7832019,0.25);
   gre->SetPointError(4,0,0.1767767);
   gre->SetPoint(5,0.9054664,0.625);
   gre->SetPointError(5,0,0.2795085);
   gre->SetPoint(6,1.107958,0.375);
   gre->SetPointError(6,0,0.2165063);
   gre->SetPoint(7,nan,0);
   gre->SetPointError(7,0,0);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6","NormalWoA_ptSpect_not_rec_ChargeInt__bin0_0",100,0.02591629,1.206325);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6->SetMaximum(1.488909);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06__6->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin0_06);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_0");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1365568,0.3763441);
   gre->SetPointError(0,0,0.04498172);
   gre->SetPoint(1,0.2780109,0.7258065);
   gre->SetPointError(1,0,0.06246747);
   gre->SetPoint(2,0.4220426,0.983871);
   gre->SetPointError(2,0,0.07272983);
   gre->SetPoint(3,0.5748953,1);
   gre->SetPointError(3,0,0.07332356);
   gre->SetPoint(4,0.7183709,0.5913978);
   gre->SetPointError(4,0,0.05638757);
   gre->SetPoint(5,0.8930637,0.4408602);
   gre->SetPointError(5,0,0.04868487);
   gre->SetPoint(6,1.181148,0.3172043);
   gre->SetPointError(6,0,0.04129648);
   gre->SetPoint(7,1.543743,0.005376344);
   gre->SetPointError(7,0,0.005376344);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7","NormalWoA_ptSpect_not_rec_ChargeInt__bin1_0",100,0,1.684462);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7->SetMaximum(1.180656);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07__7->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin1_07);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_0");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.131845,0.6073665);
   gre->SetPointError(0,0,0.0009789132);
   gre->SetPoint(1,0.2773381,0.9905019);
   gre->SetPointError(1,0,0.001250105);
   gre->SetPoint(2,0.4234827,1);
   gre->SetPointError(2,0,0.001256084);
   gre->SetPoint(3,0.5675178,0.7678202);
   gre->SetPointError(3,0,0.001100648);
   gre->SetPoint(4,0.7164583,0.3962371);
   gre->SetPointError(4,0,0.0007906719);
   gre->SetPoint(5,0.8852307,0.2411303);
   gre->SetPointError(5,0,0.0006168004);
   gre->SetPoint(6,1.158405,0.1454604);
   gre->SetPointError(6,0,0.0004790614);
   gre->SetPoint(7,1.53162,0.00240922);
   gre->SetPointError(7,0,6.165339e-05);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8","NormalWoA_ptSpect_not_rec_ChargeInt__bin2_0",100,0,1.671597);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8->SetMaximum(1.101147);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08__8->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin2_08);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_0");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1318,0.6024363);
   gre->SetPointError(0,0,0.000540684);
   gre->SetPoint(1,0.2774462,0.9846522);
   gre->SetPointError(1,0,0.0006912409);
   gre->SetPoint(2,0.4234792,1);
   gre->SetPointError(2,0,0.0006966072);
   gre->SetPoint(3,0.5674084,0.7536399);
   gre->SetPointError(3,0,0.0006047417);
   gre->SetPoint(4,0.7164506,0.3877954);
   gre->SetPointError(4,0,0.0004337997);
   gre->SetPoint(5,0.8847618,0.2340174);
   gre->SetPointError(5,0,0.0003369862);
   gre->SetPoint(6,1.152647,0.1321081);
   gre->SetPointError(6,0,0.0002531936);
   gre->SetPoint(7,1.5311,0.001784307);
   gre->SetPointError(7,0,2.942543e-05);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9","NormalWoA_ptSpect_not_rec_ChargeInt__bin3_0",100,0,1.67103);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9->SetMaximum(1.100591);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09__9->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin3_09);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_0");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_0");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1319506,0.6092606);
   gre->SetPointError(0,0,0.001056187);
   gre->SetPoint(1,0.2772387,0.9868647);
   gre->SetPointError(1,0,0.001344213);
   gre->SetPoint(2,0.4237089,1);
   gre->SetPointError(2,0,0.001353129);
   gre->SetPoint(3,0.5673719,0.7655329);
   gre->SetPointError(3,0,0.001183917);
   gre->SetPoint(4,0.7165129,0.3920357);
   gre->SetPointError(4,0,0.0008472314);
   gre->SetPoint(5,0.8849803,0.2376969);
   gre->SetPointError(5,0,0.0006597069);
   gre->SetPoint(6,1.156512,0.1402148);
   gre->SetPointError(6,0,0.0005066828);
   gre->SetPoint(7,1.531491,0.001940816);
   gre->SetPointError(7,0,5.961169e-05);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10","NormalWoA_ptSpect_not_rec_ChargeInt__bin4_0",100,0,1.671445);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10->SetMaximum(1.1013);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010__10->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin4_010);
   
   gre->Draw(" p");
   
   TPaveText *pt = new TPaveText(0.01,0.9384156,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(0);
   TText *text = pt->AddText("Normal_ptSpect_not_rec_ChargeInt__bin0_0");
   pt->Draw();
   c_1->Modified();
   c->cd();
  
// ------------>Primitives in pad: c_2
   c_2 = new TPad("c_2", "c_2",0.51,0.51,0.99,0.99);
   c_2->Draw();
   c_2->cd();
   c_2->Range(nan,-0.1539344,nan,1.38541);
   c_2->SetFillColor(0);
   c_2->SetBorderMode(0);
   c_2->SetBorderSize(2);
   c_2->SetFrameBorderMode(0);
   c_2->SetFrameBorderMode(0);
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin0_1");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin0_1");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,nan,0);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,0.2866004,0.01428571);
   gre->SetPointError(1,0,0.01428571);
   gre->SetPoint(2,0.4007961,0.01428571);
   gre->SetPointError(2,0,0.01428571);
   gre->SetPoint(3,0.5764613,0.02857143);
   gre->SetPointError(3,0,0.02020305);
   gre->SetPoint(4,0.7391571,0.02857143);
   gre->SetPointError(4,0,0.02020305);
   gre->SetPoint(5,0.858524,0.04285714);
   gre->SetPointError(5,0,0.02474358);
   gre->SetPoint(6,1.299786,0.5571429);
   gre->SetPointError(6,0,0.08921426);
   gre->SetPoint(7,1.869553,1);
   gre->SetPointError(7,0,0.1195229);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin0_111__11 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin0_111__11","Normal_ptSpect_not_rec_ChargeInt__bin0_1",100,nan,nan);
   Normal_ptSpect_not_rec_ChargeInt__bin0_111__11->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_111__11->SetMaximum(1.231475);
   Normal_ptSpect_not_rec_ChargeInt__bin0_111__11->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_111__11->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_111__11->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin0_111__11->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin0_111);
   
   gre->Draw("ap");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin1_1");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin1_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1077749,0.1489362);
   gre->SetPointError(0,0,0.05629258);
   gre->SetPoint(1,0.284043,0.2553191);
   gre->SetPointError(1,0,0.07370429);
   gre->SetPoint(2,0.4326739,0.4255319);
   gre->SetPointError(2,0,0.09515183);
   gre->SetPoint(3,0.5659135,0.4042553);
   gre->SetPointError(3,0,0.09274253);
   gre->SetPoint(4,0.7286654,0.4468085);
   gre->SetPointError(4,0,0.09750161);
   gre->SetPoint(5,0.8984593,0.2765957);
   gre->SetPointError(5,0,0.07671385);
   gre->SetPoint(6,1.211881,0.9574468);
   gre->SetPointError(6,0,0.1427277);
   gre->SetPoint(7,1.984107,1);
   gre->SetPointError(7,0,0.145865);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin1_112__12 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin1_112__12","Normal_ptSpect_not_rec_ChargeInt__bin1_1",100,0,2.17174);
   Normal_ptSpect_not_rec_ChargeInt__bin1_112__12->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_112__12->SetMaximum(1.251187);
   Normal_ptSpect_not_rec_ChargeInt__bin1_112__12->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_112__12->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_112__12->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin1_112__12->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin1_112);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin2_1");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin2_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1319395,0.2071928);
   gre->SetPointError(0,0,0.002071204);
   gre->SetPoint(1,0.2801622,0.3898298);
   gre->SetPointError(1,0,0.002841011);
   gre->SetPoint(2,0.4273224,0.5065013);
   gre->SetPointError(2,0,0.003238364);
   gre->SetPoint(3,0.5752143,0.5494223);
   gre->SetPointError(3,0,0.003372784);
   gre->SetPoint(4,0.7241936,0.5332519);
   gre->SetPointError(4,0,0.00332278);
   gre->SetPoint(5,0.8961132,0.6238561);
   gre->SetPointError(5,0,0.003593996);
   gre->SetPoint(6,1.219695,1);
   gre->SetPointError(6,0,0.004550252);
   gre->SetPoint(7,1.876572,0.6215578);
   gre->SetPointError(7,0,0.00358737);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin2_113__13 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin2_113__13","Normal_ptSpect_not_rec_ChargeInt__bin2_1",100,0,2.051036);
   Normal_ptSpect_not_rec_ChargeInt__bin2_113__13->SetMinimum(0.1251788);
   Normal_ptSpect_not_rec_ChargeInt__bin2_113__13->SetMaximum(1.084493);
   Normal_ptSpect_not_rec_ChargeInt__bin2_113__13->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_113__13->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_113__13->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin2_113__13->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin2_113);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin3_1");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin3_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1321023,0.1968996);
   gre->SetPointError(0,0,0.0009397931);
   gre->SetPoint(1,0.2802324,0.3707118);
   gre->SetPointError(1,0,0.00128952);
   gre->SetPoint(2,0.4267718,0.4867361);
   gre->SetPointError(2,0,0.001477599);
   gre->SetPoint(3,0.5751499,0.5315875);
   gre->SetPointError(3,0,0.001544178);
   gre->SetPoint(4,0.7239994,0.5159956);
   gre->SetPointError(4,0,0.001521363);
   gre->SetPoint(5,0.8968878,0.606017);
   gre->SetPointError(5,0,0.00164874);
   gre->SetPoint(6,1.222327,1);
   gre->SetPointError(6,0,0.002117922);
   gre->SetPoint(7,1.865443,0.6076497);
   gre->SetPointError(7,0,0.00165096);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin3_114__14 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin3_114__14","Normal_ptSpect_not_rec_ChargeInt__bin3_1",100,0,2.038777);
   Normal_ptSpect_not_rec_ChargeInt__bin3_114__14->SetMinimum(0.115344);
   Normal_ptSpect_not_rec_ChargeInt__bin3_114__14->SetMaximum(1.082734);
   Normal_ptSpect_not_rec_ChargeInt__bin3_114__14->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_114__14->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_114__14->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin3_114__14->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin3_114);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin4_1");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin4_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1324635,0.2055607);
   gre->SetPointError(0,0,0.001416112);
   gre->SetPoint(1,0.2804911,0.3804692);
   gre->SetPointError(1,0,0.001926581);
   gre->SetPoint(2,0.4267519,0.4969221);
   gre->SetPointError(2,0,0.002201768);
   gre->SetPoint(3,0.5753928,0.5402176);
   gre->SetPointError(3,0,0.002295682);
   gre->SetPoint(4,0.7237109,0.5227257);
   gre->SetPointError(4,0,0.00225821);
   gre->SetPoint(5,0.896905,0.6051022);
   gre->SetPointError(5,0,0.002429639);
   gre->SetPoint(6,1.221667,1);
   gre->SetPointError(6,0,0.003123399);
   gre->SetPoint(7,1.88098,0.6360763);
   gre->SetPointError(7,0,0.002491048);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin4_115__15 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin4_115__15","Normal_ptSpect_not_rec_ChargeInt__bin4_1",100,0,2.055832);
   Normal_ptSpect_not_rec_ChargeInt__bin4_115__15->SetMinimum(0.1242467);
   Normal_ptSpect_not_rec_ChargeInt__bin4_115__15->SetMaximum(1.083021);
   Normal_ptSpect_not_rec_ChargeInt__bin4_115__15->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_115__15->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_115__15->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin4_115__15->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin4_115);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_1");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_1");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,nan,0);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,0.325306,0.1428571);
   gre->SetPointError(1,0,0.1428571);
   gre->SetPoint(2,0.4353878,0.2857143);
   gre->SetPointError(2,0,0.2020305);
   gre->SetPoint(3,nan,0);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,nan,0);
   gre->SetPointError(4,0,0);
   gre->SetPoint(5,0.9375295,0.2857143);
   gre->SetPointError(5,0,0.2020305);
   gre->SetPoint(6,1.258918,1);
   gre->SetPointError(6,0,0.3779645);
   gre->SetPoint(7,1.664533,0.4285714);
   gre->SetPointError(7,0,0.2474358);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16","NormalWoA_ptSpect_not_rec_ChargeInt__bin0_1",100,nan,nan);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16->SetMaximum(1.515761);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116__16->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin0_116);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_1");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1183827,0.1176471);
   gre->SetPointError(0,0,0.04159452);
   gre->SetPoint(1,0.2948781,0.04411765);
   gre->SetPointError(1,0,0.02547134);
   gre->SetPoint(2,0.4201128,0.1764706);
   gre->SetPointError(2,0,0.05094267);
   gre->SetPoint(3,0.592391,0.2647059);
   gre->SetPointError(3,0,0.06239177);
   gre->SetPoint(4,0.733906,0.2647059);
   gre->SetPointError(4,0,0.06239177);
   gre->SetPoint(5,0.8944535,0.5294118);
   gre->SetPointError(5,0,0.0882353);
   gre->SetPoint(6,1.262164,0.8970588);
   gre->SetPointError(6,0,0.1148566);
   gre->SetPoint(7,1.88686,1);
   gre->SetPointError(7,0,0.1212678);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17","NormalWoA_ptSpect_not_rec_ChargeInt__bin1_1",100,0,2.063708);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17->SetMaximum(1.23153);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117__17->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin1_117);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_1");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1323045,0.2302836);
   gre->SetPointError(0,0,0.001331564);
   gre->SetPoint(1,0.2801465,0.4297692);
   gre->SetPointError(1,0,0.001819065);
   gre->SetPoint(2,0.4266675,0.555494);
   gre->SetPointError(2,0,0.002068094);
   gre->SetPoint(3,0.5749549,0.59587);
   gre->SetPointError(3,0,0.002141935);
   gre->SetPoint(4,0.7239273,0.565511);
   gre->SetPointError(4,0,0.002086657);
   gre->SetPoint(5,0.8961359,0.6440918);
   gre->SetPointError(5,0,0.002226919);
   gre->SetPoint(6,1.218917,1);
   gre->SetPointError(6,0,0.002774793);
   gre->SetPoint(7,1.890698,0.6580278);
   gre->SetPointError(7,0,0.002250882);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18","NormalWoA_ptSpect_not_rec_ChargeInt__bin2_1",100,0,2.066538);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18->SetMinimum(0.1515697);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18->SetMaximum(1.080157);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118__18->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin2_118);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_1");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1322784,0.2325845);
   gre->SetPointError(0,0,0.0007396355);
   gre->SetPoint(1,0.2801518,0.4306003);
   gre->SetPointError(1,0,0.001006386);
   gre->SetPoint(2,0.4267338,0.557851);
   gre->SetPointError(2,0,0.001145477);
   gre->SetPoint(3,0.5750536,0.5961807);
   gre->SetPointError(3,0,0.001184176);
   gre->SetPoint(4,0.7238197,0.5675557);
   gre->SetPointError(4,0,0.001155398);
   gre->SetPoint(5,0.8963154,0.6379562);
   gre->SetPointError(5,0,0.001224963);
   gre->SetPoint(6,1.219352,1);
   gre->SetPointError(6,0,0.001533654);
   gre->SetPoint(7,1.87186,0.603865);
   gre->SetPointError(7,0,0.001191783);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19","NormalWoA_ptSpect_not_rec_ChargeInt__bin3_1",100,0,2.045818);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19->SetMinimum(0.154876);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19->SetMaximum(1.078503);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119__19->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin3_119);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_1");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_1");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1324055,0.2432129);
   gre->SetPointError(0,0,0.001426892);
   gre->SetPoint(1,0.2799266,0.4439496);
   gre->SetPointError(1,0,0.001927812);
   gre->SetPoint(2,0.4267337,0.5651835);
   gre->SetPointError(2,0,0.002175167);
   gre->SetPoint(3,0.5751363,0.5998744);
   gre->SetPointError(3,0,0.002240929);
   gre->SetPoint(4,0.723424,0.555523);
   gre->SetPointError(4,0,0.002156497);
   gre->SetPoint(5,0.8962561,0.6222427);
   gre->SetPointError(5,0,0.002282326);
   gre->SetPoint(6,1.221605,1);
   gre->SetPointError(6,0,0.002893329);
   gre->SetPoint(7,1.892024,0.6641748);
   gre->SetPointError(7,0,0.002357974);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20","NormalWoA_ptSpect_not_rec_ChargeInt__bin4_1",100,0,2.067985);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20->SetMinimum(0.1656753);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20->SetMaximum(1.079004);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120__20->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin4_120);
   
   gre->Draw(" p");
   
   pt = new TPaveText(nan,0.9435597,0.01,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(0);
   text = pt->AddText("Normal_ptSpect_not_rec_ChargeInt__bin0_1");
   pt->Draw();
   c_2->Modified();
   c->cd();
  
// ------------>Primitives in pad: c_3
   c_3 = new TPad("c_3", "c_3",0.01,0.01,0.49,0.49);
   c_3->Draw();
   c_3->cd();
   c_3->Range(-0.3152667,-0.1385715,2.8374,1.247144);
   c_3->SetFillColor(0);
   c_3->SetBorderMode(0);
   c_3->SetBorderSize(2);
   c_3->SetFrameBorderMode(0);
   c_3->SetFrameBorderMode(0);
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin2_2");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin2_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.1275951,0.055293);
   gre->SetPointError(0,0,0.002951322);
   gre->SetPoint(1,0.2836142,0.09814115);
   gre->SetPointError(1,0,0.003931942);
   gre->SetPoint(2,0.4283535,0.1446125);
   gre->SetPointError(2,0,0.004772923);
   gre->SetPoint(3,0.5781447,0.173598);
   gre->SetPointError(3,0,0.005229424);
   gre->SetPoint(4,0.7260581,0.1874606);
   gre->SetPointError(4,0,0.005434212);
   gre->SetPoint(5,0.9005331,0.2605545);
   gre->SetPointError(5,0,0.006406648);
   gre->SetPoint(6,1.238544,0.5521424);
   gre->SetPointError(6,0,0.009326251);
   gre->SetPoint(7,2.304448,1);
   gre->SetPointError(7,0,0.01255109);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin2_221__21 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin2_221__21","Normal_ptSpect_not_rec_ChargeInt__bin2_2",100,0,2.522134);
   Normal_ptSpect_not_rec_ChargeInt__bin2_221__21->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_221__21->SetMaximum(1.108572);
   Normal_ptSpect_not_rec_ChargeInt__bin2_221__21->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_221__21->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin2_221__21->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin2_221__21->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin2_221);
   
   gre->Draw("ap");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin0_2");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin0_2");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,nan,0);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,nan,0);
   gre->SetPointError(1,0,0);
   gre->SetPoint(2,nan,0);
   gre->SetPointError(2,0,0);
   gre->SetPoint(3,nan,0);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,nan,0);
   gre->SetPointError(4,0,0);
   gre->SetPoint(5,nan,0);
   gre->SetPointError(5,0,0);
   gre->SetPoint(6,nan,0);
   gre->SetPointError(6,0,0);
   gre->SetPoint(7,2.942596,1);
   gre->SetPointError(7,0,0.3162278);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin0_222__22 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin0_222__22","Normal_ptSpect_not_rec_ChargeInt__bin0_2",100,nan,nan);
   Normal_ptSpect_not_rec_ChargeInt__bin0_222__22->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_222__22->SetMaximum(1.447851);
   Normal_ptSpect_not_rec_ChargeInt__bin0_222__22->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_222__22->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin0_222__22->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin0_222__22->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin0_222);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin1_2");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin1_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,nan,0);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,nan,0);
   gre->SetPointError(1,0,0);
   gre->SetPoint(2,0.3807934,0.25);
   gre->SetPointError(2,0,0.25);
   gre->SetPoint(3,nan,0);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,nan,0);
   gre->SetPointError(4,0,0);
   gre->SetPoint(5,nan,0);
   gre->SetPointError(5,0,0);
   gre->SetPoint(6,1.333806,0.5);
   gre->SetPointError(6,0,0.3535534);
   gre->SetPoint(7,2.559893,1);
   gre->SetPointError(7,0,0.5);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin1_223__23 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin1_223__23","Normal_ptSpect_not_rec_ChargeInt__bin1_2",100,nan,nan);
   Normal_ptSpect_not_rec_ChargeInt__bin1_223__23->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_223__23->SetMaximum(1.65);
   Normal_ptSpect_not_rec_ChargeInt__bin1_223__23->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_223__23->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin1_223__23->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin1_223__23->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin1_223);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin3_2");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin3_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.131622,0.04785335);
   gre->SetPointError(0,0,0.001198207);
   gre->SetPoint(1,0.2804222,0.08988629);
   gre->SetPointError(1,0,0.001642187);
   gre->SetPoint(2,0.4283551,0.1313192);
   gre->SetPointError(2,0,0.001984906);
   gre->SetPoint(3,0.5772507,0.1579911);
   gre->SetPointError(3,0,0.002177169);
   gre->SetPoint(4,0.7264562,0.1724821);
   gre->SetPointError(4,0,0.002274824);
   gre->SetPoint(5,0.9000484,0.2356965);
   gre->SetPointError(5,0,0.002659209);
   gre->SetPoint(6,1.240548,0.5339174);
   gre->SetPointError(6,0,0.004002329);
   gre->SetPoint(7,2.311318,1);
   gre->SetPointError(7,0,0.005477417);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin3_224__24 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin3_224__24","Normal_ptSpect_not_rec_ChargeInt__bin3_2",100,0,2.529287);
   Normal_ptSpect_not_rec_ChargeInt__bin3_224__24->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_224__24->SetMaximum(1.10136);
   Normal_ptSpect_not_rec_ChargeInt__bin3_224__24->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_224__24->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin3_224__24->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin3_224__24->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin3_224);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("Normal_ptSpect_not_rec_ChargeInt__bin4_2");
   gre->SetTitle("Normal_ptSpect_not_rec_ChargeInt__bin4_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,0.13549,0.04411431);
   gre->SetPointError(0,0,0.001581575);
   gre->SetPoint(1,0.2811017,0.08879565);
   gre->SetPointError(1,0,0.00224386);
   gre->SetPoint(2,0.4274351,0.1321161);
   gre->SetPointError(2,0,0.002737019);
   gre->SetPoint(3,0.5761198,0.1522454);
   gre->SetPointError(3,0,0.002938137);
   gre->SetPoint(4,0.7260172,0.1701066);
   gre->SetPointError(4,0,0.003105707);
   gre->SetPoint(5,0.900407,0.2266954);
   gre->SetPointError(5,0,0.003585265);
   gre->SetPoint(6,1.240177,0.5108301);
   gre->SetPointError(6,0,0.005381932);
   gre->SetPoint(7,2.342766,1);
   gre->SetPointError(7,0,0.007530086);
   
   TH1F *Normal_ptSpect_not_rec_ChargeInt__bin4_225__25 = new TH1F("Normal_ptSpect_not_rec_ChargeInt__bin4_225__25","Normal_ptSpect_not_rec_ChargeInt__bin4_2",100,0,2.563494);
   Normal_ptSpect_not_rec_ChargeInt__bin4_225__25->SetMinimum(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_225__25->SetMaximum(1.10403);
   Normal_ptSpect_not_rec_ChargeInt__bin4_225__25->SetDirectory(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_225__25->SetStats(0);
   Normal_ptSpect_not_rec_ChargeInt__bin4_225__25->GetXaxis()->SetTitle("kT [GeV]");
   Normal_ptSpect_not_rec_ChargeInt__bin4_225__25->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(Normal_ptSpect_not_rec_ChargeInt__bin4_225);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_2");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_2");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,nan,nan);
   gre->SetPointError(0,0,nan);
   gre->SetPoint(1,nan,nan);
   gre->SetPointError(1,0,nan);
   gre->SetPoint(2,nan,nan);
   gre->SetPointError(2,0,nan);
   gre->SetPoint(3,nan,nan);
   gre->SetPointError(3,0,nan);
   gre->SetPoint(4,nan,nan);
   gre->SetPointError(4,0,nan);
   gre->SetPoint(5,nan,nan);
   gre->SetPointError(5,0,nan);
   gre->SetPoint(6,nan,nan);
   gre->SetPointError(6,0,nan);
   gre->SetPoint(7,nan,nan);
   gre->SetPointError(7,0,nan);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26","NormalWoA_ptSpect_not_rec_ChargeInt__bin0_2",100,nan,nan);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26->SetMinimum(nan);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26->SetMaximum(nan);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226__26->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin0_226);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_2");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,nan,0);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,nan,0);
   gre->SetPointError(1,0,0);
   gre->SetPoint(2,0.4885305,0.1428571);
   gre->SetPointError(2,0,0.1428571);
   gre->SetPoint(3,nan,0);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,0.7614236,0.2857143);
   gre->SetPointError(4,0,0.2020305);
   gre->SetPoint(5,nan,0);
   gre->SetPointError(5,0,0);
   gre->SetPoint(6,1.344912,0.8571429);
   gre->SetPointError(6,0,0.3499271);
   gre->SetPoint(7,2.909281,1);
   gre->SetPointError(7,0,0.3779645);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27","NormalWoA_ptSpect_not_rec_ChargeInt__bin1_2",100,nan,nan);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27->SetMaximum(1.515761);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227__27->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin1_227);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_2");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1311,0.05248894);
   gre->SetPointError(0,0,0.001703863);
   gre->SetPoint(1,0.281326,0.09568584);
   gre->SetPointError(1,0,0.002300513);
   gre->SetPoint(2,0.4288277,0.1431416);
   gre->SetPointError(2,0,0.002813739);
   gre->SetPoint(3,0.5763156,0.1693584);
   gre->SetPointError(3,0,0.003060583);
   gre->SetPoint(4,0.7244214,0.1911504);
   gre->SetPointError(4,0,0.003251535);
   gre->SetPoint(5,0.8990018,0.2422566);
   gre->SetPointError(5,0,0.003660485);
   gre->SetPoint(6,1.239966,0.5438606);
   gre->SetPointError(6,0,0.005484595);
   gre->SetPoint(7,2.364018,1);
   gre->SetPointError(7,0,0.007437051);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28","NormalWoA_ptSpect_not_rec_ChargeInt__bin2_2",100,0,2.58731);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28->SetMaximum(1.103102);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228__28->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin2_228);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_2");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#00ff00");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1316645,0.0523349);
   gre->SetPointError(0,0,0.0009115864);
   gre->SetPoint(1,0.2810261,0.09908064);
   gre->SetPointError(1,0,0.001254286);
   gre->SetPoint(2,0.428311,0.1399355);
   gre->SetPointError(2,0,0.001490617);
   gre->SetPoint(3,0.5767416,0.1663729);
   gre->SetPointError(3,0,0.001625337);
   gre->SetPoint(4,0.7261307,0.1850617);
   gre->SetPointError(4,0,0.001714196);
   gre->SetPoint(5,0.8997017,0.2481462);
   gre->SetPointError(5,0,0.001984979);
   gre->SetPoint(6,1.238921,0.5441655);
   gre->SetPointError(6,0,0.002939461);
   gre->SetPoint(7,2.327682,1);
   gre->SetPointError(7,0,0.00398476);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29","NormalWoA_ptSpect_not_rec_ChargeInt__bin3_2",100,0,2.547283);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29->SetMaximum(1.099241);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229__29->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin3_229);
   
   gre->Draw(" p");
   
   gre = new TGraphErrors(8);
   gre->SetName("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_2");
   gre->SetTitle("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_2");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff00ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(28);
   gre->SetPoint(0,0.1339816,0.04678643);
   gre->SetPointError(0,0,0.001428302);
   gre->SetPoint(1,0.2822353,0.09732275);
   gre->SetPointError(1,0,0.00206);
   gre->SetPoint(2,0.4292481,0.1371762);
   gre->SetPointError(2,0,0.00244568);
   gre->SetPoint(3,0.5772521,0.1629895);
   gre->SetPointError(3,0,0.002665875);
   gre->SetPoint(4,0.7259321,0.1743263);
   gre->SetPointError(4,0,0.002757031);
   gre->SetPoint(5,0.8997614,0.2343246);
   gre->SetPointError(5,0,0.003196458);
   gre->SetPoint(6,1.240527,0.5068021);
   gre->SetPointError(6,0,0.004700882);
   gre->SetPoint(7,2.375082,1);
   gre->SetPointError(7,0,0.006603286);
   
   TH1F *NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30 = new TH1F("NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30","NormalWoA_ptSpect_not_rec_ChargeInt__bin4_2",100,0,2.599192);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30->SetMinimum(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30->SetMaximum(1.102728);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30->SetDirectory(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30->SetStats(0);
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30->GetXaxis()->SetTitle("kT [GeV]");
   NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230__30->GetYaxis()->SetTitle("normalized counts");
   gre->SetHistogram(NormalWoA_ptSpect_not_rec_ChargeInt__bin4_230);
   
   gre->Draw(" p");
   
   pt = new TPaveText(0.01,0.9384156,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(0);
   text = pt->AddText("Normal_ptSpect_not_rec_ChargeInt__bin2_2");
   pt->Draw();
   c_3->Modified();
   c->cd();
  
// ------------>Primitives in pad: c_4
   c_4 = new TPad("c_4", "c_4",0.51,0.01,0.99,0.49);
   c_4->Draw();
   c_4->cd();
   c_4->Range(0,0,1,1);
   c_4->SetFillColor(0);
   c_4->SetBorderMode(0);
   c_4->SetBorderSize(2);
   c_4->SetFrameBorderMode(0);
   c_4->Modified();
   c->cd();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
