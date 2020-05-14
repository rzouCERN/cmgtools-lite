void makePlot(TString TVet, TString BVet, TString Year)
{
// *** Constants ***
  const int nBins = 100;
  const double Offset = 1.2;
  const TString fileName = "2lss_diff_Top-tagged/priority/TTHnobb_fxfx_Friend_DefaultModule_"+TVet+BVet+"Veto_"+Year+".root";

// *** Files ***
  TFile f(fileName);
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_DefaultModule_NoTopMediumBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_NoTopNoBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_NoTopMediumBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_NoTopLooseBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_TopNoBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_TopMediumBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_TopLooseBottomVeto.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_NoTopNoBottomVeto_NoPTRes.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_NoTopMediumBottomVeto_NoPTRes.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_NoTopLooseBottomVeto_NoPTRes.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_TopNoBottomVeto_NoPTRes.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_TopMediumBottomVeto_NoPTRes.root");
//  TFile f("2lss_diff_Top-tagged/TTHnobb_fxfx_Friend_TopLooseBottomVeto_NoPTRes.root");
  TTree *ft = (TTree*) f.Get("Friends");

// *** Declare variables ***
  // reconstruction loop dependent vars
  float Hreco_pTVisPlusNu;
  float Hreco_pTHvis;
  float Hreco_pTTrueGen;
  float Hreco_pTTrueGenPlusNu;
  float Hreco_pTHgen;
  float Hreco_nQFromWFromH;

  // reconstruction loop dependent vars
  ft->SetBranchAddress("Hreco_pTVisPlusNu",&Hreco_pTVisPlusNu);
  ft->SetBranchAddress("Hreco_pTHvis",&Hreco_pTHvis);
  ft->SetBranchAddress("Hreco_pTTrueGen",&Hreco_pTTrueGen);
  ft->SetBranchAddress("Hreco_pTTrueGenPlusNu",&Hreco_pTTrueGenPlusNu);
  ft->SetBranchAddress("Hreco_pTHgen",&Hreco_pTHgen);
  ft->SetBranchAddress("Hreco_nQFromWFromH",&Hreco_nQFromWFromH);

// *** Declare histograms *** 
  TH1D *hst_pTHvis = new TH1D("hst_pTHvis","hst_pTHvis",nBins,-0.5,599.5);
  TH1D *hst_pTVisPlusNu = new TH1D("pTVisPlusNu","pTVisPlusNu",nBins,-0.5,599.5);
  TH1D *hst_pTHgen = new TH1D("hst_pTHgen","hst_pTHgen",nBins,-0.5,599.5);
  TH1D *hst_pTTrueGen = new TH1D("hst_pTTrueGen","hst_pTTrueGen",nBins,-0.5,599.5);
  TH1D *hst_pTTrueGenPlusNu = new TH1D("pTTrueGenPlusNu","pTTrueGenPlusNu",nBins,-0.5,599.5);
  
// *** Declare counts ***
  int nEvents = 0;
  int nQFromWFromH = 0;
  for (int i=0; i<ft->GetEntries(); i++)
  {
    ft->GetEntry(i);
    nEvents ++;
    if (Hreco_nQFromWFromH != 2) nQFromWFromH ++;

    // Plot pT dists if entry != -99
    if (Hreco_pTHvis > 0) hst_pTHvis->Fill(Hreco_pTHvis);
    if (Hreco_pTVisPlusNu > 0) hst_pTVisPlusNu->Fill(Hreco_pTVisPlusNu);
    if (Hreco_pTHgen > 0) hst_pTHgen->Fill(Hreco_pTHgen);
    if (Hreco_pTTrueGen > 0) hst_pTTrueGen->Fill(Hreco_pTTrueGen);
    if (Hreco_pTTrueGenPlusNu > 0) hst_pTTrueGenPlusNu->Fill(Hreco_pTTrueGenPlusNu);
  }

  // Scale all distributions
  hst_pTHvis->Scale(1/hst_pTHvis->GetEntries());
  hst_pTVisPlusNu->Scale(1/hst_pTVisPlusNu->GetEntries());
  hst_pTHgen->Scale(1/hst_pTHgen->GetEntries());
  hst_pTTrueGen->Scale(1/hst_pTTrueGen->GetEntries());
  hst_pTTrueGenPlusNu->Scale(1/hst_pTTrueGenPlusNu->GetEntries());

  // Print counts
  cout << "N Events:    " << nEvents << endl;
  cout << "Passed reco: " << hst_pTHvis->GetEntries() << endl;

  gSystem->Exec("mkdir 1DDistPlots");
  TString VetoString = TVet+BVet+"Veto";

  TCanvas *can2 = new TCanvas();
  hst_pTHvis->SetLineColor(kBlue);
  hst_pTHgen->SetLineColor(kRed);
  hst_pTVisPlusNu->SetLineColor(kGreen);
  hst_pTTrueGen->SetLineColor(kMagenta);
  hst_pTTrueGenPlusNu->SetLineColor(kBlack);
  hst_pTHvis->GetXaxis()->SetTitle("pT");
  hst_pTHvis->GetYaxis()->SetTitle("a.u.");
  hst_pTHvis->GetXaxis()->CenterTitle();
  hst_pTHvis->GetYaxis()->CenterTitle();
  hst_pTHvis->GetXaxis()->SetTitleOffset(Offset);
  hst_pTHvis->GetYaxis()->SetTitleOffset(Offset);
  hst_pTHvis->GetYaxis()->SetRangeUser(0,0.07);
  hst_pTHvis->SetTitle("pT Distributions ("+Year+")");
  hst_pTHvis->SetStats(0);
  hst_pTHvis->Draw("HIST");
  hst_pTHgen->Draw("HIST SAME");
  hst_pTVisPlusNu->Draw("HIST SAME");
  hst_pTTrueGen->Draw("HIST SAME");
  hst_pTTrueGenPlusNu->Draw("HIST SAME");
  TLatex tex2(.4,.85,VetoString);
  tex2.SetTextSize(0.03);
  tex2.SetNDC(kTRUE);
  tex2.Draw();

  TLegend *leg2 = new TLegend(0.50,0.55,0.85,0.80);
  leg2->AddEntry(hst_pTHvis,"Reco Higgs pT");
  leg2->AddEntry(hst_pTHgen,"Gen Higgs pT");
  leg2->AddEntry(hst_pTVisPlusNu,"Reco + Gen(Nu)");
  leg2->AddEntry(hst_pTTrueGen,"Gen(q1) + Gen(q2) + Gen(l)");
  leg2->AddEntry(hst_pTTrueGenPlusNu,"Gen(q1) + Gen(q2) + Gen(l) + Gen(nu)");
  leg2->Draw();
  can2->SaveAs("1DDistPlots/plot_PTDists_"+TVet+BVet+Year+".png");
}

void plot1DDists()
{
  makePlot("NoTop", "MediumBottom", "2016");
  makePlot("NoTop", "MediumBottom", "2017");
  makePlot("NoTop", "MediumBottom", "2018");
}
