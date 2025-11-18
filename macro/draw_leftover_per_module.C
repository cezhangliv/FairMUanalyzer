// File: draw_leftover_per_module.C
// Usage (batch):
//   root -l -b -q 'draw_leftover_per_module.C("input.root","leftoverhits")'
// leftoverhits_stat0.{png,pdf}, leftoverhits_stat1.{png,pdf}, leftoverhits_stat2.{png,pdf}

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <iostream>
#include <sstream>

void draw_one_stat(TFile* fin, int stat, const TString& outprefix, bool logy=false) {
  
  TString cname = Form("c_stat%d", stat);
  TCanvas* c = new TCanvas(cname, cname, 1800, 1100);
  c->Divide(3,2,0.005,0.005); // 3列2行

  
  //gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFontSize(0.045);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.045);

  
  for (int m=0; m<6; ++m) {
    c->cd(m+1);
    gPad->SetTicks(1,1);
    gPad->SetGrid(0,0);
    if (logy) gPad->SetLogy();

    TString hname = Form("case_h1d_LeftoverhitsPerModule_stat%d_Module_%d_t1mmm", stat, m);
    TH1* h = nullptr;
    fin->GetObject(hname, h);

    if (!h) {
    
      TPaveText *pt = new TPaveText(0.15,0.40,0.85,0.60,"NDC");
      pt->SetFillColor(0);
      pt->SetBorderSize(0);
      pt->AddText(Form("Missing histogram:"));
      pt->AddText(hname);
      pt->Draw();
      continue;
    }

    
    h->SetLineWidth(2);
    h->SetTitle(Form("Leftover hits per module  |  stat%d  module%d", stat, m));
    h->GetXaxis()->SetTitle("Leftover hits");
    h->GetYaxis()->SetTitle("Entries");

    
    if (logy) {
      double maxv = h->GetMaximum();
      if (maxv <= 0) maxv = 1.0;
      h->SetMinimum(0.5);
      h->SetMaximum(maxv*10.);
    }

    h->Draw("HIST");
    cout<<stat<<" "<<m<<": "<<h->GetBinContent(1)<<endl;
    cout<<stat<<" "<<m<<": "<<h->GetBinContent(2)<<endl;
  
    //latex.DrawLatex(0.14, 0.92, Form("stat%d", stat));
    //latex.DrawLatex(0.80, 0.92, Form("Module %d", m));
  }

  
  TString out_png = Form("%s_stat%d.png", outprefix.Data(), stat);
  TString out_pdf = Form("%s_stat%d.pdf", outprefix.Data(), stat);
  c->SaveAs(out_png);
  c->SaveAs(out_pdf);
}

void draw_leftover_per_module(
  const char* infile="../result/FairMUanalyzer_01Oct25_GA1001_single_muon_interaction_1_CDbugfix11July25_run8_muedaq04-1750227094-1750228937_MF1_output.root", 
  const char* outprefix="leftoverhits") {
  
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "[ERROR] Cannot open file: " << infile << std::endl;
    return;
  }

  
  gROOT->SetBatch(kTRUE);
  //gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.04);
  //gStyle->SetTitleX(0.12);

  
  for (int stat=0; stat<3; ++stat) {
    draw_one_stat(fin, stat, outprefix, /*logy=*/false);
  }

  fin->Close();
  delete fin;
  std::cout << "[INFO] Done. Outputs: " << outprefix
            << "_stat{0,1,2}.{png,pdf}" << std::endl;
}
