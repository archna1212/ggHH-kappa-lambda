#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

struct HistSet
{
    TH1D *mbb;
    TH1D *mgg;
    TH1D *mhh;
    TH1D *ptbb;
    TH1D *ptgg;
    TH1D *drbb;
    TH1D *drgg;
};

TLorentzVector makeP4(
    double pt,
    double eta,
    double phi,
    double mass)
{
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, mass);
    return v;
}

void analyze(
    const char *file,
    HistSet &h,
    double kl)
{
    TChain chain("Delphes");
    chain.Add(file);

    ExRootTreeReader reader(&chain);

    TClonesArray *jets =
        reader.UseBranch("Jet");

    TClonesArray *photons =
        reader.UseBranch("Photon");

    Long64_t entries = reader.GetEntries();

    Long64_t passJets = 0;
    Long64_t passB = 0;
    Long64_t passPho = 0;
    Long64_t passMbb = 0;
    Long64_t passMgg = 0;
    Long64_t passFinal = 0;

    for(Long64_t entry = 0; entry < entries; ++entry)
    {
        reader.ReadEntry(entry);

        vector<Jet*> goodJets;
        vector<Jet*> bjets;
        vector<Photon*> goodPho;

        // --------------------------------------------------
        // Jets
        // --------------------------------------------------

        for(int i = 0; i < jets->GetEntriesFast(); ++i)
        {
            Jet *j = (Jet*)jets->At(i);

            if(j->PT < 25.0)
                continue;

            if(fabs(j->Eta) > 2.5)
                continue;

            goodJets.push_back(j);

            if(j->BTag == 1)
                bjets.push_back(j);
        }

        if(goodJets.size() < 2)
            continue;

        passJets++;

        if(bjets.size() < 2)
            continue;

        passB++;

        // --------------------------------------------------
        // Photons
        // --------------------------------------------------

        for(int i = 0; i < photons->GetEntriesFast(); ++i)
        {
            Photon *p = (Photon*)photons->At(i);

            if(p->PT < 25.0)
                continue;

            if(fabs(p->Eta) > 2.5)
                continue;

            goodPho.push_back(p);
        }

        if(goodPho.size() < 2)
            continue;

        passPho++;

        // --------------------------------------------------
        // Select b-jet pair closest to mH = 125 GeV
        // --------------------------------------------------

        double bestDiff = 1e9;
        int bestI = -1;
        int bestJ = -1;

        for(size_t i = 0; i < bjets.size(); ++i)
        {
            for(size_t j = i + 1; j < bjets.size(); ++j)
            {
                TLorentzVector bi = makeP4(
                    bjets[i]->PT,
                    bjets[i]->Eta,
                    bjets[i]->Phi,
                    bjets[i]->Mass
                );

                TLorentzVector bj = makeP4(
                    bjets[j]->PT,
                    bjets[j]->Eta,
                    bjets[j]->Phi,
                    bjets[j]->Mass
                );

                double mbb = (bi + bj).M();
                double diff = fabs(mbb - 125.0);

                if(diff < bestDiff)
                {
                    bestDiff = diff;
                    bestI = (int)i;
                    bestJ = (int)j;
                }
            }
        }

        if(bestI < 0 || bestJ < 0)
            continue;

        TLorentzVector b1 = makeP4(
            bjets[bestI]->PT,
            bjets[bestI]->Eta,
            bjets[bestI]->Phi,
            bjets[bestI]->Mass
        );

        TLorentzVector b2 = makeP4(
            bjets[bestJ]->PT,
            bjets[bestJ]->Eta,
            bjets[bestJ]->Phi,
            bjets[bestJ]->Mass
        );

        // --------------------------------------------------
        // Leading two photons
        // --------------------------------------------------

        sort(
            goodPho.begin(),
            goodPho.end(),
            [](Photon *a, Photon *b)
            {
                return a->PT > b->PT;
            }
        );

        // Photons are massless.
        TLorentzVector g1 = makeP4(
            goodPho[0]->PT,
            goodPho[0]->Eta,
            goodPho[0]->Phi,
            0.0
        );

        TLorentzVector g2 = makeP4(
            goodPho[1]->PT,
            goodPho[1]->Eta,
            goodPho[1]->Phi,
            0.0
        );

        TLorentzVector Hbb = b1 + b2;
        TLorentzVector Hgg = g1 + g2;
        TLorentzVector HH = Hbb + Hgg;

        // --------------------------------------------------
        // Higgs-mass selections
        // --------------------------------------------------

        double mbb = Hbb.M();
        double mgg = Hgg.M();

        if(mbb < 100.0 || mbb > 150.0)
            continue;

        passMbb++;

        if(mgg < 120.0 || mgg > 130.0)
            continue;

        passMgg++;
        passFinal++;

        // --------------------------------------------------
        // Final selected observables
        // --------------------------------------------------

        h.mbb->Fill(mbb);
        h.mgg->Fill(mgg);
        h.mhh->Fill(HH.M());
        h.ptbb->Fill(Hbb.Pt());
        h.ptgg->Fill(Hgg.Pt());
        h.drbb->Fill(b1.DeltaR(b2));
        h.drgg->Fill(g1.DeltaR(g2));
    }

    cout << endl;
    cout << "====================================" << endl;
    cout << "KL = " << kl << endl;
    cout << "Generated entries = " << entries << endl;
    cout << ">=2 good jets    = " << passJets << endl;
    cout << ">=2 b jets       = " << passB << endl;
    cout << ">=2 photons      = " << passPho << endl;
    cout << "m_bb window      = " << passMbb << endl;
    cout << "m_gg window      = " << passMgg << endl;
    cout << "FINAL SELECTED   = " << passFinal << endl;
    cout << "====================================" << endl;
}

void normalize(TH1D *h)
{
    double integral = h->Integral();

    if(integral > 0)
        h->Scale(1.0 / integral);
}

void drawFive(
    TH1D *hm15,
    TH1D *h0,
    TH1D *h1,
    TH1D *h2,
    TH1D *h5,
    const char *xtitle,
    const char *outname)
{
    TCanvas *c =
        new TCanvas("c", "c", 900, 700);

    normalize(hm15);
    normalize(h0);
    normalize(h1);
    normalize(h2);
    normalize(h5);
    hm15->SetLineColor(kBlack);   hm15->SetLineStyle(1);
    h0->SetLineColor(kBlue+1);    h0->SetLineStyle(2);
    h1->SetLineColor(kRed+1);     h1->SetLineStyle(1);
    h2->SetLineColor(kGreen+2);    h2->SetLineStyle(2);
    h5->SetLineColor(kMagenta+1); h5->SetLineStyle(1);

    hm15->SetLineWidth(3);
    h0->SetLineWidth(3);
    h1->SetLineWidth(3);
    h2->SetLineWidth(3);
    h5->SetLineWidth(3);


    hm15->SetStats(0);
    h0->SetStats(0);
    h1->SetStats(0);
    h2->SetStats(0);
    h5->SetStats(0);

    hm15->GetXaxis()->SetTitle(xtitle);
    hm15->GetYaxis()->SetTitle("Normalized entries");

    double ymax = max(
        max(hm15->GetMaximum(), h0->GetMaximum()),
        max(h1->GetMaximum(), max(h2->GetMaximum(), h5->GetMaximum()))
    );

    hm15->SetMaximum(1.25 * ymax);

    hm15->Draw("hist");
    h0->Draw("hist same");
    h1->Draw("hist same");
    h2->Draw("hist same");
    h5->Draw("hist same");

    TLegend *leg =
        new TLegend(0.60, 0.63, 0.88, 0.88);

    leg->AddEntry(hm15, "#kappa_{#lambda}=-1.5", "l");
    leg->AddEntry(h0, "#kappa_{#lambda}=0", "l");
    leg->AddEntry(h1, "#kappa_{#lambda}=1", "l");
    leg->AddEntry(h2, "#kappa_{#lambda}=2", "l");
    leg->AddEntry(h5, "#kappa_{#lambda}=5", "l");

    leg->Draw();

    c->SaveAs(outname);

    delete leg;
    delete c;
}

void drawFour(
    TH1D *h0,
    TH1D *h1,
    TH1D *h2,
    TH1D *h5,
    const char *xtitle,
    const char *outname)
{
    TCanvas *c =
        new TCanvas("c", "c", 900, 700);

    normalize(h0);
    normalize(h1);
    normalize(h2);
    normalize(h5);

    h0->SetLineWidth(3);
    h1->SetLineWidth(3);
    h2->SetLineWidth(3);
    h5->SetLineWidth(3);


    h0->SetStats(0);
    h1->SetStats(0);
    h2->SetStats(0);
    h5->SetStats(0);

    h0->GetXaxis()->SetTitle(xtitle);
    h0->GetYaxis()->SetTitle("Normalized entries");

    double ymax = max(
        max(h0->GetMaximum(), h1->GetMaximum()),
        max(h2->GetMaximum(), h5->GetMaximum())
    );

    h0->SetMaximum(1.25 * ymax);

    h0->Draw("hist");
    h1->Draw("hist same");
    h2->Draw("hist same");
    h5->Draw("hist same");

    TLegend *leg =
        new TLegend(0.65, 0.68, 0.88, 0.88);

    leg->AddEntry(h0, "#kappa_{#lambda}=0", "l");
    leg->AddEntry(h1, "#kappa_{#lambda}=1", "l");
    leg->AddEntry(h2, "#kappa_{#lambda}=2", "l");
    leg->AddEntry(h5, "#kappa_{#lambda}=5", "l");

    leg->Draw();

    c->SaveAs(outname);

    delete leg;
    delete c;
}

HistSet makeHistSet(const char *prefix)
{
    HistSet h;

    h.mbb = new TH1D(
        Form("%s_mbb", prefix),
        "",
        40, 50, 200
    );

    h.mgg = new TH1D(
        Form("%s_mgg", prefix),
        "",
        40, 80, 160
    );

    h.mhh = new TH1D(
        Form("%s_mhh", prefix),
        "",
        50, 200, 1200
    );

    h.ptbb = new TH1D(
        Form("%s_ptbb", prefix),
        "",
        40, 0, 500
    );

    h.ptgg = new TH1D(
        Form("%s_ptgg", prefix),
        "",
        40, 0, 500
    );

    h.drbb = new TH1D(
        Form("%s_drbb", prefix),
        "",
        40, 0, 4
    );

    h.drgg = new TH1D(
        Form("%s_drgg", prefix),
        "",
        40, 0, 4
    );

    return h;
}

void differential_kl_scan_m15()
{
    const int nKL = 5;

    double klValues[nKL] = {-1.5, 0, 1, 2, 5};

    const char *files[nKL] =
    {
        "Events/run_klm15_forcedbbgg_10k_delphes.root",
        "Events/run_kl0_forcedbbgg_10k_delphes.root",
        "Events/run_kl1_forcedbbgg_10k_delphes.root",
        "Events/run_kl2_forcedbbgg_10k_delphes.root",
        "Events/run_kl5_forcedbbgg_10k_delphes.root"
    };

    HistSet h[nKL];

    for(int i = 0; i < nKL; ++i)
    {
        h[i] =
            makeHistSet(Form("kl_%g", klValues[i]));
    }

    // ------------------------------------------------------
    // Analyze all four samples
    // ------------------------------------------------------

    for(int i = 0; i < nKL; ++i)
    {
        analyze(
            files[i],
            h[i],
            klValues[i]
        );
    }

    // ------------------------------------------------------
    // Draw the seven normalized detector-level shapes
    // ------------------------------------------------------

    drawFive(
        h[0].mbb,
        h[1].mbb,
        h[2].mbb,
        h[3].mbb,
        h[4].mbb,
        "m_{bb} [GeV]",
        "analysis_plots/scan_mbb_kl.png"
    );

    drawFive(
        h[0].mgg,
        h[1].mgg,
        h[2].mgg,
        h[3].mgg,
        h[4].mgg,
        "m_{#gamma#gamma} [GeV]",
        "analysis_plots/scan_mgg_kl.png"
    );

    drawFive(
        h[0].mhh,
        h[1].mhh,
        h[2].mhh,
        h[3].mhh,
        h[4].mhh,
        "m_{HH} [GeV]",
        "analysis_plots/scan_mhh_kl.png"
    );

    drawFive(
        h[0].ptbb,
        h[1].ptbb,
        h[2].ptbb,
        h[3].ptbb,
        h[4].ptbb,
        "p_{T}(H_{bb}) [GeV]",
        "analysis_plots/scan_ptbb_kl.png"
    );

    drawFive(
        h[0].ptgg,
        h[1].ptgg,
        h[2].ptgg,
        h[3].ptgg,
        h[4].ptgg,
        "p_{T}(H_{#gamma#gamma}) [GeV]",
        "analysis_plots/scan_ptgg_kl.png"
    );

    drawFive(
        h[0].drbb,
        h[1].drbb,
        h[2].drbb,
        h[3].drbb,
        h[4].drbb,
        "#DeltaR_{bb}",
        "analysis_plots/scan_drbb_kl.png"
    );

    drawFive(
        h[0].drgg,
        h[1].drgg,
        h[2].drgg,
        h[3].drgg,
        h[4].drgg,
        "#DeltaR_{#gamma#gamma}",
        "analysis_plots/scan_drgg_kl.png"
    );

    cout << endl;
    cout << "==================================================" << endl;
    cout << " FIVE-POINT DETECTOR SHAPE SCAN COMPLETE" << endl;
    cout << "==================================================" << endl;
}
