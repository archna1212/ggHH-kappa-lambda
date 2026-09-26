#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "ExRootTreeReader.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

struct Result
{
    double kl;
    Long64_t generated;
    Long64_t selected;
    TH1D *shape;
    TH1D *expected;
};

Result analyze(
    double kl,
    const char *filename,
    double physicalYield)
{
    Result r;
    r.kl = kl;

    TChain chain("Delphes");
    chain.Add(filename);

    ExRootTreeReader reader(&chain);

    TClonesArray *jets =
        reader.UseBranch("Jet");

    TClonesArray *photons =
        reader.UseBranch("Photon");

    r.generated = reader.GetEntries();
    r.selected = 0;

    r.shape = new TH1D(
        Form("mhh_shape_kl_%g", kl),
        "",
        40, 200.0, 1400.0
    );

    r.expected = new TH1D(
        Form("mhh_expected_kl_%g", kl),
        "",
        40, 200.0, 1400.0
    );

    // --------------------------------------------------
    // Exact final selection from cutflow_kl_scan.C
    // --------------------------------------------------

    for(Long64_t entry = 0;
        entry < r.generated;
        ++entry)
    {
        reader.ReadEntry(entry);

        vector<Jet*> goodJets;
        vector<Jet*> bjets;
        vector<Photon*> goodPhotons;

        // -----------------------------
        // Jets
        // -----------------------------

        for(int j = 0;
            j < jets->GetEntries();
            ++j)
        {
            Jet *jet = (Jet*)jets->At(j);

            if(jet->PT > 25.0 &&
               fabs(jet->Eta) < 2.5)
            {
                goodJets.push_back(jet);

                if(jet->BTag)
                    bjets.push_back(jet);
            }
        }

        if(goodJets.size() < 2)
            continue;

        if(bjets.size() < 2)
            continue;

        // -----------------------------
        // Photons
        // -----------------------------

        for(int j = 0;
            j < photons->GetEntries();
            ++j)
        {
            Photon *photon =
                (Photon*)photons->At(j);

            if(photon->PT > 25.0 &&
               fabs(photon->Eta) < 2.5)
            {
                goodPhotons.push_back(photon);
            }
        }

        if(goodPhotons.size() < 2)
            continue;

        // -----------------------------
        // Best bb pair:
        // closest to 125 GeV
        // -----------------------------

        double bestDifference = 1e9;

        int bestI = -1;
        int bestJ = -1;

        for(int a = 0;
            a < (int)bjets.size();
            ++a)
        {
            for(int b = a + 1;
                b < (int)bjets.size();
                ++b)
            {
                TLorentzVector b1, b2;

                b1.SetPtEtaPhiM(
                    bjets[a]->PT,
                    bjets[a]->Eta,
                    bjets[a]->Phi,
                    bjets[a]->Mass
                );

                b2.SetPtEtaPhiM(
                    bjets[b]->PT,
                    bjets[b]->Eta,
                    bjets[b]->Phi,
                    bjets[b]->Mass
                );

                double mbbPair =
                    (b1 + b2).M();

                double diff =
                    fabs(mbbPair - 125.0);

                if(diff < bestDifference)
                {
                    bestDifference = diff;
                    bestI = a;
                    bestJ = b;
                }
            }
        }

        if(bestI < 0 || bestJ < 0)
            continue;

        TLorentzVector b1, b2;

        b1.SetPtEtaPhiM(
            bjets[bestI]->PT,
            bjets[bestI]->Eta,
            bjets[bestI]->Phi,
            bjets[bestI]->Mass
        );

        b2.SetPtEtaPhiM(
            bjets[bestJ]->PT,
            bjets[bestJ]->Eta,
            bjets[bestJ]->Phi,
            bjets[bestJ]->Mass
        );

        TLorentzVector Hbb = b1 + b2;

        // -----------------------------
        // First two selected photons
        // -----------------------------

        TLorentzVector g1, g2;

        g1.SetPtEtaPhiM(
            goodPhotons[0]->PT,
            goodPhotons[0]->Eta,
            goodPhotons[0]->Phi,
            0.0
        );

        g2.SetPtEtaPhiM(
            goodPhotons[1]->PT,
            goodPhotons[1]->Eta,
            goodPhotons[1]->Phi,
            0.0
        );

        TLorentzVector Hgg = g1 + g2;

        double mbb = Hbb.M();
        double mgg = Hgg.M();

        // -----------------------------
        // Higgs mass windows
        // -----------------------------

        if(!(mbb > 100.0 && mbb < 150.0))
            continue;

        if(!(mgg > 120.0 && mgg < 130.0))
            continue;

        // -----------------------------
        // Final HH mass
        // -----------------------------

        TLorentzVector HH = Hbb + Hgg;

        r.shape->Fill(HH.M());

        r.selected++;
    }

    // --------------------------------------------------
    // Normalize shape
    // --------------------------------------------------

    double integral = r.shape->Integral();

    if(integral > 0.0)
    {
        r.shape->Scale(1.0 / integral);

        // Physical expected events/bin
        for(int bin = 1;
            bin <= r.shape->GetNbinsX();
            ++bin)
        {
            r.expected->SetBinContent(
                bin,
                r.shape->GetBinContent(bin)
                * physicalYield
            );
        }
    }

    return r;
}

void drawShape(
    Result &rm15,
    Result &r0,
    Result &r1,
    Result &r2, Result &r5)
{
    TCanvas *c =
        new TCanvas(
            "c_shape",
            "Final mHH shapes",
            1600,
            1200
        );

    rm15.shape->SetStats(0);
    r0.shape->SetStats(0);
    r1.shape->SetStats(0);
    r2.shape->SetStats(0);
    r5.shape->SetStats(0);

    rm15.shape->SetLineWidth(3);
    r0.shape->SetLineWidth(3);
    r1.shape->SetLineWidth(3);
    r2.shape->SetLineWidth(3);
    r5.shape->SetLineWidth(3);

    rm15.shape->SetLineColor(kBlack);
    r0.shape->SetLineColor(kBlue);
    r1.shape->SetLineColor(kRed);
    r2.shape->SetLineColor(kGreen+2);
    r5.shape->SetLineColor(kMagenta);

    rm15.shape->SetLineStyle(1);
    r0.shape->SetLineStyle(2);
    r1.shape->SetLineStyle(1);
    r2.shape->SetLineStyle(2);
    r5.shape->SetLineStyle(1);

    r0.shape->GetXaxis()->SetTitle(
        "m_{HH} [GeV]"
    );

    r0.shape->GetYaxis()->SetTitle(
        "Normalized entries"
    );

    double ymax = max(
        max(
            rm15.shape->GetMaximum(),
            r0.shape->GetMaximum()
        ),
        max(
            max(
                r1.shape->GetMaximum(),
                r2.shape->GetMaximum()
            ),
            r5.shape->GetMaximum()
        )
    );

    r0.shape->SetMaximum(1.25 * ymax);

    rm15.shape->Draw("hist");
    r0.shape->Draw("hist same");
    r1.shape->Draw("hist same");
    r2.shape->Draw("hist same");
    r5.shape->Draw("hist same");

    TLegend *leg =
        new TLegend(
            0.64, 0.67,
            0.88, 0.88
        );

    leg->AddEntry(
        rm15.shape,
        "#kappa_{#lambda}=-1.5",
        "l"
    );

    leg->AddEntry(
        r0.shape,
        "#kappa_{#lambda}=0",
        "l"
    );

    leg->AddEntry(
        r1.shape,
        "#kappa_{#lambda}=1",
        "l"
    );

    leg->AddEntry(
        r2.shape,
        "#kappa_{#lambda}=2",
        "l"
    );

    leg->AddEntry(
        r5.shape,
        "#kappa_{#lambda}=5",
        "l"
    );

    leg->Draw();

    c->SaveAs(
        "analysis_plots/final_mhh_shapes_kl_m15.png"
    );

    delete leg;
    delete c;
}

void drawExpected(
    Result &rm15,
    Result &r0,
    Result &r1,
    Result &r2, Result &r5)
{
    TCanvas *c =
        new TCanvas(
            "c_expected",
            "Expected mHH distribution",
            1600,
            1200
        );

    rm15.expected->SetStats(0);
    r0.expected->SetStats(0);
    r1.expected->SetStats(0);
    r2.expected->SetStats(0);
    r5.expected->SetStats(0);

    rm15.expected->SetLineWidth(3);
    r0.expected->SetLineWidth(3);
    r1.expected->SetLineWidth(3);
    r2.expected->SetLineWidth(3);
    r5.expected->SetLineWidth(3);

    rm15.expected->SetLineColor(kBlack);
    r0.expected->SetLineColor(kBlue);
    r1.expected->SetLineColor(kRed);
    r2.expected->SetLineColor(kGreen+2);
    r5.expected->SetLineColor(kMagenta);

    rm15.expected->SetLineStyle(1);
    r0.expected->SetLineStyle(2);
    r1.expected->SetLineStyle(1);
    r2.expected->SetLineStyle(2);
    r5.expected->SetLineStyle(1);

    r0.expected->GetXaxis()->SetTitle(
        "m_{HH} [GeV]"
    );

    r0.expected->GetYaxis()->SetTitle(
        "Expected events / bin"
    );

    double ymax = max(
        max(
            rm15.expected->GetMaximum(),
            r0.expected->GetMaximum()
        ),
        max(
            max(
                r1.expected->GetMaximum(),
                r2.expected->GetMaximum()
            ),
            r5.expected->GetMaximum()
        )
    );

    r0.expected->SetMaximum(1.25 * ymax);

    rm15.expected->Draw("hist");
    r0.expected->Draw("hist same");
    r1.expected->Draw("hist same");
    r2.expected->Draw("hist same");
    r5.expected->Draw("hist same");

    TLegend *leg =
        new TLegend(
            0.64, 0.67,
            0.88, 0.88
        );

    leg->AddEntry(
        rm15.expected,
        "#kappa_{#lambda}=-1.5",
        "l"
    );

    leg->AddEntry(
        r0.expected,
        "#kappa_{#lambda}=0",
        "l"
    );

    leg->AddEntry(
        r1.expected,
        "#kappa_{#lambda}=1",
        "l"
    );

    leg->AddEntry(
        r2.expected,
        "#kappa_{#lambda}=2",
        "l"
    );

    leg->AddEntry(
        r5.expected,
        "#kappa_{#lambda}=5",
        "l"
    );

    leg->Draw();

    c->SaveAs(
        "analysis_plots/final_mhh_expected_kl_m15.png"
    );

    delete leg;
    delete c;
}

void final_mhh_templates_kl_m15_pub()
{
    // ----------------------------------------------
    // Four benchmark points
    // ----------------------------------------------

    double klValues[5] = {-1.5, 0, 1, 2, 5};

    const char *files[5] =
    {
        "Events/run_klm15_forcedbbgg_10k_delphes.root",
        "Events/run_kl0_forcedbbgg_10k_delphes.root",
        "Events/run_kl1_forcedbbgg_10k_delphes.root",
        "Events/run_kl2_forcedbbgg_10k_delphes.root",
        "Events/run_kl5_forcedbbgg_10k_delphes.root"
    };

    // ----------------------------------------------
    // Inclusive HH cross sections [pb]
    // ----------------------------------------------

    double xs[5] =
    {
        0.070073,
        0.030540,
        0.014504,
        0.006786,
        0.033068
    };

    // ----------------------------------------------
    // Detector efficiencies
    // FINAL / 10000
    // ----------------------------------------------

    double eff[5] =
    {
        0.08760,
        0.0942,
        0.1043,
        0.1034,
        0.0703
    };

    // ----------------------------------------------
    // Branching ratio
    // ----------------------------------------------

    double BRbb = 0.5824;
    double BRgg = 0.00227;

    double BRbbgg =
        2.0 * BRbb * BRgg;

    // ----------------------------------------------
    // Luminosity
    // ----------------------------------------------

    double lumi = 3000.0;

    double yields[5];

    for(int i = 0; i < 5; ++i)
    {
        // pb * fb^-1 * 1000 = events
        yields[i] =
            xs[i]
            * BRbbgg
            * eff[i]
            * lumi
            * 1000.0;
    }

    // ----------------------------------------------
    // Analyze
    // ----------------------------------------------

    Result rm15 =
        analyze(
            -1.5,
            files[0],
            yields[0]
        );

    Result r0 =
        analyze(
            0,
            files[1],
            yields[1]
        );

    Result r1 =
        analyze(
            1,
            files[2],
            yields[2]
        );

    Result r2 =
        analyze(
            2,
            files[3],
            yields[3]
        );

    Result r5 =
        analyze(
            5,
            files[4],
            yields[4]
        );

    // ----------------------------------------------
    // Print summary
    // ----------------------------------------------

    cout << endl;
    cout << "============================================================"
         << endl;

    cout << " FINAL-SELECTED mHH TEMPLATES"
         << endl;

    cout << "============================================================"
         << endl;

    cout << left
         << setw(8)  << "KL"
         << setw(12) << "Generated"
         << setw(12) << "Selected"
         << setw(15) << "Efficiency"
         << setw(18) << "Physical yield"
         << endl;

    cout << string(65, '-') << endl;

    Result results[5] =
    {
        rm15, r0, r1, r2, r5
    };

    for(int i = 0; i < 5; ++i)
    {
        double efficiency =
            100.0
            * results[i].selected
            / (double)results[i].generated;

        cout << left
             << setw(8)  << results[i].kl
             << setw(12) << results[i].generated
             << setw(12) << results[i].selected
             << setw(15) << efficiency
             << setw(18) << yields[i]
             << endl;
    }

    cout << endl;
    cout << "BR(HH -> bb gamma gamma) = "
         << BRbbgg << endl;

    cout << "Luminosity = "
         << lumi << " fb^-1" << endl;

    cout << "============================================================"
         << endl;

    // ----------------------------------------------
    // Draw
    // ----------------------------------------------

    drawShape(rm15, r0, r1, r2, r5);
    drawExpected(rm15, r0, r1, r2, r5);

    // ----------------------------------------------
    // Save templates
    // ----------------------------------------------

    TFile out(
        "analysis_plots/final_mhh_templates_kl_m15.root",
        "RECREATE"
    );

    rm15.shape->Write();
    r0.shape->Write();
    r1.shape->Write();
    r2.shape->Write();
    r5.shape->Write();

    rm15.expected->Write();
    r0.expected->Write();
    r1.expected->Write();
    r2.expected->Write();
    r5.expected->Write();

    out.Close();

    cout << endl;
    cout << "Created:" << endl;
    cout << "  analysis_plots/final_mhh_shapes_kl_m15.png" << endl;
    cout << "  analysis_plots/final_mhh_expected_kl_m15.png" << endl;
    cout << "  analysis_plots/final_mhh_templates_kl_m15.root" << endl;
    cout << endl;
}
