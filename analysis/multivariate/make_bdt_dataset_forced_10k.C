#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

using namespace std;

void FillBDTTree(
    const char *filename,
    int label,
    TTree *outTree,
    double &mHH,
    double &pTHbb,
    double &pTHgg,
    double &dRbb,
    int &target
)
{
    TFile *file = TFile::Open(filename);

    if(!file || file->IsZombie())
    {
        cerr << "ERROR: Cannot open " << filename << endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("Delphes");

    if(!tree)
    {
        cerr << "ERROR: Delphes tree not found in "
             << filename << endl;
        file->Close();
        return;
    }

    ExRootTreeReader reader(tree);

    TClonesArray *branchJet =
        reader.UseBranch("Jet");

    TClonesArray *branchPhoton =
        reader.UseBranch("Photon");

    Long64_t nEntries = reader.GetEntries();

    target = label;

    for(Long64_t entry = 0; entry < nEntries; entry++)
    {
        reader.ReadEntry(entry);

        vector<Jet*> bjets;
        vector<Photon*> photons;

        // b-tagged jets
        for(int i = 0; i < branchJet->GetEntries(); i++)
        {
            Jet *jet = (Jet*)branchJet->At(i);

            if(!jet) continue;

            if(jet->BTag != 0 &&
               jet->PT > 20.0 &&
               fabs(jet->Eta) < 2.5)
            {
                bjets.push_back(jet);
            }
        }

        // photons
        for(int i = 0; i < branchPhoton->GetEntries(); i++)
        {
            Photon *ph = (Photon*)branchPhoton->At(i);

            if(!ph) continue;

            if(ph->PT > 25.0 &&
               fabs(ph->Eta) < 2.5)
            {
                photons.push_back(ph);
            }
        }

        if(bjets.size() < 2 || photons.size() < 2)
            continue;

        // Sort by pT
        sort(bjets.begin(), bjets.end(),
             [](Jet *a, Jet *b)
             {
                 return a->PT > b->PT;
             });

        sort(photons.begin(), photons.end(),
             [](Photon *a, Photon *b)
             {
                 return a->PT > b->PT;
             });

        Jet *b1 = bjets[0];
        Jet *b2 = bjets[1];

        Photon *g1 = photons[0];
        Photon *g2 = photons[1];

        TLorentzVector pb1, pb2, pg1, pg2;

        pb1.SetPtEtaPhiM(
            b1->PT, b1->Eta, b1->Phi, b1->Mass
        );

        pb2.SetPtEtaPhiM(
            b2->PT, b2->Eta, b2->Phi, b2->Mass
        );

        pg1.SetPtEtaPhiM(
            g1->PT, g1->Eta, g1->Phi, 0.0
        );

        pg2.SetPtEtaPhiM(
            g2->PT, g2->Eta, g2->Phi, 0.0
        );

        TLorentzVector Hbb = pb1 + pb2;
        TLorentzVector Hgg = pg1 + pg2;
        TLorentzVector HH  = Hbb + Hgg;

        mHH   = HH.M();
        pTHbb = Hbb.Pt();
        pTHgg = Hgg.Pt();
        dRbb  = pb1.DeltaR(pb2);

        outTree->Fill();
    }

    file->Close();
}


void make_bdt_dataset_forced_10k()
{
    cout << endl;
    cout << "==========================================" << endl;
    cout << "CREATING BDT DATASET" << endl;
    cout << "==========================================" << endl;

    TFile *outFile =
        new TFile(
            "analysis_plots/bdt_dataset_forced_10k.root",
            "RECREATE"
        );

    TTree *tree =
        new TTree("BDTTree","KL1 vs KL2");

    double mHH   = 0.0;
    double pTHbb = 0.0;
    double pTHgg = 0.0;
    double dRbb  = 0.0;
    int target   = -1;

    tree->Branch("mHH",   &mHH,   "mHH/D");
    tree->Branch("pTHbb", &pTHbb, "pTHbb/D");
    tree->Branch("pTHgg", &pTHgg, "pTHgg/D");
    tree->Branch("dRbb",  &dRbb,  "dRbb/D");
    tree->Branch("target",&target,"target/I");

    // KL = 1
    cout << "Reading KL = 1..." << endl;

    FillBDTTree(
        "Events/run_kl1_forcedbbgg_10k_delphes.root",
        0,
        tree,
        mHH,
        pTHbb,
        pTHgg,
        dRbb,
        target
    );

    Long64_t nKL1 = tree->GetEntries();

    cout << "KL = 1 entries: "
         << nKL1 << endl;

    // KL = 2
    cout << "Reading KL = 2..." << endl;

    FillBDTTree(
        "Events/run_kl2_forcedbbgg_10k_delphes.root",
        1,
        tree,
        mHH,
        pTHbb,
        pTHgg,
        dRbb,
        target
    );

    Long64_t nTotal = tree->GetEntries();
    Long64_t nKL2 = nTotal - nKL1;

    cout << "KL = 2 entries: "
         << nKL2 << endl;

    outFile->cd();
    tree->Write();
    outFile->Close();

    cout << endl;
    cout << "==========================================" << endl;
    cout << "BDT DATASET CREATED" << endl;
    cout << "Total entries = "
         << nTotal << endl;
    cout << "Output: analysis_plots/bdt_dataset_forced_10k.root" << endl;
    cout << "==========================================" << endl;
}
