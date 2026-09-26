# Higgs Self-Coupling in HH → bb̄γγ

Detector-level phenomenology study of the dependence of Higgs-pair production and
reconstructed HH → bb̄γγ observables on the Higgs self-coupling modifier κλ.

## Scope

The study traces κλ information from:

gg → HH production rate
→ production kinematics
→ forced H → bb̄ and H → γγ topology
→ Pythia8 showering
→ HepMC
→ Delphes fast detector simulation
→ reconstructed observables.

This is a phenomenological detector-level study, not an ATLAS or CMS
experimental analysis.

## Detector-level benchmark points

κλ = -1.5, 0, 1, 2, 5

Each detector benchmark contains 10,000 generated events.

## Final reconstruction selection

- jets: pT > 25 GeV, |η| < 2.5
- at least two jets
- at least two b-tagged jets
- photons: pT > 25 GeV, |η| < 2.5
- at least two photons
- b-jet pair chosen as the pair closest to mH = 125 GeV
- two leading photons
- 100 < mbb < 150 GeV
- 120 < mγγ < 130 GeV

The corrected detector-shape analysis reproduces the established final selected
event counts:

| κλ | Selected | Efficiency |
|---:|---:|---:|
| -1.5 | 876 | 8.76% |
| 0 | 942 | 9.42% |
| 1 | 1043 | 10.43% |
| 2 | 1034 | 10.34% |
| 5 | 703 | 7.03% |

## Reconstructed observables

The detector-level shape study considers:

- mbb
- mγγ
- mHH
- pT(Hbb)
- pT(Hγγ)
- ΔRbb
- ΔRγγ

The six non-mHH publication figures are included in `results/figures/`.
The mHH result is handled separately by the dedicated mHH analysis.

## Multivariate study

The compact BDT is a controlled signal-vs-signal benchmark.
It must not be interpreted as a signal-background experimental classifier
or as evidence for an experimental sensitivity improvement.

## Limitations

The current study does not include:

- a full bb̄γγ background model;
- detector/systematic nuisance parameters;
- a full experimental likelihood;
- collision data;
- a realistic ATLAS/CMS statistical model.

The selected yields use fast detector simulation and finite Monte Carlo samples.

## Reproducibility

The repository contains analysis code and methodology, but not the large generated
ROOT/LHE/HepMC event samples or complete MG5/Delphes installations.

Users should install compatible versions of MG5_aMC@NLO, Pythia8, HepMC,
Delphes and ROOT and follow the documented workflow.
