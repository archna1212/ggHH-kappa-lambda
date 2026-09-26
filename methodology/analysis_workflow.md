# Analysis workflow

The detector-level analysis follows:

1. MadGraph5_aMC@NLO generates loop-induced gg -> HH benchmark samples.
2. The rare H -> bbbar and H -> gammagamma topology is constructed using the validated forced-decay procedure.
3. Events are showered and hadronised with Pythia8.
4. Events are transported through HepMC.
5. Delphes fast detector simulation is applied using the ATLAS detector card.
6. ROOT/ExRootAnalysis performs reconstruction and observable extraction.

## Detector benchmarks

The detector-level samples contain 10,000 events for each:

- kappa_lambda = -1.5
- kappa_lambda = 0
- kappa_lambda = 1
- kappa_lambda = 2
- kappa_lambda = 5

## Final selection

Jets:
- pT > 25 GeV
- |eta| < 2.5
- at least two selected jets
- at least two b-tagged jets

Photons:
- pT > 25 GeV
- |eta| < 2.5
- at least two selected photons

Higgs reconstruction:
- choose the b-tagged jet pair closest to 125 GeV
- choose the two leading photons
- require 100 < m_bb < 150 GeV
- require 120 < m_gammagamma < 130 GeV

The reconstructed Higgs candidates are

H_bb = b1 + b2

H_gammagamma = gamma1 + gamma2

H_HH = H_bb + H_gammagamma

## Validated selected-event counts

| kappa_lambda | selected events | efficiency |
|---:|---:|---:|
| -1.5 | 876 | 8.76% |
| 0 | 942 | 9.42% |
| 1 | 1043 | 10.43% |
| 2 | 1034 | 10.34% |
| 5 | 703 | 7.03% |

These counts were reproduced independently by the corrected final-selection shape analysis.

## Reconstructed observables

The detector-level shape analysis considers:

- m_bb
- m_gammagamma
- m_HH
- pT(H_bb)
- pT(H_gammagamma)
- DeltaR_bb
- DeltaR_gammagamma

The publication figures supplied in `results/figures/` correspond to the six non-m_HH observables. The dedicated m_HH analysis is kept separate.

## Statistical scope

The present repository does not provide a background-aware experimental likelihood.

The BDT study is signal-vs-signal and is therefore an information-content cross-check rather than a signal-background experimental classifier.

The signal-only rate-plus-shape statistical study should not be interpreted as an experimental exclusion or discovery sensitivity.
