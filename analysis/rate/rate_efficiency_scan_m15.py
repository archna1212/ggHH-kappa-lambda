import csv

kl = [-1.5, 0, 1, 2, 5]
xs = [0.070073, 0.030540, 0.014504, 0.006786, 0.033068]
generated = [10000, 10000, 10000, 10000, 10000]
selected = [876, 942, 1043, 1034, 703]

BRbb = 0.5824
BRgg = 0.00227
BRbbgg = 2.0 * BRbb * BRgg
lumi = 3000.0

rows = []

for k, x, gen, sel in zip(kl, xs, generated, selected):
    eff = sel / gen
    yield_events = x * BRbbgg * eff * lumi * 1000.0

    rows.append({
        "kappa_lambda": k,
        "cross_section_pb": x,
        "generated": gen,
        "selected": sel,
        "efficiency": eff,
        "efficiency_percent": 100.0 * eff,
        "expected_selected_events_3000fb": yield_events
    })

with open("analysis_plots/rate_efficiency_scan_m15.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("Created analysis_plots/rate_efficiency_scan_m15.csv")
for r in rows:
    print(
        f"kappa={r['kappa_lambda']:>4}: "
        f"sigma={r['cross_section_pb']:.6f} pb, "
        f"eff={r['efficiency_percent']:.2f}%, "
        f"N={r['expected_selected_events_3000fb']:.5f}"
    )
