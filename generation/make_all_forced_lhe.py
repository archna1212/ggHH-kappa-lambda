#!/usr/bin/env python3

import gzip
import math
import os
import random
from pathlib import Path

# ============================================================
# Configuration
# ============================================================

BASE = Path.home() / "Desktop/MG5_aMC_v3_6_3/ggHH_bbgg"
EVENTS = BASE / "Events"

KAPPAS = [0, 1, 2, 5]
N_EVENTS = 10000

# Existing KL1/KL2 files will be kept.
SKIP_EXISTING = True

# Higgs/daughter masses in GeV
MB = 4.8
MGAMMA = 0.0

# Reproducible random seeds
SEED_BASE = 10000


# ============================================================
# Basic four-vector utilities
# ============================================================

def invariant_mass(px, py, pz, E):
    m2 = E * E - px * px - py * py - pz * pz
    return math.sqrt(max(m2, 0.0))


def boost(px, py, pz, E, bx, by, bz):
    """
    Boost a four-vector from the parent's rest frame
    to the lab frame using beta = p/E.
    """

    beta2 = bx * bx + by * by + bz * bz

    if beta2 < 1e-16:
        return px, py, pz, E

    gamma = 1.0 / math.sqrt(1.0 - beta2)
    bp = bx * px + by * py + bz * pz
    gamma2 = (gamma - 1.0) / beta2

    px_lab = px + gamma2 * bp * bx + gamma * bx * E
    py_lab = py + gamma2 * bp * by + gamma * by * E
    pz_lab = pz + gamma2 * bp * bz + gamma * bz * E
    E_lab = gamma * (E + bp)

    return px_lab, py_lab, pz_lab, E_lab


def two_body_decay(parent, m1, m2, rng):
    """
    Isotropic two-body decay of a parent four-vector.

    parent = (px, py, pz, E)

    Returns:
        daughter1 four-vector
        daughter2 four-vector
    """

    px, py, pz, E = parent

    # Parent invariant mass
    M = invariant_mass(px, py, pz, E)

    if M <= m1 + m2:
        raise RuntimeError(
            f"Cannot decay parent with M={M:.6f} GeV "
            f"into masses {m1} and {m2}"
        )

    # Two-body momentum in parent rest frame
    term1 = M * M - (m1 + m2) ** 2
    term2 = M * M - (m1 - m2) ** 2

    pstar = math.sqrt(max(term1 * term2, 0.0)) / (2.0 * M)

    E1star = math.sqrt(m1 * m1 + pstar * pstar)
    E2star = math.sqrt(m2 * m2 + pstar * pstar)

    # Isotropic direction
    cos_theta = rng.uniform(-1.0, 1.0)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    phi = rng.uniform(0.0, 2.0 * math.pi)

    px_star = pstar * sin_theta * math.cos(phi)
    py_star = pstar * sin_theta * math.sin(phi)
    pz_star = pstar * cos_theta

    # Daughter 1
    d1_rest = (
        px_star,
        py_star,
        pz_star,
        E1star,
    )

    # Daughter 2: opposite momentum
    d2_rest = (
        -px_star,
        -py_star,
        -pz_star,
        E2star,
    )

    # Parent velocity
    bx = px / E
    by = py / E
    bz = pz / E

    d1_lab = boost(*d1_rest, bx, by, bz)
    d2_lab = boost(*d2_rest, bx, by, bz)

    return d1_lab, d2_lab


# ============================================================
# LHE particle formatting
# ============================================================

def format_particle(
    pid,
    status,
    mother1,
    mother2,
    color1,
    color2,
    px,
    py,
    pz,
    E,
    mass,
    lifetime=0.0,
    spin=0.0,
):
    return (
        f"{pid:8d} "
        f"{status:2d} "
        f"{mother1:4d} "
        f"{mother2:4d} "
        f"{color1:4d} "
        f"{color2:4d} "
        f"{px: .11e} "
        f"{py: .11e} "
        f"{pz: .11e} "
        f"{E: .11e} "
        f"{mass: .11e} "
        f"{lifetime: .11e} "
        f"{spin: .11e}\n"
    )


def parse_particle(line):
    fields = line.split()

    if len(fields) < 13:
        return None

    return {
        "id": int(fields[0]),
        "status": int(fields[1]),
        "m1": int(fields[2]),
        "m2": int(fields[3]),
        "c1": int(fields[4]),
        "c2": int(fields[5]),
        "px": float(fields[6]),
        "py": float(fields[7]),
        "pz": float(fields[8]),
        "E": float(fields[9]),
        "M": float(fields[10]),
        "lifetime": float(fields[11]),
        "spin": float(fields[12]),
    }


# ============================================================
# Transform one event
# ============================================================

def transform_event(event_lines, rng):
    """
    Convert:

        gg -> HH

    into an explicit forced decay topology:

        H1 -> b bbar
        H2 -> gamma gamma
    """

    # Find event-header line
    event_header_index = None

    for i, line in enumerate(event_lines):
        stripped = line.strip()

        if stripped.startswith("<event>"):
            event_header_index = i
            break

    if event_header_index is None:
        raise RuntimeError("Missing <event>")

    # Locate closing event tag
    end_index = None

    for i in range(event_header_index + 1, len(event_lines)):
        if event_lines[i].strip() == "</event>":
            end_index = i
            break

    if end_index is None:
        raise RuntimeError("Missing </event>")

    # Particle lines begin after the event header.
    header_line = event_lines[event_header_index + 1].strip()
    header_fields = header_line.split()

    if not header_fields:
        raise RuntimeError("Empty LHE event header")

    original_nparticles = int(header_fields[0])

    particle_lines = event_lines[
        event_header_index + 2:
        event_header_index + 2 + original_nparticles
    ]

    if len(particle_lines) != original_nparticles:
        raise RuntimeError(
            f"Expected {original_nparticles} particle lines, "
            f"found {len(particle_lines)}"
        )

    particles = []

    for line in particle_lines:
        particle = parse_particle(line)

        if particle is None:
            raise RuntimeError(f"Could not parse particle line:\n{line}")

        particles.append(particle)

    # Find the two Higgs bosons
    higgs_indices = [
        i for i, p in enumerate(particles)
        if p["id"] == 25
    ]

    if len(higgs_indices) != 2:
        raise RuntimeError(
            f"Expected exactly 2 Higgs bosons, "
            f"found {len(higgs_indices)}"
        )

    h1_index = higgs_indices[0]
    h2_index = higgs_indices[1]

    h1 = particles[h1_index]
    h2 = particles[h2_index]

    # Preserve Higgs four-momenta, but mark Higgs as decayed.
    h1["status"] = 2
    h2["status"] = 2

    # Parent particle indices in LHE are 1-based.
    h1_number = h1_index + 1
    h2_number = h2_index + 1

    # --------------------------------------------------------
    # H1 -> b bbar
    # --------------------------------------------------------

    h1_parent = (
        h1["px"],
        h1["py"],
        h1["pz"],
        h1["E"],
    )

    b, bbar = two_body_decay(
        h1_parent,
        MB,
        MB,
        rng,
    )

    b_px, b_py, b_pz, b_E = b
    bb_px, bb_py, bb_pz, bb_E = bbar

    b_mass = invariant_mass(
        b_px, b_py, b_pz, b_E
    )

    bbar_mass = invariant_mass(
        bb_px, bb_py, bb_pz, bb_E
    )

    # --------------------------------------------------------
    # H2 -> gamma gamma
    # --------------------------------------------------------

    h2_parent = (
        h2["px"],
        h2["py"],
        h2["pz"],
        h2["E"],
    )

    gamma1, gamma2 = two_body_decay(
        h2_parent,
        MGAMMA,
        MGAMMA,
        rng,
    )

    g1_px, g1_py, g1_pz, g1_E = gamma1
    g2_px, g2_py, g2_pz, g2_E = gamma2

    # --------------------------------------------------------
    # Build output particle list
    # --------------------------------------------------------

    output_particles = particles.copy()

    output_particles.extend([
        {
            "id": 5,
            "status": 1,
            "m1": h1_number,
            "m2": h1_number,
            "c1": 503,
            "c2": 0,
            "px": b_px,
            "py": b_py,
            "pz": b_pz,
            "E": b_E,
            "M": b_mass,
            "lifetime": 0.0,
            "spin": 0.0,
        },
        {
            "id": -5,
            "status": 1,
            "m1": h1_number,
            "m2": h1_number,
            "c1": 0,
            "c2": 503,
            "px": bb_px,
            "py": bb_py,
            "pz": bb_pz,
            "E": bb_E,
            "M": bbar_mass,
            "lifetime": 0.0,
            "spin": 0.0,
        },
        {
            "id": 22,
            "status": 1,
            "m1": h2_number,
            "m2": h2_number,
            "c1": 0,
            "c2": 0,
            "px": g1_px,
            "py": g1_py,
            "pz": g1_pz,
            "E": g1_E,
            "M": 0.0,
            "lifetime": 0.0,
            "spin": 0.0,
        },
        {
            "id": 22,
            "status": 1,
            "m1": h2_number,
            "m2": h2_number,
            "c1": 0,
            "c2": 0,
            "px": g2_px,
            "py": g2_py,
            "pz": g2_pz,
            "E": g2_E,
            "M": 0.0,
            "lifetime": 0.0,
            "spin": 0.0,
        },
    ])

    # --------------------------------------------------------
    # Rewrite event header: 4 particles -> 8 particles
    # --------------------------------------------------------

    header_fields[0] = str(len(output_particles))
    new_header = " ".join(header_fields) + "\n"

    output = []

    # Keep <event>
    output.append(event_lines[event_header_index])

    # New event header
    output.append(new_header)

    # Original particles
    for p in output_particles:
        output.append(
            format_particle(
                p["id"],
                p["status"],
                p["m1"],
                p["m2"],
                p["c1"],
                p["c2"],
                p["px"],
                p["py"],
                p["pz"],
                p["E"],
                p["M"],
                p["lifetime"],
                p["spin"],
            )
        )

    # Preserve anything after the particle block
    # (mgrwt, weights, etc.)
    trailer_start = event_header_index + 2 + original_nparticles

    for line in event_lines[trailer_start:end_index]:
        output.append(line)

    output.append("</event>\n")

    return output


# ============================================================
# Transform one LHE file
# ============================================================

def make_forced_sample(input_file, output_file, n_events, seed):
    rng = random.Random(seed)

    print()
    print("=" * 70)
    print(f"Input : {input_file}")
    print(f"Output: {output_file}")
    print(f"Events: {n_events}")
    print(f"Seed  : {seed}")
    print("=" * 70)

    processed = 0

    with gzip.open(input_file, "rt") as fin, \
         gzip.open(output_file, "wt") as fout:

        # Copy file header until first event.
        # Do NOT write the first <event> here;
        # transform_event() will write it.
        for line in fin:
            if line.strip() == "<event>":
                current_event = [line]
                break

            fout.write(line)
        else:
            raise RuntimeError("No <event> found")

        for line in fin:
            current_event.append(line)

            if line.strip() == "</event>":
                if processed < n_events:
                    transformed = transform_event(
                        current_event,
                        rng,
                    )

                    for out_line in transformed:
                        fout.write(out_line)

                    processed += 1

                    if processed % 1000 == 0:
                        print(
                            f"  processed {processed}/{n_events}"
                        )

                current_event = []

                if processed >= n_events:
                    # Skip remaining original events, but copy
                    # the closing LesHouchesEvents tag.
                    for remaining in fin:
                        if "</LesHouchesEvents>" in remaining:
                            fout.write(remaining)
                            break
                    break

    if processed != n_events:
        raise RuntimeError(
            f"Only processed {processed} events, "
            f"requested {n_events}"
        )

    print(f"Done: {processed} events written.")


# ============================================================
# Main batch driver
# ============================================================

def main():

    print()
    print("FORCED HH -> bb gamma gamma BATCH GENERATOR")
    print("=" * 70)
    print(f"Required kappa values: {KAPPAS}")
    print(f"Events per sample     : {N_EVENTS}")
    print("=" * 70)

    for kl in KAPPAS:

        input_file = (
            EVENTS
            / f"run_kl{kl}_30k"
            / "unweighted_events.lhe.gz"
        )

        output_file = (
            EVENTS
            / f"run_kl{kl}_forcedbbgg_10k.lhe.gz"
        )

        if not input_file.exists():
            print()
            print(f"[SKIP] Input does not exist:")
            print(f"       {input_file}")
            continue

        if output_file.exists() and SKIP_EXISTING:
            print()
            print(f"[SKIP] Output already exists:")
            print(f"       {output_file}")
            continue

        seed = SEED_BASE + kl + 100

        make_forced_sample(
            input_file,
            output_file,
            N_EVENTS,
            seed,
        )

    print()
    print("=" * 70)
    print("ALL REQUESTED FORCED SAMPLES COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
