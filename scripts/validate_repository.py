#!/usr/bin/env python3
"""Static consistency checks for the EeFF-water release package.

No external Python packages are required.
"""

from __future__ import annotations

from pathlib import Path
import csv
import re
import sys

ROOT = Path(__file__).resolve().parents[1]


def fail(msg: str) -> None:
    raise AssertionError(msg)


def numeric_key(path: Path):
    m = re.search(r"_(\d+)$", path.name)
    return int(m.group(1)) if m else path.name


def parse_data_file(path: Path):
    """Return (n_atoms, n_nuclei, spins) from a supplied LAMMPS data file."""
    lines = path.read_text(errors="replace").splitlines()
    n_atoms = None
    for line in lines[:30]:
        m = re.match(r"\s*(\d+)\s+atoms\s*$", line)
        if m:
            n_atoms = int(m.group(1))
            break
    if n_atoms is None:
        fail(f"Could not read atom count from {path}")

    in_atoms = False
    parsed = 0
    n_nuclei = 0
    spins = set()
    for line in lines:
        s = line.strip()
        if s.startswith("Atoms"):
            in_atoms = True
            continue
        if not in_atoms:
            continue
        if not s:
            continue
        if re.match(r"^[A-Za-z]", s):
            if parsed:
                break
            continue
        parts = s.split()
        if len(parts) < 8 or not re.match(r"^[+-]?\d+$", parts[0]):
            continue
        try:
            atom_type = int(parts[1])
            spin = int(float(parts[6]))
        except ValueError:
            continue
        parsed += 1
        spins.add(spin)
        if atom_type in (1, 2):
            n_nuclei += 1
        if parsed == n_atoms:
            break

    if parsed != n_atoms:
        fail(f"Parsed {parsed} atoms but header says {n_atoms}: {path}")
    return n_atoms, n_nuclei, spins


def check_parameters():
    order = ROOT / "parameters" / "parameter_order.csv"
    with order.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 59:
        fail(f"parameter_order.csv has {len(rows)} rows; expected 59")
    for i, row in enumerate(rows):
        if int(row["argument_index_zero_based"]) != i:
            fail("parameter_order.csv indices are not contiguous 0-58")

    inc = (ROOT / "parameters" / "eeff_water_pair_style.inc").read_text()
    line = next((x for x in inc.splitlines() if x.strip().startswith("pair_style")), None)
    if not line:
        fail("No pair_style command found in parameter include")
    tokens = line.split()
    if tokens[:2] != ["pair_style", "eff/cut"]:
        fail("Unexpected pair_style header")
    values = tokens[2:]
    if len(values) != 59:
        fail(f"pair_style include has {len(values)} numeric arguments; expected 59")
    for row, token in zip(rows, values):
        if abs(float(row["value"]) - float(token)) > 1e-12:
            fail(f"Parameter mismatch at argument {row['argument_index_zero_based']}")

    mapping = ROOT / "parameters" / "manuscript_parameter_mapping.csv"
    with mapping.open(encoding="utf-8", newline="") as f:
        mrows = list(csv.DictReader(f))
    if len(mrows) != 59:
        fail("manuscript_parameter_mapping.csv must contain 59 rows")
    if any(r["used_by_eeff_water"] != "no" for r in mrows[57:]):
        fail("Compatibility arguments 57-58 must be marked inactive")

    cpp = (ROOT / "src" / "pair_eff_cut.cpp").read_text()
    if "if (narg != 59)" not in cpp:
        fail("PairEffCut::settings() does not enforce exactly 59 arguments")


def check_training_data():
    base = ROOT / "training_data"
    dirs = sorted(p for p in base.iterdir() if p.is_dir())
    total_configs = 0
    all_spins = set()
    for d in dirs:
        data_files = sorted(d.glob("data.*"), key=numeric_key)
        if not data_files:
            continue
        energies = [x for x in (d / "energy.txt").read_text().splitlines() if x.strip()]
        force_lines = [x for x in (d / "force.txt").read_text().splitlines() if x.strip()]
        if len(energies) != len(data_files):
            fail(f"Energy count mismatch in {d.name}: {len(energies)} vs {len(data_files)}")
        expected_force_lines = 0
        for df in data_files:
            _, n_nuclei, spins = parse_data_file(df)
            expected_force_lines += 3 * n_nuclei
            all_spins.update(spins)
        if len(force_lines) != expected_force_lines:
            fail(f"Force-label count mismatch in {d.name}: {len(force_lines)} vs {expected_force_lines}")
        total_configs += len(data_files)
    if total_configs != 1772:
        fail(f"Training-data total is {total_configs}; expected 1772")
    if 3 in {abs(x) for x in all_spins}:
        fail("Spin-3 particles found in training data")


def check_examples():
    expected = [
        ("ambient_water_512h2o/in.ambient", "ambient_water_512h2o/data.512h2o_300k_1atm"),
        ("compressed_liquid_water/in.compressed", "compressed_liquid_water/data.512h2o"),
        ("proton_defect_water/in.h3o", "proton_defect_water/data.511h2o_1h3o"),
        ("proton_defect_water/in.oh", "proton_defect_water/data.511h2o_1oh"),
        ("superionic_water/in.bcc", "superionic_water/data.bcc_432h2o"),
        ("superionic_water/in.fcc", "superionic_water/data.fcc_500h2o"),
        ("superionic_water/in.hcp", "superionic_water/data.hcp_512h2o"),
        ("water_cluster_6/in.minimize", "water_cluster_6/data.perturbed_6_1"),
    ]
    all_spins = set()
    for inp, data in expected:
        ip = ROOT / "examples" / inp
        dp = ROOT / "examples" / data
        if not ip.exists(): fail(f"Missing example input: {ip}")
        if not dp.exists(): fail(f"Missing example data: {dp}")
        txt = ip.read_text()
        if "include ../../parameters/eeff_water_pair_style.inc" not in txt:
            fail(f"Example does not use the shared parameter include: {ip}")
        if re.search(r"^\s*pair_style\s+eff/cut", txt, flags=re.M):
            fail(f"Example duplicates pair_style parameters instead of using the shared include: {ip}")
        dump_ids = set(re.findall(r"^\s*dump\s+(\S+)", txt, flags=re.M))
        for dump_id in re.findall(r"^\s*dump_modify\s+(\S+)", txt, flags=re.M):
            if dump_id not in dump_ids:
                fail(f"dump_modify refers to undefined dump ID '{dump_id}' in {ip}")
        _, _, spins = parse_data_file(dp)
        all_spins.update(spins)
    if 3 in {abs(x) for x in all_spins}:
        fail("Spin-3 particles found in representative example data")

    ambient = (ROOT / "examples/ambient_water_512h2o/in.ambient").read_text()
    if not re.search(r"^run\s+5000000\s*$", ambient, flags=re.M):
        fail("Ambient-water production segment must be 5,000,000 steps (1 ns at 0.2 fs)")

    # Broken slash-delimited index forms from the draft package must not remain.
    for p in (ROOT / "examples").rglob("in.*"):
        for line in p.read_text().splitlines():
            if re.search(r"^\s*variable\s+\S+\s+index\s+\S+/\S+", line):
                fail(f"Slash-delimited variable index remains in {p}: {line}")


def check_csv_files():
    for p in (ROOT / "Source_Data").rglob("*.csv"):
        with p.open(encoding="utf-8", newline="") as f:
            rows = list(csv.reader(f))
        if not rows or rows[0] != ["r_A", "g_r"]:
            fail(f"Unexpected CSV header in {p}")
        for i, row in enumerate(rows[1:], start=2):
            if len(row) != 2:
                fail(f"CSV row {i} in {p} does not have two columns")
            float(row[0]); float(row[1])


def check_public_docs():
    files = [
        ROOT / "README.md",
        ROOT / "CITATION.cff",
        ROOT / "CODE_AVAILABILITY.md",
        ROOT / "LICENSE.md",
        ROOT / "parameters/README.md",
        ROOT / "lammps_patch/BASE_VERSION.txt",
        ROOT / "lammps_patch/INSTALL.md",
    ]
    bad_patterns = ["<USER_OR_ORG>", "YYYY-MM-DD", "10.xxxx", "Add all manuscript authors", "TODO before"]
    for p in files:
        txt = p.read_text(encoding="utf-8")
        for pat in bad_patterns:
            if pat in txt:
                fail(f"Placeholder '{pat}' remains in {p}")


def main():
    checks = [
        ("parameters", check_parameters),
        ("training data", check_training_data),
        ("examples", check_examples),
        ("source-data CSV format", check_csv_files),
        ("public documentation", check_public_docs),
    ]
    for name, fn in checks:
        fn()
        print(f"[OK] {name}")
    print("[OK] EeFF-water repository static validation completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        sys.exit(1)
