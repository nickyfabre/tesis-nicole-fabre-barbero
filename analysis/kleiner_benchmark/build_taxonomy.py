"""Build taxonomy.json (the tree used by thesis_benchmark.py) from the NCBI taxonomy.

Needs ete3 with its local NCBI database:
    pip install ete3
    python build_taxonomy.py --recipe ../../files/amounts_species.tsv --out taxonomy.json

Names in the recipe are matched to NCBI automatically. Names that are ambiguous or
spelled differently from NCBI are resolved through strain_names.tsv (next to this
script): a small, documented table, one row per strain, with the NCBI name or taxid
and where that information comes from. It is data, not code, so it can be checked
and cited in the Methods.

Output format (what taxoshare.realdata.load_taxonomy reads):
    {"label2taxid": {"K12": 562, ...},
     "lineage": {"562": [[1, "root", "no rank"], ..., [562, "Escherichia coli", "species"]], ...}}

Each strain label becomes its own leaf under its taxid, so the two S. aureus strains
(137, 259) stay separate sibling leaves instead of merging. Unresolved names stop the
script with an error instead of being silently mocked.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_overrides(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    # skip comment lines, also when a spreadsheet program has wrapped them in quotes
    lines = [l for l in open(path, encoding="utf-8-sig") if l.strip() and not l.lstrip('"').startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if rows and ("Label" not in rows[0] or "NCBI_name_or_taxid" not in rows[0]):
        raise SystemExit(f"{path.name}: expected columns Label and NCBI_name_or_taxid, found {list(rows[0])}")
    out = {r["Label"].strip(): r["NCBI_name_or_taxid"].strip().strip('"') for r in rows
           if (r.get("NCBI_name_or_taxid") or "").strip()}
    print(f"name overrides from {path.name}: " + (", ".join(f"{k} -> {v}" for k, v in out.items()) or "none"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recipe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--names", type=Path, default=HERE / "strain_names.tsv")
    a = ap.parse_args()

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    overrides = read_overrides(a.names)

    rows = list(csv.DictReader(open(a.recipe, encoding="utf-8"), delimiter="\t"))
    label2taxid, lineage, missing = {}, {}, []
    for r in rows:
        label, name = r["Label"].strip(), r["Species"].strip()
        query = overrides.get(label, name)
        if query.isdigit():
            tid = int(query)
        else:
            hit = ncbi.get_name_translator([query]).get(query)
            if not hit:
                missing.append(f"{label}: {name!r} (queried {query!r})")
                continue
            tid = int(hit[0])
        label2taxid[label] = tid
        chain = ncbi.get_lineage(tid)
        names = ncbi.get_taxid_translator(chain)
        ranks = ncbi.get_rank(chain)
        lineage[str(tid)] = [[int(t), names[t], ranks[t]] for t in chain]

    if missing:
        raise SystemExit("Could not resolve:\n  " + "\n  ".join(missing) +
                         f"\nFill in their row in {a.names.name} (NCBI name or taxid) and re-run.")
    a.out.write_text(json.dumps({"label2taxid": label2taxid, "lineage": lineage}, indent=1))
    print(f"wrote {a.out}: {len(label2taxid)} strains, {len(set(label2taxid.values()))} distinct taxids")
    for label, tid in label2taxid.items():
        print(f"  {label:<6} -> {tid:<8} {lineage[str(tid)][-1][1]}")


if __name__ == "__main__":
    main()
