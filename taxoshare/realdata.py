"""Loading the quantMS/Triqler peptide table and the mock-community recipe."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .tree import Tree

CONDITION_TO_RECIPE = {"Constant": "C", "Protein": "P", "Uneven": "U"}


def load_taxonomy(path: Path, labels: list[str]) -> Tree:
    """Tree whose leaves are the strain *labels*, built from NCBI lineages.

    Each label is attached as a leaf beneath its taxid node, so two labels that
    resolve to the same taxid (the two Staphylococcus aureus strains here) become
    sibling leaves under a shared parent instead of silently merging into one.
    Labels NCBI could not resolve hang off a synthetic 'unresolved' node.
    """
    data = json.loads(Path(path).read_text())
    label2taxid, lineage = data["label2taxid"], data["lineage"]

    parent: dict[str, str | None] = {"root": None, "unresolved": "root"}
    for taxid, chain in lineage.items():
        prev = "root"
        for tid, _name, _rank in chain:
            node = f"tax{tid}"
            if node not in parent:
                parent[node] = prev
            prev = node

    for label in labels:
        tid = label2taxid.get(label)
        parent[label] = f"tax{tid}" if tid else "unresolved"
    return Tree(parent)


def load_recipe(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    labels = [r["Label"].strip() for r in rows]
    truth: dict[str, np.ndarray] = {}
    for cond, key in CONDITION_TO_RECIPE.items():
        col = f"Input Protein Amount in {key} (ug) normalized"
        vals = np.array([float(r[col]) for r in rows], dtype=float)
        truth[cond] = vals / vals.sum()
    return labels, truth


def load_peptides(path: Path, labels: set[str]) -> pd.DataFrame:
    """One row per peptide observation, with its candidate species set.

    Decoys and rows whose accessions name no community member (contaminants) are
    dropped. `pep` is 1 - searchScore: the quantMS table already reports a PSM
    posterior probability.
    """
    df = pd.read_csv(path, sep="\t")
    df = df[~df["proteins"].astype(str).str.contains("DECOY", na=False)].copy()

    def species_of(accessions: str) -> frozenset[str]:
        out = set()
        for acc in str(accessions).split(";"):
            prefix = acc.strip().split("_")[0]
            if prefix in labels:
                out.add(prefix)
        return frozenset(out)

    df["candidates"] = df["proteins"].map(species_of)
    df = df[df["candidates"].map(len) > 0].copy()
    df["pep"] = (1.0 - df["searchScore"].astype(float)).clip(lower=0.0, upper=1.0)
    df["intensity"] = df["intensity"].astype(float)
    return df
