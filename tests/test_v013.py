"""Tests for the v0.1.3 additions.

Covers the top-level naming API (``name_reaction`` / ``batch_name_reaction``),
the RDChiral-based ``Reaction.explain`` template, the ``Compound`` class, and a
guard that importing the package stays lazy (no eager torch/rxnmapper import).
"""
import math
import subprocess
import sys

import numpy as np
import pytest
import requests

import rxn_insight as ri
from rxn_insight import Compound, Reaction, batch_name_reaction, name_reaction
from rxn_insight.utils import get_similarity

SUZUKI = "BrC1=CC=CC=C1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


def test_public_api():
    """The stable public surface is exported."""
    for symbol in (
        "Reaction",
        "Molecule",
        "Compound",
        "ORDDatabase",
        "Database",
        "name_reaction",
        "batch_name_reaction",
    ):
        assert symbol in ri.__all__
        assert hasattr(ri, symbol)


def test_name_reaction_top_level():
    """``ri.name_reaction`` maps a reaction SMILES straight to a name."""
    assert name_reaction(SUZUKI) == "Suzuki coupling with boronic acids"


def test_batch_name_reaction():
    """``ri.batch_name_reaction`` names a list and honours n_jobs."""
    names = batch_name_reaction([SUZUKI, SUZUKI], n_jobs=1)
    assert names == ["Suzuki coupling with boronic acids"] * 2


def test_explain_uses_rdchiral_template():
    """``explain`` returns an RDChiral SMIRKS template and analysis keys."""
    explanation = Reaction(SUZUKI).explain()
    template = explanation["template"]
    assert isinstance(template, str) and ">>" in template
    for key in ("reaction_center", "bond_changes", "classification"):
        assert key in explanation


def test_compound_from_name():
    """``Compound`` resolves a chemical name to a structure (needs network)."""
    try:
        compound = Compound("benzene")
    except Exception as exc:  # OPSIN/PubChem endpoint may be unreachable
        pytest.skip(f"name-resolution service unavailable: {exc}")
    assert compound.smiles == "c1ccccc1"
    assert compound.mol is not None


@pytest.mark.parametrize("kwargs", [{}, {"use_opsin": False}], ids=["opsin", "pubchem"])
def test_compound_unresolvable_name_raises(kwargs):
    """Both resolution routes report an unresolvable name as a ValueError.

    They used to disagree: OPSIN raised KeyError('cml') while PubChem returned a
    half-built object with no ``smiles`` attribute at all.
    """
    try:
        with pytest.raises(ValueError, match="could not resolve"):
            Compound("definitely-not-a-real-chemical-xyzzy", **kwargs)
    except requests.RequestException as exc:  # endpoint unreachable
        pytest.skip(f"name-resolution service unavailable: {exc}")


@pytest.mark.parametrize(
    "metric",
    ["jaccard", "dice", "kulczynski1", "rogerstanimoto",
     "russellrao", "sokalmichener", "sokalsneath"],
)
def test_binary_similarity_metrics(metric):
    """Every advertised binary metric is dispatchable and returns a float."""
    v1 = np.array([1, 0, 1, 1])
    v2 = np.array([1, 1, 1, 0])
    assert isinstance(get_similarity(v1, v2, metric=metric), float)


def test_kulczynski1_handles_identical_vectors():
    """Vectors that never mismatch must not raise (SciPy returned inf here).

    ``kulczynski1`` divides by the number of mismatches, so comparing a
    fingerprint with itself would otherwise be a ZeroDivisionError.
    """
    v = np.array([1, 0, 1, 1])
    assert math.isinf(get_similarity(v, v, metric="kulczynski1"))
    assert get_similarity(v, v, metric="sokalmichener") == 1.0


def test_unknown_metric_raises():
    """An unsupported metric name is rejected."""
    v = np.array([1, 0, 1, 1])
    with pytest.raises(ValueError):
        get_similarity(v, v, metric="not-a-metric")


def test_import_is_lazy():
    """Importing rxn_insight does not eagerly import torch/rxnmapper.

    Runs in a fresh interpreter so it is independent of other tests that may
    have already triggered the lazy rxnmapper/torch import.
    """
    code = (
        "import sys, rxn_insight; "
        "heavy = [m for m in ('torch', 'rxnmapper') if m in sys.modules]; "
        "sys.exit('eagerly imported: ' + ','.join(heavy) if heavy else 0)"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
