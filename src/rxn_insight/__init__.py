from rxn_insight.reaction import Reaction
from rxn_insight.molecule import Molecule, Compound
from rxn_insight.ord import ORDDatabase
from rxn_insight.database import Database
from rxn_insight import naming as _naming

__all__ = [
    "Reaction",
    "Molecule",
    "Compound",
    "ORDDatabase",
    "Database",
    "name_reaction",
    "batch_name_reaction",
]

_DEFAULT_SMIRKS_DB = None


def _default_smirks_db():
    """Load and cache the bundled, curated SMIRKS database."""
    global _DEFAULT_SMIRKS_DB
    if _DEFAULT_SMIRKS_DB is None:
        _DEFAULT_SMIRKS_DB = _naming._load_smirks_db()
    return _DEFAULT_SMIRKS_DB


def name_reaction(reaction_smiles, smirks_db=None):
    """Return the named reaction for a reaction SMILES.

    Bare-minimum helper: a reaction SMILES goes in, the named reaction
    (e.g. ``"Suzuki coupling"``) comes out — or ``"OtherReaction"`` when
    nothing matches. The bundled SMIRKS database is loaded and cached on
    first use unless an explicit ``smirks_db`` is supplied.
    """
    if smirks_db is None:
        smirks_db = _default_smirks_db()
    return _naming.name_reaction(reaction_smiles, smirks_db)


def batch_name_reaction(reactions, n_jobs=None, smirks_db=None):
    """Name a list of reaction SMILES in parallel.

    Thin wrapper around :func:`rxn_insight.naming.name_reactions_batch`
    (joblib-parallel). ``n_jobs`` defaults to half the available CPU cores;
    pass ``1`` to disable parallelism. Returns a list of names aligned with
    the input.
    """
    return _naming.name_reactions_batch(reactions, n_jobs=n_jobs, smirks_db=smirks_db)
