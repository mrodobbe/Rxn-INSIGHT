import itertools
import multiprocessing as mp
import os
import re
from importlib import resources
from typing import Dict, List, Optional, Set, Tuple, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd

_ATOM_MAP_RE = re.compile(r":\d+")


# Module-level cache for compiled SMIRKS reaction objects.
# Each worker process in joblib builds its own cache on first use.
_compiled_smirks: dict = {}


def _get_compiled_rxn(smirks: str):
    """Return a compiled ChemicalReaction for the given SMIRKS, caching the result."""
    try:
        return _compiled_smirks[smirks]
    except KeyError:
        rxn_obj = AllChem.ReactionFromSmarts(smirks)
        _compiled_smirks[smirks] = rxn_obj
        return rxn_obj


def name_reaction(rxn: str, smirks_db: pd.DataFrame) -> str:
    """Determines the name of the reaction from a database based on SMIRKS transformations.

    Args:
        rxn (str): The name of the reaction in Reaction SMILES format.
        smirks_db (pd.DataFrame): DataFrame containing SMIRKS patterns and corresponding reaction names.

    Returns:
        str: The name of the reaction, or 'OtherReaction' if no specific name can be determined.
    """
    parts = rxn.split(">>")
    if len(parts) != 2:
        return "OtherReaction"
    reactants_smiles, products_smiles = parts
    reactants = reactants_smiles.split(".")
    products = products_smiles.split(".")

    if (
            len(reactants) > 4 or len(products) > 4
    ):  # There are no templates for reactions with more than four reactants.
        return "OtherReaction"

    new_products = []  # Try to canonicalize SMILES

    for product in products:
        try:
            new_products.append(
                Chem.MolToSmiles(Chem.MolFromSmiles(product), isomericSmiles=False)
            )
        except:
            new_products.append(product)

    num_reactants = len(reactants)

    # Allow templates with nreact <= num_reactants so that rc_only
    # templates (which drop reagent components) still match.
    selected_rxns = smirks_db[smirks_db["nreact"] <= num_reactants]
    react_mols = [Chem.MolFromSmiles(reactant) for reactant in reactants]

    # TODO: Further refine reactions by superclass

    for i in selected_rxns.index:  # Iterate over all reactants to find a match
        smirks = selected_rxns["smirks"][i]
        rxn = _get_compiled_rxn(smirks)
        nreact_template = int(selected_rxns["nreact"][i])

        # Build reactant tuples: subsets of the right size × permutations.
        # When nreact_template == num_reactants this is equivalent to the
        # original full-permutation logic.
        if nreact_template == num_reactants:
            if num_reactants == 1:
                all_tuples = [tuple(react_mols)]
            else:
                all_tuples = list(itertools.permutations(react_mols))
        else:
            all_tuples = []
            for subset in itertools.combinations(react_mols, nreact_template):
                if nreact_template == 1:
                    all_tuples.append(subset)
                else:
                    all_tuples.extend(itertools.permutations(subset))

        matched = False
        for tup in all_tuples:
            try:
                outcomes = rxn.RunReactants(tup)
            except Exception:
                continue
            for prods in outcomes:
                try:
                    prod = Chem.MolToSmiles(prods[0], isomericSmiles=False)
                except Exception:
                    continue
                if prod in new_products:
                    matched = True
                    break
            if matched:
                break

        if matched:
            return selected_rxns["name"][i].strip("{}")

    return "OtherReaction"


def test_smirks(rxn: str, smirks: str) -> dict:
    """Test whether a SMIRKS template fires on a reaction.

    Tries all reactant subsets and permutations so that rc_only templates
    (fewer reactant components than the actual reaction) are handled.

    Args:
        rxn: Reaction SMILES (``reactants>>products``).
        smirks: SMIRKS template string.

    Returns:
        Dict with ``applicable`` (bool), ``correct`` (bool), and
        ``products`` (set of predicted SMILES).
    """
    out = {"applicable": False, "correct": False, "products": set()}
    parts = rxn.split(">>")
    if len(parts) != 2:
        return out
    reactants = [Chem.MolFromSmiles(s) for s in parts[0].split(".")]
    expected = {
        Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False)
        for s in parts[1].split(".")
    }

    try:
        rxn_obj = AllChem.ReactionFromSmarts(smirks)
    except Exception:
        return out
    if rxn_obj is None:
        return out

    nreact = rxn_obj.GetNumReactantTemplates()
    for subset in itertools.combinations(reactants, nreact):
        for perm in itertools.permutations(subset):
            try:
                outcomes = rxn_obj.RunReactants(perm)
            except Exception:
                continue
            for prods in outcomes:
                try:
                    smi = Chem.MolToSmiles(prods[0], isomericSmiles=False)
                except Exception:
                    continue
                out["applicable"] = True
                out["products"].add(smi)
                if smi in expected:
                    out["correct"] = True

    return out


def _external_data_file(filename: str):
    """Return the path to *filename* under the ``RXN_INSIGHT_DATA`` directory.

    Users who hold the full (unpublished) reaction-naming datasets can point the
    ``RXN_INSIGHT_DATA`` environment variable at a folder containing them; the
    package then loads those transparently. Returns ``None`` when the variable is
    unset or the file is absent, so the bundled public data (or a graceful
    fallback) is used instead.
    """
    root = os.environ.get("RXN_INSIGHT_DATA")
    if root:
        candidate = os.path.join(root, filename)
        if os.path.isfile(candidate):
            return candidate
    return None


def _load_smirks_db() -> pd.DataFrame:
    """Load and curate the SMIRKS database.

    Prefers a fuller external ``final_ordered_smirks_db.json`` (via
    ``RXN_INSIGHT_DATA``) when present, otherwise the bundled ``smirks.json``.
    """
    from rxn_insight.utils import curate_smirks

    external = _external_data_file("final_ordered_smirks_db.json")
    if external is not None:
        smirks = pd.read_json(external, orient="records", lines=True)
    else:
        smirks = pd.read_json(
            resources.files(f"{__package__}.data").joinpath("smirks.json"),
            orient="records", lines=True,
        )
    return curate_smirks(smirks)


# ── Class-name lookup ───────────────────────────────────────────────────

_class_mapping: Optional[dict] = None


_SUPERCLASS_NAMES = {
    "1": "Heteroatom Alkylation and Arylation",
    "2": "Acylation",
    "3": "C-C Coupling",
    "4": "Aromatic Heterocycle Formation",
    "5": "Deprotection",
    "6": "Reduction",
    "7": "Oxidation",
    "8": "Functional Group Interconversion",
    "9": "Functional Group Addition",
    "10": "Protection",
    "11": "Miscellaneous",
}


def _load_class_mapping() -> dict:
    """Load the reaction-class code→name mapping.

    Starts from the eleven top-level superclass names. When an external
    ``structured_mapping.json`` is available via ``RXN_INSIGHT_DATA``, the
    detailed tier names are layered on top; without it only the superclass
    (tier-1) names are known and deeper tiers fall back to the raw code.
    """
    global _class_mapping
    if _class_mapping is None:
        mapping = dict(_SUPERCLASS_NAMES)
        external = _external_data_file("structured_mapping.json")
        if external is not None:
            import json
            with open(external, encoding="utf-8") as f:
                mapping.update(json.load(f))
        _class_mapping = mapping
    return _class_mapping


def get_class_name(
    code: str,
    tier: Optional[int] = None,
    class_mapping: Optional[Union[dict, str]] = None,
) -> Union[str, Dict[str, str]]:
    """Look up human-readable names for a reaction class code.

    Args:
        code: Dot-separated class code (e.g. ``"3.1.1.2.1"``).
        tier: If given, return only the name at that tier depth (1–5).
            If ``None``, return a dict with all tier levels.
        class_mapping: Optional class code→name lookup, given either as a ``dict``
            or as a path to a JSON file of that form (e.g. a private
            ``class_names.json`` / ``structured_mapping.json``). When ``None``, the
            bundled mapping is used — which, without the private data, knows only
            the eleven superclass (tier-1) names.

    Returns:
        A single name string when *tier* is specified, or a dict like
        ``{"tier_1": "C-C Coupling", "tier_2": "Suzuki-Miyaura Coupling", ...}``
        when *tier* is ``None``.

    Examples:
        >>> get_class_name("3.1.1.2.1", tier=1)
        'C-C Coupling'
        >>> get_class_name("3.1.1.2.1", tier=2, class_mapping="class_names.json")
        'Suzuki-Miyaura Coupling'
    """
    if class_mapping is None:
        mapping = _load_class_mapping()
    elif isinstance(class_mapping, dict):
        mapping = {**_SUPERCLASS_NAMES, **class_mapping}
    else:
        import json
        with open(class_mapping, encoding="utf-8") as f:
            mapping = {**_SUPERCLASS_NAMES, **json.load(f)}
    parts = code.split(".")
    result: Dict[str, str] = {}
    for depth in range(1, len(parts) + 1):
        prefix = ".".join(parts[:depth])
        result[f"tier_{depth}"] = mapping.get(prefix, prefix)

    if tier is not None:
        return result.get(f"tier_{tier}", code)
    return result


def _sanitize_reaction(rxn: str) -> str:
    """Strip atom maps and drop the agent section from a reaction SMILES.

    Accepts ``reactants>>products``, ``reactants>agents>products``, or
    mapped variants thereof and returns a clean ``reactants>>products``
    string suitable for ``name_reaction``.
    """
    parts = rxn.split(">")
    if len(parts) == 3:
        reactants, _, products = parts
    elif len(parts) == 2:
        reactants, products = parts
    else:
        return rxn
    reactants = _ATOM_MAP_RE.sub("", reactants)
    products = _ATOM_MAP_RE.sub("", products)
    return f"{reactants}>>{products}"


def _name_reaction_from_records(
    rxn: str,
    smirks_records: List[tuple],
) -> str:
    """Name a reaction using a lightweight list of (name, smirks, nreact) tuples."""
    parts = rxn.split(">>")
    if len(parts) != 2:
        return "OtherReaction"
    reactants_smiles, products_smiles = parts
    reactants = reactants_smiles.split(".")
    products = products_smiles.split(".")

    if len(reactants) > 4 or len(products) > 4:
        return "OtherReaction"

    new_products = []
    for product in products:
        try:
            new_products.append(
                Chem.MolToSmiles(Chem.MolFromSmiles(product), isomericSmiles=False)
            )
        except Exception:
            new_products.append(product)

    react_mols = [Chem.MolFromSmiles(r) for r in reactants]
    if any(m is None for m in react_mols):
        return "OtherReaction"

    num_reactants = len(react_mols)

    for name, smirks, nreact in smirks_records:
        if nreact > num_reactants:
            continue

        rxn_obj = _get_compiled_rxn(smirks)
        if rxn_obj is None:
            continue

        if nreact == num_reactants:
            if num_reactants == 1:
                all_tuples = [tuple(react_mols)]
            else:
                all_tuples = list(itertools.permutations(react_mols))
        else:
            all_tuples = []
            for subset in itertools.combinations(react_mols, nreact):
                if nreact == 1:
                    all_tuples.append(subset)
                else:
                    all_tuples.extend(itertools.permutations(subset))

        matched = False
        for tup in all_tuples:
            try:
                outcomes = rxn_obj.RunReactants(tup)
            except Exception:
                continue
            for prods in outcomes:
                try:
                    prod = Chem.MolToSmiles(prods[0], isomericSmiles=False)
                except Exception:
                    continue
                if prod in new_products:
                    matched = True
                    break
            if matched:
                break

        if matched:
            return name.strip("{}")

    return "OtherReaction"


def _name_chunk(chunk: List[str], smirks_db) -> List[str]:
    """Worker: name a chunk of reactions.

    ``smirks_db`` is either a DataFrame (single-process path) or a list of
    ``(name, smirks, nreact)`` tuples (parallel path, much cheaper to pickle).
    """
    if isinstance(smirks_db, list):
        return [_name_reaction_from_records(rxn, smirks_db) for rxn in chunk]
    return [name_reaction(rxn, smirks_db) for rxn in chunk]


def name_reactions_batch(
    reactions: List[str],
    n_jobs: Optional[int] = None,
    smirks_db: Optional[pd.DataFrame] = None,
    progress: bool = True,
) -> List[str]:
    """Name a list of reactions by matching against the SMIRKS pattern database.

    Args:
        reactions: List of reaction SMILES strings.  Accepts mapped or
            unmapped SMILES in ``reactants>>products`` or
            ``reactants>agents>products`` format — atom maps and agents
            are stripped automatically before matching.
        n_jobs: Number of parallel workers. Defaults to half the available CPUs.
            Use ``1`` to disable parallelism.
        smirks_db: Pre-loaded SMIRKS database. If ``None``, the bundled
            ``smirks.json`` is loaded automatically.
        progress: Show a tqdm progress bar (default: True).

    Returns:
        List of reaction names, one per input reaction. Reactions that cannot
        be matched return ``'OtherReaction'``.
    """
    if smirks_db is None:
        smirks_db = _load_smirks_db()

    if n_jobs is None:
        n_jobs = max(1, int(0.5 * mp.cpu_count()))

    # Sanitize once upfront, before dispatching to workers
    reactions = [_sanitize_reaction(rxn) for rxn in reactions]

    if n_jobs == 1:
        return [
            name_reaction(rxn, smirks_db)
            for rxn in tqdm(reactions, desc="Naming reactions", disable=not progress)
        ]

    # Convert DataFrame to lightweight tuples for cheaper IPC serialization
    records = list(zip(
        smirks_db["name"].tolist(),
        smirks_db["smirks"].tolist(),
        smirks_db["nreact"].astype(int).tolist(),
    ))

    # Chunk to reduce IPC overhead and give tqdm sensible granularity
    n = len(reactions)
    effective_jobs = n_jobs if n_jobs > 0 else mp.cpu_count()
    chunk_size = max(50, n // (effective_jobs * 8))
    chunks = [reactions[i: i + chunk_size] for i in range(0, n, chunk_size)]

    raw_chunks: List[List[str]] = Parallel(n_jobs=n_jobs)(
        delayed(_name_chunk)(chunk, records)
        for chunk in tqdm(chunks, desc="Naming reactions", disable=not progress)
    )

    return [name for chunk_res in raw_chunks for name in chunk_res]


