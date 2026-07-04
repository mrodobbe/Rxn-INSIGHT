import itertools
import multiprocessing as mp
import os
import re
from collections import Counter, defaultdict, deque
from importlib import resources
from typing import Dict, List, Optional, Set, Tuple, Union

from rdkit import Chem
from rdkit.Chem import AllChem
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd

_ATOM_MAP_RE = re.compile(r":\d+")

# --- SMIRKS generalisation constants ---
_HALIDES = frozenset({"F", "Cl", "Br", "I"})
_ATOM_BLOCK_PAT = re.compile(r"\[([^\]]+)\]")
_MAP_END_RE = re.compile(r":(\d+)$")


def _get_mapped_neighbors(smirks_side: str) -> Dict[int, Set[int]]:
    """Return adjacency of mapped atoms on one side of a SMIRKS.

    Uses RDKit SMARTS parser to build the molecular graph, then for each
    mapped atom returns the set of map numbers of its directly bonded
    mapped neighbours.  Unmapped neighbours are ignored.

    Returns:
        ``{map_num: {neighbour_map_num, ...}}`` for every mapped atom.
    """
    mol = Chem.MolFromSmarts(smirks_side)
    if mol is None:
        return {}

    map_to_idx: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn > 0:
            map_to_idx[mn] = atom.GetIdx()

    adj: Dict[int, Set[int]] = {mn: set() for mn in map_to_idx}
    for mn, idx in map_to_idx.items():
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            nbr_mn = nbr.GetAtomMapNum()
            if nbr_mn > 0:
                adj[mn].add(nbr_mn)
    return adj


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


# ===================================================================
# SMIRKS Generalisation
# ===================================================================


def _parse_atom_block(content: str) -> dict:
    """Parse SMARTS atom block content (text between ``[`` and ``]``).

    Handles formats produced by ``get_atomic_smarts``::

        C;H0;D3;+0:1   →  element C, H=0, D=3, charge=0, map=1
        c;H1;D2;+0:3   →  aromatic c, H=1, D=2, charge=0, map=3
        Cl;H0;D1;+0     →  unmapped Cl leaving group
        *:3              →  wildcard with map 3

    Args:
        content: Text between ``[`` and ``]``.

    Returns:
        Dict with keys: elements, aromatic, H, D, charge, map_num,
        stereo, is_wildcard, raw.
    """
    result: dict = {
        "elements": [],
        "aromatic": False,
        "H": None,
        "D": None,
        "charge": None,
        "map_num": None,
        "stereo": None,
        "is_wildcard": False,
        "raw": content,
    }

    # Extract map number at the very end  (:N)
    m = _MAP_END_RE.search(content)
    if m:
        result["map_num"] = int(m.group(1))
        content = content[: m.start()]

    # Wildcard
    if content.strip() == "*":
        result["is_wildcard"] = True
        result["elements"] = ["*"]
        return result

    # Split by semicolons
    parts = [p.strip() for p in content.split(";") if p.strip()]
    if not parts:
        return result

    # First part: element(s), possibly comma-separated
    elements = [e.strip() for e in parts[0].split(",")]
    result["elements"] = elements
    result["aromatic"] = elements[0][0].islower() if elements[0] else False

    # Remaining parts: H, D, charge, stereo, ring (ring ignored)
    for part in parts[1:]:
        hm = re.match(r"^H(\d+)$", part)
        if hm:
            result["H"] = int(hm.group(1))
            continue
        dm = re.match(r"^D(\d+)$", part)
        if dm:
            result["D"] = int(dm.group(1))
            continue
        cm = re.match(r"^([+-]\d+)$", part)
        if cm:
            result["charge"] = int(cm.group(1))
            continue
        if part in ("@", "@@"):
            result["stereo"] = part
            continue
        # Ring info (!R, R1, r6) — ignored for generalisation

    return result


def _get_atomic_number(elem: str) -> int:
    """Get atomic number from an element symbol, handling aromatic lowercase."""
    canonical = elem.capitalize()
    try:
        return Chem.GetPeriodicTable().GetAtomicNumber(canonical)
    except Exception:
        return 0


_MAP_IN_BRACKET_RE = re.compile(r":(\d+)(?=\])")


def _canonical_renumber(smirks: str) -> str:
    """Renumber atom maps to a canonical order based on atom roles.

    Templates that represent the same transformation but use different
    map number assignments (e.g. N at ``:1`` vs N at ``:2``) get
    identical numbering after this step, enabling correct grouping.

    Sort key per map number (in priority order):
    RC atoms first (detected via H/charge *changes* AND connectivity
    changes between R and P sides), then by element (atomic number).
    RC atoms additionally sort by R/P H counts; context atoms sort
    only by element and original map number to ensure ring/non-ring
    variants align.
    """
    parts = smirks.split(">>")
    if len(parts) != 2:
        return smirks

    r_atoms: Dict[int, dict] = {}
    p_atoms: Dict[int, dict] = {}
    for m in _ATOM_BLOCK_PAT.finditer(parts[0]):
        parsed = _parse_atom_block(m.group(1))
        if parsed["map_num"] is not None:
            r_atoms[parsed["map_num"]] = parsed
    for m in _ATOM_BLOCK_PAT.finditer(parts[1]):
        parsed = _parse_atom_block(m.group(1))
        if parsed["map_num"] is not None:
            p_atoms[parsed["map_num"]] = parsed

    all_maps = sorted(set(r_atoms) | set(p_atoms))
    if not all_maps:
        return smirks

    # RC detection for canonical renumbering: use R/P property changes
    # (H change, charge change) AND connectivity changes (different
    # mapped neighbours between R and P sides — detects bond-forming
    # atoms like the alpha C in N-alkylation or Br in halogenation).
    rc: Set[int] = set()
    for mn in all_maps:
        ra, pa = r_atoms.get(mn), p_atoms.get(mn)
        if ra is None or pa is None:
            rc.add(mn)
            continue
        if ra["H"] != pa["H"] or ra["charge"] != pa["charge"]:
            rc.add(mn)

    # Connectivity-based RC: atoms whose bonded mapped neighbours change
    r_adj = _get_mapped_neighbors(parts[0])
    p_adj = _get_mapped_neighbors(parts[1])
    for mn in all_maps:
        if mn not in rc:
            if r_adj.get(mn, set()) != p_adj.get(mn, set()):
                rc.add(mn)

    def _sort_key(mn: int):
        ra, pa = r_atoms.get(mn), p_atoms.get(mn)
        ref = ra or pa
        anum = _get_atomic_number(ref["elements"][0]) if ref else 0
        if mn in rc:
            # RC atoms: sort by element, then R-side H, P-side H
            return (
                0,
                anum,
                ra["H"] if ra and ra["H"] is not None else 99,
                pa["H"] if pa and pa["H"] is not None else 99,
                mn,
            )
        # Context atoms: sort by element and original map only —
        # no H/D/charge so that ring vs non-ring variants align.
        return (1, anum, mn)

    sorted_maps = sorted(all_maps, key=_sort_key)
    old_to_new = {old: idx + 1 for idx, old in enumerate(sorted_maps)}

    if all(old_to_new[m] == m for m in all_maps):
        return smirks

    def _replace_map(match: re.Match) -> str:
        return f":{old_to_new.get(int(match.group(1)), int(match.group(1)))}"

    return _MAP_IN_BRACKET_RE.sub(_replace_map, smirks)


def _mapped_atom_signature(smirks: str) -> str:
    """Compute a topology signature based on mapped atoms and halide LGs only.

    Unmapped non-halide context atoms (ring tails, branches) are ignored so
    that templates differing only in context neighbourhood are grouped
    together.  The signature is a sorted atom inventory per side — sufficient
    for same-name templates that always share the same connectivity.
    """
    sides = smirks.split(">>")
    if len(sides) != 2:
        return smirks

    sig_parts: List[str] = []
    for side in sides:
        atoms: List[str] = []
        has_halide = False
        for m in _ATOM_BLOCK_PAT.finditer(side):
            parsed = _parse_atom_block(m.group(1))
            if parsed["is_wildcard"]:
                if parsed["map_num"] is not None:
                    atoms.append(f"*:{parsed['map_num']}")
            elif parsed["map_num"] is not None:
                anum = _get_atomic_number(parsed["elements"][0])
                atoms.append(f"#{anum}:{parsed['map_num']}")
            elif parsed.get("elements") and set(parsed["elements"]) <= _HALIDES:
                has_halide = True
        atoms.sort()
        sig_parts.append(",".join(atoms) + ("|X" if has_halide else ""))

    return ">>".join(sig_parts)


def _rc_signature(smirks: str, rc_maps: Set[int]) -> str:
    """Compute a grouping signature based on RC atoms only.

    Unlike ``_mapped_atom_signature`` which includes ALL mapped atoms,
    this signature only considers atoms whose properties change between
    R and P sides (reaction-centre atoms).  This allows templates with
    different numbers of context atoms (e.g. ring vs non-ring variants)
    to be grouped together, since their RC atoms are identical.

    The signature encodes per-side: element, H count, and canonical map
    number of each RC atom, plus component counts.  Unmapped atoms
    (including halide LGs) are NOT included — they are reagent details,
    not part of the core transformation.
    """
    sides = smirks.split(">>")
    if len(sides) != 2:
        return smirks

    nreact = _count_components(sides[0])
    nprod = _count_components(sides[1])

    sig_parts: List[str] = []
    for side in sides:
        atoms: List[str] = []
        for m in _ATOM_BLOCK_PAT.finditer(side):
            parsed = _parse_atom_block(m.group(1))
            if parsed["map_num"] is not None and parsed["map_num"] in rc_maps:
                if parsed["is_wildcard"]:
                    atoms.append(f"*:{parsed['map_num']}")
                else:
                    anum = _get_atomic_number(parsed["elements"][0])
                    h = parsed["H"] if parsed["H"] is not None else ""
                    atoms.append(f"#{anum}H{h}:{parsed['map_num']}")
        atoms.sort()
        sig_parts.append(",".join(atoms))

    return f"{nreact}>{nprod}|" + ">>".join(sig_parts)


def _tokenize_smirks(smirks: str) -> List[Tuple[str, str]]:
    """Split SMIRKS into ``(type, value)`` token pairs.

    Types: ``"atom"`` for ``[…]`` blocks, ``"sep"`` for everything else
    (bonds, parentheses, ``>>``).
    """
    tokens: List[Tuple[str, str]] = []
    last = 0
    for m in _ATOM_BLOCK_PAT.finditer(smirks):
        if m.start() > last:
            tokens.append(("sep", smirks[last: m.start()]))
        tokens.append(("atom", m.group(0)))
        last = m.end()
    if last < len(smirks):
        tokens.append(("sep", smirks[last:]))
    return tokens


def _identify_rc_maps(smirks: str) -> Set[int]:
    """Identify reaction-centre map numbers.

    Uses two signals (checked in order):

    1. **D (degree) presence** — ``relax_context=True`` omits D for context
       atoms, so any atom with D in its SMARTS is RC.
    2. **R/P property changes** — fallback for non-relaxed templates:
       atoms whose H, D, or charge differ between reactant and product sides.
    """
    parts = smirks.split(">>")
    if len(parts) != 2:
        return set()

    r_atoms: Dict[int, dict] = {}
    p_atoms: Dict[int, dict] = {}

    for m in _ATOM_BLOCK_PAT.finditer(parts[0]):
        parsed = _parse_atom_block(m.group(1))
        if parsed["map_num"] is not None:
            r_atoms[parsed["map_num"]] = parsed

    for m in _ATOM_BLOCK_PAT.finditer(parts[1]):
        parsed = _parse_atom_block(m.group(1))
        if parsed["map_num"] is not None:
            p_atoms[parsed["map_num"]] = parsed

    rc: Set[int] = set()
    for map_num in set(r_atoms) | set(p_atoms):
        ra = r_atoms.get(map_num)
        pa = p_atoms.get(map_num)

        # Only on one side → RC
        if ra is None or pa is None:
            rc.add(map_num)
            continue

        # D present on either side → RC (relax_context discriminator)
        if ra["D"] is not None or pa["D"] is not None:
            rc.add(map_num)
            continue

        # Property changes → RC (fallback for non-relaxed templates)
        if ra["H"] != pa["H"] or ra["charge"] != pa["charge"]:
            rc.add(map_num)

    return rc


def _generalize_leaving_group(blocks: List[dict]) -> str:
    """Generalise unmapped halide leaving-group atoms.

    Produces ``[F,Cl,Br,I]`` list notation sorted by atomic number.
    """
    if not blocks:
        return "[*]"

    all_elements: Set[str] = set()
    for b in blocks:
        all_elements.update(b.get("elements", []))

    if all_elements and all_elements <= _HALIDES:
        order = ["F", "Cl", "Br", "I"]
        sorted_elems = [h for h in order if h in all_elements]
        # Keep H and charge if present
        h_counts = Counter(b["H"] for b in blocks if b["H"] is not None)
        charge_counts = Counter(
            b["charge"] for b in blocks if b["charge"] is not None
        )
        parts_list = [",".join(sorted_elems)]
        if h_counts:
            parts_list.append(f"H{h_counts.most_common(1)[0][0]}")
        if charge_counts:
            chg = charge_counts.most_common(1)[0][0]
            parts_list.append(f"+{chg}" if chg >= 0 else str(chg))
        return "[" + ";".join(parts_list) + "]"

    return "[*]"


def _generalize_atom_block(
    blocks: List[dict],
    is_rc: bool,
    context_level: str,
) -> str:
    """Create a generalised SMARTS atom block from multiple instances.

    Args:
        blocks: Parsed atom blocks at the same position across templates.
        is_rc: Whether this atom is in the reaction centre.
        context_level: ``"wildcard"`` replaces context atoms with ``[*:N]``;
            ``"element"`` keeps element + H (if uniform) + charge, using
            ``#N`` notation when aromaticity varies.
    """
    if not blocks:
        return "[*]"

    map_num = blocks[0].get("map_num")

    # Wildcards stay as-is
    if all(b.get("is_wildcard") for b in blocks):
        return f"[*:{map_num}]" if map_num is not None else "[*]"

    # --- Context atoms (not RC) ---
    if not is_rc:
        if context_level == "wildcard":
            return f"[*:{map_num}]" if map_num is not None else "[*]"

        # "element" level: smart generalisation
        aromaticities = set(b["aromatic"] for b in blocks)
        mixed_arom = len(aromaticities) > 1
        if mixed_arom:
            # Mix of aromatic and aliphatic → use #N
            anum = _get_atomic_number(blocks[0]["elements"][0])
            elem_str = f"#{anum}"
        else:
            elem_str = blocks[0]["elements"][0]

        # H: majority vote — but drop when aromaticity varies (#N mode)
        # because H counts are unreliable across aromatic/aliphatic
        if mixed_arom:
            h_part = ""
        else:
            h_counts_ctx = Counter(
                b["H"] for b in blocks if b["H"] is not None
            )
            if h_counts_ctx:
                h_part = f";H{h_counts_ctx.most_common(1)[0][0]}"
            else:
                h_part = ""

        # Charge: majority vote
        charge_counts = Counter(
            b["charge"] for b in blocks if b["charge"] is not None
        )
        chg = charge_counts.most_common(1)[0][0] if charge_counts else 0
        chg_part = f";+{chg}" if chg >= 0 else f";{chg}"

        smarts = f"[{elem_str}{h_part}{chg_part}"
        if map_num is not None:
            smarts += f":{map_num}"
        smarts += "]"
        return smarts

    # --- Unmapped atoms → halide LG or wildcard ---
    if map_num is None:
        all_elements: Set[str] = set()
        for b in blocks:
            all_elements.update(b.get("elements", []))
        if all_elements and all_elements <= _HALIDES:
            return _generalize_leaving_group(blocks)
        return "[*]"

    # --- RC atom (mapped): keep element, H, D, charge (majority vote) ---
    elements_rc: Set[str] = set()
    for b in blocks:
        elements_rc.update(b.get("elements", []))

    aromaticities = set(b["aromatic"] for b in blocks)
    if len(aromaticities) > 1:
        anum = _get_atomic_number(blocks[0]["elements"][0])
        elem_str = f"#{anum}"
    else:
        elem_str = ",".join(sorted(elements_rc))

    h_counts = Counter(b["H"] for b in blocks if b["H"] is not None)
    d_counts = Counter(b["D"] for b in blocks if b["D"] is not None)
    charge_counts = Counter(b["charge"] for b in blocks if b["charge"] is not None)

    parts_list = [elem_str]
    if h_counts:
        parts_list.append(f"H{h_counts.most_common(1)[0][0]}")
    if d_counts:
        parts_list.append(f"D{d_counts.most_common(1)[0][0]}")
    if charge_counts:
        chg = charge_counts.most_common(1)[0][0]
        parts_list.append(f"+{chg}" if chg >= 0 else str(chg))

    smarts = "[" + ";".join(parts_list)
    smarts += f":{map_num}]"
    return smarts


def _count_components(smirks_side: str) -> int:
    """Count disconnected components in one side of a SMIRKS.

    Counts top-level ``.`` separators (outside brackets and parentheses).
    """
    depth = 0
    in_bracket = False
    count = 1
    for ch in smirks_side:
        if ch == "[":
            in_bracket = True
        elif ch == "]":
            in_bracket = False
        elif not in_bracket:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "." and depth == 0:
                count += 1
    return count


def _build_rc_only_smirks(
    rc_maps: Set[int],
    gen_r: Dict[int, str],
    gen_p: Dict[int, str],
    skeleton: str,
    gen_halide: Optional[str] = None,
) -> Optional[str]:
    """Build SMIRKS keeping RC atoms and unmapped atoms, dropping context.

    Unmapped atoms reachable from RC atoms (leaving groups, functional
    handles like B(OH)2) are preserved.  Only mapped non-RC (context)
    atoms are stripped.  Unmapped single-atom halides are replaced with
    *gen_halide* when provided.

    Returns ``None`` if parsing fails or the result is empty.
    """
    skel_parts = skeleton.split(">>")
    if len(skel_parts) != 2:
        return None

    _BOND_SYM = {
        Chem.BondType.SINGLE: "-",
        Chem.BondType.DOUBLE: "=",
        Chem.BondType.TRIPLE: "#",
        Chem.BondType.AROMATIC: ":",
    }

    # Original component counts from skeleton (needed to pad with [*]).
    orig_ncomp = [_count_components(s) for s in skel_parts]

    side_strs: List[str] = []
    for side_idx, (side, gen_map) in enumerate(
        zip(skel_parts, [gen_r, gen_p])
    ):
        mol = Chem.MolFromSmarts(side)
        if mol is None:
            return None

        n_atoms = mol.GetNumAtoms()

        # Classify each atom: "rc" / "context" / "unmapped"
        role: List[str] = [""] * n_atoms
        idx_to_map: Dict[int, int] = {}
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            mn = atom.GetAtomMapNum()
            if mn > 0:
                idx_to_map[i] = mn
                role[i] = "rc" if mn in rc_maps else "context"
            else:
                role[i] = "unmapped"

        # Original atom blocks from skeleton string (positional — index i
        # corresponds to the i-th [...] block in the SMARTS string, which
        # matches atom index i from MolFromSmarts).
        orig_blocks: List[str] = [
            m.group(0) for m in _ATOM_BLOCK_PAT.finditer(side)
        ]

        # Kept atoms: RC + unmapped reachable from RC via BFS (never
        # crossing through context atoms).
        kept: Set[int] = set()
        bfs_queue: deque = deque()
        for i in range(n_atoms):
            if role[i] == "rc":
                kept.add(i)
                bfs_queue.append(i)
        while bfs_queue:
            i = bfs_queue.popleft()
            for nbr in mol.GetAtomWithIdx(i).GetNeighbors():
                j = nbr.GetIdx()
                if j not in kept and role[j] == "unmapped":
                    kept.add(j)
                    bfs_queue.append(j)

        if not any(role[i] == "rc" for i in kept):
            return None

        # Adjacency among kept atoms
        adj: Dict[int, List[Tuple[int, str]]] = defaultdict(list)
        for bond in mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a1 in kept and a2 in kept:
                bsym = _BOND_SYM.get(bond.GetBondType(), "~")
                adj[a1].append((a2, bsym))
                adj[a2].append((a1, bsym))

        # Resolve atom block string for a kept atom index.
        def _get_block(idx: int) -> str:
            if role[idx] == "rc":
                return gen_map.get(idx_to_map[idx], f"[*:{idx_to_map[idx]}]")
            # Unmapped: replace halide LGs with merged notation
            if gen_halide and idx < len(orig_blocks):
                inner = orig_blocks[idx][1:-1]
                parsed = _parse_atom_block(inner)
                if (parsed.get("elements")
                        and set(parsed["elements"]) <= _HALIDES):
                    return gen_halide
            return orig_blocks[idx] if idx < len(orig_blocks) else "[*]"

        # Recursive DFS producing SMIRKS with proper branch notation.
        visited: Set[int] = set()

        def _dfs(idx: int) -> str:
            visited.add(idx)
            result = _get_block(idx)
            unvisited = [
                (nbr, bsym) for nbr, bsym in adj[idx]
                if nbr not in visited
            ]
            if not unvisited:
                return result
            # Last neighbour continues the main chain; the rest become
            # parenthesised branches.
            for i, (nbr, bsym) in enumerate(unvisited):
                if nbr in visited:
                    continue
                sub = _dfs(nbr)
                if i < len(unvisited) - 1:
                    result += f"({bsym}{sub})"
                else:
                    result += f"{bsym}{sub}"
            return result

        # Build components starting from RC atoms sorted by map number.
        rc_starts = sorted(
            [i for i in kept if role[i] == "rc"],
            key=lambda i: idx_to_map[i],
        )
        components: List[str] = []
        for start in rc_starts:
            if start in visited:
                continue
            components.append(_dfs(start))
        # Pick up any remaining unmapped-only components (unlikely).
        for i in sorted(kept):
            if i not in visited:
                components.append(_dfs(i))

        # Pad with [*] wildcards to preserve the original component count.
        # Dropped components (e.g. NBS in halogenation) had no RC atoms
        # and no reachable unmapped atoms, but RunReactants still needs the
        # right number of reactant templates for nreact filtering to work.
        while len(components) < orig_ncomp[side_idx]:
            components.append("[*]")

        side_strs.append(".".join(components))

    return f"{side_strs[0]}>>{side_strs[1]}"


def _generalize_smirks_core(
    template_counts: Dict[str, int],
    max_output: int,
    context_level: str,
) -> List[Tuple[str, int, int]]:
    """Core SMIRKS generalisation returning full statistics.

    Returns:
        List of ``(generalised_smirks, group_total_count,
        group_n_templates)`` tuples sorted by descending count.
    """
    if not template_counts:
        return []

    # Phase 0: Canonical renumbering — normalise map assignments so that
    # templates representing the same transformation but with different
    # original map numbering get identical map numbers.
    renumbered_counts: Dict[str, int] = {}
    for template, count in template_counts.items():
        canon = _canonical_renumber(template)
        renumbered_counts[canon] = renumbered_counts.get(canon, 0) + count

    # Phase 0.5: Compute RC maps GLOBALLY across ALL (renumbered) templates.
    # An atom is RC if any of its properties CHANGE between R and P sides
    # (H, D, or charge) in the MAJORITY of templates.  Using property
    # *changes* (not D-presence) is critical because full-specificity
    # templates have D on every atom — D-presence would mark everything RC.
    change_votes: Dict[int, int] = defaultdict(int)
    no_change_votes: Dict[int, int] = defaultdict(int)
    one_side_only: Set[int] = set()

    for template in renumbered_counts:
        parts = template.split(">>")
        if len(parts) != 2:
            continue
        r_atoms: Dict[int, dict] = {}
        p_atoms: Dict[int, dict] = {}
        for m in _ATOM_BLOCK_PAT.finditer(parts[0]):
            parsed = _parse_atom_block(m.group(1))
            if parsed["map_num"] is not None:
                r_atoms[parsed["map_num"]] = parsed
        for m in _ATOM_BLOCK_PAT.finditer(parts[1]):
            parsed = _parse_atom_block(m.group(1))
            if parsed["map_num"] is not None:
                p_atoms[parsed["map_num"]] = parsed

        for map_num in set(r_atoms) | set(p_atoms):
            ra = r_atoms.get(map_num)
            pa = p_atoms.get(map_num)
            if ra is None or pa is None:
                one_side_only.add(map_num)
                continue
            changed = (
                ra["H"] != pa["H"]
                or ra["D"] != pa["D"]
                or ra["charge"] != pa["charge"]
            )
            if changed:
                change_votes[map_num] += 1
            else:
                no_change_votes[map_num] += 1

    rc_maps: Set[int] = set(one_side_only)
    for map_num in set(change_votes) | set(no_change_votes):
        if change_votes.get(map_num, 0) > no_change_votes.get(map_num, 0):
            rc_maps.add(map_num)

    # Phase 0.6: Connectivity-based RC detection — atoms whose set of
    # bonded mapped neighbours changes between R and P are RC even if
    # their own properties (H, D, charge) are unchanged (e.g. Br in
    # aromatic bromination: Br-Br → Ar-Br, properties same, bonds differ).
    conn_change_votes: Dict[int, int] = defaultdict(int)
    conn_no_change_votes: Dict[int, int] = defaultdict(int)
    for template in renumbered_counts:
        parts = template.split(">>")
        if len(parts) != 2:
            continue
        r_adj = _get_mapped_neighbors(parts[0])
        p_adj = _get_mapped_neighbors(parts[1])
        all_maps_t = set(r_adj) | set(p_adj)
        for mn in all_maps_t:
            if mn in rc_maps:
                continue  # already RC, skip
            r_nbrs = r_adj.get(mn, set())
            p_nbrs = p_adj.get(mn, set())
            if r_nbrs != p_nbrs:
                conn_change_votes[mn] += 1
            else:
                conn_no_change_votes[mn] += 1

    for mn in set(conn_change_votes) | set(conn_no_change_votes):
        if conn_change_votes.get(mn, 0) > conn_no_change_votes.get(mn, 0):
            rc_maps.add(mn)

    # Phase 1: Group by RC-only signature — only RC atoms (whose properties
    # change between R and P) are included.  This allows ring and non-ring
    # variants of the same transformation to be grouped together.
    groups: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    for template, count in renumbered_counts.items():
        sig = _rc_signature(template, rc_maps)
        groups[sig].append((template, count))

    # Phase 2: Generalise each group using map-number-based matching
    results: List[Tuple[str, int, int]] = []
    for _sig, group_templates in groups.items():
        total_count = sum(c for _, c in group_templates)

        # Collect parsed atoms by map number and side across ALL templates
        r_by_map: Dict[int, List[dict]] = defaultdict(list)
        p_by_map: Dict[int, List[dict]] = defaultdict(list)
        all_halide_blocks: List[dict] = []

        for template, _count in group_templates:
            tparts = template.split(">>")
            if len(tparts) != 2:
                continue
            for m_obj in _ATOM_BLOCK_PAT.finditer(tparts[0]):
                parsed = _parse_atom_block(m_obj.group(1))
                if parsed["map_num"] is not None:
                    r_by_map[parsed["map_num"]].append(parsed)
                elif (parsed.get("elements")
                      and set(parsed["elements"]) <= _HALIDES):
                    all_halide_blocks.append(parsed)
            for m_obj in _ATOM_BLOCK_PAT.finditer(tparts[1]):
                parsed = _parse_atom_block(m_obj.group(1))
                if parsed["map_num"] is not None:
                    p_by_map[parsed["map_num"]].append(parsed)

        # Generalise each mapped position (separately for R and P sides)
        gen_r: Dict[int, str] = {}
        gen_p: Dict[int, str] = {}
        for map_num in set(r_by_map) | set(p_by_map):
            is_rc = map_num in rc_maps
            if map_num in r_by_map:
                gen_r[map_num] = _generalize_atom_block(
                    r_by_map[map_num], is_rc, context_level,
                )
            if map_num in p_by_map:
                gen_p[map_num] = _generalize_atom_block(
                    p_by_map[map_num], is_rc, context_level,
                )

        # Generalise halide leaving groups (across all templates in group)
        gen_halide = (
            _generalize_leaving_group(all_halide_blocks)
            if all_halide_blocks else None
        )

        # Find the simplest template (fewest atom blocks) as skeleton
        skeleton = min(
            group_templates, key=lambda x: len(_ATOM_BLOCK_PAT.findall(x[0]))
        )[0]

        if context_level == "rc_only":
            # Build minimal SMIRKS: RC atoms + unmapped atoms, no context
            result_smirks = _build_rc_only_smirks(
                rc_maps, gen_r, gen_p, skeleton, gen_halide,
            )
            if result_smirks is None:
                continue
        else:
            skel_parts = skeleton.split(">>")
            if len(skel_parts) != 2:
                continue

            # Replace atoms in the skeleton by map number
            def _replace_side(side: str, gen_map: Dict[int, str]) -> str:
                def _replacer(match: re.Match) -> str:
                    parsed = _parse_atom_block(match.group(1))
                    if parsed["map_num"] is not None:
                        return gen_map.get(
                            parsed["map_num"], match.group(0),
                        )
                    # Unmapped atom — keep as-is unless it's a halide LG
                    elems = set(parsed.get("elements", []))
                    if elems and elems <= _HALIDES and gen_halide:
                        return gen_halide
                    return match.group(0)
                return _ATOM_BLOCK_PAT.sub(_replacer, side)

            gen_r_side = _replace_side(skel_parts[0], gen_r)
            gen_p_side = _replace_side(skel_parts[1], gen_p)
            result_smirks = f"{gen_r_side}>>{gen_p_side}"

        # Validate with RDKit
        try:
            rxn = AllChem.ReactionFromSmarts(result_smirks)
            if rxn is None:
                continue
        except Exception:
            continue

        results.append((result_smirks, total_count, len(group_templates)))

    # Phase 3: Sort by count descending, return top N
    results.sort(key=lambda x: x[1], reverse=True)
    return results[:max_output]


def generalize_smirks(
    templates: Union[List[str], Dict[str, int]],
    max_output: int = 5,
    context_level: str = "wildcard",
) -> List[str]:
    """Generalise many specific SMIRKS into fewer general patterns.

    Takes a collection of specific SMIRKS templates (e.g. extracted at
    ``rr0rp1_ring0``) and reduces them to a small number of general
    patterns by:

    1. Grouping templates by topology signature (same connectivity).
    2. Merging atom properties within each group (majority vote for
       reaction-centre atoms, wildcards for context atoms).
    3. Ranking by frequency and returning the top patterns.

    Args:
        templates: SMIRKS template strings, either as a list (each counted
            once) or a ``{template: count}`` dict.
        max_output: Maximum number of generalised SMIRKS to return.
        context_level: How to treat context (non-RC) atoms.
            ``"wildcard"`` → ``[*:N]``; ``"element"`` → keep element + charge.

    Returns:
        List of generalised SMIRKS strings, sorted by descending frequency.
    """
    if isinstance(templates, dict):
        template_counts = templates
    else:
        template_counts = {t: 1 for t in templates}

    return [
        s for s, _, _ in _generalize_smirks_core(
            template_counts, max_output, context_level
        )
    ]


def _invert_template_to_name(
    template_to_name: Dict[str, Counter],
    ambiguity: str = "majority",
) -> Tuple[Dict[str, Counter], Dict[str, List[str]]]:
    """Invert a template → Counter(names) mapping to name → Counter(templates).

    Args:
        template_to_name: Maps each SMIRKS template to a Counter of names
            it is associated with.
        ambiguity: How to handle templates mapping to multiple names.
            ``"majority"`` — assign to the most frequent name.
            ``"all"`` — assign to every associated name.
            ``"unambiguous"`` — skip ambiguous templates entirely.

    Returns:
        ``(name_to_templates, ambiguous_report)`` where *ambiguous_report*
        maps each ambiguous template to its list of candidate names.
    """
    name_to_templates: Dict[str, Counter] = defaultdict(Counter)
    ambiguous_report: Dict[str, List[str]] = {}

    for template, name_counter in template_to_name.items():
        names = list(name_counter.keys())
        if len(names) > 1:
            ambiguous_report[template] = names

        if ambiguity == "majority":
            winner = name_counter.most_common(1)[0][0]
            total = sum(name_counter.values())
            name_to_templates[winner][template] += total
        elif ambiguity == "all":
            for name, count in name_counter.items():
                name_to_templates[name][template] += count
        elif ambiguity == "unambiguous":
            if len(names) == 1:
                name = names[0]
                name_to_templates[name][template] += name_counter[name]
        else:
            raise ValueError(f"Unknown ambiguity mode: {ambiguity!r}")

    return dict(name_to_templates), ambiguous_report


def _generalize_for_name(
    name: str,
    template_counts: Counter,
    max_output: int,
    context_level: str,
    n_ambiguous: int,
) -> List[dict]:
    """Worker: generalise templates for a single reaction name."""
    results = _generalize_smirks_core(
        dict(template_counts), max_output, context_level
    )
    total_input = sum(template_counts.values())
    n_input_templates = len(template_counts)

    rows = []
    for smirks, group_count, _group_n_templates in results:
        parts = smirks.split(">>")
        nreact = _count_components(parts[0]) if len(parts) == 2 else 1
        nproduct = _count_components(parts[1]) if len(parts) == 2 else 1
        coverage = group_count / total_input if total_input > 0 else 0.0

        rows.append({
            "name": name,
            "smirks": smirks,
            "nreact": nreact,
            "nproduct": nproduct,
            "coverage": round(coverage, 4),
            "n_input_templates": n_input_templates,
            "n_ambiguous": n_ambiguous,
        })

    return rows


def build_smirks_db(
    template_to_name: Dict[str, Counter],
    max_output: int = 5,
    context_level: str = "wildcard",
    ambiguity: str = "majority",
    n_jobs: Optional[int] = None,
    progress: bool = True,
) -> pd.DataFrame:
    """Build a generalised SMIRKS naming database from template → name mapping.

    This is the main entry point for creating a compact reaction-naming
    database from a large set of specific SMIRKS templates.

    Args:
        template_to_name: Maps each SMIRKS template to a ``Counter`` of
            reaction names it is associated with.  Most templates map to a
            single name, but ambiguous templates are handled according to
            *ambiguity*.
        max_output: Maximum number of generalised SMIRKS per reaction name.
        context_level: ``"wildcard"`` or ``"element"`` — how to treat
            context atoms (non-RC atoms added by radius expansion).
        ambiguity: How to handle templates mapping to multiple names.
            ``"majority"`` — assign to the most frequent name.
            ``"all"`` — assign to every associated name.
            ``"unambiguous"`` — skip ambiguous templates entirely.
        n_jobs: Number of parallel workers. ``None`` defaults to half the
            available CPUs.
        progress: Show a tqdm progress bar.

    Returns:
        DataFrame with columns ``[name, smirks, nreact, nproduct, coverage,
        n_input_templates, n_ambiguous]``.
    """
    # Step 1: Invert mapping
    name_to_templates, ambiguous_report = _invert_template_to_name(
        template_to_name, ambiguity=ambiguity
    )

    # Count ambiguous templates per name
    ambig_per_name: Dict[str, int] = defaultdict(int)
    for template, names in ambiguous_report.items():
        for name in names:
            ambig_per_name[name] += 1

    names = list(name_to_templates.keys())

    if n_jobs is None:
        n_jobs = max(1, int(0.5 * mp.cpu_count()))

    # Step 2: Generalise per name
    if n_jobs == 1 or len(names) <= 4:
        all_rows: List[dict] = []
        for name in tqdm(names, desc="Generalising SMIRKS", disable=not progress):
            rows = _generalize_for_name(
                name,
                name_to_templates[name],
                max_output,
                context_level,
                ambig_per_name.get(name, 0),
            )
            all_rows.extend(rows)
    else:
        raw_results = Parallel(n_jobs=n_jobs)(
            delayed(_generalize_for_name)(
                name,
                name_to_templates[name],
                max_output,
                context_level,
                ambig_per_name.get(name, 0),
            )
            for name in tqdm(
                names, desc="Generalising SMIRKS", disable=not progress
            )
        )
        all_rows = [row for rows in raw_results for row in rows]

    if not all_rows:
        return pd.DataFrame(
            columns=[
                "name", "smirks", "nreact", "nproduct", "coverage",
                "n_input_templates", "n_ambiguous",
            ]
        )

    return pd.DataFrame(all_rows)
