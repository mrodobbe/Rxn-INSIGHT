"""
Module for integration with the Open Reaction Database (ORD).

This module provides functionality to load, process, and analyze chemical reactions
from the Open Reaction Database (ORD) protocol buffer files. It extends Rxn-INSIGHT's
Database class to handle the specialized ORD format.

Dependencies:
    - ord_schema: The Open Reaction Database schema package
    - protoc-wheel-0: Protocol buffers compiler
    - google.protobuf: Protocol buffers library
"""

try:
    from ord_schema.message_helpers import load_message
    from ord_schema.proto import dataset_pb2
    from google.protobuf.json_format import MessageToJson
except ImportError:
    print(
        """
        Open-Reaction-Database modules are missing. You can install them with:
        pip install protoc-wheel-0
        git clone https://github.com/Open-Reaction-Database/ord-schema.git
        cd ord_schema
        python setup.py install
        """
    )

import json
import pandas as pd
from rxn_insight.database import Database


class ORDDatabase(Database):
    """
    A class to manage and analyze datasets from the Open Reaction Database.

    This class extends Rxn-INSIGHT's Database class to provide specialized
    functionality for loading, processing, and analyzing ORD protocol buffer files.
    It handles the conversion from ORD's complex structure to Rxn-INSIGHT's format.

    Attributes:
        ord_file (str): Path to the ORD protocol buffer file.
        ord_dataset: Loaded ORD dataset as a protocol buffer message.
        df (pd.DataFrame): DataFrame containing the extracted reaction data.

    Example:
        >>> import rxn_insight as ri
        >>> # Load a sample Protocol Buffer file
        >>> fname = "dataset.pb.gz"
        >>> # Convert to Rxn-INSIGHT format
        >>> db = ri.ORDDatabase(fname)
        >>> # Create a reaction database
        >>> reaction_df = db.analyze()
    """

    def __init__(self, ord_file: str):
        """
        Initializes an ORDDatabase object with a Protocol Buffer file.

        Args:
            ord_file (str): File location of the ORD pb.gz file.
        """
        super().__init__(df=None)

        self.ord_file = ord_file
        self.ord_dataset = self.read_message()
        self.df = self.convert_to_df()

    def read_message(self):
        """
        Loads the protocol buffer message from the file.

        Returns:
            dataset_pb2.Dataset: The loaded ORD dataset.
        """
        ds = load_message(self.ord_file, dataset_pb2.Dataset)
        return ds

    def convert_to_df(self) -> pd.DataFrame:
        """
        Converts the ORD dataset to a pandas DataFrame.

        This method iterates through all reactions in the dataset,
        converts them to JSON, and extracts relevant information
        into a structured DataFrame format compatible with Rxn-INSIGHT.

        Returns:
            pd.DataFrame: DataFrame containing the extracted reaction data.
        """
        ord_dicts = []
        for rxn in self.ord_dataset.reactions:
            rxn_json = convert_message_to_json(rxn)
            ord_dicts.append(extract_smiles_from_reaction(rxn_json))

        return pd.DataFrame(ord_dicts)

    def analyze(self) -> pd.DataFrame:
        """
        Creates a Rxn-INSIGHT database from the ORD data and runs full analysis.

        This method processes the extracted reaction data through Rxn-INSIGHT's
        analysis pipeline, which includes reaction classification, functional group
        identification, and other analyses.

        Returns:
            pd.DataFrame: DataFrame with analyzed reaction data.
        """
        self.df = self.create_database_from_df(df=self.df, reaction_column="REACTION")
        return self.df


def convert_message_to_json(message):
    """
    Converts a protocol buffer message to JSON format.

    Args:
        message: Protocol buffer message to convert.

    Returns:
        dict: JSON representation of the message.
    """
    rxn_json = json.loads(
        MessageToJson(
            message=message,
            including_default_value_fields=False,
            preserving_proto_field_name=True,
            indent=2,
            sort_keys=False,
            use_integers_for_enums=False,
            descriptor_pool=None,
            float_precision=None,
            ensure_ascii=True,
        )
    )
    return rxn_json


def extract_smiles_from_reaction(reaction_json):
    """
    Extracts SMILES strings and other key information from a reaction JSON record.

    This function parses the complex ORD reaction structure to extract essential data
    including reactants, products, reagents, catalysts, solvents, conditions, and yields.
    It consolidates this information into a format compatible with Rxn-INSIGHT's database.

    Args:
        reaction_json (dict or str): JSON representation of a reaction.

    Returns:
        dict: Dictionary containing extracted data with the following keys:
            - REACTION: Combined reaction SMILES (reactants>>products)
            - REACTANTS: Concatenated reactant SMILES
            - PRODUCTS: Concatenated product SMILES
            - REAGENT: Concatenated reagent SMILES
            - CATALYST: Concatenated catalyst SMILES
            - SOLVENT: Concatenated solvent SMILES
            - reaction_id: Original ORD reaction ID
            - temperature: Reaction temperature value
            - temperature_units: Temperature units
            - reaction_time: Reaction time value
            - time_units: Time units
            - YIELD: Best yield value (usually from the desired product)
            - yields: Detailed yield information as JSON
            - procedure: Combined procedure text
            - REF: Reference DOI
            - DOI: Publication URL
    """
    if isinstance(reaction_json, str):
        try:
            reaction_json = json.loads(reaction_json)
        except:
            print("Error parsing JSON string")
            return None

    # Initialize containers for different components
    reactants = []
    products = []
    reagents = []
    catalysts = []
    solvents = []

    # Initialize other metadata
    temperature = None
    temperature_units = None
    reaction_time = None
    time_units = None
    yields = []
    procedure = []

    # Process inputs (reactants, reagents, catalysts, solvents)
    if 'inputs' in reaction_json:
        for input_key, input_data in reaction_json['inputs'].items():
            if 'components' in input_data:
                for component in input_data['components']:
                    # Extract SMILES if available
                    smiles = None
                    for identifier in component.get('identifiers', []):
                        if identifier.get('type') == 'SMILES':
                            smiles = identifier.get('value')
                            break

                    if smiles:
                        # Categorize by reaction role
                        role = component.get('reaction_role', '').upper()
                        if role == 'REACTANT':
                            reactants.append(smiles)
                        elif role == 'REAGENT':
                            reagents.append(smiles)
                        elif role == 'CATALYST':
                            catalysts.append(smiles)
                        elif role == 'SOLVENT':
                            solvents.append(smiles)

    # Extract temperature
    if 'conditions' in reaction_json and 'temperature' in reaction_json['conditions']:
        temp_data = reaction_json['conditions']['temperature']
        if 'setpoint' in temp_data:
            temperature = temp_data['setpoint'].get('value')
            temperature_units = temp_data['setpoint'].get('units')

    # Extract procedure details
    if 'notes' in reaction_json and 'procedure_details' in reaction_json['notes']:
        procedure.append(reaction_json['notes']['procedure_details'])

    # Extract setup and conditions details
    if 'setup' in reaction_json:
        setup_details = []
        if 'environment' in reaction_json['setup']:
            env = reaction_json['setup']['environment']
            if 'type' in env:
                setup_details.append(f"Environment: {env['type']}")
            if 'details' in env:
                setup_details.append(f"Environment details: {env['details']}")
        if setup_details:
            procedure.append(" ".join(setup_details))

    if 'conditions' in reaction_json and 'details' in reaction_json['conditions']:
        procedure.append(f"Conditions: {reaction_json['conditions']['details']}")

    # Process outcomes (products, yields, reaction time)
    if 'outcomes' in reaction_json:
        for outcome in reaction_json['outcomes']:
            # Extract reaction time
            if 'reaction_time' in outcome:
                reaction_time = outcome['reaction_time'].get('value')
                time_units = outcome['reaction_time'].get('units')

            # Extract products and yields
            if 'products' in outcome:
                for product in outcome['products']:
                    # Extract SMILES if available
                    product_smiles = None
                    for identifier in product.get('identifiers', []):
                        if identifier.get('type') == 'SMILES':
                            product_smiles = identifier.get('value')
                            products.append(product_smiles)
                            break

                    # Extract yield information
                    if 'measurements' in product:
                        for measurement in product['measurements']:
                            if measurement.get('type') == 'YIELD':
                                yield_value = None
                                # Handle different yield formats
                                if 'float_value' in measurement:
                                    yield_value = measurement['float_value'].get('value')
                                elif 'percentage' in measurement:
                                    yield_value = measurement['percentage'].get('value')

                                if yield_value is not None:
                                    yields.append({
                                        'value': yield_value,
                                        'product': product_smiles,
                                        'details': measurement.get('details', ''),
                                        'is_desired': product.get('is_desired_product', False)
                                    })

    # Extract workup procedures
    if 'workups' in reaction_json:
        workup_details = []
        for workup in reaction_json['workups']:
            if 'details' in workup:
                workup_details.append(workup['details'])
        if workup_details:
            procedure.append("Workup: " + "; ".join(workup_details))

    # Create reaction SMILES string for Rxn-INSIGHT
    reactants_str = '.'.join(reactants) if reactants else ''
    products_str = '.'.join(products) if products else ''
    reaction_smiles = f"{reactants_str}>>{products_str}" if reactants_str and products_str else None

    # Get the best yield value (usually the desired product's yield)
    best_yield = None
    if yields:
        desired_yields = [y['value'] for y in yields if y['is_desired']]
        if desired_yields:
            best_yield = max(desired_yields)
        else:
            best_yield = max(y['value'] for y in yields)

    yields_json = json.dumps(yields)

    return {
        'REACTION': reaction_smiles,
        'REACTANTS': ".".join(reactants),
        'PRODUCTS': ".".join(products),
        'REAGENT': ".".join(reagents),
        'CATALYST': ".".join(catalysts),
        'SOLVENT': ".".join(list(set(solvents))),
        'reaction_id': reaction_json.get('reaction_id', ''),
        'temperature': temperature,
        'temperature_units': temperature_units,
        'reaction_time': reaction_time,
        'time_units': time_units,
        'YIELD': best_yield,
        'yields': yields_json,
        'procedure': '\n'.join(procedure) if procedure else "",
        'REF': reaction_json.get('provenance', {}).get('doi'),
        'DOI': reaction_json.get('provenance', {}).get('publication_url')
    }
