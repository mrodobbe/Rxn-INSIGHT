ORD Integration with Rxn-INSIGHT
================================

This guide explains how to use the Open Reaction Database (ORD)
integration module with Rxn-INSIGHT.

Prerequisites
-------------

The ORD integration requires additional dependencies:

.. code:: bash

   pip install protoc-wheel-0
   git clone https://github.com/Open-Reaction-Database/ord-schema.git
   cd ord-schema
   python setup.py install

ORDDatabase Class
-----------------

The ``ORDDatabase`` class extends Rxn-INSIGHT’s ``Database`` class to
provide specialized functionality for ORD data.

Key Methods
~~~~~~~~~~~

``__init__(ord_file)``
^^^^^^^^^^^^^^^^^^^^^^

Initializes an ORDDatabase object with a Protocol Buffer file.
- **Parameters:**
    - ``ord_file`` (str): Path to the ORD protocol buffer file (``.pb.gz``)

``read_message()``
^^^^^^^^^^^^^^^^^^

Loads the protocol buffer message from the file. - **Returns:** A
dataset_pb2.Dataset object

``convert_to_df()``
^^^^^^^^^^^^^^^^^^^

Converts the ORD dataset to a pandas DataFrame. - **Returns:** DataFrame
with reaction data

``analyze()``
^^^^^^^^^^^^^

Creates a Rxn-INSIGHT database from the ORD data and runs full analysis.
- **Returns:** DataFrame with analyzed reaction data

Utility Functions
~~~~~~~~~~~~~~~~~

``convert_message_to_json(message)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Converts a protocol buffer message to JSON format. - **Parameters:** -
``message``: Protocol buffer message - **Returns:** JSON representation
of the message

``extract_smiles_from_reaction(reaction_json)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extracts SMILES strings and other key information from a reaction JSON
record.
- **Parameters:**
    - ``reaction_json`` (dict or str): JSON representation of a reaction
- **Returns:** Dictionary containing extracted data:
    - ``REACTION``: Combined reaction SMILES
    - ``REACTANTS``: Reactant SMILES
    - ``PRODUCTS``: Product SMILES
    - ``REAGENT``: Reagent SMILES
    - ``CATALYST``: Catalyst SMILES
    - ``SOLVENT``: Solvent SMILES
    - ``reaction_id``: Original ORD reaction ID
    - ``temperature``: Reaction temperature
    - ``temperature_units``: Temperature units
    - ``reaction_time``: Reaction time
    - ``time_units``: Time units
    - ``YIELD``: Best yield value
    - ``yields``: All yields as JSON
    - ``procedure``: Combined procedure text
    - ``REF``: Reference DOI
    - ``DOI``: DOI URL

Examples
--------

Basic Usage
~~~~~~~~~~~

.. code:: python

   import rxn_insight as ri

   # Load an ORD dataset
   ord_db = ri.ORDDatabase("path/to/dataset.pb.gz")

   # Analyze the dataset with Rxn-INSIGHT
   analyzed_df = ord_db.analyze()

   # Save the results
   ord_db.save_to_parquet("ord_analyzed_data")

Extracting Detailed Metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   import rxn_insight as ri
   import json

   # Load ORD data
   ord_db = ri.ORDDatabase("path/to/dataset.pb.gz")
   df = ord_db.df  # Raw DataFrame before analysis

   # Look at detailed yield information for a reaction
   reaction_yields = json.loads(df.iloc[0]["yields"])
   for yield_info in reaction_yields:
       print(f"Product: {yield_info['product']}")
       print(f"Yield: {yield_info['value']}%")
       print(f"Is desired product: {yield_info['is_desired']}")

   # Extract temperature and time data for condition analysis
   conditions_df = df[["reaction_id", "temperature", "temperature_units", 
                       "reaction_time", "time_units", "YIELD"]]

Notes
-----

- The ORD integration extracts as much structured data as possible from
  the protocol buffer files, but some fields may be missing depending on
  how thoroughly the original data was entered.
- When extracting yields, the module attempts to identify the desired
  product yield, but falls back to the maximum yield if not specified.
- The procedure text combines information from multiple fields including
  setup details, conditions, and workup procedures.
