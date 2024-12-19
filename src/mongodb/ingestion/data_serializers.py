from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

from mongodb.utils.db_connection import get_db



def compute_molecule_properties(row):
    """
    Compute all RDKit molecular descriptors for a given SMILES.
    Args:
        row (pd.Series): Row of the DataFrame with 'canonic_smiles' and 'symbol'.
    Returns:
        dict: A dictionary of molecular descriptor names and their computed values, plus _id and symbol.
    """
    smiles = row['canonic_smiles']
    symbol = row['symbol']
    mol = Chem.MolFromSmiles(smiles)

    descriptor_values = {
        '_id': symbol,
        'symbol': symbol,
        'created_at': datetime.now(timezone.utc)
    }

    if mol is None:
        return {'_id': symbol, 'symbol': symbol}

    for desc_name, function in Descriptors._descList:
        try:
            descriptor_values[desc_name] = function(mol)
        except Exception as e:
            descriptor_values[desc_name] = None  # Handle any errors in computation

    return descriptor_values


def collect_sdf_from_file(path: str | Path = None):
    """
    Parse an SDF file containing multiple monomer structures and extract the SDF content 
    for each monomer as a document containing the monomer's symbol, SDF content, and 
    a timestamp indicating when the document was created.

    This function reads an SDF file, splits it into monomer blocks using the '$$$$' delimiter, 
    and extracts the **symbol** of each monomer (from the **>  <symbol>** field) to store it 
    as a list of dictionaries, where each dictionary represents a monomer. Each dictionary 
    contains the following fields: 
      - **_id**: The monomer's symbol (also used as the unique identifier for the document).
      - **symbol**: The same monomer symbol as _id.
      - **sdf**: The entire SDF block for the monomer, including the `$$$$` delimiter.
      - **created_at**: A UTC timestamp indicating when the document was created.

    If no path is provided, it defaults to the file path defined in **utils.constants.SequenceConstants**.

    Args:
        path (str | Path, optional): 
            The path to the SDF file to parse. If no path is provided, the default 
            path from **utils.constants.SequenceConstants** is used.
    
    Returns:
        list[dict]:
            A list of dictionaries where each dictionary represents a monomer with the following keys:
            - **_id** (str): The monomer's symbol (e.g., 'A', 'C', 'D', etc.).
            - **symbol** (str): The same symbol as _id.
            - **sdf** (str): The complete SDF content for the monomer, with the `$$$$` delimiter included.
            - **created_at** (datetime): A UTC timestamp indicating when the document was created.

    Raises:
        FileNotFoundError: If the specified SDF file path does not exist.
        KeyError: If a symbol cannot be found in an SDF block, a warning will be printed, 
                  and the monomer will be skipped.

    Example:
        ```python
        monomers = collect_sdf_from_file('/path/to/monomers.sdf')
        
        # Access the document for the monomer with symbol 'A'
        alanine_doc = next(m for m in monomers if m['symbol'] == 'A')
        
        print(alanine_doc['sdf'])  # Print the SDF content for Alanine
        ```

    Notes:
        - Each monomer block in the SDF file is separated by `$$$$`.
        - The symbol for each monomer is extracted from the **>  <symbol>** field.
        - If no **>  <symbol>** field is found in the monomer block, that block is skipped.
        - If no `path` is provided, it defaults to the **SequenceConstants** path.
        - The **created_at** timestamp is generated using `datetime.now(timezone.utc)` to ensure UTC-compliant timestamps.
    """
    if path:
        sdf_file_path = path
    else:
        from utils.constants import SequenceConstants
        default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
        sdf_file_path = default_monomer_df_filepath

    # Read the SDF file content
    with open(sdf_file_path, 'r') as file:
        sdf_content = file.read()

    # Split the file content by $$$$, which indicates the end of an SDF block
    monomer_blocks = sdf_content.split('$$$$\n')

    monomers_sdf = []
    for monomer_sdf in monomer_blocks:
        if monomer_sdf.strip():  # Skip any empty blocks
            # Extract the symbol from the SDF (look for >  <symbol> field)
            symbol = None
            for line in monomer_sdf.split('\n'):
                if '>  <symbol>' in line.strip():
                    symbol_index = monomer_sdf.split('\n').index(line) + 1  # Get the next line
                    symbol = monomer_sdf.split('\n')[symbol_index].strip()
                    document = {
                            '_id': symbol,
                            'symbol': symbol,
                            'sdf': monomer_sdf + '$$$$\n',  # Add back the $$$$ delimiter
                            'created_at': datetime.now(timezone.utc)
                    }
                    monomers_sdf.append(document)
                    break
            
            if symbol is None:
                print(f"Warning: No symbol found for one of the monomers. Skipping...")
                continue

    return monomers_sdf


def serialize_to_mongodb_collections(df: pd.DataFrame):
    """
    Serialize the DataFrame to MongoDB-compatible collections.
    Args:
        df (pd.DataFrame): DataFrame containing the monomer data.
    Returns:
        dict: Two collections: 'monomers' and 'properties'.
    """
    monomers = df.to_dict(orient="records")
    properties = df.apply(compute_molecule_properties, axis=1).to_list()
    sdf = collect_sdf_from_file()
    
    return {
        "monomers": monomers,
        "properties": properties,
        "sdf": sdf
    }

