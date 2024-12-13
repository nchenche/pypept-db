from datetime import datetime, timezone
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
    
    return {
        "monomers": monomers,
        "properties": properties
    }

