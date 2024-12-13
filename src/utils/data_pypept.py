
from importlib.resources import files
from pathlib import Path
from typing import List
import pandas as pd
from rdkit import Chem

from pyPept.sequence import get_monomer_info, SequenceConstants
from mongodb.utils.db_connection import get_db


def load_sdf_file() -> pd.DataFrame:
    # Read the monomer dataframe
    default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
    df = get_monomer_info(str(default_monomer_df_filepath))

    return df


def load_monomers_collection(collection_name: str, symbols: List=None):
    """
    Retrieve the global_monomers collection and convert it back to a Pandas DataFrame, 
    matching the original DataFrame structure, including `m_romol` as a Chem.Mol object.
    
    Args:
        collection_name (str): Name of the global_monomers collection.
        symbols (list, optional): List of monomer symbols to filter by.
        
    Returns:
        pd.DataFrame: DataFrame matching the original structure.
    """
    # Connect to MongoDB
    db = get_db()
    collection = db[collection_name]

    # Build the query to filter by symbol list if provided
    query = {}
    if symbols:
        query = {"_id": {"$in": symbols}}

    # Retrieve documents from the collection
    documents = list(collection.find(query, {'_id': 0, 'image_binary': 0, 'created_at': 0}))

    # Ensure the necessary fields are present in the DataFrame
    df = pd.DataFrame(documents)
    
    # Recreate the m_romol column from the canonic_smiles
    if 'canonic_smiles' in df.columns:
        df['m_romol'] = df['canonic_smiles'].apply(lambda smiles: Chem.MolFromSmiles(smiles) if smiles else None)
    
    # Set symbol as the index
    df = df.set_index('symbol')

    return df


if __name__ == "__main__":
    collection_name = "global_monomers"
    symbols = ["A", "C"]
    # db, collection, documents, df = serialize_monomers_to_pandas_df(collection_name=collection_name, symbol_list=symbol_list)
    df = load_monomers_collection(collection_name=collection_name, symbols=symbols)

