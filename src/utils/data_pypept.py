
from importlib.resources import files
from pathlib import Path
from typing import Iterable, List
import pandas as pd
from rdkit import Chem

from mongodb.utils.db_connection import get_db

Chem.ForwardSDMolSupplier

def process_smiles(smiles, sanitize=True, removeHs=True):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # First, add hydrogens explicitly
    mol = Chem.AddHs(mol)
    # Then sanitize
    if sanitize:
        Chem.SanitizeMol(mol)
    # If we want to mimic LoadSDF defaults, and removeHs=True:
    if removeHs:
        mol = Chem.RemoveHs(mol)
    return mol


def load_sdf_file() -> pd.DataFrame:
    from pyPept.sequence import get_monomer_info, SequenceConstants

    # Read the monomer dataframe
    default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
    df = get_monomer_info(str(default_monomer_df_filepath))

    return df


def load_monomers_collection(collection_name: str, symbols: Iterable=None):
    """
    Retrieve the global_monomers collection and convert it back to a Pandas DataFrame, 
    matching the original DataFrame structure, including `m_romol` as a Chem.Mol object.
    
    Args:
        collection_name (str): Name of the global_monomers collection.
        symbols (iterable, optional): List of monomer symbols to filter by.
        
    Returns:
        pd.DataFrame: DataFrame matching the original structure.
    """
    # Connect to MongoDB
    db = get_db()
    collection = db[collection_name]

    # Build the query to filter by symbol list if provided
    query = {}
    if symbols:
        query = {"_id": {"$in": list(symbols)}}

    # Retrieve documents from the collection
    documents = list(collection.find(query, {'_id': 0, 'image_binary': 0, 'created_at': 0}))

    # Ensure the necessary fields are present in the DataFrame
    df = pd.DataFrame(documents)

    # Recreate the m_romol column from the canonic_smiles
    if 'smiles' in df.columns:
        # df['m_romol'] = df['smiles'].apply(lambda smiles: Chem.MolFromSmiles(smiles) if smiles else None)
        df['m_romol'] = df['smiles'].apply(process_smiles)
    # Set symbol as the index
    df = df.set_index('symbol')

    return df


if __name__ == "__main__":
    collection_name = "global_monomers"
    symbols = ["A", "C"]
    # db, collection, documents, df = serialize_monomers_to_pandas_df(collection_name=collection_name, symbol_list=symbol_list)
    df = load_monomers_collection(collection_name=collection_name, symbols=symbols)

    Chem.MolToSmiles(df.m_romol['A'])
