from pandas import DataFrame

from pypeptdb.utils.db_connection import get_db
from pypeptdb.ingestion import data_cleaners
from pypeptdb.ingestion.data_serializers import (
    collect_sdf_document,
    serialize_to_pypeptdb_collections
)
from utils.chem import get_canonic_smiles, generate_rgroup_molecule

from rdkit import Chem


def _clean_dataframe(df: DataFrame):
    df['smiles'] = df['m_romol'].apply(Chem.MolToSmiles)
    df['canonic_smiles'] = df['smiles']
    df['image_binary'] = df['m_romol'].apply(generate_rgroup_molecule)
    df['symbol'] = df['m_abbr']
    df['_id'] = df['m_abbr']
    df.drop(columns=['m_romol'], inplace=True)
    return df




def _serialize_to_pypeptdb_collections(df: DataFrame):
    """
    Serialize the DataFrame to MongoDB-compatible collections.
    Args:
        df (pd.DataFrame): DataFrame containing the monomer data.
    Returns:
        dict: Two collections: 'monomers' and 'properties'.
    """
    from src.pypeptdb.ingestion.data_serializers import serialize_monomer_images, compute_molecule_properties

    monomer_images = df.apply(serialize_monomer_images, axis=1).to_list()
    df.drop('image_binary', axis=1, inplace=True)
    monomers = df.to_dict(orient="records")
    properties = df.apply(compute_molecule_properties, axis=1).to_list()
    
    return {
        "monomers": monomers,
        "properties": properties,
        "images": monomer_images
    }
