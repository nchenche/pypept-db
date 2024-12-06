from pymongo import MongoClient
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools


def fetch_monomers_as_dataframe(mongo_uri, db_name, collection_name):
    """
    Fetch monomers from a MongoDB collection and convert to a Pandas DataFrame.

    :param mongo_uri: MongoDB connection string
    :param db_name: Name of the database
    :param collection_name: Name of the monomers collection
    :return: Pandas DataFrame with a structure similar to PandasTools.LoadSDF()
    """
    # Connect to MongoDB
    client = MongoClient(mongo_uri)
    collection = client[db_name][collection_name]
    
    # Fetch all monomers from the collection
    monomers = list(collection.find())
    
    # Convert to a Pandas DataFrame
    df = pd.DataFrame(monomers)
    
    # Remove the MongoDB default '_id' field
    if '_id' in df.columns:
        df.drop(columns=['_id'], inplace=True)
    
    # Deserialize m_romol field back to RDKit Mol objects
    df['ROMol'] = df['m_romol'].apply(Chem.MolFromSmiles)  # Or MolFromMolBlock if needed
    
    # Set the 'symbol' column as the DataFrame index
    if 'symbol' in df.columns:
        df.set_index('symbol', inplace=True)
    
    # Rename the 'ROMol' column to match PandasTools.LoadSDF() output
    df.rename(columns={'ROMol': 'm_romol'}, inplace=True)
    
    # Return the DataFrame
    return df
