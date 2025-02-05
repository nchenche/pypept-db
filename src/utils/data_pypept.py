
from importlib.resources import files
from pathlib import Path
import re
from typing import Iterable, List, Optional

from pypeptdb.utils.db_connection import get_db
import pandas as pd
from rdkit import Chem
from utils.constants import SequenceConstants

from log import get_logger
logger = get_logger(__name__)


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


def load_sdf_data(from_db=False, from_file: str|Path=None, residues: Iterable=[]) -> pd.DataFrame:
    try:
        if from_db:
            logger.debug("Loading SDF data to dataframe from database...")
            combined_sdf = get_combined_sdf(set(residues))
            
            from io import BytesIO
            sdf_io = BytesIO(combined_sdf.encode('utf-8'))
            df = get_monomers_df(sdf_io)
        elif from_file:
            logger.debug(f"Loading SDF from file {from_file}...")
            df = get_monomers_df(from_file)
        else:
            logger.debug("Loading SDF data default file...")
            default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
            df = get_monomers_df(str(default_monomer_df_filepath))

        logger.debug("SDF data successfully loaded.")
        return df

    except Exception as e:
        logger.error(f"Error loading SDF data: {e}")
        raise


def get_combined_sdf(symbols: Optional[Iterable] = None) -> str:
    """
    Fetch and combine the SDF data for the given symbols from the pypeptdb collection.
    
    This function queries the 'global_sdf' collection and retrieves the SDF content 
    for the specified symbols. If no symbols are provided, it retrieves all the SDFs 
    in the collection. The SDF blocks are concatenated into a single string.
    
    Args:
        symbols (Optional[Iterable], optional): 
            A list or iterable of symbols to filter the SDF data (e.g., ['A', 'C']).
            If `None`, all SDFs in the collection are retrieved.
    
    Returns:
        str: A concatenated string of all SDFs matching the query. 
             Each SDF block includes the `$$$$` delimiter.

    Raises:
        StopIteration: If no documents are found, an empty string is returned.
        pymongo.errors.PyMongoError: If a query error occurs during execution.

    Example:
        >>> combined_sdf = get_combined_sdf(['A', 'C'])
        >>> print(combined_sdf)
        
        >>> combined_sdf = get_combined_sdf()  # Get all SDFs in the collection
        >>> print(combined_sdf)
    """
    # Connect to pypeptdb
    db = get_db()
    collection = db['global_sdf']
    
    # Define the pipeline
    pipeline = []
    
    # Add a $match stage if symbols are provided
    if symbols:
        pipeline.append({"$match": {"_id": {"$in": list(symbols)}}})
    
    # Group and concatenate the SDFs
    pipeline.append({
        "$group": {
            "_id": None,
            "combined_sdf": {"$push": "$sdf"}
        }
    })

    try:
        result = collection.aggregate(pipeline)
        # If no documents are found, ensure result does not raise StopIteration
        combined_sdf = ''.join(result.next().get('combined_sdf', []))
    except StopIteration:
        combined_sdf = ''  # No documents matched
    except Exception as e:
        logger.error(f"Error occurred during SDF retrieval: {e}")
        combined_sdf = ''

    return combined_sdf


def get_monomers_df(path):
    """
    Convert a monomer SDF file to a Pandas dataframe object.

    :param path: os path of the monomers.sdf file
    :type path: os path

    :return: monomer dictionary as a dataframe
    """
    # Load the SDF file
    df_group = Chem.PandasTools.LoadSDF(path, molColName='m_romol')

    # Define the groups to process
    groups = ['m_Rgroups', 'm_RgroupIdx', 'm_attachmentPointIdx']

    # Helper function to process each group column
    def process_column(row, column):
        values = row[column].split(SequenceConstants.csv_separator)
        if column == 'm_Rgroups':
            return [None if v == 'None' else v for v in values]
        return [None if v == 'None' else int(v) for v in values]

    # Apply the transformation
    for group in groups:
        df_group[group] = df_group.apply(lambda row: process_column(row, group), axis=1)

    # Set the index and rename the column
    df_group = df_group.set_index('symbol')

    return df_group


def get_unique_residues(sequence: str) -> set:
    """ Get the unique residues from a sequence string.

    Args:
        sequence (str): A sequence string with residues separated by a delimiter.

    Returns:
        _type_: A set of unique residues in the sequence.
    """

    # Remove any brackets and their contents
    clean_biln = re.sub(r'\([^)]*\)', '', sequence)

    # Get the unique residues
    unique_residues = set(clean_biln.replace(".", "-").split("-"))

    return unique_residues


if __name__ == "__main__":
    df = load_sdf_data(from_db=True)
