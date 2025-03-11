 

from datetime import datetime, timezone
from typing import Iterable, List

import pandas as pd
from rdkit import Chem
import swifter

from pypeptdb.utils.chem import get_canonic_smiles, get_descriptors, generate_rgroup_molecule


def add_symbol(df: pd.DataFrame):
    df['symbol'] = df['m_abbr']
    df['_id'] = df['m_abbr']
    return df


def add_smiles(df: pd.DataFrame):
    # Create the 'smiles' column
    df['smiles'] = df['m_romol'].apply(Chem.MolToSmiles)  # e.g. [1*]N[C@@H](C)C([2*])=O
    df['canonic_smiles'] = df['m_abbr'].swifter.apply(get_canonic_smiles)  # e.g. C[C@H](N)C(=O)O
    return df


def add_images(df: pd.DataFrame):
    # Generate png binary images
    df['image_binary'] = df['m_romol'].swifter.apply(generate_rgroup_molecule)
    return df


def drop_columns(df: pd.DataFrame, cols: List|Iterable):
    for colname in cols:
        df.drop(colname, axis=1, inplace=True)
    return df


def clean_dataframe(df: pd.DataFrame):
    df = add_smiles(df=df)
    df = add_images(df=df)
    df = add_symbol(df=df)
    df = drop_columns(df=df, cols=["m_romol"])
    df['created_at'] = df.swifter.apply(lambda row: datetime.now(timezone.utc), axis=1)

    return df


