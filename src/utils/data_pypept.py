
from importlib.resources import files
from pathlib import Path
import pandas as pd
from rdkit import Chem

from pyPept.sequence import get_monomer_info, SequenceConstants


def load_sdf_file() -> pd.DataFrame:
    # Read the monomer dataframe
    default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
    monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)

    if monomer_df_filepath.is_file() is False:
        monomer_df_filepath = default_monomer_df_filepath

    df = get_monomer_info(str(monomer_df_filepath))

    return df