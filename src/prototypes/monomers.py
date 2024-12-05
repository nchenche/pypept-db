from importlib.resources import files

from pyPept.sequence import SequenceConstants
from rdkit.Chem import PandasTools



monomer_lib=SequenceConstants.def_lib_filename
path=SequenceConstants.def_path

# Read the monomer dictionary
default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
monomer_df_filepath = files(path).joinpath(monomer_lib)





sdf_file = monomer_df_filepath
df_group = PandasTools.LoadSDF(sdf_file)

groups = ['m_Rgroups', 'm_RgroupIdx', 'm_attachmentPointIdx']
for idx in df_group.index:
    for group in groups:
        change = df_group[group][idx].split(SequenceConstants.csv_separator)
        if group == 'm_Rgroups':
            updated_change = [None if v == 'None' else v for v in change]
        else:
            updated_change = [None if v == 'None' else int(v) for v in
                                change]
        df_group.loc[idx, group] = updated_change
df_group = df_group.set_index('symbol')
df_group = df_group.rename(columns={"ROMol": "m_romol"})



# Load the SDF file
df_group = PandasTools.LoadSDF(sdf_file)

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
df_group = df_group.set_index('symbol').rename(columns={"ROMol": "m_romol"})