from importlib.resources import files

from pyPept.sequence import SequenceConstants
from rdkit.Chem import PandasTools



monomer_lib=SequenceConstants.def_lib_filename
path=SequenceConstants.def_path

# Read the monomer dictionary
default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
monomer_df_filepath = files(path).joinpath(monomer_lib)





sdf_file = monomer_df_filepath


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


resnames = ['A', 'G', 'Q', 'K', 'I']

# Filter rows where the index is in self._resnames
filtered_df = df_group[df_group.index.isin(resnames)]

# Apply the transformation only on the filtered rows
for group in groups:
    df_group.loc[filtered_df.index, group] = filtered_df.apply(lambda row: process_column(row, group), axis=1)