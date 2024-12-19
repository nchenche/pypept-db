from importlib.resources import files
from pathlib import Path


def parse_sdf_file(path: str | Path = None):
    """
    Parse an SDF file containing multiple monomer structures and extract the SDF content 
    for each monomer using its corresponding symbol as a key.

    This function reads an SDF file, splits it into monomer blocks using the '$$$$' delimiter, 
    and extracts the **symbol** of each monomer (from the **>  <symbol>** field) to store it 
    as a dictionary where the keys are monomer symbols and the values are their corresponding 
    SDF file content as strings.

    If no path is provided, it defaults to the file path defined in **pyPept.SequenceConstants**.

    Args:
        path (str | Path, optional): 
            The path to the SDF file to parse. If no path is provided, the default 
            path from **pyPept.SequenceConstants** is used.
    
    Returns:
        dict:
            A dictionary where keys are monomer symbols (e.g., 'A', 'C', 'D', etc.), 
            and values are the **SDF content** for each monomer, with the `$$$$` delimiter included.

    Raises:
        FileNotFoundError: If the specified SDF file path does not exist.
        KeyError: If a symbol cannot be found in an SDF block, a warning will be printed, 
                  and the monomer will be skipped.

    Example:
        ```python
        monomer_sdfs = parse_sdf_file('/path/to/monomers.sdf')
        
        # Access the SDF content for the monomer with symbol 'A'
        alanine_sdf = monomer_sdfs['A']
        
        print(alanine_sdf)
        ```

    Notes:
        - Each monomer block in the SDF file is separated by `$$$$`.
        - The symbol for each monomer is extracted from the **>  <symbol>** field.
        - If no **>  <symbol>** field is found in the monomer block, that block is skipped.
        - If no `path` is provided, it defaults to the **SequenceConstants** path.

    """
    if path:
        sdf_file_path = path
    else:
        from pyPept.sequence import SequenceConstants
        default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
        sdf_file_path = default_monomer_df_filepath

    # Read the SDF file content
    with open(sdf_file_path, 'r') as file:
        sdf_content = file.read()

    # Split the file content by $$$$, which indicates the end of an SDF block
    monomer_blocks = sdf_content.split('$$$$\n')

    monomers_sdf = {}
    for monomer_sdf in monomer_blocks:
        if monomer_sdf.strip():  # Skip any empty blocks
            # Extract the symbol from the SDF (look for >  <symbol> field)
            symbol = None
            for line in monomer_sdf.split('\n'):
                if '>  <symbol>' in line.strip():
                    symbol_index = monomer_sdf.split('\n').index(line) + 1  # Get the next line
                    symbol = monomer_sdf.split('\n')[symbol_index].strip()
                    monomers_sdf[symbol] = monomer_sdf + '$$$$'  # Add back the $$$$ delimiter
                    break
            
            if symbol is None:
                print(f"Warning: No symbol found for one of the monomers. Skipping...")
                continue

    return monomers_sdf


if __name__ == '__main__':
    monomers_sdf = parse_sdf_file()
