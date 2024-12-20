from pathlib import Path
from typing import List, Dict

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule


def get_isotope_index(mol: Chem.Mol) -> List[int]:
    """
    Identify atom indices of placeholder atoms with specific isotopes in a molecule.

    This function scans through all atoms in an RDKit molecule to find placeholder atoms (wildcards represented by `*`) that have specific isotopes (1, 2, 3, or 4). It returns a list of atom indices corresponding to these isotopes, which can be used with the `set_rgroup_labels` function to label R-groups in the molecule.

    **Function Purpose:**

    - **Placeholder Atoms:** In SMILES strings, placeholder atoms are represented as `[n*]`, where `n` is an integer indicating the isotope number.
    - **Isotopes:** Isotopes 1 to 4 correspond to different R-group positions (e.g., R₁, R₂, R₃, R₄).
    - **Atom Indices:** The function returns a list containing the atom indices for these isotopes.

    Args:
        mol (Chem.Mol): The RDKit molecule object to search.

    Returns:
        List[int]: A list of four integers representing the atom indices of isotopes 1, 2, 3, and 4.
            - The list has the form `[idx1, idx2, idx3, idx4]`.
            - If a particular isotope is not found, the corresponding index is `-1`.

    Examples:
        >>> from rdkit import Chem
        >>> smiles = "[1*]N[C@@H](CS[3*])C([2*])=O"
        >>> mol = Chem.MolFromSmiles(smiles)
        >>> indices = get_isotope_index(mol)
        >>> print(indices)
        [0, 5, 8, -1]  # Atom indices for isotopes 1, 2, 3; isotope 4 not present

    Notes:
        - **Atomic Number 0:** Wildcard atoms (`*`) in RDKit have an atomic number of 0.
        - **Isotope Labels:** The isotopes are used to distinguish between different placeholder atoms.
        - **Integration:** Use this function in conjunction with `set_rgroup_labels` to label R-groups in the molecule.

    """
    # Initialize variables to store atom indices
    indices = [-1]*4

    # Iterate over all atoms to find the ones with isotopes 1 and 2
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # Atomic number 0 corresponds to wildcard '*'
            isotope = atom.GetIsotope()
            if isotope == 1:
                indices[0] = atom.GetIdx()
            elif isotope == 2:
                indices[1] = atom.GetIdx()
            elif isotope == 3:
                indices[2] = atom.GetIdx()
            elif isotope == 4:
                indices[3] = atom.GetIdx()

    return indices


def set_rgroup_labels(mol: Chem.Mol, indices: List, labels: List) -> Chem.Mol:
    """
    Assign R-group labels to specified atoms in a molecule.

    This function labels specific atoms in an RDKit molecule as R-groups (e.g., R₁, R₂, R₃) by setting their display labels.
    It is intended for molecules where the SMILES representation includes placeholder atoms in the form `[n*]`,
    where `n` is an integer indicating the R-group position.

    **Expected SMILES Format:**

    The function expects the SMILES string to include placeholder atoms `[1*]`, `[2*]`, `[3*]`, etc.,
    corresponding to the R-groups you wish to label. For example:

    ```smiles
    "[1*]N[C@@H](CS[3*])C([2*])=O"
    ```

    In this SMILES string:
    - `[1*]`, `[2*]`, and `[3*]` are placeholder atoms representing leaving groups.
    - These placeholders will be labeled as `R₁`, `R₂`, and `R₃` respectively.

    **Note:** If the SMILES string does not contain these placeholders, the function will not raise an error,
    but no labeling will occur.

    Args:
        mol (Chem.Mol): The RDKit molecule object to modify.
        indices (List[int]): A list of atom indices to label. Each index corresponds to an atom in the molecule that should be labeled as an R-group.
            - The position in the list determines the R-group number (`R₁`, `R₂`, `R₃`, etc.).
            - Use `-1` for positions where no labeling is needed.

    Returns:
        Chem.Mol: The modified molecule with R-group labels assigned to the specified atoms.

    Examples:
        >>> from rdkit import Chem
        >>> smiles = "[1*]N[C@@H](CS[3*])C([2*])=O"
        >>> mol = Chem.MolFromSmiles(smiles)
        >>> indices = [0, 5, 8]  # Atom indices for [1*], [2*], and [3*]
        >>> labeled_mol = set_rgroup_labels(mol, indices)
        >>> # The atoms at indices 0, 5, and 8 are now labeled as R₁, R₂, and R₃

    """
    for n, idx in enumerate(indices):
        if idx == -1:
            continue
        mol.GetAtomWithIdx(idx).SetProp("_displayLabel", labels[n])
        # mol.GetAtomWithIdx(idx).SetProp("_displayLabel",f"R<sub>{n+1}</sub>")
    
    return mol


def generate_molecule_image_binary(mol: Chem.Mol, image_format: str="png", width: int=300, height: int=300, is_alpha_background: bool=False):
    # set drawing function as png or svg compatible
    drawing_function = rdMolDraw2D.MolDraw2DCairo if image_format == 'png' else rdMolDraw2D.MolDraw2DSVG

    # Instantiate a drawer object
    drawer = drawing_function(width=width, height=height)
    drawer_options = drawer.drawOptions()

    # Set background transparency
    if is_alpha_background:
        drawer_options.setBackgroundColour((0, 0, 0, 0))

    ## Color specific atoms
    # colors = {
    #     0: (153/255,255/255,51/255, 0.5),
    #     5: (0.2,0.8,0.5, 0.5)
    # }
    # rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=[0, 5], highlightAtomColors=colors, highlightAtomRadii={0: 0.75})

    # Prepare and draw the molecule
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)

    # Finish drawing and get the drawing text
    drawer.FinishDrawing()
    drawing_text = drawer.GetDrawingText()

    return drawing_text


def write_molecule_image(mol: Chem.Mol, outdir: str | Path, out_basename: str, image_format: str = 'png', is_alpha_background=False, width: int = 300, height: int = 300):
    """
    Generate and save an image of a chemical structure from a SMILES string.

    This function takes a SMILES (Simplified Molecular Input Line Entry System) string representing a chemical compound and generates an image depicting its molecular structure. The image is saved to the specified output directory with the given name and format.

    Args:
        smiles (str): The SMILES string representing the chemical structure.
        outdir (str or Path): The directory where the output image will be saved.
        outname (str): The name of the output image file (without the file extension).
        image_format (str, optional): The format of the output image. Supported formats are `'png'` and `'svg'`. Defaults to `'png'`.
        is_alpha_background (bool, optional): If `True`, the background of the image will be transparent (alpha channel enabled). Applicable only when `image_format` is `'png'`. Defaults to `False`.
        width (int, optional): The width of the output image in pixels. Defaults to `250`.
        height (int, optional): The height of the output image in pixels. Defaults to `250`.

    Returns:
        None

    Raises:
        ValueError: If the SMILES string is invalid or cannot be parsed.
        IOError: If the image cannot be saved to the specified directory.
        OSError: If the output directory does not exist or is not writable.

    Examples:
        >>> draw_from_smiles('C1=CC=CC=C1', 'output_images', 'benzene', image_format='svg')
        >>> draw_from_smiles('CCO', '/images', 'ethanol', is_alpha_background=True, width=300, height=300)
    """

    clean_basename = Path(out_basename).stem
    image_output = str( Path(outdir) / f'{clean_basename}.{image_format}')

    drawing_text = generate_molecule_image_binary(
        mol=mol,
        image_format=image_format,
        width=width, height=height,
        is_alpha_background=is_alpha_background
    )

    # Write the image
    with open(image_output, "wb") as f:
        f.write(drawing_text)


def generate_rgroup_molecule(mol: Chem.Mol):
    # Get "isotope" indices
    indices = get_isotope_index(mol=mol)

    # Set R group label(s)
    if sum(indices) != -4:
        rgroups = [f"R<sub>{i+1}</sub>" if x != -1 else None for i,x in enumerate(indices) ]
        mol = set_rgroup_labels(mol=mol, indices=indices, labels=rgroups)

    drawing_text = generate_molecule_image_binary(
        mol=mol,
        image_format="png",
        is_alpha_background=True
    )

    return drawing_text


def get_canonic_smiles(abbr: str) -> str:
    """Computes and returns a complete SMILE of the abbr monomer.

    Args:
        abbr (str): Valid abbreviation name of a pyPept monomer

    Returns:
        str: SMILES of the monomer abbr

    Examples:
        >>> from rdkit import Chem
        >>> smiles = "[1*]N[C@@H](CS[3*])C([2*])=O"
        >>> mol = Chem.MolFromSmiles(smiles)
        >>> indices = [0, 5, 8]  # Atom indices for [1*], [2*], and [3*]
        >>> labeled_mol = set_rgroup_labels(mol, indices)
        >>> # The atoms at indices 0, 5, and 8 are now labeled as R₁, R₂, and R₃

    """
    # Start the Sequence object
    seq = Sequence(abbr)

    # Generate the complete SMILES
    mol = Molecule(seq)
    smiles = mol.get_molecule(fmt='Smiles')

    return smiles


def get_descriptors(smiles: str) -> Dict:
    """Computes all RDkit descriptors

    Args:
        smiles (str): Smiles format of the molecule

    Returns:
        Dict: Dictionary with rdkit descriptors as keys and their values
    """
    mol = Chem.MolFromSmiles(smiles)

    # Create a dictionary to hold descriptor values
    descriptor_values = {}

    # Loop over all descriptors in the _descList
    for desc_name, function in Descriptors._descList:  # from rdkit.Chem import Descriptors
        try:
            descriptor_values[desc_name] = function(mol)
        except:
            descriptor_values[desc_name] = None

    return descriptor_values


def convert_rgroup_smiles(rgroup_smiles, rgroup_values):
    """
    Replace R-group placeholders ([1*], [2*], [3*], [4*]) in a SMILES with specific R-group values.
    
    Args:
        rgroup_smiles (str): The SMILES with R-groups (e.g., "[1*]N[C@@H](C)C([2*])=O")
        rgroup_values (list): The substituents for each R-group (e.g., ['H', 'OH', None, None])
        
    Returns:
        str: The resulting SMILES with R-groups replaced (e.g., "C[C@H](N)C(=O)O")
    """
    # Loop over each possible R-group [1*], [2*], [3*], [4*]
    for i, rgroup in enumerate(rgroup_values, start=1):
        placeholder = f'[{i}*]'
        
        if rgroup is None:
            continue  # Skip if no R-group is defined (keep [i*] as is)
        
        if rgroup == 'H':
            # If R-group is H, remove the placeholder completely (implicit hydrogen)
            rgroup_smiles = rgroup_smiles.replace(placeholder, '')
        else:
            # Otherwise, replace the placeholder with the explicit R-group (like OH, CH3, etc.)
            rgroup_smiles = rgroup_smiles.replace(placeholder, rgroup)
    
    return rgroup_smiles
