from rdkit import Chem

def handle_protonation(mol, ph):
    """
    Handle N-/C-terminal and side-chain protonation/deprotonation 
    for zwitterion generation at a specific pH.
    Args:
        mol (rdkit.Chem.Mol): RDKit molecule.
        ph (float): pH value.
    Returns:
        rdkit.Chem.Mol: Modified molecule.
    """
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    
    for atom in mol.GetAtoms():
        atom_symbol = atom.GetSymbol()
        
        # Handle N-terminal (NH3+ ↔ NH2)
        if atom_symbol == 'N' and atom.GetFormalCharge() == 1:  # NH3+
            if ph > 9.0:
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(2)  # NH2
        
        # Handle C-terminal (COOH ↔ COO-)
        if atom_symbol == 'C' and atom.GetIsAromatic() == False:  # Check non-aromatic carbons
            bonded_oxygen_atoms = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
            if len(bonded_oxygen_atoms) >= 2:  # C bonded to two oxygens (likely COOH)
                for neighbor in bonded_oxygen_atoms:
                    if neighbor.GetFormalCharge() == 0 and neighbor.GetTotalNumHs() > 0:  # Check for COOH
                        if ph > 3.0:
                            neighbor.SetFormalCharge(-1)
                            neighbor.SetNumExplicitHs(0)  # Remove hydrogen from COO-
                            break  # Only deprotonate one oxygen

        # Handle Aspartate / Glutamate side-chain (COOH ↔ COO-)
        if atom_symbol == 'C':
            bonded_oxygen_atoms = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
            if len(bonded_oxygen_atoms) >= 2:
                for neighbor in bonded_oxygen_atoms:
                    if neighbor.GetFormalCharge() == 0 and neighbor.GetTotalNumHs() > 0:  # Check for side-chain COOH
                        if ph > 4.0:
                            neighbor.SetFormalCharge(-1)
                            neighbor.SetNumExplicitHs(0)  # Remove hydrogen from COO-
                            break  # Only deprotonate one oxygen

        # Handle Lysine side-chain (NH2 ↔ NH3+)
        if atom_symbol == 'N' and atom.GetFormalCharge() == 0:
            if any(n.GetSymbol() == 'C' for n in atom.GetNeighbors()):  # Side-chain nitrogen
                if ph < 10.5:
                    atom.SetFormalCharge(1)
                    atom.SetNumExplicitHs(3)  # NH3+

        # Handle Arginine side-chain (guanidinium)
        if atom_symbol == 'N' and atom.GetFormalCharge() == 0:
            if any(n.GetSymbol() == 'C' for n in atom.GetNeighbors()):  # Side-chain nitrogen
                if ph < 12.5:
                    atom.SetFormalCharge(1)
                    atom.SetNumExplicitHs(2)  # Guanidinium (NH2+)

        # Handle Histidine side-chain (imidazole)
        if atom_symbol == 'N' and atom.GetFormalCharge() == 0:
            if any(n.GetIsAromatic() for n in atom.GetNeighbors()):  # Part of imidazole
                if ph < 6.0:
                    atom.SetFormalCharge(1)
                    atom.SetNumExplicitHs(1)  # Protonate one N in imidazole

    # Remove explicit hydrogens and sanitize
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    
    return mol


if __name__ == "__main__":
    # Example usage
    smiles = 'C[C@H](N)C(=O)O'  # Alanine
    mol = Chem.MolFromSmiles(smiles)
    mol = handle_protonation(mol, ph=7.4)  # Apply protonation logic for pH 7.4

    # Get the new SMILES
    new_smiles = Chem.MolToSmiles(mol)
    print(f"Zwitterion SMILES at pH 7.4: {new_smiles}")
