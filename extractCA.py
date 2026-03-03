import os
import numpy as np
from Bio.PDB import PDBParser

def extract_calpha_coords(pdb_file):
    """
    Return an (N,3) numpy array of CA atom coordinates from a PDB file.

    Example:
    [12.4,  8.1, -3.5],   # Cα of residue 1
    [13.0,  8.7, -2.9],   # Cα of residue 2
    [14.1,  9.2, -1.8],   # Cα of residue 3
    ...

    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("str", pdb_file)
    
    ca_coords = list()

    for model in structure:
        chains = list(model)
        first_chain = chains[0]
        for residue in first_chain:
            if "CA" in residue:
                ca_coords.append(residue["CA"].get_coord())
    
    return np.array(ca_coords)

def load_all_calpha(folder):
    """Load all PDBs in a directory and return padded CA arrays (M, L, 3)."""

    pdb_files = sorted([
        os.path.join(folder, f)
        for f in os.listdir(folder)
        if f.lower().endswith(".pdb")
    ])

    all_cords = list()

    # Extract all C-alpha coordinates of all PDB files
    for pdb in pdb_files:
        coords = extract_calpha_coords(pdb)
        all_cords.append(coords)
    
    ca_matrix = np.stack(all_cords, axis=0) # (M, 99, 3)
        
    return ca_matrix

if __name__ == "__main__":
    protein_name = "alpha_synuclein_1-99"
    folder = "Ribbons/TestRibbons/"
    ca_matrix = load_all_calpha(folder)
    print("Final array shape:", ca_matrix.shape)
    np.save(f"Calpha/{protein_name}_calpha_coordinates.npy", ca_matrix)
