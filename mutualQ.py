import numpy as np

def align_and_mask_coords(struct_coord, ref_coord):
    """
    Align struct_coord and ref_coord by a joint NaN mask so that
    row i in both arrays always corresponds to the same residue index.
    """
    if struct_coord.shape != ref_coord.shape:
        raise ValueError("struct_coord and ref_coord must have same shape")

    joint_mask = (~np.isnan(struct_coord).any(axis=1) &
                  ~np.isnan(ref_coord).any(axis=1))

    struct_clean = struct_coord[joint_mask]
    ref_clean = ref_coord[joint_mask]

    if struct_clean.shape[0] != ref_clean.shape[0]:
        raise RuntimeError("Masked coordinate arrays have inconsistent length")

    # optional: keep track of which residue indices remain
    residue_indices = np.nonzero(joint_mask)[0]

    return struct_clean, ref_clean, residue_indices

def compute_Q_score(struct_coord, ref_coord):
    """
    Compute the Q-score (equation 5 of report) between a structure and a reference.
    struct_coords: (N_residues, 3) numpy array
    ref_coords: (N_residues, 3) numpy array
    """
    struct_coord, ref_coord, residue_indices = align_and_mask_coords(
        struct_coord, ref_coord
    )
    
    # mask out NaN rows, in case some centroids are NaN
    mask = ~np.isnan(struct_coord).any(axis=1) & ~np.isnan(ref_coord).any(axis=1)
    struct_coord = struct_coord[mask]
    ref_coord = ref_coord[mask]

    # Remove NaN rows (padded residues)
    N = struct_coord.shape[0]

    if N < 4:
        return np.nan  # formula requires at least 4 residues

    # Compute pairwise distance

    Q_sum = 0.0

    for i in range(N-2):
        for j in range(i+3,N):
            rij = np.linalg.norm(struct_coord[i]-struct_coord[j])
            rij_N = np.linalg.norm(ref_coord[i]-ref_coord[j])
            diff = (rij - rij_N)**2
            sigma0 = 1.0
            #sigma_ij = sigma0 * (1.0 + abs(i-j)**0.15)
            sigma_ij = abs(i-j)**0.15 #https://pmc.ncbi.nlm.nih.gov/articles/PMC3406225/#S2 (see supllementary page 14)
            Q_sum += np.exp(-diff / (2.0 * sigma_ij**2)) 
    
    Q = (2.0 / ((N-2)*(N-3))) * Q_sum
    return Q

def compute_mutual_Q(struct_A, struct_B):
    q_AB = compute_Q_score(struct_A, struct_B)
    q_BA = compute_Q_score(struct_B, struct_A)

    if np.isnan(q_AB) or np.isnan(q_BA):
        return np.nan
    return 0.5 * (q_AB + q_BA)  # mutual-Q

if __name__ == "__main__":
    protein_name = "alpha_synuclein_1-99"

    # Load your saved Cα coordinates: shape (M, L, 3), M ≈ 100
    centroid_matrix = np.load(f"Calpha/{protein_name}_calpha_coordinates.npy")
    M = centroid_matrix.shape[0]

    mutual_Q_matrix = np.zeros((M, M), dtype=float) # Create a 100*100 mutual-Q value matrix
    for i in range(M):
        mutual_Q_matrix[i, i] = 1.0  # self-similarity
        for j in range(i + 1, M):
            mq = compute_mutual_Q(centroid_matrix[i], centroid_matrix[j])
            mutual_Q_matrix[i, j] = mq
            mutual_Q_matrix[j, i] = mq  # symmetry

    np.save(f"Qmatrix/{protein_name}_mutual_Q_matrix.npy", mutual_Q_matrix)
    print("mutual_Q_matrix shape:", mutual_Q_matrix.shape)

