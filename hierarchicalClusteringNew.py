import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# --- CONFIG ---
protein_name = "alpha_synuclein_1-99"
pdb_folder = "Ribbons/TestRibbons/"
mutual_q_file = f"Qmatrix/{protein_name}_mutual_Q_matrix.npy"

# choose how many clusters (polymorphs) you want
# e.g., 5–6, similar to the paper
target_n_clusters = 6

# --- 1. Load mutual-Q matrix ---
mutual_Q = np.load(mutual_q_file)
M = mutual_Q.shape[0]

# --- 2. Distance matrix: d = 1 - Q ---
dist_matrix = 1.0 - mutual_Q
np.fill_diagonal(dist_matrix, 0.0)
dist_matrix = 0.5 * (dist_matrix + dist_matrix.T)
condensed_dist = squareform(dist_matrix, checks=False)

# --- 3. Hierarchical clustering (average linkage) ---
Z = linkage(condensed_dist, method="average")

# --- 4. Flat clusters: request target_n_clusters ---
cluster_ids = fcluster(Z, t=target_n_clusters, criterion="maxclust")
unique_clusters = np.unique(cluster_ids)
n_clusters = len(unique_clusters)
print(f"Number of clusters: {n_clusters}")

# --- 5. Map indices to PDB filenames ---
pdb_files = sorted(
    f for f in os.listdir(pdb_folder) if f.lower().endswith(".pdb")
)
assert len(pdb_files) == M, "Number of PDBs must match matrix size"

# full paths (if needed later)
pdb_paths = [os.path.join(pdb_folder, f) for f in pdb_files]

# --- 6. Build cluster -> member indices ---
clusters = {}
for idx, cid in enumerate(cluster_ids):
    clusters.setdefault(cid, []).append(idx)

# --- 7. For each cluster: mean Q and centroid (medoid) ---
cluster_info = []

for cid, members in clusters.items():
    members = np.array(members, dtype=int)
    sub_Q = mutual_Q[np.ix_(members, members)]

    # mean mutual-Q inside cluster (including diagonal)
    mean_q = sub_Q.mean()

    # centroid = structure with highest avg mutual-Q to others
    mean_q_per_structure = sub_Q.mean(axis=1)
    local_best = np.argmax(mean_q_per_structure)
    centroid_idx = members[local_best]

    centroid_pdb = pdb_files[centroid_idx]
    centroid_path = pdb_paths[centroid_idx]

    cluster_info.append(
        {
            "cluster_id": cid,
            "size": len(members),
            "mean_mutual_Q": float(mean_q),
            "centroid_index": int(centroid_idx),
            "centroid_pdb": centroid_pdb,
            "centroid_path": centroid_path,
        }
    )

# --- 8. Print summary & save text report ---
out_base = f"Clusters/final_results/{protein_name}_{target_n_clusters}_new_clusters"
os.makedirs(os.path.dirname(out_base + ".txt"), exist_ok=True)

with open(out_base + ".txt", "w") as f:
    for info in sorted(cluster_info, key=lambda x: x["cluster_id"]):
        line = (
            f"Cluster {info['cluster_id']:>3} | "
            f"size: {info['size']:>3} | "
            f"mean Q: {info['mean_mutual_Q']:.3f} | "
            f"centroid: idx {info['centroid_index']} "
            f"({info['centroid_pdb']})"
        )
        f.write(line + "\n")
        print(line)

# --- 9. Clustered mutual-Q heatmap with boundaries ---
sorted_indices = np.argsort(cluster_ids)
sorted_Q = mutual_Q[np.ix_(sorted_indices, sorted_indices)]
sorted_cluster_ids = cluster_ids[sorted_indices]

plt.figure(figsize=(6, 5))
im = plt.imshow(sorted_Q, vmin=0.20, vmax=0.8, cmap="turbo")
plt.colorbar(im, label="Mutual-Q")

current = 0
for cid in sorted(unique_clusters):
    size = np.sum(sorted_cluster_ids == cid)
    current += size
    plt.axhline(current - 0.5, color="white", linewidth=0.8)
    plt.axvline(current - 0.5, color="white", linewidth=0.8)

plt.xlabel("Structure index")
plt.ylabel("Structure index")

plt.tight_layout()
plt.savefig(out_base + ".png", dpi=300)
plt.close()

print(f"Clustered heatmap saved to: {out_base}.png")
print(f"Cluster summary saved to: {out_base}.txt")

# --- 10. Dendrogram animation (build process, stop at target_n_clusters) ---
fig, ax = plt.subplots(figsize=(8, 4))

heights = np.unique(Z[:, 2])
heights.sort()

def n_clusters_at_height(h):
    labels = fcluster(Z, t=h, criterion="distance")
    return len(np.unique(labels))

target_height = heights[-1]
for h in heights:
    k = n_clusters_at_height(h)
    if k <= target_n_clusters:
        target_height = h
        break

print(f"Distance threshold for ≈{target_n_clusters} clusters: {target_height:.4f}")

# heights used in animation (do not go beyond target_height)
anim_heights = [h for h in heights if h <= target_height]
if len(anim_heights) == 0:
    anim_heights = [target_height]

# extend with copies of the last height so final frame is held ~5 s
interval_ms = 300                 # frame interval
extra_time_s = 5.0
extra_frames = int(extra_time_s * 1000 / interval_ms)  # ~17
anim_heights = anim_heights + [anim_heights[-1]] * extra_frames

n_frames = len(anim_heights)

def update(i):
    ax.clear()
    h = anim_heights[i]
    k = n_clusters_at_height(h)

    dendrogram(
        Z,
        no_labels=True,
        color_threshold=h,
        ax=ax,
    )
    ax.set_ylabel("Distance (1 - Q)")
    ax.set_xlabel("Structure index")
    ax.set_title(f"Dendrogram build (cut ≤ {h:.3f})  |  clusters = {k}")
    ax.axhline(h, color="red", linestyle="--", linewidth=1)
    return ax

anim = FuncAnimation(fig, update, frames=n_frames, interval=interval_ms, blit=False)

gif_path = f"{out_base}_dendrogram_build.gif"
anim.save(gif_path, writer="pillow", dpi=150)
plt.close(fig)

print(f"Dendrogram animation saved to: {gif_path}")