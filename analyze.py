"""Detect feigned space-time avalanches in an IID time series and plot them.

A feigned correlation is a link between a source site (r, c, t) and a
descendant (r', c', t+1) inside the Moore neighbourhood of the source
(Chebyshev distance <= 1, including the same cell) whose values satisfy
|v_source - v_descendant| <= z. There are no within-frame links: two cells
observed in the same iteration are never connected, however similar they are.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MOORE_OFFSETS = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, node):
        parent = self.parent
        root = node
        while parent.setdefault(root, root) != root:
            root = parent[root]
        while parent[node] != root:
            parent[node], node = root, parent[node]
        return root

    def union(self, a, b):
        root_a, root_b = self.find(a), self.find(b)
        if root_a != root_b:
            self.parent[root_b] = root_a


def load_frames(csv_path):
    """Reconstruct the long-form CSV into an array of shape (T, n, n)."""
    table = pd.read_csv(csv_path)
    iterations = int(table["iteration"].max()) + 1
    n = int(table["row"].max()) + 1
    if len(table) != iterations * n * n:
        raise ValueError(f"{csv_path} does not contain a complete T x n x n grid")
    ordered = table.sort_values(["iteration", "row", "col"])
    return ordered["value"].to_numpy().reshape(iterations, n, n)


def slice_bounds(offset, n):
    """Source and destination index ranges along one axis for a Moore offset."""
    src_start = max(0, -offset)
    src_stop = n - max(0, offset)
    return src_start, src_stop, src_start + offset, src_stop + offset


def offset_deltas(frames, dr, dc):
    """Absolute differences and node ids for one Moore offset across all t -> t+1."""
    _, n, _ = frames.shape
    r0, r1, rd0, rd1 = slice_bounds(dr, n)
    c0, c1, cd0, cd1 = slice_bounds(dc, n)
    source = frames[:-1, r0:r1, c0:c1]
    descendant = frames[1:, rd0:rd1, cd0:cd1]
    return np.abs(source - descendant), (r0, c0, rd0, cd0)


def build_graph(frames, z):
    """Union every space-time site pair that satisfies the feigned-correlation rule."""
    _, n, _ = frames.shape
    uf = UnionFind()
    edges = []
    deltas = []

    for dr, dc in MOORE_OFFSETS:
        delta, (r0, c0, rd0, cd0) = offset_deltas(frames, dr, dc)
        deltas.append(delta.ravel())

        for t, i, j in zip(*np.nonzero(delta <= z)):
            src_node = int((t * n + (r0 + i)) * n + (c0 + j))
            dst_node = int(((t + 1) * n + (rd0 + i)) * n + (cd0 + j))
            uf.union(src_node, dst_node)
            edges.append((src_node, dst_node))

    return uf, edges, np.concatenate(deltas)


AVALANCHE_COLUMNS = [
    "size",
    "duration",
    "spatial_extent",
    "start_iteration",
    "end_iteration",
    "edge_count",
]


def collect_avalanches(uf, edges, n):
    """Group connected space-time sites into avalanches and measure each one."""
    components = {}
    for node in list(uf.parent):
        components.setdefault(uf.find(node), []).append(node)

    edges_per_root = {}
    for src_node, _ in edges:
        root = uf.find(src_node)
        edges_per_root[root] = edges_per_root.get(root, 0) + 1

    records = []
    for root, nodes in components.items():
        times = [node // (n * n) for node in nodes]
        cells = {node % (n * n) for node in nodes}
        records.append(
            {
                "size": len(nodes),
                "duration": max(times) - min(times) + 1,
                "spatial_extent": len(cells),
                "start_iteration": min(times),
                "end_iteration": max(times),
                "edge_count": edges_per_root.get(root, 0),
            }
        )
    frame = pd.DataFrame(records, columns=AVALANCHE_COLUMNS)
    return frame.sort_values("size", ascending=False)


def plot_distribution(values, xlabel, title, path):
    counts = pd.Series(values).value_counts().sort_index()
    figure, axes = plt.subplots(figsize=(6, 4.5))
    axes.scatter(counts.index, counts.to_numpy(), s=28, color="#2b6cb0")
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlabel(xlabel)
    axes.set_ylabel("number of avalanches")
    axes.set_title(title)
    axes.grid(True, which="both", linewidth=0.3, alpha=0.5)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def plot_edge_rate(sweep, path):
    figure, axes = plt.subplots(figsize=(6, 4.5))
    axes.plot(sweep["z"], sweep["edge_rate"], marker="o", color="#b7791f")
    axes.set_xscale("log")
    axes.set_xlabel("tolerance z")
    axes.set_ylabel("feigned-correlation edge rate R(z)")
    axes.set_title("Fraction of candidate temporal pairs within z")
    axes.grid(True, which="both", linewidth=0.3, alpha=0.5)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--z", type=float, required=True)
    parser.add_argument("--plot-dir", required=True)
    parser.add_argument(
        "--z-sweep",
        default="0.01,0.03,0.1,0.3,1,3,10,30",
        help="comma-separated tolerances for the edge-rate diagnostic",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)

    frames = load_frames(args.csv)
    iterations, n, _ = frames.shape

    uf, edges, deltas = build_graph(frames, args.z)
    avalanches = collect_avalanches(uf, edges, n)

    total_candidates = deltas.size
    edge_rate = float(np.mean(deltas <= args.z))

    sweep_values = [float(item) for item in args.z_sweep.split(",") if item.strip()]
    sweep = pd.DataFrame(
        {
            "z": sweep_values,
            "edge_rate": [float(np.mean(deltas <= value)) for value in sweep_values],
        }
    )

    avalanche_path = os.path.join(args.plot_dir, "avalanches.csv")
    sweep_path = os.path.join(args.plot_dir, "edge_rate_vs_z.csv")
    avalanches.to_csv(avalanche_path, index=False)
    sweep.to_csv(sweep_path, index=False)

    if not avalanches.empty:
        plot_distribution(
            avalanches["size"],
            "avalanche size (space-time sites)",
            "Avalanche size distribution",
            os.path.join(args.plot_dir, "avalanches_vs_size.png"),
        )
        plot_distribution(
            avalanches["duration"],
            "avalanche duration (iterations)",
            "Avalanche duration distribution",
            os.path.join(args.plot_dir, "avalanches_vs_duration.png"),
        )
        plot_distribution(
            avalanches["spatial_extent"],
            "spatial extent (unique cells)",
            "Avalanche spatial extent distribution",
            os.path.join(args.plot_dir, "avalanches_vs_spatial_extent.png"),
        )
    plot_edge_rate(sweep, os.path.join(args.plot_dir, "edge_rate_vs_z.png"))

    print(f"source: {args.csv}  frames={iterations}  grid={n}x{n}  z={args.z}")
    print(f"candidate temporal pairs: {total_candidates}")
    print(f"feigned-correlation edges: {len(edges)}  R(z)={edge_rate:.6f}")
    print(f"avalanches detected: {len(avalanches)}")
    if not avalanches.empty:
        print(
            f"largest: size={avalanches['size'].iloc[0]} "
            f"duration={avalanches['duration'].iloc[0]} "
            f"spatial_extent={avalanches['spatial_extent'].iloc[0]}"
        )
        print(
            f"median size={avalanches['size'].median()} "
            f"median duration={avalanches['duration'].median()}"
        )
    print(f"wrote metrics to {avalanche_path} and {sweep_path}")


if __name__ == "__main__":
    main()
