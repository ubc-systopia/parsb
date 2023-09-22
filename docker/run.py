import docker
import argparse
import multiprocessing
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Docker import Docker

# don't plot adjacency matrices with more than this many vertices
MAX_N_VERTICES_PLOT = 10_000


def plot_edgelists(el_path, mp_path, plt_path):
    """Read a graph's edgelist into a 2D numpy array representing the
    adjacency matrix of the graph
    Read a vertex ordering mapping
    Relabel the edges of the graph using the vertex ordering
    Plot both (isomorphic) adjacency matrices

    Args:
        el_path (str): path to edge list
        mp_path (str): path to vertex ordering

    Returns:
        np.array: graph's adjacency matrix
    """
    arr = pd.read_csv(
        el_path, sep=" ", header=None, names=["src", "dest"], index_col=False
    ).values

    mp = pd.read_csv(
        mp_path,
        sep=" ",
        skiprows=2,
        header=None,
        names=["old", "new"],
        index_col=False,
    ).values
    n = int(np.max(arr) + 1)
    if n > MAX_N_VERTICES_PLOT:
        print(f"Skipping plots of adjacency matrices for verification.")
        print(f"Graph has {n} vertices > {MAX_N_VERTICES_PLOT}.")
        return

    mat = np.zeros((n, n))
    map_mat = np.zeros((n, n))
    m = len(arr)
    for i in range(len(arr)):
        src = int(arr[i][0])
        dest = int(arr[i][1])
        mat[src, dest] = 1
        mp_u = mp[src][1]
        mp_v = mp[dest][1]
        map_mat[mp_v, mp_u] = 1

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
    ax1.spy(mat, markersize=1)
    ax2.spy(map_mat, markersize=1)
    ax1.set_title("Original")
    ax2.set_title("SlashBurn")
    for ax in [ax1, ax2]:
        ax.set_xlim(0, n)
        ax.set_ylim(0, n)
        ax.set_aspect("equal", "datalim")
        ax.invert_yaxis()
    plt.tight_layout()

    print(f"Plotting Original and SlashBurn adjacency matrices at:\n\t {plt_path}")
    plt.savefig(plt_path, dpi=200)
    plt.cla()
    plt.close()

    return


def main(args):
    client = docker.from_env()
    parsb_root = str(Path(os.path.realpath(__file__)).parent.parent.absolute())
    d = Docker(client, args, parsb_root)
    d.info()
    d.pull()
    if not d.container_running():
        d.start_container()
    else:
        print(f"{d.container_name} already running..")
    d.run()

    if args.plot_verify:
        plot_path = os.path.join(d.graph_dir, f"{d.graph_fname}.plot.png")
        plot_edgelists(d.graph_path, d.output_path, plot_path)
        return

    return


def add_arguments(parser):
    parser.add_argument("-f", "--graph-path", required=True, help="path to edgelist")

    parser.add_argument(
        "-o", "--output-path", required=False, help="path to store slashburn ordering"
    )

    parser.add_argument(
        "-p",
        "--percent",
        required=False,
        help="SlashBurn hyperparameter, default=0.005; k=pn hubs will be removed in every SlashBurn iteration, where n is the number of vertices in the graph",
    )

    parser.add_argument(
        "-t",
        "--num-threads",
        required=False,
        help="Default: number of available threads.",
    )

    parser.add_argument(
        "--plot-verify",
        action=argparse.BooleanOptionalAction,
        help="if set plot the original adjacency matrix and "
        "slashburn adjacency matrix",
    )
    return parser


if __name__ == "__main__":
    argparse_desc = """
    Compute the SlashBurn ordering of an input graph
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter
    )
    parser = add_arguments(parser)

    args = parser.parse_args()

    main(args)
