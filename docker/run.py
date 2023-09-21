import docker
import argparse
import multiprocessing
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# docker constants
DOCKER_USER_NAME = "atrostan"
DOCKER_REPO = "parsb"
DOCKER_TAG = "latest"

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


def get_parent_dir(file_path):
    p = Path(file_path)
    return str(p.parent.absolute())


def get_filename(file_path):
    return os.path.basename(file_path)


def get_num_threads():
    return multiprocessing.cpu_count()


class Docker:
    client = None
    input_args = None
    graph_path = ""
    output_path = ""
    graph_dir = ""
    graph_fname = ""
    order_fname = ""
    mounts = []
    docker_env = {}
    container_name = "parsb-container"
    docker_repo = f"{DOCKER_USER_NAME}/{DOCKER_REPO}"
    docker_tag = f"{DOCKER_TAG}"
    docker_data_dir = "/data"
    percent = 0.005
    num_threads = get_num_threads()
    block_width = 65_536
    block_size = 1024  # unused constant, but required for spray compilation
    time_execution = False

    def __init__(self, _client, _input_args):
        self.client = _client
        self.input_args = _input_args
        self.set_default_params()
        self.parse_fnames()
        # get directory that stores edgelist (to mount into docker)
        self.graph_dir = get_parent_dir(self.graph_path)
        self.mounts.append(
            docker.types.Mount(
                target=self.docker_data_dir, source=self.graph_dir, type="bind"
            )
        )
        self.docker_env = {
            # compile time
            "BSIZE": self.block_size,
            "BWIDTH": self.block_width,
            "TIME": int(self.time_execution),
            # runtime
            "GRAPH_PATH": os.path.join(
                self.docker_data_dir,
                self.graph_fname,
            ),
            "OUTPUT_PATH": os.path.join(
                self.docker_data_dir,
                self.order_fname,
            ),
            "PERCENT": self.percent,
            "NUM_THREADS": self.num_threads,
        }

    def set_default_params(self):
        if self.input_args.percent:
            self.percent = self.input_args.percent
        else:
            self.percent = 0.005

        if self.input_args.num_threads:
            self.num_threads = self.input_args.num_threads
        else:
            self.num_threads = get_num_threads()

        if self.input_args.block_width:
            self.block_width = self.input_args.block_width
        else:
            self.block_width = 65_536

        if self.input_args.time_execution:
            self.time_execution = self.input_args.time_execution
        else:
            self.time_execution = False

        self.graph_path = self.input_args.graph_path
        self.output_path = self.input_args.output_path
        return

    def print_sep(self):
        print(f"-" * 86)

    def info(self):
        """Display the parameters for this run of parsb"""
        desc_strs = {
            "Graph Edgelist Path": self.graph_path,
            "SlashBurn Output Path": self.output_path,
            "SlashBurn Percent": self.percent,
            "SlashBurn Number of Threads": self.num_threads,
            "Spray Block Width": self.block_width,
            "Timing Subroutine Execution?": self.time_execution,
        }

        def print_aligned_info(desc_str, val):
            print(f"{desc_str : <30}{val}")

        print(f"parsb run Parameters: \n")

        for k, v in desc_strs.items():
            print_aligned_info(k, v)
        self.print_sep()

        return

    def parse_fnames(self):
        self.graph_fname = get_filename(self.graph_path)
        self.order_fname = get_filename(self.output_path)
        return

    def pull(self):
        ret = self.client.images.pull(repository=self.docker_repo, tag=self.docker_tag)
        return

    def container_running(self):
        """Check if the parsb container is already running"""
        clist = self.client.containers.list()
        for c in clist:
            if c.name == self.container_name:
                return True

        return False

    def get_parsb_container(self):
        clist = self.client.containers.list()
        for c in clist:
            if c.name == self.container_name:
                return c

        return None

    def stop_container(self):
        container = self.get_parsb_container()
        if container:
            container.stop()
            container.remove()
        return

    def start_container(self):
        """if parsb container not already running,
        start container, mounting directory that stores input graph
        """
        img_name = f"{self.docker_repo}:{self.docker_tag}"
        if not self.container_running():
            print(f"Starting {self.container_name} from {img_name} image..")
            container = self.client.containers.run(
                image=img_name,
                detach=True,
                tty=True,
                name=self.container_name,
                mounts=self.mounts,
            )
            return container
        return None

    def run(self):
        container = self.get_parsb_container()
        if not container:
            return
        run_cmd = ["/bin/bash", "/root/parsb/docker/run.sh"]
        ret = container.exec_run(cmd=run_cmd, environment=self.docker_env)
        output = ret.output.decode("utf-8")
        [print(l) for l in output.split("\n")]

        return


def main(args):
    client = docker.from_env()
    d = Docker(client, args)
    d.info()
    d.pull()
    # d.stop_container()
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


if __name__ == "__main__":
    argparse_desc = """
    Compute the SlashBurn ordering of an input graph
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-f", "--graph-path", required=True, help="path to edgelist")

    parser.add_argument(
        "-o", "--output-path", required=True, help="path to store slashburn ordering"
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
        "-b",
        "--block-width",
        required=False,
        help="Spray Hyperparameter; Default: 65536; Recommended: experiment by setting to few multiples of size of L2 cache",
    )

    parser.add_argument(
        "--time-execution",
        action=argparse.BooleanOptionalAction,
        help="if set, each subroutine of parsb will be timed for debugging",
    )

    parser.add_argument(
        "--plot-verify",
        action=argparse.BooleanOptionalAction,
        help="if set plot the original adjacency matrix and "
        "slashburn adjacency matrix",
    )

    args = parser.parse_args()

    main(args)
