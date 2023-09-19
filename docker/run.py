import docker
import argparse
import multiprocessing
from pathlib import Path
import os

# docker constants
DOCKER_USER_NAME = "atrostan"
DOCKER_REPO = "parsb"
DOCKER_TAG = "latest"


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
        # get directory that stores edgelist
        # (to mount into docker)
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

    def set_default_params(
        self,
    ):
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

    def info(
        self,
    ):
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

        print(f"-" * 86)
        print(f"parsb run Parameters: \n")

        for k, v in desc_strs.items():
            print_aligned_info(k, v)
        print(f"-" * 86)

        return

    def parse_fnames(
        self,
    ):
        self.graph_fname = get_filename(self.graph_path)
        self.order_fname = get_filename(self.output_path)
        return

    def pull(
        self,
    ):
        ret = self.client.images.pull(repository=self.docker_repo, tag=self.docker_tag)
        print(f"{ret = }")
        return

    def container_running(
        self,
    ):
        """Check if the parsb container is already running"""
        clist = self.client.containers.list()
        for c in clist:
            if c.name == self.container_name:
                return True

        return False

    def get_parsb_container(
        self,
    ):
        clist = self.client.containers.list()
        print(f"{clist = }")
        for c in clist:
            if c.name == self.container_name:
                return c

        return None

    def stop_container(
        self,
    ):
        container = self.get_parsb_container()
        if container:
            container.stop()
            container.remove()
        return

    def run_container(
        self,
    ):
        """if already running, stop parsb docker container,
        and remove it
        (re)run container, mounting
        """

        if not self.container_running():
            container = self.client.containers.run(
                image=f"{self.docker_repo}:{self.docker_tag}",
                detach=True,
                tty=True,
                name=self.container_name,
                mounts=self.mounts,
            )
            print(f"{container = }")
            return container
        return None

    def run(
        self,
    ):
        container = self.get_parsb_container()
        print(f"{container = }")
        if not container:
            return
        run_cmd = ["/bin/bash", "/root/parsb/docker/run.sh"]
        ret = container.exec_run(cmd=run_cmd, environment=self.docker_env)
        print(f"{ret = }")

        return

    def configure(
        self,
    ):
        return

    def build(
        self,
    ):
        return

    def compile(
        self,
    ):
        self.configure()
        self.build()
        return


def main(args):
    client = docker.from_env()
    d = Docker(client, args)
    d.info()
    d.pull()
    d.stop_container()
    d.run_container()
    d.run()
    return


if __name__ == "__main__":
    argparse_desc = """
    Compute the SlashBurn ordering of an input graph
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-g", "--graph-path", required=True, help="path to edgelist")

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
    args = parser.parse_args()

    main(args)
