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
    parsb_root = ""
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
    docker_parsb_dir = "/root/parsb"
    percent = 0.005
    num_threads = get_num_threads()
    block_width = 65_536
    block_size = 1024  # unused constant, but required for spray compilation
    time_execution = False

    def __init__(self, _client, _input_args, _parsb_root):
        self.client = _client
        self.input_args = _input_args
        self.parsb_root = _parsb_root
        self.output_path = os.path.join(self.parsb_root, "data", "graphs", "slashburn_order")
        self.set_default_params()
        if self.graph_path and self.output_path:
            self.parse_fnames()
        # get directory that stores edgelist (to mount into docker)
        if self.graph_path:
            self.graph_dir = get_parent_dir(self.graph_path)
            self.mounts += [
                docker.types.Mount(
                    target=self.docker_data_dir, source=self.graph_dir, type="bind"
                )
            ]

        # mount directory of this repository
        self.mounts += [
            docker.types.Mount(
                target=self.docker_parsb_dir, source=self.parsb_root, type="bind"
            ),
        ]

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

        for name, value in vars(self.input_args).items():
            if value:
                match name:
                    case "graph_path":
                        self.graph_path = value
                    case "ouput_path":
                        self.output_path = value
                    case "percent":
                        self.percent = value
                    case "num_threads":
                        self.num_threads = value
                    case "block_width":
                        self.block_width = value
                    case "time_execution":
                        self.time_execution = value
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
            print(f"Stopping {container}..")
            container.stop()
            print(f"Removing {container}..")
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

    def build(self):
        container = self.get_parsb_container()
        if not container:
            return

        run_cmd = ["/bin/bash", "/root/parsb/docker/build.sh"]
        ret = container.exec_run(cmd=run_cmd, environment=self.docker_env)
        output = ret.output.decode("utf-8")
        [print(l) for l in output.split("\n")]
        return

    def run(self):
        container = self.get_parsb_container()
        if not container:
            return
        run_cmd = ["/bin/bash", "/root/parsb/docker/run.sh"]
        ret = container.exec_run(cmd=run_cmd, environment=self.docker_env)
        output = ret.output.decode("utf-8")
        [print(l) for l in output.split("\n")]

        return
