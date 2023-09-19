import docker
import argparse
import multiprocessing

# docker constants
DOCKER_USER_NAME = "atrostan"
DOCKER_REPO = "parsb"
DOCKER_TAG = "latest"


def get_num_threads():
    return multiprocessing.cpu_count()


class Docker():
    client = None
    input_args = None
    graph_path = ""
    output_path = ""
    percent = 0.005
    num_threads = get_num_threads()
    block_width = 65_536
    block_size = 1024  # unused constant, but required for spray compilation
    time_execution = False

    def __init__(self, _client, _input_args):
        self.client = _client
        self.input_args = _input_args
        self.set_default_params()

    def set_default_params(self,):
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

    def info(self,):
        """Display the parameters for this run of parsb
        """
        desc_strs = {
            "Graph Edgelist Path": self.graph_path,
            "SlashBurn Output Path": self.output_path,
            "SlashBurn Percent": self.percent,
            "SlashBurn Number of Threads": self.num_threads,
            "Spray Block Width": self.block_width,
            "Timing Subroutine Execution?": self.time_execution
        }

        def print_aligned_info(desc_str, val):
            print(f"{desc_str : <30}{val}")

        print(f"-"*86)
        print(f"parsb run Parameters: \n")

        for k, v in desc_strs.items():
            print_aligned_info(k, v)
        print(f"-"*86)

        return

    def pull(self,):
        self.client.images.pull(
            repository=f"{DOCKER_USER_NAME}/{DOCKER_REPO}",
            tag=f"{DOCKER_TAG}"
        )
        return

    # def res

    def configure(self,):
        return

    def build(self,):
        return

    def compile(self,):
        self.configure()
        self.build()
        return


def main(args):

    client = docker.from_env()

    d = Docker(client, args)
    d.info()

    return


if __name__ == '__main__':
    argparse_desc = """
    Compute the SlashBurn ordering of an input graph
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-g', '--graph-path', required=True,
                        help='path to edgelist')

    parser.add_argument('-o', '--output-path', required=True,
                        help='path to store slashburn ordering')

    parser.add_argument('-p', '--percent', required=False,
                        help='SlashBurn hyperparameter, default=0.005; k=pn hubs will be removed in every SlashBurn iteration, where n is the number of vertices in the graph')

    parser.add_argument('-t', '--num-threads', required=False,
                        help='Default: number of available threads.')

    parser.add_argument('-b', '--block-width', required=False,
                        help='Spray Hyperparameter; Default: 65536; Recommended: experiment by setting to few multiples of size of L2 cache')

    parser.add_argument('--time-execution', action=argparse.BooleanOptionalAction,
                        help='if set, each subroutine of parsb will be timed for debugging')
    args = parser.parse_args()

    main(args)
