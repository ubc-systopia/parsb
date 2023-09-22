import docker
import argparse
import multiprocessing
from pathlib import Path
import os
from Docker import Docker


def main(args):
    client = docker.from_env()
    parsb_root = str(Path(os.path.realpath(__file__)).parent.parent.absolute())
    d = Docker(client, args, parsb_root)
    d.pull()
    if not d.container_running():
        d.start_container()
    else:
        print(f"{d.container_name} already running..")

    d.build()
    d.stop_container()
    return


if __name__ == "__main__":
    argparse_desc = """
    Configure and Build parsb in Docker
    """
    parser = argparse.ArgumentParser(
        description=argparse_desc, formatter_class=argparse.RawTextHelpFormatter
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
