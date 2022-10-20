import uuid
import argparse

name = "minimal_polymer"
job_id = uuid.uuid4()

parser = argparse.ArgumentParser()
parser.add_argument("--n_polymers", nargs="?", default=1, const=1,  type=int, help="Number of polymers to simulate.")
parser.add_argument("--n_bead_per_chain", nargs="?", default=50, const=50, type=int, help="Number of beads per one chain.")
args = parser.parse_args()
print(args)