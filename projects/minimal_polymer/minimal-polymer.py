#
# Copyright (C) 2013-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
This sample sets up a polymer.
"""
import espressomd
espressomd.assert_features(["WCA"])
from espressomd import thermostat
from espressomd import interactions
from espressomd import polymer
from espressomd.io.writer import vtf  # pylint: disable=import-error
import numpy as np
import uuid
import argparse

# Default simulation values
#############################################################
name = "minimal-polymer"
job_id = uuid.uuid4()

parser = argparse.ArgumentParser()
parser.add_argument("--n_polymers", nargs="?", default=1, const=1,  type=int, help="Number of polymers to simulate.")
parser.add_argument("--n_beads_per_chain", nargs="?", default=5, const=5, type=int, help="Number of beads per one chain.")
parser.add_argument("--box_l", nargs="?", default=100, const=100, type=int, help="Size of the simulation box.")
args = parser.parse_args()

# System parameters
#############################################################

system = espressomd.System(box_l=[args.box_l, args.box_l, args.box_l])
system.set_random_state_PRNG()
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
system.cell_system.set_n_square(use_verlet_lists=False)

system.non_bonded_inter[0, 0].wca.set_params(
    epsilon=1, sigma=1)

fene = interactions.FeneBond(k=10, d_r_max=2)
system.bonded_inter.add(fene)


positions = polymer.positions(n_polymers=args.n_polymers,
                              beads_per_chain=args.n_beads_per_chain,
                              bond_length=1.0,
                              seed=3210)
for i, pos in enumerate(positions[0]):
    id = len(system.part)
    system.part.add(id=id, pos=pos)
    if i > 0:
        system.part[id].add_bond((fene, id - 1))


#############################################################
# Files
#############################################################
logfile = open(f"{name}_{args.n_beads_per_chain}_{job_id}.log", 'w')
obsfile = open(f"{name}_{args.n_beads_per_chain}_{job_id}.obs", 'w')

#############################################################
#      Warmup                                               #
#############################################################

warm_steps = 10
wca_cap = 1
system.force_cap = wca_cap
i = 0
act_min_dist = system.analysis.min_dist()

# warmup with zero temperature to remove overlaps
system.thermostat.set_langevin(kT=0.0, gamma=1.0)

# slowly ramp un up the cap
while (act_min_dist < 0.95):
    logfile.write("min_dist: {} \t force cap: {}\n".format(act_min_dist, wca_cap))
    system.integrator.run(warm_steps)
    system.part[:].v = [0, 0, 0]
    # Warmup criterion
    act_min_dist = system.analysis.min_dist()
    wca_cap = wca_cap * 1.01
    system.force_cap = wca_cap

# remove force cap
wca_cap = 0
system.force_cap = wca_cap
system.integrator.run(warm_steps * 10)

# restore simulation temperature
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.integrator.run(warm_steps * 10)
logfile.write("Finished warmup\n")


#############################################################
#      Integration                                          #
#############################################################

logfile.write("simulating...\n")
t_steps = 1000
for t in range(t_steps):
    logfile.write("step {} of {}\n".format(t, t_steps))
    system.integrator.run(warm_steps)
    rg = espressomd.analyze.Analysis.calc_rg(chain_start=0, number_of_chains=1, chain_lenght=len(system.part[:]))
    obsfile.write(f"{rg[0]}\n")


