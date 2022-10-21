#
# Copyright (C) 2013-2022 The ESPResSo project
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
Set up a linear polymer.
"""
import espressomd
import espressomd.thermostat
import espressomd.interactions
import espressomd.polymer
import espressomd.io.writer.vtf  # pylint: disable=import-error
import uuid
import argparse

espressomd.assert_features(["WCA"])

import numpy as np

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
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4
system.cell_system.set_n_square(use_verlet_lists=False)
#outfile = open('polymer.vtf', 'w')

system.non_bonded_inter[0, 0].wca.set_params(epsilon=1, sigma=1)

fene = espressomd.interactions.FeneBond(k=10, d_r_max=2)
system.bonded_inter.add(fene)
#espressomd.polymer.create_polymer(N_P=args.n_polymers, MPC=args.n_beads_per_chain, bond_length=1.0, bond=fene)

positions = espressomd.polymer.linear_polymer_positions(
    n_polymers=args.n_polymers, beads_per_chain=args.n_beads_per_chain, bond_length=1.0, seed=3210)
previous_part = None
for pos in positions[0]:
    part = system.part.add(pos=pos)
    if previous_part:
        part.add_bond((fene, previous_part))
    previous_part = part

#espressomd.io.writer.vtf.writevsf(system, outfile)

#############################################################
# Files
#############################################################
logfile = open(f"{name}_{args.n_beads_per_chain}_{job_id}.log", 'w')
xyzfile = open(f"{name}_{args.n_beads_per_chain}_{job_id}.xyz", 'w')
obsfile = open(f"{name}_{args.n_beads_per_chain}_{job_id}.obs", 'w')


#############################################################
#      Warmup                                               #
#############################################################

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=0.01)
while system.analysis.min_dist() < 0.95:
    logfile.write(f"minimization: {system.analysis.energy()['total']:+.2e}\n")
    system.integrator.run(20)

logfile.write(f"minimization: {system.analysis.energy()['total']:+.2e}\n")
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)


#############################################################
#      Integration                                          #
#############################################################

logfile.write("simulating...\n")
t_steps = 1000
for t in range(t_steps):
    logfile.write(f"step {t + 1} of {t_steps}\n", end='\r', flush=True)
    system.integrator.run(10)
    gyration_tensor = system.analysis.gyration_tensor()
    # xyzfile.write(f"%\n")
    # for p in system.part:
    #     xyzfile.write(f"part {p.pos[0]} {p.pos[1]} {p.pos[2]}")
    #     xyzfile.write(f"\n")
    obsfile.write(f"{gyration_tensor['Rg^2']}\n")
