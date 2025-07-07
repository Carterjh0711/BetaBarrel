#oval_constriction
#Use /net/software/containers/pyrosetta.sif
#ls -1 | wc -1 *pdb
#qlogin --mem=24g -c 1 -p cpu

import os
import sys
import glob
import pyrosetta
import argparse
from pyrosetta import rosetta
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring.constraints import DihedralConstraint
from pyrosetta.rosetta.core.scoring.func import CircularHarmonicFunc
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.relax import FastRelax
import faulthandler
faulthandler.enable()

def min_torsion(residues, pose):
	movemap = rosetta.core.kinematics.MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(False)
	for res in residues:
		movemap.set_bb(res, True)
		movemap.set_chi(res, True)
	scorefxn = pyrosetta.get_fa_scorefxn()
	scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1.25) #default is 1.0
	scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, 0.75) #default is 0.5

	min_mover = rosetta.protocols.minimization_packing.MinMover()
	min_mover.movemap(movemap)
	min_mover.score_function(scorefxn)
	min_mover.min_type("dfpmin_armijo_nonmonotone")

	min_mover.apply(pose)


def fast_relax(residues, pose):
	#Alternative to min_torsion mover
	scorefxn = get_fa_scorefxn()
	scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1.0)
	scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, 0.5)
	relax = FastRelax()
	relax.set_scorefxn(scorefxn)
	movemap = rosetta.core.kinematics.MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(False)
	for res in residues:
		movemap.set_bb(res, True)
		movemap.set_chi(res, True)
	relax.set_movemap(movemap)
	relax.apply(pose)


def dihedral_constraint(residues, pose):
	for res in residues:
		# AtomIDs for phi: C(i-1), N(i), CA(i), C(i)
		phi_atoms = (
   			AtomID(pose.residue(res - 1).atom_index("C"), res - 1),
    		AtomID(pose.residue(res).atom_index("N"), res),
    		AtomID(pose.residue(res).atom_index("CA"), res),
    		AtomID(pose.residue(res).atom_index("C"), res)
		)

		# AtomIDs for psi: N(i), CA(i), C(i), N(i+1)
		psi_atoms = (
    		AtomID(pose.residue(res).atom_index("N"), res),
    		AtomID(pose.residue(res).atom_index("CA"), res),
    		AtomID(pose.residue(res).atom_index("C"), res),
    		AtomID(pose.residue(res + 1).atom_index("N"), res + 1)
		)

		# Define beta strand target angles in radians
		phi_target = -2.36   # ~-135 degrees
		psi_target =  2.36   # ~135 degrees

		phi_func = CircularHarmonicFunc(phi_target, 0.3)
		psi_func = CircularHarmonicFunc(psi_target, 0.3)

		# Create the constraints
		phi_constraint = DihedralConstraint(*phi_atoms, phi_func)
		psi_constraint = DihedralConstraint(*psi_atoms, psi_func)

		pose.add_constraint(phi_constraint)
		pose.add_constraint(psi_constraint)


def sameAA_packer(pose):
	#Seems to work without an interactive node. Does not allow Amino Acids to change
	scorefxn = pyrosetta.create_score_function('ref2015')
	task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
	task.restrict_to_repacking()
	packer = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task)
	packer.apply(pose)


def changeAA_packer(pose):
	#Allows Amino Acids to change
	scorefxn = pyrosetta.create_score_function('ref2015')
	task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
	packer = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task)
	packer.apply(pose)


def slice_chain(pose, chain_id):
    chain_pose = rosetta.core.pose.Pose()
    residue_indices = rosetta.utility.vector1_unsigned_long()
    for i in range(1, pose.total_residue() + 1):
        if pose.pdb_info().chain(i) == chain_id:
            residue_indices.append(i)
    rosetta.core.pose.pdbslice(chain_pose, pose, residue_indices)
    return chain_pose


def h_bond_constraints(loop_residues, pose):
	hbond_pairs = []
	for i in loop_residues:
		for j in loop_residues:
			dist = pose.residue(i).xyz("O").distance(pose.residue(j).xyz("N"))
			if dist < 3.2:
				hbond_pairs.append((i,j))
				#print(f"Possible H-bond: A{i} O <--> B{j} N ({dist:.2f} Å)")
				#print(hbond_pairs)

	for res_i, res_j in hbond_pairs:
		try:
			atom_id_i = AtomID(pose.residue(res_i).atom_index("O"), res_i)
			atom_id_j = AtomID(pose.residue(res_j).atom_index("N"), res_j)
			func = HarmonicFunc(2.8, 0.3)
			constraint = AtomPairConstraint(atom_id_i, atom_id_j, func)
			pose.add_constraint(constraint)

		except Exception as e:
			print(f"Skipping {res_i}-{res_j}: {e}")


def select_residues(start_res, end_res, pose, pdb_length):
	dssp = rosetta.core.scoring.dssp.Dssp(pose)
	dssp.insert_ss_into_pose(pose)
	#selected_residues = [i for i in range(1, pose.total_residue() + 1)
	#						if pose.secstruct(i) == "E"]

	selected = []
	for i in range(9):
		chain = list(range(start_res + (pdb_length * i), end_res + (pdb_length * i) + 1))
		for res in chain:
			selected.append(res)
	selected = [i for i in selected if pose.secstruct(i) != 'H']
	print("Selected residues:")
	for i in selected:
		print(f"{i}: {pose.residue(i).name()}")

	return selected


def main():
	pyrosetta.init()
	#pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]
	#for i in range(len(pdb_files)):
	parser = argparse.ArgumentParser(description="Asymmetric constriction")
	parser.add_argument("-pdb", help="pdb to be modified")
	parser.add_argument("-output", help="output directory")
	parser.add_argument("-move_through", action='store_true')
	parser.add_argument("-fast", action ='store_true')
	args = parser.parse_args()
	pose = pyrosetta.pose_from_pdb(args.pdb)
	print("pose_loaded")

	chain_id = "A"
	pdb_length = sum(1 for i in range(1, pose.total_residue() + 1) if pose.pdb_info().chain(i) == chain_id)
	print(f"Monomer has {pdb_length} residues")

	selected_residues = select_residues(70, 90, pose, pdb_length)
	reversed_residues = list(reversed(selected_residues))

	num_selected = sum(1 for i in selected_residues if pose.pdb_info().chain(i) == chain_id)
	print(f'Selected {num_selected} residues per chain')

	for iteration in range(5000):
		if (iteration + 1) % 10 == 0:
			if args.move_through:
				selected_residues = selected_residues[num_selected:] + selected_residues[0:num_selected]
		pose.remove_constraints()
		if iteration % 2 == 0:
			#Residues in forward direction, odd iterations
			dihedral_constraint(selected_residues, pose)
			print("Dihedral Complete")
			h_bond_constraints(selected_residues, pose) 
			print("H Bond Complete")
			if args.fast:
				print("Fast Started")
				fast_relax(selected_residues, pose)
				print("Fast Relax Complete")
			else:
				print("Min started")
				min_torsion(selected_residues, pose)
				print("Minimization Complete")
		else:
			#Residues in reverse direction, even iterations
			#if iteration == 1:
				#print(reversed_residues)
			dihedral_constraint(reversed_residues, pose)
			print("Dihedral Complete, reversed")
			h_bond_constraints(reversed_residues, pose)
			print("H Bond Complete, reversed") 
			if args.fast:
				print("Fast Started")
				fast_relax(selected_residues, pose)
				print("Fast Relax Complete")
			else:
				print('Min started')
				min_torsion(reversed_residues, pose)
				print("Minimization Complete, reversed")

		changeAA_packer(pose)
		print(f"Packing Complete {iteration+1}")
			
		if (iteration+1) % 100 == 0:
			pose.dump_pdb(f"./output/{args.output}/changeAA_5000xMin_{iteration+1}.pdb")
			print(f"Iteration {iteration+1} Complete")
	print("Encourage Beta Sheet Complete")	


if __name__ == "__main__":
	main()