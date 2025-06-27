#refine beta barrel
#Use /net/software/containers/pyrosetta.sif
#ls -1 | wc -1 *pdb
#qlogin --mem=24g -c 1 -p cpu

import os
import sys
import glob
import pyrosetta
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
	scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint, 1.0)
	scorefxn.set_weight(rosetta.core.scoring.dihedral_constraint, 0.5)

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


def packer(pose):
	#Requires the use of an interactive node. Uses a lot of memory. Do not use this version.
	tf = rosetta.core.pack.task.TaskFactory()
	tf.push_back(rosetta.core.pack.task.operation.InitializeFromCommandline())
	packer = rosetta.protocols.minimization_packing.PackRotamersMover()
	packer.task_factory(tf)
	packer.score_function(pyrosetta.get_fa_scorefxn())
	packer.apply(pose)
	print("Packing done")


def new_packer(pose):
	#Seems to work without an interactive node.
	scorefxn = pyrosetta.create_score_function('ref2015')
	task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
	task.restrict_to_repacking()
	packer = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task)
	packer.apply(pose)


def initialize_pyrosetta(symdef_file):
    pyrosetta.init(extra_options=f"-ex1 -ex2aro -symmetry:symmetry_definition {symdef_file}")


def get_symdef_file(symmetry_type):
    symdef_dir = "./symm"
    symdef_file = f'{symmetry_type}.sym'
    for file in os.listdir(symdef_dir):
        if symmetry_type in file:
            symdef_file = os.path.join(symdef_dir, file)
            break
    if symdef_file is None:
        print(f"Error: No symmetry definition file found for {symmetry_type}")
        sys.exit(1)
    return symdef_file


def slice_chain(pose, chain_id):
    chain_pose = rosetta.core.pose.Pose()
    residue_indices = rosetta.utility.vector1_unsigned_long()
    for i in range(1, pose.total_residue() + 1):
        if pose.pdb_info().chain(i) == chain_id:
            residue_indices.append(i)
    rosetta.core.pose.pdbslice(chain_pose, pose, residue_indices)
    return chain_pose


def apply_symmetry(pdb_path, symmetry_type, symmetrized_output_dir, repack):
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    symmetrized_output_path = os.path.join(symmetrized_output_dir, f"sym_{pdb_name}_{symmetry_type}.pdb")

    # Get the appropriate symmetry definition file
    symdef_file = get_symdef_file(symmetry_type)

    # Initialize PyRosetta with the specified symmetry definition file
    initialize_pyrosetta(symdef_file)

    # Load the input structure in PyRosetta
    pose = pyrosetta.pose_from_file(pdb_path)

    # Ensure only chain A is used
    chain_a_pose = slice_chain(pose, "A")

    # Apply symmetry setup
    rosetta.protocols.symmetry.SetupForSymmetryMover(symdef_file).apply(chain_a_pose)

    if repack:
        # Repack residues
        scorefxn = pyrosetta.create_score_function('ref2015')
        task = rosetta.core.pack.task.TaskFactory.create_packer_task(chain_a_pose)
        task.restrict_to_repacking()
        packer = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task)
        packer.apply(chain_a_pose)

    # Save the symmetric pose
    chain_a_pose.dump_pdb(symmetrized_output_path)

    return symmetrized_output_path


def get_sym_pdb(pdb_files, index):
	pdb_path = f'./{str(pdb_files[index])}'
	symmetry_type = 'C9'
	symmetrized_output_dir = './sym_pdb/'
	repack = True 
	sym_output_file = apply_symmetry(pdb_path, symmetry_type, symmetrized_output_dir, repack)

	return sym_output_file


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
	selected = [i for i in selected if pose.secstruct(i) != 'H' or 'E']
	print("Selected residues:")
	for i in selected:
		print(f"{i}: {pose.residue(i).name()}")

	#selected_chain_residues = {}
	#pose_pdbinfo = pose.pdb_info()
	#for i in range(1, pose.total_residue() + 1):
	#	chain = pose_pdbinfo.chain(i)
	#	resnum = pose_pdbinfo.number(i)
	#	if resnum in selected:
	#		selected_chain_residues.setdefault(chain, []).append(i)

	return selected


def main():
	#pyrosetta.init()

	pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]
	for i in range(len(pdb_files)):
		sym_pdb = get_sym_pdb(pdb_files, i)
		pose = pyrosetta.pose_from_file(sym_pdb)
		print("pose_loaded")

		chain_id = "A"
		pdb_length = sum(1 for i in range(1, pose.total_residue() + 1) if pose.pdb_info().chain(i) == chain_id)
		print(f"Monomer has {pdb_length} residues")
		original_selected_residues = select_residues(70, 90, pose, pdb_length)
		original_reversed_residues = list(reversed(original_selected_residues))

		for iteration in range(100):
			pose.remove_constraints()
			selected_residues = select_residues(70, 90, pose, pdb_length)
			reversed_residues = list(reversed(selected_residues))
			if iteration % 2 == 0:
				#Residues in forward direction, odd iterations
				dihedral_constraint(selected_residues, pose)
				print("Dihedral Complete")
				h_bond_constraints(original_selected_residues, pose) #allow hydrogen bonding to residues already a Beta Sheet
				print("H Bond Complete")
				min_torsion(selected_residues, pose)
				print("Minimization Complete")
			else:
				#Residues in reverse direction, even iterations
				if iteration == 1:
					print(reversed_residues)
				dihedral_constraint(reversed_residues, pose)
				print("Dihedral Complete, reversed")
				h_bond_constraints(original_reversed_residues, pose) #allow hydrogen bonding to residues already a Beta Sheet
				print("H Bond Complete, reversed") 
				min_torsion(reversed_residues, pose)
				print("Minimization Complete, reversed")

			#fast_relax(selected_residues, pose)
			#print("Fast Relax Complete")
			new_packer(pose)
			print(f"Packing Complete {iteration+1}")
			
			if (iteration+1) % 50 == 0:
				pose.dump_pdb(f"Update_residues_100xMin_{iteration+1}.pdb")
				print(f"Iteration {iteration+1} Complete")
		print("Encourage Beta Sheet Complete")	


if __name__ == "__main__":
	main()