"""SpaceDock Launcher

Requirements : 
Python3.9+

This scripts allows to launch a SpaceDocking given a list of GOLD docking pose path 
or a folder containing GOLD docking poses.

Typical command line : 
python3.9 SpaceDock_launcher.py -o DRD3 -r "1" -pbb results_docking/

"""
import sys
sys.path.append('SpaceDock/')
import os
import argparse
import reagent as rct
import reaction as rnx

def select_reaction(args_):

	"""Function that return a list of selected reactions according to -r or --reaction argument
		
		Parameters
		----------
		args_ : str
			argument -r or --reaction

		Return
		------
		selected_reactions : list
			selected reactions e.g. if -r "1 10" ["Amide, Sulfonamide"]

	"""

	selected_reactions = []

	if len(args_) == 1 and args_[0] == "all":
		return ["Amide", "Sulfonamide", "Benzoxazole"]
	else:
		for arg_rnx in args_:
			if arg_rnx == "all":
				print("Bad arguments in -r or --reaction. 'all' is present but some reactions are choosen.")
				sys.exit(2)
			if arg_rnx == "1":
				selected_reactions.append("Amide")
			if arg_rnx == "10":
				selected_reactions.append("Sulfonamide")
			if arg_rnx == "27":
				selected_reactions.append("Benzoxazole")


		return selected_reactions

def calculate_packet(id_packet_, n_packets_, n_poses_):

	"""Function that calculate the number of reagent docking pose to annote THE FIRST REAGENT OF A REACTION 
		according to -ip and -np argument

		e.g. these is 100 amines and 100 acids for Amide reaction
		if -ip 1 and -np 1 (default case), start = 0 and end = 100
		bcause -ip (packet id) = 1 ajd -np (packet number) = 1
		All amine and acids will be annoted

		e.g. -ip 1 and -np 10
		There is 10 packets and we selected packet id 1
		start = 0 and end = 9
		10 amines (from 0 to 9) will be annoted, all acids will be annoted

		Parameters
		----------
		id_packet_ : int
			argument -ip, packet id
		n_packets_ : int
			argument -np, packet number
		n_packets_ : int
			docking pose number
				
		Return
		------
		start, end : int
			starting and end point to annote first reagent of a reaction

	"""

	lengh_packets = int(n_poses_ / n_packets_)

	start = (id_packet_ - 1) * lengh_packets

	if id_packet_ != n_packets_:
		end = id_packet_ * lengh_packets - 1
	else:
		end = n_poses_

	return start, end

def search_path_docking_poses(path_docking_):

	"""Function that search for GOLD docking poses into a given path folder

		Parameters
		----------
		path_docking_ : path containing mol2 docking poses or can be GOLD results path

		Return
		------
		start, end : int
			starting and end point to annote first reagent of a reaction

	"""

	list_path_docking_poses = []
	for (folder, sub_folder, files) in os.walk(path_docking_):
		for file in files:
			if file.find("gold_soln_") != -1:
				with open('{}/{}'.format(folder, file), 'r') as file_pose_hdl:
					mol2_str = file_pose_hdl.read()
					index_cut_begin = mol2_str.find("@<TRIPOS>MOLECULE")
					name_bb = mol2_str[index_cut_begin:].split("\n")[1].split("|")[0]
				list_path_docking_poses += ['{}\t{}/{}'.format(name_bb, folder, file)]
	return list_path_docking_poses



if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-lbb', '--listbbposes', type=str, required=False,
            help='GOLD docking pose path stored in tsv file : name_BB \\t path_BBs_pose')
	parser.add_argument('-pbb', '--pathbbposes', type=str, required=False,
            help='GOLD docking pose folder')
	parser.add_argument('-t', '--tag', type=str, required=False, default="data/tag_table_REAL.tsv",
            help='Building block tag table. Default datas/tag_table_REAL.tsv')
	parser.add_argument('-o','--output', type=str, required=True,
            help='Output name .mol2 and .tsv')
	parser.add_argument('-ip','--idpacket', type=int, required=False, default=1,
            help='Packet id. Default 1')
	parser.add_argument('-np','--npacket', type=int, required=False, default=1,
            help='Packet number. Default 1')
	parser.add_argument('-r', '--reaction', type=str, required=True, choices=["all", "1", "10", "27"], nargs="+",
            help='Reaction selection. all : every implented reactions. 1 : Amide, 10 : Sulfonamide, 27 : Benzoxazole.')

	### --- IChem --- ###
	parser.add_argument('--ichem', type=str, required=False,
            help='enable IChem IFP filtering. \n ichem.conf is needed')

	args = parser.parse_args()

	selected_reactions = select_reaction(args.reaction)
	print("Selected reactions : ", selected_reactions, "\n")

	if args.ichem != None:
		ichem_path = None
		thresold_FULL = None
		thresold_POLAR = None
		ligand = None
		protein = None
		with open(args.ichem, 'r') as file_ichem_conf_hdl:
			lines_ichem_conf = file_ichem_conf_hdl.read().split("\n")
			for line in lines_ichem_conf:
				if line != "":
					if line.strip()[0] != "#":
						if line.find("ICHEM_PATH") != -1:
							ichem_path = line.split()[1]
						if line.find("THRESOLD_FULL")  != -1:
							thresold_FULL = float(line.split()[1])
						if line.find("THRESOLD_POLAR")  != -1:
							thresold_POLAR = float(line.split()[1])
						if line.find("LIGAND")  != -1:
							ligand = line.split()[1]
						if line.find("PROTEIN")  != -1:
							protein = line.split()[1]
		ichem_conf = (ichem_path, thresold_FULL, thresold_POLAR, ligand, protein)
		for parameter in ichem_conf:
			if parameter == None:
				print("Missing parameters in ichem conf file")
				print("Please check the ichem conf file")
				sys.exit(2)

		print("IChem Parameters : ")
		print("PATH", ichem_path)
		print("Thresold", "FULL :", thresold_FULL, "POLAR :", thresold_POLAR)
		print("Rerefence ligand", ligand)
		print("Protein", protein)
		print("\n")
	else:
		ichem_conf = None

	path_tag_table = args.tag
	dico_tag_table = rct.Loader.load_tag_table(path_tag_table, selected_reactions)

	if args.listbbposes != None and args.pathbbposes == None:
		with open(args.bbposes, 'r') as file_poses_hdl:
			list_path_docking_poses = file_poses_hdl.read().split("\n")
		if len(list_path_docking_poses) == 0:
			print("No docking poses founded in", args.bbposes)
			sys.exit(2)
	elif args.listbbposes == None and args.pathbbposes != None:
		list_path_docking_poses = search_path_docking_poses(args.pathbbposes)
		if len(list_path_docking_poses) == 0:
			print("No docking poses founded in", args.pathbbposes)
			sys.exit(2)
	else:
		print("Wrong arguments, list or a path for building block docking poses")
		sys.exit(2)



	start, end = calculate_packet(args.idpacket, args.npacket, len(list_path_docking_poses))

	loaded_poses = rct.Loader.load_poses(list_path_docking_poses, dico_tag_table, selected_reactions, start, end)
	
	name_output = args.output

	print("\n")
	for reaction in selected_reactions:
		if reaction == "Amide":
			amines, carboxylic_acids = loaded_poses["Amide"]
			if amines != [] and carboxylic_acids != []:
				print("Starting", reaction, "reaction...")
				rnx.Synthesis.amide(amines, carboxylic_acids, name_output, ichem_conf)
				print(reaction, "done \n")
			else:
				print(reaction, "is skipped due to lack of BBs capable to do that reaction \n")

		elif reaction == "Sulfonamide":
			amines, sulfonyl_chloride = loaded_poses["Sulfonamide"]
			if amines != [] and sulfonyl_chloride != []:
				print("Starting", reaction, "reaction...")
				rnx.Synthesis.sulfonamide(amines, sulfonyl_chloride, name_output, ichem_conf)
				print(reaction, "done \n")
			else:
				print(reaction, "is skipped due to lack of BBs capable to do that reaction \n")

		elif reaction == "Benzoxazole":
			aminophenols, benzaldehydes = loaded_poses["Benzoxazole"]
			if aminophenols != [] and benzaldehydes != []:
				print("Starting", reaction, "reaction...")
				rnx.Synthesis.benzoxazole(aminophenols, benzaldehydes, name_output, ichem_conf)
				print(reaction, "done \n")
			else:
				print(reaction, "is skipped due to lack of BBs capable to do that reaction \n")















