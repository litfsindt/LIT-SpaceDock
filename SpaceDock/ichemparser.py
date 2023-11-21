"""IChem Parser for SpaceDock

"""


import subprocess


def gen_Tanimoto_IFP(ichem_path_, site_, pose_, ligand_, polar_ = False):

	"""Call IChem to compute IFPs and tanimoto between a pose and a reference ligand

		Parameters
		----------
		ichem_path_ : str
			IChem path executable file
		site_ : str
			protein or site path as mol2
		pose_ : str
			docking pose or molecule path as mol2
		ligand_ : str
			reference ligand path as mol2
		polar_ : boolean
			Compute or not polar IFP (Default False)
		
		Return
		------
		output_IChem : str
			IChem IFP output as string (containing bitstring and tanimoto between pose and ref)
	"""

	if polar_ == False:
		command = subprocess.run([ichem_path_, "IFP", site_, pose_, ligand_], stdout = subprocess.PIPE)
	else:
		command = subprocess.run([ichem_path_, "--polar", "IFP", site_, pose_, ligand_], stdout = subprocess.PIPE)
	output_IChem = command.stdout.decode()

	return output_IChem

def get_TC_IFP(ichem_path_, site_, pose_, ligand_, seuil_tc_full_, seuil_tc_polar_):

	"""Extract tanimoto coefficient from IChem IFP output.
		It filtered the poses with tanimoto coefficient according seuil_tc_full_ and seuil_tc_polar_
		
		For SpaceDock "pose_" is TEMP_ICHEM.mol2
		recombination that 

		Parameters
		----------
		ichem_path_ : str
			IChem path executable file
		site_ : str
			protein or site path as mol2
		pose_ : str
			docking pose or molecule path as mol2
		ligand_ : str
			reference ligand path as mol2
		seuil_tc_full_ : float
			Thresold to discard poses than have tanimoto coeff on FULL IFP less than seuil_tc_full_
		seuil_tc_polar_
			Thresold to discard poses than have tanimoto coeff on POLAR IFP less than seuil_tc_polar_
		
		Return
		------
		filtered_poses : dict
			poses that passed both full and polar filters
			key : molecule name and values : tc_full and tc_polar
	"""

	### This can be optimised by extracting polar bit from IFP FULL ###
	### Has to be done ###
	line_output_IChem_IFP = gen_Tanimoto_IFP(ichem_path_, site_, pose_, ligand_).split("\n")
	line_output_IChem_IFP_polar = gen_Tanimoto_IFP(ichem_path_, site_, pose_, ligand_, polar_ = True).split("\n")

	tcs_poses = {}
	for line in line_output_IChem_IFP:
		splited_line = line.split("\t")
		
		### --- This IF is to avoid extraction the bitstring --- ###
		### We only want to extract line that contains Tanimoto Coefficient ###
		if len(splited_line) == 3 and splited_line[0] != splited_line[1] and (splited_line[2].find("0") != -1 or splited_line[2].find("1") != -1):
			nom_ref, pose, tc = line.split("\t")
			tcs_poses[pose] = tc

	for line in line_output_IChem_IFP_polar:
		splited_line = line.split("\t")

		### --- This IF is to avoid extraction the bitstring --- ###
		### We only want to extract line that contains Tanimoto Coefficient ###
		if len(splited_line) == 3 and splited_line[0] != splited_line[1] and (splited_line[2].find("0") != -1 or splited_line[2].find("1") != -1):
			nom_ref, pose, tc_polar = line.split("\t")
			tcs_poses[pose] = [tcs_poses[pose], tc_polar]	

	filtered_poses = {}

	### --- Pose filtering after tanimoto coefficient extraction ---###

	for pose, tcs in tcs_poses.items():
		tc, tc_polar = tcs
		if float(tc) >= seuil_tc_full_ and float(tc_polar) >= seuil_tc_polar_:
			filtered_poses[pose] = [tc, tc_polar]

	return filtered_poses

def write_temp_file(dico_mol2_):

	"""Write TEMPORARY molecule as mol2 from dico_mol2_ who contains molecule loaded as Mol2 object with mol2parser.py
		These molecules will be submitted for IFP calculation and filetered

		Molecule are saved with TEMPORARY name. Because IChem IFP limits the number of char in the name.
		It's annoying so we have to do an index to match the temporary name and the real name of the molecule

		Parameters
		----------
		dico_mol2_ : dict
			dictionary containing as key name of the molecule and as values mol2 object and others outpout informations
		
		Return
		------
		dico_index : dict
			dictionary to match the temporary name and the true name
			
	"""

	dico_index = {}
	compteur = 1

	with open("TEMP_ICHEM.mol2", 'w') as file_temp_mol2:
		for nom_recombinaison, mol2_and_values in dico_mol2_.items():
			mol2, values = mol2_and_values
			temp_name = 'TEMP_{}'.format(compteur)
			dico_index[temp_name] = nom_recombinaison
			mol2.molecule.set_name(temp_name)
			file_temp_mol2.write('{}\n'.format(mol2.build()))
			compteur += 1

	return dico_index