import math
import os
import mol2parser as m2r
import ichemparser as IPS

class Reaction:

	"""Parent class that defines general topological rules to connect two docking poses of chemical reagent
	"""

	def __init__(self):
		self.distance_between_center_mass_down = 5
		self.distance_between_center_mass_up = 7.5
		self.distance_between_connectable_atoms_down = 1 
		self.distance_between_connectable_atoms_up = 3
		self.clash_max = 4

class Amide(Reaction):

	"""Class that defines local topological rules for Amide reaction
	"""

	def __init__(self):
		Reaction.__init__(self)
		self.angle_amide_bond_down = 55
		self.angle_amide_bond_up = 155
		self.sum_angle_amide_bond_down = 145
		self.sum_angle_amide_bond_up = 265
		self.dihedral_amide_bond = 120

class Sulfonamide(Reaction):

	"""Class that defines local topological rules for Sulfonamide reaction
	"""

	def __init__(self):
		Reaction.__init__(self)

class Benzoxazole(Reaction):

	"""Class that defines local topological rules for Benzoxazole reaction
	"""

	def __init__(self):
		Reaction.__init__(self)
		self.angle_cycle_down = 55
		self.angle_cycle_up = 155


class Topological:

	"""Class used to compute the topological filters for each SpaceDock reaction

		Methods
		-------
		angle(coord_1_, coord_2_, coord_3_)
			return an angle (in degree) computed with law of cosines method

		dihedral_angle(coord_1_, coord_2_, coord_3_, coord_4_)
			return a dihedral angle (in degree) computed with law of cosines method

		check_clashes(coords_reactant1_, coords_reactant2_, reaction_, thresold_distance_ = 2.5)
			return True or False if the number of allowed clashes is exceded

		IFP(dico_merged_mol2_, ligand_, protein_)
			return recombination filtered by IFPs

	"""

	def angle(coord_1_, coord_2_, coord_3_):

		"""compute an angle (in degree) computed with law of cosines method

			  B	
			 /-\
			A   C

			Parameters
			----------
			coord_1_ : list
				3D coordinates of atom A eg [1.0, 2.0, 1.0]
			coord_2_ : list
				3D coordinates of atom B eg [1.0, 2.0, 1.0]
			coord_3_ : list
				3D coordinates of atom C eg [1.0, 2.0, 1.0]

			returns
			-------
			float
				angle in degree


		"""

		a = math.dist(coord_2_, coord_1_)
		b = math.dist(coord_2_, coord_3_)
		c = math.dist(coord_1_, coord_3_)

		loc = (-(c*c) + (a*a) + (b*b)) / (2*a*b)
		if loc < 1 and loc > -1:
			angle = math.acos(loc)
		else:
			angle = 360

		return angle * 180/3.14

	def dihedral_angle(coord_1_, coord_2_, coord_3_, coord_4_):
		
		"""compute a dihedral angle in degree
			
			  B	  D
			 /-\-/
			A   C
			
			Parameters
			----------
			coord_1_ : list
				3D coordinates of atom A eg [1.0, 2.0, 1.0]
			coord_2_ : list
				3D coordinates of atom B eg [1.0, 2.0, 1.0]
			coord_3_ : list
				3D coordinates of atom C eg [1.0, 2.0, 1.0]
			coord_4_ : list
				3D coordinates of atom D eg [1.0, 2.0, 1.0]

			returns
			-------
			float
				dihedral angle in degree

		"""

		def norme(vecteur_):
			"""
				return the norm of a vector
				
				Parameters
				----------
				vecteur_  : list
					3D coordinates of vector

			"""

			return math.dist(vecteur_, (0,0,0))

		def scalar_product(vecteur_1_, vecteur_2_):
			"""
				return the scalar product of two vectors
				
				Parameters
				----------
				vecteur_1_  : list
					3D coordinates of the 1st vector
				vecteur_2_  : list
					3D coordinates of the 2nd vector

			"""

			x_vect1, y_vect1, z_vect1 = vecteur_1_
			x_vect2, y_vect2, z_vect2 = vecteur_2_
			return (x_vect1 * x_vect2 + y_vect1 * y_vect2 + z_vect1 * z_vect2)

		def vectorial_product(vecteur_1_, vecteur_2_):
			"""
				return the vectorial product of two vectors
				
				Parameters
				----------
				vecteur_1_  : list
					3D coordinates of the 1st vector
				vecteur_2_  : list
					3D coordinates of the 2nd vector

			"""
			x_vect1, y_vect1, z_vect1 = vecteur_1_
			x_vect2, y_vect2, z_vect2 = vecteur_2_

			return (y_vect1 * z_vect2 - z_vect1 * y_vect2,
					z_vect1 * x_vect2 - x_vect1 * z_vect2,
					x_vect1 * y_vect2 - y_vect1 * x_vect2)

		x_coord1, y_coord1, z_coord1 = coord_1_
		x_coord2, y_coord2, z_coord2 = coord_2_
		x_coord3, y_coord3, z_coord3 = coord_3_
		x_coord4, y_coord4, z_coord4 = coord_4_

		vecteur_12 = (float(x_coord1) - float(x_coord2), float(y_coord1) - float(y_coord2), float(z_coord1) - float(z_coord2))
		vecteur_23 = (float(x_coord2) - float(x_coord3), float(y_coord2) - float(y_coord3), float(z_coord2) - float(z_coord3))
		vecteur_34 = (float(x_coord3) - float(x_coord4), float(y_coord3) - float(y_coord4), float(z_coord3) - float(z_coord4))
	 
		vecteur_normal1 = vectorial_product(vecteur_12, vecteur_23)
		vecteur_normal2 = vectorial_product(vecteur_23, vecteur_34)
	 
		if scalar_product(vecteur_23 , vectorial_product(vecteur_normal1, vecteur_normal2)) < 0:
			
			try:
				calc = round(scalar_product(vecteur_normal1,vecteur_normal2)/(norme(vecteur_normal1)*norme(vecteur_normal2)), 2)
				dihedral = math.acos(calc)*180/math.pi
			except:
				dihedral = 0
			return dihedral
		else:
			try:
				calc = round(scalar_product(vecteur_normal1,vecteur_normal2)/(norme(vecteur_normal1)*norme(vecteur_normal2)), 2)
				dihedral = math.acos(calc)*180/math.pi
			except:
				dihedral = 0
			return dihedral

	def check_clashes(coords_reactant1_, coords_reactant2_, max_clash_, thresold_distance_ = 2.5):


		"""return True or False if the number of allowed clashes is exceded
			clash : Distance between two atoms < 2.5 (default value)
			number of allowed clash : defined in the class Reaction (by default 4)

			Parameters
			----------
			coords_reactant1_ : list
				list of list of 3D coordinates of reactant 1 [[1.0, 2.0, 1.0], [1.0, 3.0, 1.0], ...]
			coords_reactant2_ : list
				list of list of 3D coordinates of reactant 2 [[1.0, 2.0, 1.0], [1.0, 3.0, 1.0], ...]
			max_clash_ : <Reaction>
				number of allowed clash (defined inside of <Reaction> object)
			thresold_distance_ : float
				thresold that define if there is a clash of not (default 2.5)
			
			returns
			-------
			bolean
				True or False if the variable clashes is upper the limit
		"""

		clashes = 0	
		
		for coord_reactant1 in coords_reactant1_:
			x1, y1, z1 = coord_reactant1
			for coord_reactant2 in coords_reactant2_:
				x2, y2, z2 = coord_reactant2

				dist_clash = math.dist([x1, y1, z1], [x2, y2, z2])

				if dist_clash < thresold_distance_:
					clashes += 1

				if clashes > max_clash_:
					return clashes > max_clash_

		return clashes > max_clash_

	def IFP(dico_merged_mol2_, ichem_conf_):

		"""Compute IFPs calculation with IChem and return only SpaceDock recombination that

			Parameters
			----------
			dico_merged_mol2_ : dict
				dictionary of SpaceDock recombination
			ichem_conf_ : ichem_conf file given as SpaceDock_launcher parameter
				conf file containing IChem path exe, ligand path, ...
			
			returns
			-------
			dico_merged_mol2_filtered
				dictionary of SpaceDock recombination filtered by IFPs
		"""

		ichem_path, thresold_FULL, thresold_POLAR, ligand, protein = ichem_conf_

		dico_merged_mol2_filtered = {}
		
		dico_index = IPS.write_temp_file(dico_merged_mol2_)
		dico_name_merged_mol2_filtered = IPS.get_TC_IFP(ichem_path, protein, "TEMP_ICHEM.mol2", ligand, thresold_FULL, thresold_POLAR)
		os.system('rm TEMP_ICHEM.mol2')
		
		for name_temp, tcs in dico_name_merged_mol2_filtered.items():
			name_merged = dico_index[name_temp]
			tc, tc_polar = tcs
			mol2, values = dico_merged_mol2_[name_merged]
			values.append(tc)
			values.append(tc_polar)
			mol2.molecule.set_name(name_merged)
			dico_merged_mol2_filtered[name_merged] = [mol2, values]

		return dico_merged_mol2_filtered

class Output:

	"""Class that contains fonction to generate SpaceDock outputs

		Methods
		-------
		name_and_num_pose_from_path(path_)
			split the path of a docking pose to retrive the file name and id
		write_output(dico_merged_mol2_filtered_, name_output_, num_reaction_)
			Write outputs files (mol2 and tsv)
	"""

	def name_and_num_pose_from_path(path_):
		"""Split the path of a docking pose to retrive the file name and id 
			
			It is use to generate a unique name for SpaceDock poses

			Parameters
			----------
			path_ : str
				docking pose path
				
			Return
			------
			tuple
				file name and id of a pose
		"""

		path_splited = path_.split("/")

		name_pose = path_splited[len(path_splited) - 1]
		name_splited = name_pose.split("_")

		num_pose = name_splited[len(name_splited) - 1][:-5]

		return name_pose, num_pose

	def write_output(dico_merged_mol2_filtered_, name_output_, num_reaction_):

		"""Write outputs
			
			Give the same prefix for .tsv and .mol2 file

			Parameters
			----------
			name_output_ : str
				name of the output for the .mol2 and .tsv files
			num_reaction_ : str
				number of the reaction to seperate .mol2 and .tsv for each reaction.
				In case multiple reaction are selected
				
		"""

		with open('{}_{}.mol2'.format(name_output_, num_reaction_), 'a') as file_output_mol2_hdl, open('{}_{}.tsv'.format(name_output_, num_reaction_), 'a') as file_output_values_hdl:
			for name_merged, outputs in dico_merged_mol2_filtered_.items():
				mol2, values = outputs
				file_output_mol2_hdl.write(mol2.build())

				file_output_values_hdl.write('{}'.format(name_merged))
				for value in values:
					file_output_values_hdl.write('\t{}'.format(value))
				file_output_values_hdl.write('\n')

class Synthesis:

	"""Class that contains fonction to apply topological rules for each reaction

		Methods
		-------
		amide(amines_, carboxylic_acids_, ligand_, protein_, name_output_)
			Amide workflow : cmass -> connectable atoms -> angles amide bond -> clashes

		sulfonamide(amines_, sulfonyl_chlorides_, ligand_, protein_, name_output_)
			Sulfonamide workflow : cmass -> connectable atoms -> clashes

		benzoxazole(aminophenols_, benzaldehydes_, ligand_, protein_, name_output_)
			Benzoxazole workflow : cmass -> connectable atoms -> angles cycle -> clashes

		Variables
		---------
		thresold_write : int
			When the number of recombination reach this number, SpaceDock poses are writed as .mol2


	"""

	thresold_write = 100000

	def amide(amines_, carboxylic_acids_, name_output_, ichem_conf_):

		"""Amide workflow : cmass -> connectable atoms -> angles amide bond -> clashes 	

			Series of IF statement
			
			topological_value = SOMETHING
			if topological_value > "or" < THRESOLD:
				new_topological_value = SOMETHING
				if ...
					...

			IChem can be activated to filter out combination by IFPs similarity
			ichem.conf has to be filled

			Parameters
			----------
			amines_ : list
				list of building block amine docking poses (object : Reactant::Amine)
			carboxylic_acids_ : list
				list of building block carboxylic acid docking poses (object : Reactant::CarboxylicAcid)
			name_output_ : str
				.tsv and .mol2 output name
			ichem_ : tuple 
				By default None is any ichem.conf file is given as paramter.
				tuple containing (path_ichem, thresold_FULL, thresold_POLAR, ligand and protein as mol2)

		"""

		dico_merged_mol2 = {} # This dictionary will contains SpaceDock poses (Two merged building block)

		reaction = Amide() # Setting the rules for Amide

		for amine in amines_:
			for carboxylic_acid in carboxylic_acids_:

				try:
					### --- 1st filter centre of mass --- ###
					distance_between_center_mass = math.dist(amine.coord_center_of_mass, carboxylic_acid.coord_center_of_mass)
					if reaction.distance_between_center_mass_down <= distance_between_center_mass <= reaction.distance_between_center_mass_up:

						### --- 2nd filter distance between connectable atoms "C.2-N.am" --- ###
						distance_between_N_C2 = math.dist(amine.coord_N, carboxylic_acid.coord_C2)
						if reaction.distance_between_connectable_atoms_down <= distance_between_N_C2 <= reaction.distance_between_connectable_atoms_up:

							### --- 3rd filter amide bond angles "C-N.am-C.2 and C-C.2-N.am" --- ###
							### * /!\ amine.coord_neighbors_N is a list, in case of SECONDARY amine ###
							### It as not one neighbor (e.g. cyclic amine)
							### We decide to only compute one angle in case of SECONDARY amine ###
							### We assume if one is very bad (e.g. 20 deg.), the other will be bad too ###

							### * Amide bond geometry does not allow us to be permissive here ###
							### That is why we add a sum filter of the two computed angles ###
							### To avoid strained case that will be discarded during energy minimisation step ###
							### e.g. two low angles 55 deg. will be discarded ###
							### e.g. one angle 55 deg. and another one 120 deg are fine ###

							angle_neighborN_N_C2 = Topological.angle(amine.coord_neighbors_N[0], amine.coord_N, carboxylic_acid.coord_C2)
							angle_neighborC2_C2_N = Topological.angle(carboxylic_acid.coord_neighbor_C2, carboxylic_acid.coord_C2, amine.coord_N)

							if reaction.angle_amide_bond_down <= angle_neighborN_N_C2 <= reaction.angle_amide_bond_up\
								and reaction.angle_amide_bond_down <= angle_neighborC2_C2_N <= reaction.angle_amide_bond_up\
								and reaction.sum_angle_amide_bond_down <= angle_neighborN_N_C2 + angle_neighborC2_C2_N <= reaction.sum_angle_amide_bond_up:

								### --- 4th filter clashes --- ###
								if Topological.check_clashes(amine.coord_atoms, carboxylic_acid.coord_atoms, reaction.clash_max) == False:

									kept_H = None # Will constain amine H atom ids that will will be not removed
									kept_Oco2 = None # Will constain carboxylic acid O.co2 atom ids that will be not removed
									best_dihedral_Oco2_C2_Nam_H = 0

									### --- Not a filter : SECONDARY and PRIMARY amine case --- ###
									### It is used to save the atoms that show the best dihedral angle OCNH ###
									### It helps the minimizor to not start from very strained position ###

									### /!\ Secondary : There is no H to save, dihedral is not calculated ###
									### The only one present will be removed to create the C-N bond ###
									### We only go threw carboxylic oxygens to keep the only who does not point toward the future C-N bond ###

									if amine.type_amine == "SECONDARY":
										best_dihedral_Oco2_C2_Nam_H = 180 # is set to something different than 0, necessary to pass next IF statement
										best_angle_Oco2_C2_Nam = 0
										
										for id_Oco2, coord_Oco2 in carboxylic_acid.ids_and_coords_Oco2:

											angle_Oco2_C2_Nam = Topological.angle(carboxylic_acid.coord_C2, coord_Oco2, amine.coord_N)
											if angle_Oco2_C2_Nam > best_angle_Oco2_C2_Nam:
												best_angle_Oco2_C2_Nam = angle_Oco2_C2_Nam
												kept_Oco2 = id_Oco2 

									elif amine.type_amine == "PRIMARY":
										for id_H, coord_H in amine.ids_and_coords_H:

											for id_Oco2, coord_Oco2 in carboxylic_acid.ids_and_coords_Oco2:

												angle_H_Nam_C2 = Topological.angle(coord_H, amine.coord_N, carboxylic_acid.coord_C2) # Get best dihedral angle to avoid H pointing toward amide bond
												angle_Oco2_C2_Nam = Topological.angle(coord_Oco2, carboxylic_acid.coord_C2, amine.coord_N) # Get best dihedral angle to avoid O pointing toward amide bond

												if angle_H_Nam_C2 >= 45 and angle_Oco2_C2_Nam >= 45:
												
													dihedral_angle_Oco2_C2_Nam_H = abs(Topological.dihedral_angle(coord_H, amine.coord_N, carboxylic_acid.coord_C2, coord_Oco2))
													if dihedral_angle_Oco2_C2_Nam_H >= reaction.dihedral_amide_bond and dihedral_angle_Oco2_C2_Nam_H > best_dihedral_Oco2_C2_Nam_H:
														best_dihedral_Oco2_C2_Nam_H = dihedral_angle_Oco2_C2_Nam_H
														kept_H = id_H 
														kept_Oco2 = id_Oco2 

									### --- Not a filter : Just a check to be sure that the dihedral angle has been correctly computed --- ###

									if best_dihedral_Oco2_C2_Nam_H != 0:

										### --- Congratz ! --- ###
										### If it reaches this point, it means that two chemical reagent docking poses can be combined ###

										### --- Important informations for the merger --- ###

										atoms_to_delete_amine = [id_H for id_H, coord_H in amine.ids_and_coords_H if id_H != kept_H]
										atoms_to_modify_amine = {amine.id_N : "N.am"}
										bonds_to_modify_amine = {}

										atoms_to_delete_carboxylic_acid = [id_Oco2 for id_Oco2, coord_Oco2 in carboxylic_acid.ids_and_coords_Oco2 if id_Oco2 != kept_Oco2]
										atoms_to_modify_carboxylic_acid = {kept_Oco2 : "O.2"}
										bonds_to_modify_carboxylic_acid = {(carboxylic_acid.id_C2, kept_Oco2) : "2"}

										connexion_rules = [(("1", amine.id_N), ("2", carboxylic_acid.id_C2), "am")]

										dico_information_mol2 = {"1" : [m2r.Loader.load_mol2(amine.path), atoms_to_delete_amine, atoms_to_modify_amine, bonds_to_modify_amine], 
																	"2" : [m2r.Loader.load_mol2(carboxylic_acid.path), atoms_to_delete_carboxylic_acid, atoms_to_modify_carboxylic_acid, bonds_to_modify_carboxylic_acid]}

										### --- Output side informations --- ###

										name_pose_amide, num_pose_amide = Output.name_and_num_pose_from_path(amine.path)
										name_pose_carboxylic_acid, num_pose_carboxylic_acid = Output.name_and_num_pose_from_path(carboxylic_acid.path)

										name_merged_ = '{}_{}_{}_{}_1'.format(amine.name, carboxylic_acid.name, num_pose_amide, num_pose_carboxylic_acid)
										tsv_output = [name_pose_amide, name_pose_carboxylic_acid]

										### --- mol2 merging is launched --- ###

										dico_merged_mol2[name_merged_] = (m2r.Merger.merge_mol2(dico_information_mol2, connexion_rules, name_merged_), tsv_output)

				except:
					print("Error during recombination for Amide")

			### --- Outputs --- ###

			if len(dico_merged_mol2) > Synthesis.thresold_write:

				if ichem_conf_ != None:
					dico_merged_mol2_IFP_filtered = Topological.IFP(dico_merged_mol2, ichem_conf_)
					Output.write_output(dico_merged_mol2_IFP_filtered, name_output_, "1")
				else:
					Output.write_output(dico_merged_mol2, name_output_, "1")
				dico_merged_mol2 = {}

		if ichem_conf_ != None:
			dico_merged_mol2_IFP_filtered = Topological.IFP(dico_merged_mol2, ichem_conf_)
			Output.write_output(dico_merged_mol2_IFP_filtered, name_output_, "1")
		else:
			Output.write_output(dico_merged_mol2, name_output_, "1")

	def sulfonamide(amines_, sulfonyl_chlorides_, name_output_, ichem_conf_):

		"""Sulfonamide workflow : cmass -> connectable atoms -> clashes

			It is series of IF statement
			
			topological_value = SOMETHING
			if topological_value > "or" < THRESOLD:
				new_topological_value = SOMETHING
				if ...
					...

			IChem can be activated to filter out combination by IFPs similarity
			ichem.conf has to be filled

			Parameters
			----------
			amines_ : list
				list of building block amine docking poses (object : Reactant::Amine)
			sulfonyl_chlorides_ : list
				list of building block sulfonyl chloride docking poses (object : Reactant::SulfonylChloride)
			name_output_ : str
				.tsv and .mol2 output name
			ichem_ : tuple 
				By default None is any ichem.conf file is given as paramter.
				tuple containing (path_ichem, thresold_FULL, thresold_POLAR, ligand and protein as mol2)
		"""

		dico_merged_mol2 = {}

		reaction = Sulfonamide() # Setting the rules for Sulfonamide

		for amine in amines_:
			for sulfonyl_chloride in sulfonyl_chlorides_:

				try:
					### --- 1st filter centre of mass --- ###
					distance_between_center_mass = math.dist(amine.coord_center_of_mass, sulfonyl_chloride.coord_center_of_mass)
					if reaction.distance_between_center_mass_down <= distance_between_center_mass <= reaction.distance_between_center_mass_up:
						
						### --- 2nd filter distance between connectable atoms "S-N" --- ###
						### For sulfonam. there is not angle check due its particular geometry ###

						distance_between_S_N = math.dist(amine.coord_N, sulfonyl_chloride.coord_S)
						if reaction.distance_between_connectable_atoms_down <= distance_between_S_N <= reaction.distance_between_connectable_atoms_up:

							### --- 3rd filter clashes --- ###
							### This filter is important because there is no angle check ###
							### It avoids mainly overlaps with sulfonyle oxygens O=S=O ###

							if Topological.check_clashes(amine.coord_atoms, sulfonyl_chloride.coord_atoms, reaction.clash_max) == False:

								kept_H = None

								### --- Not a filter : PRIMARY amine case --- ###
								### Used to know how many hydrogens attached to N as to be deleted ###
								if amine.type_amine == "PRIMARY":

									best_angle_H_N_S = 0

									for id_H, coord_H in amine.ids_and_coords_H:

											### Getting the best H_N_S angle to avoid H pointing toward S-N bond ###
											angle_H_N_S = Topological.angle(coord_H, amine.coord_N, sulfonyl_chloride.coord_S) 

											if angle_H_N_S > best_angle_H_N_S:
												best_angle_H_N_S = angle_H_N_S
												kept_H = id_H

								### --- Congratz ! --- ###
								### If it reaches this point, it means that two chemical reagent docking poses can be combined ###

								### --- Important informations for the merger --- ###

								atoms_to_delete_amine = [id_H for id_H, coord_H in amine.ids_and_coords_H if id_H != kept_H]
								atoms_to_modify_amine = {amine.id_N : "N.pl3"}
								bonds_to_modify_amine = {}

								atoms_to_delete_sulfonyl_chloride = [sulfonyl_chloride.id_Cl]
								atoms_to_modify_sulfonyl_chloride = {}
								bonds_to_modify_sulfonyl_chloride = {}

								connexion_rules = [(("1", amine.id_N), ("2", sulfonyl_chloride.id_S), "1")]							

								nom_pose_amine = amine.path.split("/")[len(amine.path.split("/")) - 1]
								nom_pose_sulfonyl_chloride = sulfonyl_chloride.path.split("/")[len(sulfonyl_chloride.path.split("/")) - 1]

								dico_information_mol2 = {"1" : [m2r.Loader.load_mol2(amine.path), atoms_to_delete_amine, atoms_to_modify_amine, bonds_to_modify_amine], 
															"2" : [m2r.Loader.load_mol2(sulfonyl_chloride.path), atoms_to_delete_sulfonyl_chloride, atoms_to_modify_sulfonyl_chloride, bonds_to_modify_sulfonyl_chloride]}
								
								### --- Output side informations --- ###

								name_pose_amide, num_pose_amide = Output.name_and_num_pose_from_path(amine.path)
								name_pose_sulfonyl_chloride, num_pose_sulfonyl_chloride = Output.name_and_num_pose_from_path(sulfonyl_chloride.path)

								name_merged_ = '{}_{}_{}_{}_10'.format(amine.name, sulfonyl_chloride.name, num_pose_amide, num_pose_sulfonyl_chloride)
								tsv_output = [name_pose_amide, name_pose_sulfonyl_chloride]

								### --- mol2 merging is launched --- ###

								dico_merged_mol2[name_merged_] = (m2r.Merger.merge_mol2(dico_information_mol2, connexion_rules, name_merged_), tsv_output)


				except:
					print("Error during recombination for Sulfonamide")

			### --- Outputs --- ###

			if len(dico_merged_mol2) > Synthesis.thresold_write:

				if ichem_conf_ != None:
					dico_merged_mol2_IFP_filtered = Topological.IFP(dico_merged_mol2, ichem_conf_)
					Output.write_output(dico_merged_mol2_IFP_filtered, name_output_, "10")
				else:
					Output.write_output(dico_merged_mol2, name_output_, "10")
				dico_merged_mol2 = {}

		if ichem_conf_ != None:
			dico_merged_mol2_IFP_filtered = Topological.IFP(dico_merged_mol2, ichem_conf_)
			Output.write_output(dico_merged_mol2_IFP_filtered, name_output_, "10")
		else:
			Output.write_output(dico_merged_mol2, name_output_, "10")

	def benzoxazole(aminophenols_, benzaldehydes_, name_output_, ichem_conf_):

		"""Benzoxazole workflow : cmass -> connectable atoms -> angles oxazole cycle -> clashes

			It is series of IF statement
			
			topological_value = SOMETHING
			if topological_value > "or" < THRESOLD:
				new_topological_value = SOMETHING
				if ...
					...

			IChem can be activated to filter out combination by IFPs similarity
			ichem.conf has to be filled

			Parameters
			----------
			aminophenols_ : list
				list of building block aminophenol docking poses (object : Reactant::Aminophenol)
			benzaldehydes_ : list
				list of building block benzaldehyde docking poses (object : Reactant::Benzaldehyde)
			name_output_ : str
				.tsv and .mol2 output name
			ichem_ : tuple 
				By default None is any ichem.conf file is given as paramter.
				tuple containing (path_ichem, thresold_FULL, thresold_POLAR, ligand and protein as mol2)

		"""

		dico_merged_mol2 = {}

		reaction = Benzoxazole() # Setting the rules for Benzoxazole

		for aminophenol in aminophenols_:
			for benzaldehyde in benzaldehydes_:

				try:
					### --- 1st filter centre of mass --- ###
					distance_between_center_mass = math.dist(aminophenol.coord_center_of_mass, benzaldehyde.coord_center_of_mass)
					if reaction.distance_between_center_mass_down <= distance_between_center_mass <= reaction.distance_between_center_mass_up:

						### --- 2nd filter distance between connectable atoms "O-C.2 and N-C.2" --- ###
						### For benzoxazole there is a double distance check because we are doing a cycle ###

						distance_between_N_C2 = math.dist(aminophenol.coord_N, benzaldehyde.coord_C2)
						distance_between_O_C2 = math.dist(aminophenol.coord_O, benzaldehyde.coord_C2)

						if reaction.distance_between_connectable_atoms_down <= distance_between_N_C2 <= reaction.distance_between_connectable_atoms_up\
							and reaction.distance_between_connectable_atoms_down <= distance_between_O_C2 <= reaction.distance_between_connectable_atoms_up:

							### --- 3rd filter oxazole cycle angles "C.ar-N-C2 and C.ar-O-C2" --- ###
							angle_neighborN_N_C2 = Topological.angle(aminophenol.coord_neighbor_N, aminophenol.coord_N, benzaldehyde.coord_C2)
							angle_neighborO_O_C2 = Topological.angle(aminophenol.coord_neighbor_O, aminophenol.coord_O, benzaldehyde.coord_C2)

							if reaction.angle_cycle_down <= angle_neighborN_N_C2 <= reaction.angle_cycle_up\
								and reaction.angle_cycle_down <= angle_neighborO_O_C2 <= reaction.angle_cycle_up:

								### --- 4th filter clashes --- ###
								if Topological.check_clashes(aminophenol.coord_atoms, benzaldehyde.coord_atoms, reaction.clash_max) == False:

									### --- Congratz ! --- ###
									### If it reaches this point, it means that two chemical reagent docking poses can be combined ###

									### --- Important informations for the merger --- ###

									atoms_to_delete_aminophenol = aminophenol.ids_H
									atoms_to_modify_aminophenol = {aminophenol.id_N : "N.2"}
									bonds_to_modify_aminophenol = {}

									atoms_to_delete_benzaldehyde = benzaldehyde.ids_H + [benzaldehyde.id_O2]
									atoms_to_modify_benzaldehyde = {}
									bonds_to_modify_benzaldehyde = {}

									connexion_rules = [(("1", aminophenol.id_N), ("2", benzaldehyde.id_C2), "2"), (("1", aminophenol.id_O), ("2", benzaldehyde.id_C2), "1")]

									dico_information_mol2 = {"1" : [m2r.Loader.load_mol2(aminophenol.path), atoms_to_delete_aminophenol, atoms_to_modify_aminophenol, bonds_to_modify_aminophenol], 
																"2" : [m2r.Loader.load_mol2(benzaldehyde.path), atoms_to_delete_benzaldehyde, atoms_to_modify_benzaldehyde, bonds_to_modify_benzaldehyde]}

									### --- Output side informations --- ###

									name_pose_aminophenol, num_pose_aminophenol = Output.name_and_num_pose_from_path(aminophenol.path)
									name_pose_benzaldehyde, num_pose_benzaldehyde = Output.name_and_num_pose_from_path(benzaldehyde.path)

									name_merged_ = '{}_{}_{}_{}_27'.format(aminophenol.name, benzaldehyde.name, num_pose_aminophenol, num_pose_benzaldehyde)
									tsv_output = [name_pose_aminophenol, name_pose_benzaldehyde]

									### --- mol2 merging is launched --- ###

									dico_merged_mol2[name_merged_] = (m2r.Merger.merge_mol2(dico_information_mol2, connexion_rules, name_merged_), tsv_output)

				except:
					print("Error during recombination for Benzoxazole")

			### --- Outputs --- ###

			if len(dico_merged_mol2) > Synthesis.thresold_write:

				if ichem_conf_ != None:
					dico_merged_mol2_IFP_filtered = Topological.IFP(dico_merged_mol2, ichem_conf_)
					Output.write_output(dico_merged_mol2_IFP_filtered, name_output_, "27")
				else:
					Output.write_output(dico_merged_mol2, name_output_, "27")
				dico_merged_mol2 = {}

		if ichem_conf_ != None:
			dico_merged_mol2_IFP_filtered = Topological.IFP(dico_merged_mol2, ichem_conf_)
			Output.write_output(dico_merged_mol2_IFP_filtered, name_output_, "27")
		else:
			Output.write_output(dico_merged_mol2, name_output_, "27")
