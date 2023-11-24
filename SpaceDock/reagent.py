""" SpaceDock Reagent Version 1.0.0
This script allows to load and annote building block docking pose to be used by SpaceDock methology 

requirement : 
Python 3.9+ (it used math.dist)
mol2reader.py and reaction.py

"""

import mol2parser as m2r
from reaction import Topological

class Loader:

	"""Class that contain fonction to load and annote building block docking pose

		Methods
		-------
		load_tag_table(path_tag_table_, selected_reactions_)
			Open and put into a dictionary the tag table
		load_poses(list_path_docking_poses_, dico_tag_table_, selected_reactions_, start_, end_)
			Open and annote reagent docking pose if they correspond to at least one selection reactions
			
	"""

	def load_tag_table(path_tag_table_, selected_reactions_):

		"""Open and put into a dictionary the tag table 
			BUT ONLY for selected reactions
			This table contains all informations to annote building block docking pose
			It makes the link between the name of chemical reagent and its reactivity

			Parameters
			----------
			path_tag_table_ : str
				tag table path
			selected_reactions_ : list
				selected reactions list. e.g. ["Amide", "Sulfonamide"]
				
			Return
			------
			dico_tag_table : dict
				dictionary containing as key the name of the building block and values are chemotype, tag ids, ...

		"""

		dico_tag_table = {}

		with open(path_tag_table_, 'r') as file_tag_table_hdl:
			line_tag_table = file_tag_table_hdl.read().split("\n")
			for line in line_tag_table:
				if line != "":
					try:
						### The order of this change in the new version of Tag Table ! If update, it has to be change !! ###
						name, reaction, chemotype, tag_ids = line.split("\t")

						for selected_reaction in selected_reactions_:
							if selected_reaction.find(reaction) != -1:
								if not name in dico_tag_table:
									dico_tag_table[name] = [[reaction, chemotype, tag_ids.split()]]
								else:
									dico_tag_table[name].append([reaction, chemotype, tag_ids.split()])
					except:
						print("Incomplete tag table", line)

		return dico_tag_table

	def load_poses(list_path_docking_poses_, dico_tag_table_, selected_reactions_, start_, end_):

		"""Open and annote reagent docking pose if they correspond to at least one selection reactions

			Parameters
			----------
			list_path_docking_poses_ : list
				list of reagent docking pose path to annote
			dico_tag_table_ : dict
				tag table loaded with Loader.load_tag_table() function
			selected_reactions_ : list
				selected reactions list. e.g. ["Amide", "Sulfonamide"]
			start_ : int
				starting point to start annotation of poses
				Can be different to 1 in case packet are made to distribute on a calculation center (or few cores on local computer)
			end_ : int
				ending point to start annotation of poses
				Can be different to 1 in case packet are made to distribute on a calculation center (or few cores on local computer)
				
			Return
			------
			dico_list_poses : dict
				dictionary containing as key the reaction name and as value a tuple of corresponding annoted reagent 

		"""
		
		dico_list_poses = {"Amide" : ([], []), "Sulfonamide" : ([], []), "Benzoxazole" : ([], [])}

		compteur_passed = 0
		total = len(list_path_docking_poses_)

		for pose in list_path_docking_poses_:
			if pose != "":
				name, path = pose.split("\t")

				### We generate 3D structure with corina, it adds _i00X for each isomers ###
				name_without_stereoisomer = name[:name.find("_i")]
				
				### In case I would like to take into account docking score ###
				### not yet implemented or never maybe, compounds are too small to be correctly scored ###
				score = 0

				try:

					### reaction, chemotype, tag ids,  are extracted according the tag table ###
					reagent_informations = dico_tag_table_[name_without_stereoisomer]

					for information in reagent_informations:
						reaction, chemotype, tag_ids = information

						mol2 = m2r.Loader.load_mol2(path)

						### Annotation here ###
						### To add a new reaction ###
						### copy/paste here a new elif and add the reaction into dico_list_poses ###

						if reaction == "Amide" and "Amide" in selected_reactions_:
							if chemotype == "Amine":
								if compteur_passed >= start_ and compteur_passed <= end_:
									dico_list_poses["Amide"][0].append(Amine(name, chemotype, reaction, tag_ids, path, score, mol2))
							elif chemotype == "Carboxylic acid":
								dico_list_poses["Amide"][1].append(CarboxylicAcid(name, chemotype, reaction, tag_ids, path, score, mol2))
						
						elif reaction == "Sulfonamide" and "Sulfonamide":
							if chemotype == "Amine":
								if compteur_passed >= start_ and compteur_passed <= end_:
									dico_list_poses["Sulfonamide"][0].append(Amine(name, chemotype, reaction, tag_ids, path, score, mol2))
							elif chemotype == "Sulfonylchloride":
								dico_list_poses["Sulfonamide"][1].append(SulfonylChloride(name, chemotype, reaction, tag_ids, path, score, mol2))
						
						elif reaction == "Benzoxazole" and "Benzoxazole" in selected_reactions_:
							if chemotype == "Aminophenol":
								if compteur_passed >= start_ and compteur_passed <= end_:
									dico_list_poses["Benzoxazole"][0].append(Aminophenol(name, chemotype, reaction, tag_ids, path, score, mol2))
							elif chemotype == "Benzaldehyde":
								dico_list_poses["Benzoxazole"][1].append(Benzaldehyde(name, chemotype, reaction, tag_ids, path, score, mol2))

				except:
					print("Issue with", name, "not in tag IDs or bad structure")
					pass

			compteur_passed += 1

			if compteur_passed % int(total/10) == 0:
				print("Loading BBs... ", round(compteur_passed/total*100, 1))

		return dico_list_poses

class Reactant:

	"""Parent class to annote building block docking pose
		
		Attributes
		----------
		name_ : str
			Building block name. e.g. for Enamine EN300-XXXX
		chemotype_ : str
			chemotype of the building block
		reaction_ : str
			Reaction in which the building block is involved
		tag_ids_ : list
			list of reactive atom in the tag table
		path_ : str
			path of mol2 reagent docking pose
		score_ : float
			docking score (never used)
		forbidden_atoms_ids_ : list
			list forbidden atom ids that will be discarded for clashes calculation
		mol2_ : Mol2
			loaded mol2 reagent docking pose
			
	"""

	def __init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, forbidden_atoms_ids_, mol2_):
		self.name = name_
		self.chemotype = chemotype_
		self.reaction = reaction_
		self.path = path_
		self.score = score_
		self.coord_center_of_mass = mol2_.get_center_of_mass()
		self.coord_atoms = mol2_.get_coord_heavy_atoms(False, forbidden_atoms_ids_)

class Amine(Reactant):

	"""class to annote Amide building block docking pose

		/!\  Nitrogen position differ from reaction in the tag table
		Tag table : 
		Amide : 0
		Sulfomamide : 0

		Reactive number ids and coordinates are extracted when the object is created
		
		Attributes
		----------
		name_ : str
			Building block name. e.g. for Enamine EN300-XXXX
		chemotype_ : str
			chemotype of the building block
		reaction_ : str
			Reaction in which the building block is involved
		tag_ids_ : list
			list of reactive atom in the tag table
		path_ : str
			path of mol2 reagent docking pose
		score_ : float
			docking score (never used)
		mol2_ : Mol2
			loaded mol2 reagent docking pose
			
	"""

	def __init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, mol2_):

		position_N = 0

		self.id_N = int(tag_ids_[position_N])
		self.coord_N = mol2_.get_coord_from_atom_id(self.id_N)

		id_neighbors_N = mol2_.get_neighbors_id_from_atom_id(self.id_N)
		self.coord_neighbors_N = []

		### Define if amine is primary or secondary ###
		count_C = 0 

		### Search of H attached to N ###
		self.ids_and_coords_H = []

		for id_neighbor in id_neighbors_N:
			type_neighbor = mol2_.atoms[id_neighbor - 1].type_atom
			if type_neighbor == "H":
				self.ids_and_coords_H.append([id_neighbor, mol2_.get_coord_from_atom_id(id_neighbor)])
			else:
				self.coord_neighbors_N.append(mol2_.get_coord_from_atom_id(id_neighbor))
				count_C += 1

		if count_C == 2:
			self.type_amine = "SECONDARY"
		elif count_C == 1:
			self.type_amine = "PRIMARY"

		### For Amine, N is discarded for clash calculation ###
		forbidden_atoms_ids = [self.id_N]

		### Override, common of all building block ###
		Reactant.__init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, forbidden_atoms_ids, mol2_)

class CarboxylicAcid(Reactant):

	"""class to annote Carboxylic acids building block docking pose

		/!\ sp2 carboxylate carbon position differ from reaction in the tag table
		Tag table : 
		Amide : 0

		Reactive number ids and coordinates are extracted when the object is created
		
		Attributes
		----------
		name_ : str
			Building block name. e.g. for Enamine EN300-XXXX
		chemotype_ : str
			chemotype of the building block
		reaction_ : str
			Reaction in which the building block is involved
		tag_ids_ : list
			list of reactive atom in the tag table
		path_ : str
			path of mol2 reagent docking pose
		score_ : float
			docking score (never used)
		mol2_ : Mol2
			loaded mol2 reagent docking pose
			
	"""

	def __init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, mol2_):

		position_C2 = 0
		
		self.id_C2 = int(tag_ids_[position_C2])
		self.coord_C2 = mol2_.get_coord_from_atom_id(self.id_C2)

		id_neighbors_C2 = mol2_.get_neighbors_id_from_atom_id(self.id_C2)

		### For Carboxylic acids, C.2 is discarded for clash calculation ###
		### Oxygens will be added later ###
		forbidden_atoms_ids = [self.id_C2]

		### Search of O.co2 attached to C.2 ###
		ids_and_coords_Oco2 = []

		for id_neighbor in id_neighbors_C2:
			type_neighbor = mol2_.atoms[int(id_neighbor) - 1].type_atom
			if type_neighbor.find("O") != -1:
				ids_and_coords_Oco2 += [[id_neighbor, mol2_.get_coord_from_atom_id(id_neighbor)]]
				forbidden_atoms_ids += [id_neighbor]
			else:
				coord_neighbor_C2 = mol2_.get_coord_from_atom_id(id_neighbor)

		self.ids_and_coords_Oco2 = ids_and_coords_Oco2
		self.coord_neighbor_C2 = coord_neighbor_C2

		### Override, common of all building block ###
		Reactant.__init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, forbidden_atoms_ids, mol2_)

class SulfonylChloride(Reactant):

	"""class to annote Sulfonyl chloride building block docking pose
		
		Tag table : 
		Sulfone position : 0

		Reactive number ids and coordinates are extracted when the object is created
		
		Attributes
		----------
		name_ : str
			Building block name. e.g. for Enamine EN300-XXXX
		chemotype_ : str
			chemotype of the building block
		reaction_ : str
			Reaction in which the building block is involved
		tag_ids_ : list
			list of reactive atom in the tag table
		path_ : str
			path of mol2 reagent docking pose
		score_ : float
			docking score (never used)
		mol2_ : Mol2
			loaded mol2 reagent docking pose
			
	"""

	def __init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, mol2_):

		position_S = 0
		position_Cl = 1

		self.id_S = int(tag_ids_[position_S])
		self.coord_S = mol2_.get_coord_from_atom_id(self.id_S)

		self.id_Cl = int(tag_ids_[position_Cl])

		### For Sulfonyl chloride, S and Cl are discarded for clash calculation ###
		### There is no special angle check due to its particular geometry ###
		### So Oxygens are kept for clashes calculations to filter bad recombination ###
		forbidden_atoms_ids = [self.id_S, self.id_Cl]

		### Override, common of all building block ###
		Reactant.__init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, forbidden_atoms_ids, mol2_)

class Aminophenol(Reactant):

	"""class to annote Aminophenol building block docking pose
	
		Tag table : 
		N position : 0
		O postion : 1

		Reactive number ids and coordinates are extracted when the object is created
		
		Attributes
		----------
		name_ : str
			Building block name. e.g. for Enamine EN300-XXXX
		chemotype_ : str
			chemotype of the building block
		reaction_ : str
			Reaction in which the building block is involved
		tag_ids_ : list
			list of reactive atom in the tag table
		path_ : str
			path of mol2 reagent docking pose
		score_ : float
			docking score (never used)
		mol2_ : Mol2
			loaded mol2 reagent docking pose
			
	"""

	def __init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, mol2_):

		position_N = 3
		position_neighbor_N = 2 # Is it a C.ar
		position_O = 1
		position_neighbor_O = 0 # Is it a C.ar

		self.id_N = int(tag_ids_[position_N])
		self.coord_N = mol2_.get_coord_from_atom_id(self.id_N)
		self.coord_neighbor_N = mol2_.get_coord_from_atom_id(tag_ids_[position_neighbor_N])

		self.id_O = int(tag_ids_[position_O])
		self.coord_O = mol2_.get_coord_from_atom_id(self.id_O)
		self.coord_neighbor_O = mol2_.get_coord_from_atom_id(tag_ids_[position_neighbor_O])

		### Search of H attached to N and O (They will be deleted) ###
		id_neighbors_N_O = mol2_.get_neighbors_id_from_atom_id(self.id_N) +  mol2_.get_neighbors_id_from_atom_id(self.id_O)
		self.ids_H = []

		for id_neighbor in id_neighbors_N_O:
			type_neighbor = mol2_.atoms[int(id_neighbor) - 1].type_atom
			if type_neighbor == "H":
				self.ids_H.append(id_neighbor)

		### For Aminophenol, N and O are discarded for clash calculation ###
		forbidden_atoms_ids = [self.id_N, self.id_O]

		### Override, common of all building block ###
		Reactant.__init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, forbidden_atoms_ids, mol2_)

class Benzaldehyde(Reactant):
	
	"""class to annote Benzaldehyde building block docking pose
		
		Tag table : 
		sp2 oxygen position : 2

		Reactive number ids and coordinates are extracted when the object is created
		
		Attributes
		----------
		name_ : str
			Building block name. e.g. for Enamine EN300-XXXX
		chemotype_ : str
			chemotype of the building block
		reaction_ : str
			Reaction in which the building block is involved
		tag_ids_ : list
			list of reactive atom in the tag table
		path_ : str
			path of mol2 reagent docking pose
		score_ : float
			docking score (never used)
		mol2_ : Mol2
			loaded mol2 reagent docking pose
			
	"""

	def __init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, mol2_):

		position_C2 = 1
		position_O2 = 2

		self.id_C2 = int(tag_ids_[position_C2])
		self.coord_C2 = mol2_.get_coord_from_atom_id(self.id_C2)

		self.id_O2 = int(tag_ids_[position_O2])

		id_neighbors_C2 = mol2_.get_neighbors_id_from_atom_id(self.id_C2)
		self.ids_H = []

		for id_neighbor in id_neighbors_C2:
			type_neighbor = mol2_.atoms[int(id_neighbor) - 1].type_atom
			if type_neighbor == "H":
				self.ids_H.append(id_neighbor)

		### For Benzaldehyde, C.2 and O.2 are discarded for clash calculation ###
		### Because the connection will be C.2 to O and N aminophenol ###
		### So C.2 don't have to be taken into account for clashes ###
		forbidden_atoms_ids = [self.id_C2, self.id_O2]

		### Override, common of all building block ###
		Reactant.__init__(self, name_, chemotype_, reaction_, tag_ids_, path_, score_, forbidden_atoms_ids, mol2_)


