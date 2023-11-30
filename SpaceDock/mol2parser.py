class Mol2:


	"""Class used to represent a mol2 file
		
		Attributes
		----------
		name_ : str
			mol2 file name
		molecule_ : Molecule
			@<TRIPOS>MOLECULE extracted as object Molecule
		atoms_ : list
			@<TRIPOS>ATOM extracted as list of object Atom
		bonds_ : list
			@<TRIPOS>BOND extracted as list of object Bond
		substructures_ : list
			@<TRIPOS>SUBSTRUCTURE extracted as list of object Substructure


		Methods
		-------
		set_name(name_)
			set file name by modifing name attribute
		add_atom(atom_)
			add a atom to the list of atoms given as attributes
		add_bond(bond_)
			add a bond to the list of bonds given as attributes
		build()
			construct and return the mol2 file as a str

		get_coord_from_atom_id(id_atom_)
			return coordinate [x, y, z] of an atom

		get_center_of_mass(forbidden_atoms_ids = [])
			return center of mass [x, y, z] of the molecule

		get_neighbors_id_from_atom_id(id_atom_)
			return neighbors id of an atom

		get_neighbors_bond_type_from_atom_id(id_atom_)
			return bond types connected to an atom

		get_coord_atoms(xtract_type_atoms_ = False, forbidden_atom_ids_ = [])
			return list of coordinates [[x1, y1, z1], [x2, y2, z2], ...] of all atoms

		get_coord_heavy_atoms(xtract_type_atoms_ = False, forbidden_atom_ids_ = [])
			return list of coordinates [[x1, y1, z1], [x2, y2, z2], ...] of heavy atoms

	"""

	def __init__(self, name_, molecule_, atoms_, bonds_, substructures_):
		self.name = name_
		self.molecule = molecule_
		self.atoms = atoms_
		self.bonds = bonds_
		self.substructures = substructures_

	def set_name(self, name_):

		"""set file name

			Parameters
			----------
			name_ : str
				file name

		"""

		self.name = name_

	def add_atom(self, atom_):

		"""add a atom to the list of atoms given as attributes

			Parameters
			----------
			atom_ : Atom
				Atom object

		"""

		self.atoms.append(atom_)

	def add_bond(self, bond_):

		"""add a bond to the list of bonds given as attributes

			Parameters
			----------
			bond_ : Bond
				Bond object

		"""

		self.bonds.append(bond_)

	def build(self):

		"""Construct and return the mol2 file as a str
			This function call all build() fonction of each @<TRIPOS>XXXX section

			Returns
			-------
			'{}{}{}{}'.format(string_molecule, string_atoms, string_bonds, string_substructures) : str
				formated mol2 file

		"""
		string_molecule = '@<TRIPOS>MOLECULE\n'
		string_molecule += self.molecule.build()
		string_atoms = '@<TRIPOS>ATOM\n'
		for atom in self.atoms:
			string_atoms += '{}\n'.format(atom.build())
		string_bonds = '@<TRIPOS>BOND\n'
		for bond in self.bonds:
			string_bonds += '{}\n'.format(bond.build())
		string_substructures = '@<TRIPOS>SUBSTRUCTURE\n'
		for substructure in self.substructures:
			string_substructures += '{}\n'.format(substructure.build())

		return '{}{}{}{}'.format(string_molecule, string_atoms, string_bonds, string_substructures)

	def get_coord_from_atom_id(self, id_atom_):

		"""return coordinate [x, y, z] of an atom

			Parameters
			----------
			id_atom_ : int
				atom id

			Returns
			-------
			[x, y, z] : list
				atom coordinates

		"""

		x = self.atoms[int(id_atom_) - 1].x
		y = self.atoms[int(id_atom_) - 1].y
		z = self.atoms[int(id_atom_) - 1].z
		return [x, y, z]

	def get_center_of_mass(self, forbidden_atoms_ids = []):

		"""return coordinate [x, y, z] of an atom
			Sum all x, y, z coordinate of a molecule and divide them by number of atoms (substracted by number of forbidden atoms )


			Parameters
			----------
			forbidden_atoms_ids : list
				list of atom id that will be excluded from cmass calculation

			Returns
			-------
			[x_cmass/(n_atoms - count_skipped_atoms), y_cmass/(n_atoms - count_skipped_atoms), z_cmass/(n_atoms - count_skipped_atoms)] : list
				[x_cmass, y_cmass, z_cmass] center of mass 

		"""
		n_atoms = len(self.atoms)
		x_cmass = 0
		y_cmass = 0
		z_cmass = 0
		count_skipped_atoms = 0
		if forbidden_atoms_ids == []:
			for atom in self.atoms:
				x_cmass += atom.x
				y_cmass += atom.y
				z_cmass += atom.z
		else:
			for atom in self.atoms:
				if not atom.id_atom in forbidden_atoms_ids:
					x_cmass += atom.x
					y_cmass += atom.y
					z_cmass += atom.z
				else:
					count_skipped_atoms += 1

		return [x_cmass/(n_atoms - count_skipped_atoms), y_cmass/(n_atoms - count_skipped_atoms), z_cmass/(n_atoms - count_skipped_atoms)]

	def get_neighbors_id_from_atom_id(self, id_atom_):

		"""return coordinate neighbor ids of an atom

			Parameters
			----------
			id_atom_ : int
				atom id

			Returns
			-------
			list_neighbors : list
				list of neighbor atom ids

		"""

		list_neighbors = []
		for bond in self.bonds:
			if bond.id_atom1 == int(id_atom_):
				list_neighbors += [bond.id_atom2]
			elif bond.id_atom2 == int(id_atom_):
				list_neighbors += [bond.id_atom1]
		return list_neighbors

	def get_neighbors_bond_type_from_atom_id(self, id_atom_):

		"""return bond types linked to an atom

			Parameters
			----------
			id_atom_ : int
				atom id

			Returns
			-------
			list_neighbors_bond_type : list
				list neighbor bond types

		"""

		list_neighbors_bond_type = []
		for bond in self.bonds:
			if bond.id_atom1 == int(id_atom_):
				list_neighbors_bond_type.append(bond.type_bond)
			elif bond.id_atom2 == int(id_atom_):
				list_neighbors_bond_type.append(bond.type_bond)
		return list_neighbors_bond_type

	def get_coord_atoms(self, xtract_type_atoms_ = False, forbidden_atom_ids_ = []):

		"""return list of coordinates [[x1, y1, z1], [x2, y2, z2], ...] of all atoms
			Why allowing atom type ? In case I would like to adapt to calculation of vdw clashes

			Parameters
			----------
			xtract_type_atoms_ : boolean
				True : allow to extract coordinates and atom type

			Returns
			-------
			list_coord_atoms : list
				coordinates list or coordinates + atom type list

		"""

		list_coord_atoms = []
		for atom in self.atoms:
			if not atom.id_atom in forbidden_atom_ids_:
				if xtract_type_atoms_ == False:
					list_coord_atoms += [[atom.x, atom.y, atom.z]]
				else:
					list_coord_atoms += [[atom.x, atom.y, atom.z, atom.type_atom]]
		return list_coord_atoms

	def get_coord_heavy_atoms(self, xtract_type_atoms_ = False, forbidden_atom_ids_ = []):

		"""return list of coordinates [[x1, y1, z1], [x2, y2, z2], ...] of heavy atoms
			No H considered here
			Why allowing atom type ? In case I would like to adapt to calculation of vdw clashes

			Parameters
			----------
			xtract_type_atoms_ : boolean
				True : allow to extract coordinates and atom type

			Returns
			-------
			list_coord_atoms : list
				coordinates list or coordinates + atom type list

		"""

		list_coord_atoms = []
		for atom in self.atoms:
			if not atom.id_atom in forbidden_atom_ids_ and atom.type_atom != "H":
				if xtract_type_atoms_ == False:
					list_coord_atoms += [[atom.x, atom.y, atom.z]]
				else:
					list_coord_atoms += [[atom.x, atom.y, atom.z, atom.type_atom]]
		return list_coord_atoms

class Molecule:

	"""Class used to represent @<TRIPOS>MOLECULE
		num_feat and num_sets are not supported. Set by 0 by default
		
		Attributes
		----------
		name_ : str
			the name of the molecule
		n_atoms : int
			the number of atom in the molecule
		n_bonds : int
			the number of bonds in the molecule
		n_substructure : int
			the number of substructures in the molecule
		mol_type : str
			the molecule type : SMALL, BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE
		charge_type : str
			the type of charges associated with the molecule : NO_CHARGES, USER_CHARGES


		Methods
		-------
		set_name(name_)
			set molecule name by modifing name attribute
		set_n_atoms(n_atoms_)
			set the number of atoms by modifing n_atoms attribute
		set_n_bonds(n_bonds_)
			set the number of bonds by modifing n_bonds attribute
		build()
			Construct et return @<TRIPOS>MOLECULE formated for mol2 file
	"""

	def __init__(self, name_, n_atoms_, n_bonds_, n_substructure_, mol_type_, charge_type_):
		self.name = name_
		self.n_atoms = str(n_atoms_)
		self.n_bonds = str(n_bonds_)
		self.n_substructure = str(n_substructure_)
		self.mol_type = mol_type_
		self.charge_type = charge_type_

	def set_name(self, name_):

		"""set the name of the molecule

			Parameters
			----------
			name_ : str
				molecule name

		"""

		self.name = name_	

	def set_n_atoms(self, n_atoms_):

		"""set the number of atom in the molecule

			Parameters
			----------
			n_atoms_ : int
				atom number

		"""

		self.n_atoms = n_atoms_	

	def set_n_bonds(self, n_bonds_):

		"""the number of bonds in the molecule

			Parameters
			----------
			n_bonds_ : int
				bond number

		"""

		self.n_bonds = n_bonds_	

	def build(self):

		"""Construct et return @<TRIPOS>MOLECULE formated for mol2 file

			Returns
			----------
			string_molecule : str
				formated @<TRIPOS>MOLECULE section

		"""

		string_molecule = '{}\n  {}   {}   {}   {}   {}\n{}\n{}\n'.format(self.name, 
																self.n_atoms, self.n_bonds, self.n_substructure, "0", "0",
																self.mol_type,
																self.charge_type)

		return string_molecule

class Atom:

	"""Class used to represent @<TRIPOS>ATOM
		
		Attributes
		----------
		id_atom : int
			the id number of the atom
		name : str
			the name of the atom
		x : float
			the x coordinate of the atom
		y : float
			the y coordinate of the atom
		z : float
			the z coordinate of the atom
		type_atom : str
			 the SYBYL atom type for the atom
		id_sub : int
			the id number of the substructure containing the atom
		name_sub : str
			the name of the substructure containing the atom
		charge : float
			the charge associated with the atom


		Methods
		-------
		set_id_atom(new_id_atom_)
			set the id number of the atom by modifing id_atom attribute
		set_name(new_name_)
			set the name of the atom by modifing name attribute
		set_x(new_x_)
			set the x coordinate of the atom by modifing x attribute
		set_y(new_y_)
			set the y coordinate of the atom by modifing y attribute
		set_z(new_z_)
			set the z coordinate of the atom by modifing z attribute
		set_type_atom(new_type_atom_) 
			set the SYBYL atom type by modifing type_atom attribute
		set_id_sub(new_id_sub_)
			set the id number of the substructure by modifing id_sub attribute
		set_name_sub(new_name_sub_)
			set the name of the substructure by modifing name_sub attribute
		set_charge(new_charge_)
			set the charge by modifing charge attribute
		build()
			Construct et return corresponding atom for @<TRIPOS>ATOM formated for mol2 file
	"""

	def __init__(self, id_atom_, name_, x_, y_, z_, type_atom_, id_sub_, name_sub_, charge_):
		self.id_atom = int(id_atom_)
		self.name = name_
		self.x = float(x_)
		self.y = float(y_)
		self.z = float(z_)
		self.type_atom = type_atom_
		self.id_sub = int(id_sub_)
		self.name_sub = name_sub_
		self.charge = float(charge_)

	def set_id_atom(self, new_id_atom_):

		"""set the id number of the atom

			Parameters
			----------
			new_id_atom_ : int
				the id number of the atom

		"""

		self.id_atom = int(new_id_atom_)

	def set_name(self, new_name_):

		"""set the name of the atom

			Parameters
			----------
			new_name_ : str
				atom name

		"""

		self.name = str(new_name_)

	def set_x(self, new_x_):
		
		"""set the x coordinate of the atom

			Parameters
			----------
			new_x_ : float
				x coordinate

		"""

		self.x = float(new_x_)

	def set_y(self, new_y_):
		
		"""set the y coordinate of the atom

			Parameters
			----------
			new_y_ : float
				y coordinate

		"""

		self.y = float(new_y_)

	def set_z(self, new_z_):
		
		"""set the z coordinate of the atom

			Parameters
			----------
			new_z_ : float
				z coordinate

		"""

		self.z = float(new_z_)

	def set_type_atom(self, new_type_atom_):
		
		"""set the SYBYL atom type of the atom

			Parameters
			----------
			new_type_atom_ : str
				SYBYL atom type

		"""

		self.type_atom = str(new_type_atom_)

	def set_id_sub(self, new_id_sub_):
		
		"""set the id number of the substructure of the atom

			Parameters
			----------
			new_id_sub_ : int
				id number of the substructure

		"""

		self.id_sub = int(new_id_sub_)
	
	def set_name_sub(self, new_name_sub_):

		"""set the name of the substructure of the atom

			Parameters
			----------
			new_name_sub_ : str
				name of the substructure

		"""

		self.name_sub = str(new_name_sub_)

	def set_charge(self, new_charge_):

		"""set the charge of the atom

			Parameters
			----------
			new_charge_ : float
				charge of the atom

		"""

		self.charge = float(new_charge_)

	def build(self):

		"""Construct et return an atom of @<TRIPOS>ATOM formated for mol2 file

			Returns
			----------
			string_atom : str
				formated atom for @<TRIPOS>ATOM

		"""

		string_atom = '{:6d} {:6} {:10.4f} {:10.4f} {:10.4f} {:6} {:5d} {:8} {:10.4f}'.format(self.id_atom, self.name, 
																								self.x, self.y, self.z,
																								self.type_atom,
																								self.id_sub, self.name_sub,
																								self.charge)
		return string_atom

class Bond:

	"""Class used to represent @<TRIPOS>BOND
		
		Attributes
		----------
		id_bond : int
			the id number of the bond
		id_atom1 : int
			the id number of the atom at one end of the bond
		id_atom2 : int
			the id number of the atom at the other end of the bond
		type_bond : str
			 the SYBYL bond type


		Methods
		-------
		set_id_bond(new_id_atom_)
			set the id number of the bond by modifing id_bond attribute
		set_id_atom1(new_name_)
			set the id number of the atom at one end of the bond by modifing id_atom1 attribute
		set_id_atom2(new_x_)
			set the id number of the atom at the other end of the bond by modifing id_atom2 attribute
		set_type_bond(new_y_)
			set the SYBYL bond type
		build()
			Construct et return corresponding bond for @<TRIPOS>BOND formated for mol2 file
	"""

	def __init__(self, id_bond_, id_atom1_, id_atom2_, type_bond_):
		self.id_bond = int(id_bond_)
		self.id_atom1 = int(id_atom1_)
		self.id_atom2 = int(id_atom2_)
		self.type_bond = type_bond_

	def set_id_bond(self, new_id_bond_):

		"""set the id number of the bond

			Parameters
			----------
			new_id_bond_ : int
				the id number of the bond

		"""

		self.id_bond = int(new_id_bond_)

	def set_id_atom1(self, new_id_atom1_):


		"""set the id number of the atom at one end of the bond

			Parameters
			----------
			new_id_atom1_ : int
				the id number of the atom at one end of the bond

		"""

		self.id_atom1 = int(new_id_atom1_)

	def set_id_atom2(self, new_id_atom2_):

		"""set the id number of the atom at the other end of the bond

			Parameters
			----------
			new_id_atom2_ : int
				the id number of the atom at the other end of the bond

		"""

		self.id_atom2 = int(new_id_atom2_)

	def set_type_bond(self, new_type_bond_):


		"""set the SYBYL bond type

			Parameters
			----------
			new_type_bond_ : str
				the SYBYL bond type

		"""
		self.type_bond = new_type_bond_

	def build(self):

		"""Construct et return a bond of @<TRIPOS>BOND formated for mol2 file

			Returns
			----------
			string_bond : str
				formated bond for @<TRIPOS>BOND

		"""

		string_bond = '{:6d} {:6d} {:6d} {}'.format(self.id_bond, self.id_atom1, self.id_atom2, self.type_bond)
		return string_bond

class Substructure:

	"""Class used to represent @<TRIPOS>SUBSTRUCTURE
		Only these three attributes are essential to create a substructure

		
		Attributes
		----------
		id_sub : int
			the id number of the substructure
		name_sub : str
			the name of the substructure
		id_root : int
			the id number of the substructure’s root atom


		Methods
		-------
		build()
			Construct et return corresponding substructure for @<TRIPOS>SUBSTRUCTURE formated for mol2 file
	"""

	def __init__(self, id_sub_, name_sub_, id_root_):
		self.id_sub = id_sub_
		self.name_sub = name_sub_
		self.id_root = id_root_

	def build(self):

		"""Construct et return a substructure of @<TRIPOS>SUBSTRUCTURE formated for mol2 file

			Returns
			----------
			string_substructure : str
				formated substructure for @<TRIPOS>SUBSTRUCTURE

		"""

		string_substructure = '{:6d} {:8} {:6d}'.format(int(self.id_sub), self.name_sub, int(self.id_root))
		return string_substructure

class Loader:

	"""Class that contains fonction to load a mol2 file

		Methods
		-------
		cut_mol2(string_mol2_, string_section_)
			Cut a mol2 as a string to get a particular record types

		extract_molecule(string_molecule_)
			Extract informations from @<TRIPOS>MOLECULE and create Molecule object

		extract_atoms(string_atom_)
			Extract informations from @<TRIPOS>ATOM and create Bond object and store then in a list

		extract_bonds(string_bond_)
			Extract informations from @<TRIPOS>BOND and create Bond object and store then in a list

		extract_substructures(string_substructure_)
			Extract informations from @<TRIPOS>SUBSTRUCTURE and create Substructure object and store then in a list

		load_mol2(path_mol2_)
			Open, parse a mol2 file and create a Mol2 object

	"""

	def cut_mol2(string_mol2_, string_section_):

		"""Cut a mol2 as a string to get a particular record types

			e.g. string_section_ = "@<TRIPOS>ATOM" recover atoms of a mol2
			available : @<TRIPOS>MOLECULE, @<TRIPOS>BOND, @<TRIPOS>ATOM, @<TRIPOS>SUBSTRUCTURE

			Parameters
			----------
			string_mol2_ : str
				mol2 as a string
			string_section_ : str
				@<TRIPOS>XXXX corresponding to particular section of a mol2 e.g @<TRIPOS>ATOM

				
			Return
			------
			section : str
				@<TRIPOS> record types of a mol2 file
		"""

		section_temp = string_mol2_[string_mol2_.find(string_section_)+len(string_section_):]
		section = section_temp[:section_temp.find("@")]
		return section

	def extract_molecule(string_molecule_):

		"""Extract informations from @<TRIPOS>MOLECULE and create Molecule object

			Parameters
			----------
			string_molecule_ : str
				str of the @<TRIPOS>MOLECULE section

				
			Return
			------
			molecule : Molecule
				Molecule object containing molecule name, #atom, #bond, ...
		"""

		compteur = 1

		line_molecule = string_molecule_.split("\n")
		for line in line_molecule:
			if line != "":
				if compteur == 1:
					name = line
					compteur += 1
				elif compteur == 2:
					n_atoms, n_bonds, n_substructure, *r = line.split()
					compteur += 1
				elif compteur == 3:
					mol_type = line
					compteur += 1
				elif compteur == 4:
					charge_type = line

		molecule = Molecule(name, n_atoms, n_bonds, n_substructure, mol_type, charge_type)

		return molecule

	def extract_atoms(string_atom_):

		"""Extract informations from @<TRIPOS>ATOM and create Atom object and store then in a list

			Parameters
			----------
			string_atom_ : str
				str of the @<TRIPOS>ATOM section

				
			Return
			------
			atoms : list
				list of Atom objects containing informations about each atoms (id, type, coordinates, ...)
		"""
		atoms = []

		line_atoms = string_atom_.split("\n")
		for line in line_atoms:
			if line != "":
				id_atom, name, x, y, z, type_atom, id_sub, name_sub, charge, *r = line.split()
				atom = Atom(id_atom, name, x, y, z, type_atom, id_sub, name_sub, charge)
				atoms += [atom]

		return atoms

	def extract_bonds(string_bond_):

		"""Extract informations from @<TRIPOS>BOND and create Bond object and store then in a list

			Parameters
			----------
			string_bond_ : str
				str of the @<TRIPOS>BOND section

				
			Return
			------
			bonds : list
				list of Bond objects containing informations about each atoms (ids, type)
		"""

		bonds = []

		line_bonds = string_bond_.split("\n")
		for line in line_bonds:
			if line != "":
				id_bond, id_atom1, id_atom2, type_bond, *r = line.split()
				bond = Bond(id_bond, id_atom1, id_atom2, type_bond)
				bonds += [bond]

		return bonds

	def extract_substructures(string_substructure_):

		"""Extract informations from @<TRIPOS>SUBSTRUCTURE and create Substructure object and store then in a list

			Parameters
			----------
			string_substructure_ : str
				str of the @<TRIPOS>SUBSTRUCTURE section

				
			Return
			------
			substructures : list
				list of Substructure objects containing informations about each atoms (ids, name)
		"""

		substructures = []

		line_bonds = string_substructure_.split("\n")
		for line in line_bonds:
			if line != "":
				id_sub, name_sub, id_root, *r = line.split()
				substructure = Substructure(id_sub, name_sub, id_root)
				substructures += [substructure]

		return substructures

	def load_mol2(path_mol2_):

		"""Open, parse a mol2 file and create a Mol2 object

			Parameters
			----------
			path_mol2_ : str
				mol2 file path

				
			Return
			------
			mol2 : Mol2
				mol2 object containing informations from mol2 file
		"""

		with open(path_mol2_, 'r') as file_mol2_hdl:
			string_mol2 = file_mol2_hdl.read()

		string_molecule = Loader.cut_mol2(string_mol2, "@<TRIPOS>MOLECULE")
		molecule = Loader.extract_molecule(string_molecule)

		string_atom = Loader.cut_mol2(string_mol2, "@<TRIPOS>ATOM")
		atoms = Loader.extract_atoms(string_atom)

		string_bond = Loader.cut_mol2(string_mol2, "@<TRIPOS>BOND")
		bonds = Loader.extract_bonds(string_bond)

		if string_mol2.find("@<TRIPOS>SUBSTRUCTURE") != -1:
			string_substructure = Loader.cut_mol2(string_mol2, "@<TRIPOS>SUBSTRUCTURE")
			substructures = Loader.extract_substructures(string_substructure)
		else:
			substructures = []

		name_mol2 = path_mol2_.split("/")[len(path_mol2_.split("/")) - 1]

		mol2 = Mol2(name_mol2, molecule, atoms, bonds, substructures)

		return mol2

class Merger:

	"""Class that contains fonction to merge few molecule
		Concatenate mol2 files and give a new id number to atoms, bonds, ...

		Methods
		-------
		get_atoms_and_bonds_to_merge(mol2_, atoms_to_delete_, atoms_to_modify_, bonds_to_modify_, id_atom_starting_point_, id_bond_starting_point_, id_sub_)
			return atoms and bonds used to create the new mol2 merged file

		merge_mol2(dico_information_mol2_, connexion_rules_, new_name_)
			merge two loaded mol2 file

	"""

	def get_atoms_and_bonds_to_merge(mol2_, atoms_to_delete_, atoms_to_modify_, bonds_to_modify_, id_atom_starting_point_, id_bond_starting_point_, id_sub_):

		"""This fonction delete, modify atoms and bonds given as parameter
			it reassigns id number of atoms and bonds

			e.g. string_section_ = "@<TRIPOS>ATOM" recover atoms of a mol2
			available : @<TRIPOS>MOLECULE, @<TRIPOS>BOND, @<TRIPOS>ATOM, @<TRIPOS>SUBSTRUCTURE

			Parameters
			----------
			mol2_ : Mol2
				mol2 loaded as Mol2 object (with Loader.load_mol2)
			atoms_to_delete_ : list
				list of id atom to delete
			atoms_to_modify_ : dict
				dictionary of id atom to modify. key : atom id number, value : new SYBYL atom type.
				e.g {1: "N.am"}
			bonds_to_modify_ : dict
				dictionary of id bond to modify. key : tuple of bond ids number, value : new SYBYL bond type.
				e.g. {(1, 2) : "2"}
			id_atom_starting_point_ : int
				starting point to reassign atom id number. Typically this is the atom number of the previous mol2 to merge
			id_bond_starting_point_ : int
				starting point to reassign bond id number. Typically this is the bond number of the previous mol2 to merge
			id_sub_ : int
				new substructure id number for corresponding mol2
				
			Return
			------
			modified_atoms, modified_bonds, table_old_new_id_atom, count_delete_atom, count_delete_bound : tuple
				modified_atoms : list of atoms with new id, SYBYL atom type, charge, ...
				modified_bonds : list of atoms with new id, SYBYL bond type, ...
				table_old_new_id_atom : dictionary connectivity table to make the correspondance between previous new and old id number atom
				count_delete_atom : count of deleted atom, used to adjust id number atom
				count_delete_bound : count of deleted bond, used to adjust id number bond
		"""

		modified_atoms = []
		modified_bonds = []

		table_old_new_id_atom = {}
		count_delete_atom = 0

		for atom in mol2_.atoms:
			if atom.id_atom in atoms_to_delete_:
				count_delete_atom += 1
			else:
				if atom.id_atom in atoms_to_modify_:

					### --- Setting new SYBYL atom type --- ###
					new_atom_type = atoms_to_modify_[atom.id_atom]
					atom.set_type_atom(new_atom_type)

					### --- Setting new charge --- ###
					### Has to be modified to fit the new atom type ###
					### There is not a problem for the moment (Amide, Sulfoam. and Benzoxazole) ###
					atom.set_charge(0.0000)

				old_id_atom = atom.id_atom
				new_id_atom = int(atom.id_atom) - count_delete_atom + id_atom_starting_point_
				table_old_new_id_atom[old_id_atom] = new_id_atom
				atom.set_id_atom(new_id_atom)

				atom.set_id_sub(id_sub_)
				atom.set_name_sub('BB{}'.format(id_sub_))
				modified_atoms.append(atom)

		count_delete_bound = 0
		for bond in mol2_.bonds:
			if bond.id_atom1 in atoms_to_delete_ or bond.id_atom2 in atoms_to_delete_:
				count_delete_bound += 1
			else:
				if (bond.id_atom1, bond.id_atom2) in bonds_to_modify_:
					bond.set_type_bond(bonds_to_modify_[(bond.id_atom1, bond.id_atom2)])
				elif (bond.id_atom2, bond.id_atom1) in bonds_to_modify_:
					bond.set_type_bond(bonds_to_modify_[(bond.id_atom2, bond.id_atom1)])

				bond.set_id_atom1(table_old_new_id_atom[bond.id_atom1])
				bond.set_id_atom2(table_old_new_id_atom[bond.id_atom2])

				bond.set_id_bond(bond.id_bond - count_delete_bound + id_bond_starting_point_)
				modified_bonds.append(bond)

		return modified_atoms, modified_bonds, table_old_new_id_atom, count_delete_atom, count_delete_bound

	def merge_mol2(dico_information_mol2_, connexion_rules_, new_name_):
		
		"""This function merged few mol2 loaded as mol2 object
			it used the function get_atoms_and_bonds_to_merge() to reassign atom, bond of each mol2
			


			Parameters
			----------
			dico_information_mol2_ : dict
				key : str of id of loaded mol2 ("1", "2", ...), values : tuple containing (loaded mol2, atoms_to_delete_, atoms_to_modify_ and bonds_to_modify_)
				e.g. { "1" : Mol2, [1], {2 : "C.3"}, {(2,3) : "2"}, "2" : Mol2, [5], {6 : "C.3"}, {(6,7) : "1"}}
			connexion_rules_ : list
				list of tuple to know which atoms to connect between mol2 to merge
				e.g. [(("1", 1), ("2", 3), "am")]
				Loaded mol2 "1" with atom 1 connect loaded molecule "2" with atom 3 by amide bond
			new_name_ : dict
				
			Return
			------
			modified_atoms, modified_bonds, table_old_new_id_atom, count_delete_atom, count_delete_bound : tuple
				modified_atoms : list of atoms with new id, SYBYL atom type, charge, ...
				modified_bonds : list of atoms with new id, SYBYL bond type, ...
				table_old_new_id_atom : dictionary connectivity table to make the correspondance between previous new and old id number atom
				count_delete_atom : count of deleted atom, used to adjust id number atom
				count_delete_bound : count of deleted bond, used to adjust id number bond
		"""

		new_atoms = []
		new_bonds = []
		new_substructures = []

		dico_table_old_new_id_atom = {}
		id_atom_starting_point = 0
		id_bond_starting_point = 0
		for key, information in dico_information_mol2_.items():

			### Key : Loaded mol2 id number ###
			### Use to define substructure id ###
			mol2, atoms_to_delete, atoms_to_modify, bonds_to_modify = information
			modified_atoms, modified_bonds, table_old_new_id_atom, count_delete_atom, count_delete_bound = Merger.get_atoms_and_bonds_to_merge(mol2, atoms_to_delete, atoms_to_modify, bonds_to_modify, 
																																			id_atom_starting_point, id_bond_starting_point, key)
			id_atom_starting_point += len(modified_atoms)
			id_bond_starting_point += len(modified_bonds)
			new_atoms.extend(modified_atoms)
			new_bonds.extend(modified_bonds)
			dico_table_old_new_id_atom[key] = table_old_new_id_atom

		for rules in connexion_rules_:
			key1, id_atom_connectable_1 = rules[0]
			key2, id_atom_connectable_2 = rules[1]
			type_bond_connexion = rules[2]
			table_old_new_id_atom_key1 = dico_table_old_new_id_atom[key1] # We get the tag connexion table corresping to the mol2 key
			table_old_new_id_atom_key2 = dico_table_old_new_id_atom[key2]

			new_id_atom_connectable_1 = table_old_new_id_atom_key1[id_atom_connectable_1]
			new_id_atom_connectable_2 = table_old_new_id_atom_key2[id_atom_connectable_2]
			id_bond_connexion = len(new_bonds) + 1

			connexion_bond = Bond(id_bond_connexion, new_id_atom_connectable_1, new_id_atom_connectable_2, type_bond_connexion)
			new_bonds.append(connexion_bond)

		### Substructures are added not dynamically ###
		### It should be modified to connect more than 2 mol2 ###
		new_substructures.append(Substructure(key1, 'BB{}'.format(key1), new_id_atom_connectable_1))
		new_substructures.append(Substructure(key2, 'BB{}'.format(key2), new_id_atom_connectable_2))

		molecule = Molecule(new_name_, len(new_atoms), len(new_bonds), len(new_substructures), "SMALL", "USER_CHARGES")
		merged_mol2 = Mol2(new_name_, molecule, new_atoms, new_bonds, new_substructures)

		return merged_mol2