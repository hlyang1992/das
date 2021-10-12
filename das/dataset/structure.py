import numpy as np
import spglib as spg
from ase import Atoms
from ase.io import read as ase_read
from ase.io import write as ase_write

from das.dataset import atomic_masses, chemical_symbols
from das.dataset import atomic_numbers
from .box import Box


class StructureError(Exception):
    """
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    """

    pass


class Structure:
    def __init__(
        self, box, atomic_numbers, coords, coords_are_fraction=True, sort_flag=False, max_show_items=15, wrap=True,
    ):
        """
        Args:
            box (3x3 array): Each row should correspond to a boxvector.
            atomic_numbers （1*N array): atomic numbers of all atoms.
            coords (Nx3 array): list of fractional/cartesian coordinates of each species.
        """
        self.atomic_prop = {}  # global property
        self.global_prop = {}  # atomic property
        self._box = Box(box)
        self.atomic_numbers = atomic_numbers
        coords = np.array(coords)
        atomic_numbers = np.array(atomic_numbers, dtype=np.int)
        self._protected_atomic_prop = ["atomic_numbers", "frac_coords"]

        if len(atomic_numbers) != len(coords):
            raise StructureError(
                "The list of atomic species must be of the" " same length as the list of fractional" " coordinates."
            )
        if not coords_are_fraction:
            coords = self.box.cart_to_frac(coords)

        if wrap:
            coords = np.array(coords) % 1.0

        self.frac_coords = coords
        self._max_show_items = max_show_items
        if sort_flag:
            self.sort()

    def set_global_prop(self, name, value):
        self.global_prop[name] = value

    def get_global_prop(self, name):
        prop = self.global_prop.get(name, None)
        if prop is None:
            raise StructureError(f"property {name} is not set")
        else:
            return prop

    def set_atomic_prop(self, name, value):
        if name in self._protected_atomic_prop:
            raise StructureError(f"The atomic property {name} can't be changed.")
        value = np.array(value)
        if value.shape[0] != len(self):
            raise StructureError(f"The first of value{value.shape[0]} must be {len(self)}.")
        self.atomic_prop[name] = value

    def get_atomic_prop(self, name):
        prop = self.atomic_prop.get(name, None)
        if prop is None:
            raise StructureError(f"property {name} is not set")
        else:
            return prop

    def __len__(self):
        return len(self.atomic_numbers)

    @property
    def box(self):
        return self._box

    @property
    def cart_coords(self):
        return self.box.frac_to_cart(self.frac_coords)

    @cart_coords.setter
    def cart_coords(self, cart_coords):
        fcoords = self.box.cart_to_frac(cart_coords)
        self.frac_coords = fcoords

    @property
    def frac_coords(self):
        return self.atomic_prop["frac_coords"]

    @frac_coords.setter
    def frac_coords(self, frac_coords):
        self.atomic_prop["frac_coords"] = np.array(frac_coords)

    @property
    def atomic_numbers(self):
        return self.atomic_prop["atomic_numbers"]

    @property
    def atomic_symbols(self):
        return [chemical_symbols[i] for i in self.atomic_numbers]

    @property
    def structure(self):
        from pymatgen import Structure as StructureM

        return StructureM(self.box.matrix, self.atomic_numbers, self.frac_coords,)

    @classmethod
    def from_structure(cls, structure):
        atomic_numbers = structure.atomic_numbers
        coords = [site.frac_coords for site in structure]
        box = structure.lattice.matrix

        return cls(box, atomic_numbers, coords)

    @property
    def atoms(self):
        return Atoms(positions=self.cart_coords, numbers=self.atomic_numbers, cell=self.box.matrix, pbc=True,)

    @classmethod
    def from_atoms(cls, atoms):
        return cls(atoms.cell, atoms.numbers, atoms.positions, coords_are_fraction=False)

    @property
    def phonopy_atoms(self):
        from phonopy.structure.atoms import PhonopyAtoms

        return PhonopyAtoms(numbers=self.atomic_numbers, scaled_positions=self.frac_coords, cell=self.box.matrix)

    @classmethod
    def from_phonopy_atoms(cls, atoms):
        return cls(atoms.cell, atoms.numbers, atoms.scaled_positions, coords_are_fraction=True)

    @atomic_numbers.setter
    def atomic_numbers(self, numbers):
        self.atomic_prop["atomic_numbers"] = np.array(numbers, dtype=np.int)

    @staticmethod
    def mass_numbers(mass):
        return np.argmin(np.abs(mass - atomic_masses))

    @classmethod
    def from_lmp_dat_file(cls, filename, style="atomic", sort_by_id=True):
        a_in_lmp = ase_read(filename, format="lammps-data", sort_by_id=sort_by_id, style=style)
        # fix numbers using atomic mass
        numbers = np.array([cls.mass_numbers(i) for i in a_in_lmp.get_masses()], dtype=np.int)
        a_in_lmp.numbers = numbers

        return cls.from_atoms(a_in_lmp)

    @classmethod
    def from_file(cls, filename, **kwargs):
        atoms = ase_read(filename, **kwargs)
        return cls.from_atoms(atoms)

    def write(self, filename, format=None, sort_flag=False, **kwargs):
        if sort_flag:
            self.sort()
        atoms = self.atoms
        ase_write(filename, atoms, format=format, **kwargs)

    def write_poscar(self, filename, direct=True, sort_flag=False, **kwargs):
        if sort_flag:
            self.sort()
        atoms = self.atoms
        ase_write(filename, atoms, format="vasp", direct=direct, vasp5=True, **kwargs)

    @staticmethod
    def get_lmp_box(cell, positions):
        q, r = np.linalg.qr(cell.transpose())

        new_cell = cell @ q
        new_cell[0, 1] = 0
        new_cell[0, 2] = 0
        new_cell[1, 2] = 0

        new_positions = positions @ q
        return new_cell, new_positions

    def get_lmp_crytal(self):
        # cell, positions = self.get_lmp_box(self.box.matrix, self.cart_coords)
        # lmp_crytal.wrap()
        lmp_box_matrix = self.box.get_lmp_box()
        lmp_crytal = Structure(lmp_box_matrix, self.atomic_numbers, self.frac_coords)
        return lmp_crytal

    @staticmethod
    def unique_elements_to_numbers(unique_elements):
        return sorted([atomic_numbers[i] for i in unique_elements])

    def write_lmp_data(self, filename, unique_numbers=None, style="atomic"):
        lmp_structure = self.get_lmp_crytal()

        if unique_numbers is None:
            unique_numbers = list(set(lmp_structure.atomic_numbers))
        unique_numbers = sorted(unique_numbers)
        numbers = lmp_structure.atomic_numbers
        cell = lmp_structure.box.matrix
        positions = lmp_structure.cart_coords

        types = [unique_numbers.index(n) for n in numbers]
        lmp_str = "write by das.strcture.Structure. \n"
        lmp_str += f"{0.0:.15f} {cell[0, 0]:.15f} xlo xhi\n"
        lmp_str += f"{0.0:.15f} {cell[1, 1]:.15f} ylo yhi\n"
        lmp_str += f"{0.0:.15f} {cell[2, 2]:.15f} zlo zhi\n"
        lmp_str += f"{cell[1, 0]:.15f} {cell[2, 0]:.15f} {cell[2, 1]:.15f} xy xz yz\n\n"
        lmp_str += f"{len(numbers)} atoms\n"
        lmp_str += f"{len(unique_numbers)} atom types\n\n"
        lmp_str += "Masses\n\n"
        for i, ii in enumerate(unique_numbers):
            lmp_str += f"{i + 1} {atomic_masses[ii]:20.15f}\n"
        lmp_str += "\n"
        lmp_str += f"Atoms # {style}\n\n"
        for i in range(len(numbers)):
            lmp_str += f"{i + 1:10d} {types[i] + 1:4d} {positions[i, 0]:20.15f} "
            lmp_str += f"{positions[i, 1]:20.15f} {positions[i, 2]:20.15f}\n"

        with open(filename, "w") as f0:
            f0.write(lmp_str)

    @property
    def natoms(self):
        return len(self)

    @property
    def rho(self):
        return self.natoms / self.volume

    @property
    def volume(self):
        """
        Volume of the unit cell.
        """
        return self.box.volume

    @property
    def atomic_mass(self):
        return atomic_masses[self.atomic_numbers]

    def __str__(self):
        return repr(self)

    def __repr__(self):
        show_items = len(self) if len(self) <= self._max_show_items else self._max_show_items
        outs = "Sturcture Summary\n"
        outs += f"{self.box}\n"
        for i in range(show_items):
            tmp_string = "Atom {}: {} ({:.4f}, {:.4f}, {:.4f}) ".format(i, self.atomic_numbers[i], *self.cart_coords[i])
            tmp_string += "[{:.4f}, {:.4f}, {:.4f}]".format(*self.frac_coords[i])
            outs += tmp_string + "\n"
        if show_items < len(self):
            outs += f"....\n{len(self) - show_items} atoms not show, total {self.natoms} atoms\n"
        else:
            outs += f"total {self.natoms} atoms\n"
        return outs

    def index_coords(self, frac_coords, error=1e-4):
        frac_dd = self.frac_coords - frac_coords
        frac_dd %= 1.0
        cart_dd = self.box.frac_to_cart(frac_dd)
        dists = np.sum(cart_dd ** 2, axis=1)
        min_index = np.argmin(dists)
        return min_index if dists[min_index] <= error else None

    def translate_sites(self, indices, vector, coords_are_fraction=True, wrap=True):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            indices: Integer or List of site indices on which to perform the
                translation.
            vector: Translation vector for sites.
            coords_are_fraction(bool): Whether the vector corresponds to fractional or
                cartesian coordinates.
            wrap (bool): Whether new sites are transformed to unit
                cell
        """
        indices = self.indices(indices)
        if not coords_are_fraction:
            vector = self.box.cart_to_frac(vector)
        fcoords = self.frac_coords[indices] + vector
        if wrap:
            fcoords = fcoords % 1
        self.frac_coords[indices] = fcoords

    def copy(self):
        box = self.box.matrix.copy()
        numbers = self.atomic_numbers.copy()
        coords = self.frac_coords.copy()
        return self.Structure(box, numbers, coords)

    @staticmethod
    def indices(indices):
        if not hasattr(indices, "__iter__"):
            indices = [indices]
        return indices

    def get_neighbors_in_shell(self, origin, r1=0, r2=0, coords_are_fraction=True):
        if coords_are_fraction:
            origin = self.box.frac_to_cart(origin)
        s0 = self.structure
        tmp_dat = s0.get_neighbors_in_shell(origin, (r1 + r2) * 0.5, (r2 - r1) * 0.5, include_index=True)
        neighbors = {
            "index": np.array([i[-1] for i in tmp_dat], dtype=np.int),
            "dist": np.array([i[-2] for i in tmp_dat], dtype=np.float),
        }
        return neighbors

    def get_all_neighbors(self, r):
        s0 = self.structure
        neighbor_info = s0.get_all_neighbors(r, include_index=True, include_image=False)
        neigh_index = [[i.index for i in ii] for ii in neighbor_info]
        neigh_dist = [[i.nn_distance for i in ii] for ii in neighbor_info]
        return neigh_index, neigh_dist

    def get_distance(self, frac_coord1, frac_coord2, include_image=False):
        dists_image = self.structure.lattice.get_distance_and_image(frac_coord1, frac_coord2)
        return dists_image if include_image else dists_image[0]

    def sort(self):
        """sort atoms according to atomic numbers"""
        sort_index = np.argsort(self.atomic_numbers)
        for key in self.atomic_prop:
            self.atomic_prop[key] = self.atomic_prop[key][sort_index]

    def set_coords(self, coords, index=None, coords_are_fraction=True):
        coords = np.array(coords)
        if index is None:
            index = np.arange(coords.shape[0])
        if not coords_are_fraction:
            coords = self.box.cart_to_frac(coords)
        self.frac_coords[index] = coords

    def wrap(self):
        self.frac_coords %= 1.0

    def get_primitive_structure(self, symprec=1e-5):
        p_cell = spg.find_primitive(self.spg_cell, symprec=symprec)
        return Structure.from_spg_cell(p_cell)

    def get_conventional_structure(self, symprec=1e-5):
        c_cell = spg.refine_cell(self.spg_cell, symprec=symprec)
        return Structure.from_spg_cell(c_cell)

    def get_supercell(self, supercell_matrix):
        supercell_matrix = self.shape_supercell_matrix(supercell_matrix)
        from phonopy.structure import cells

        supercell = cells.get_supercell(self.phonopy_atoms, supercell_matrix)
        return Structure.from_phonopy_atoms(supercell)

    @staticmethod
    def shape_supercell_matrix(smat):
        if smat is None:
            _smat = np.eye(3, dtype="intc", order="C")
        elif len(np.ravel(smat)) == 3:
            _smat = np.diag(smat)
        elif len(np.ravel(smat)) == 9:
            _smat = np.reshape(smat, (3, 3))
        else:
            msg = "supercell_matrix shape has to be (3,) or (3, 3)"
            raise RuntimeError(msg)
        return _smat

    @property
    def spg_cell(self):
        spg_cell = (self.box.matrix, self.frac_coords, self.atomic_numbers)
        return spg_cell

    @classmethod
    def from_spg_cell(cls, spg_cell):
        return cls(box=spg_cell[0], coords=spg_cell[1], atomic_numbers=spg_cell[2])

    def get_symmetry_data(self, symprec=1e-5):
        # https://spglib.github.io/spglib/python-spglib.html#get-symmetry-dataset
        symmetry_data = spg.get_symmetry_dataset(self.spg_cell, symprec=symprec)
        return symmetry_data

    def show_symmetry_info(self, symprec=1e-5):
        sym_data = self.get_symmetry_data(symprec=symprec)
        primitive_lattice = sym_data["primitive_lattice"]
        res = ""
        res += f"space group: {sym_data['number']}({sym_data['international']})\n"
        res += "primitive lattice:\n"
        res +="      A : " + " ".join(f"{i:.6f}" for i in primitive_lattice[0]) +"\n"
        res +="      B : " + " ".join(f"{i:.6f}" for i in primitive_lattice[1]) +"\n"
        res +="      C : " + " ".join(f"{i:.6f}" for i in primitive_lattice[2]) +"\n"
        return res
