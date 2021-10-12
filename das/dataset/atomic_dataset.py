import json
import logging
import pickle
import time
import uuid
from itertools import islice

import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read as ase_read
from path import Path

from das.dataset import atomic_numbers
from das.utils import glob_filename
from .structure import Structure


# TODO: refactor read io, move io to Structure


def read(fn, **kwargs):
    if isinstance(fn, Path):
        fn = str(fn.abspath())
    return ase_read(fn, **kwargs)


def _read_uid_tag_file(fn):
    tag_file = fn.dirname() / "tag"
    module_logger.debug(f"read tag file {tag_file}")
    if tag_file.exists():
        with open(tag_file) as f0:
            line = f0.readline()
            uid = line.strip()
        return uid
    else:
        return str(uuid.uuid4())


module_logger = logging.getLogger("das.dataset.atomic_dataset")


def length_angle_to_box(a, b, c, alpha, beta, gamma):
    lx = a
    xy = b * np.cos(gamma)
    xz = c * np.cos(beta)
    ly = np.sqrt(np.abs(b ** 2 - xy ** 2))
    yz = (b * c * np.cos(alpha) - xy * xz) / ly
    lz = np.sqrt(np.abs(c ** 2 - xz ** 2 - yz ** 2))
    cell = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])
    return cell


class MyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "mlp_sampler.dataset.atomic_dataset" and name == "AtomicDataset":
            return AtomicDataset
        else:
            return pickle.Unpickler.find_class(self, module, name)


class AtomicDataset:
    def __init__(self, unique_numbers):
        self.logger = logging.getLogger("das.dataset.atomic_dataset.AtomicDataset")
        self._data = None
        self.unique_numbers = tuple(sorted(unique_numbers))

    @classmethod
    def from_unique_elements(cls, unique_elements):
        unique_numbers = Structure.unique_elements_to_numbers(unique_elements)
        return cls(unique_numbers)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

    def __len__(self):
        return len(self._data) if self._data is not None else 0

    def save_json(self, fn):
        if self.data:
            data_dict = self.data.to_dict()
            columns = list(data_dict.keys())
            n_items = len(data_dict[columns[0]])

            to_save_data = [{k: np.array(data_dict[k][i]).tolist() for k in columns} for i in range(n_items)]
        else:
            to_save_data = []

        dumps_data = {"unique_numbers": self.unique_numbers, "data": to_save_data}
        with open(fn, "w", encoding="utf-8") as f0:
            return json.dump(dumps_data, f0, ensure_ascii=False, indent=4)

    @classmethod
    def load_json(cls, fn):
        with open(fn, encoding="utf-8") as f0:
            data_dict = json.load(f0)
        unique_numbers = data_dict["unique_numbers"]
        tmp_data = data_dict["data"]
        dataset = cls(unique_numbers)

        if tmp_data:
            new_data = {}
            for i, item in enumerate(tmp_data):
                for k, v in item.items():
                    if new_data.get(k, None) is None:
                        new_data[k] = {}
                    new_data[k][i] = np.array(v)
            dataset.data = pd.DataFrame(new_data)
        return dataset

    def save(self, fn):
        pickle.dump(self, open(fn, "wb"))

    @staticmethod
    def load(fn):
        cls = MyUnpickler(open(fn, "rb")).load()
        return cls

    def _get_reader(self, filetype):
        if filetype == "ase":
            reader = self.read_ase_file
        elif filetype == "mtp":
            reader = self.read_mtp_file
        elif filetype == "dump":
            reader = self.read_dump_file
        elif filetype == "nc":
            reader = self.read_lmp_nc
        else:
            raise Exception("Unknown filetype")
        return reader

    def read_from_file(self, filename, filetype="ase"):
        read_file = self._get_reader(filetype)
        self.read_ase_atoms(read_file(filename))

    def read_from_directories(self, data_paths, filename, is_glob=True, processes=None, filetype="ase"):
        if is_glob:
            all_dat_files = list(glob_filename(data_paths, filename))
        else:
            all_dat_files = [Path(data_paths[0]) / filename]
        self.logger.info(f"reading from files: {all_dat_files}")
        read_file = self._get_reader(filetype)
        # processes = multiprocessing.cpu_count() if processes is None else processes
        # with multiprocessing.Pool(processes=processes) as pool:
        #     all_ase_atoms = pool.map(read_file, all_dat_files)
        all_ase_atoms = [read_file(ii) for ii in all_dat_files]
        all_ase_atoms = [j for i in all_ase_atoms for j in i]

        self.read_ase_atoms(all_ase_atoms)

    def read_ase_atoms(self, atoms_list):
        atoms_frame = self.ase_atoms_to_atomic_dict(atoms_list)
        self.add_to_database(atoms_frame)

    @staticmethod
    def database_to_ase_atoms(pd_data):
        atomic_dict_list = pd_data.to_dict(orient="records")
        atoms_list = [AtomicDataset.atomic_dict_to_ase_atoms(i) for i in atomic_dict_list]
        return atoms_list

    @property
    def ase_atoms(self):
        if self.data is None:
            return []
        atomic_dict_list = self.data.to_dict(orient="records")
        atoms_list = [AtomicDataset.atomic_dict_to_ase_atoms(i) for i in atomic_dict_list]
        return atoms_list

    def get_atoms_by_uid(self, uid):
        res = self.data[self.data.uid == uid].to_dict(orient="records")
        if not res:
            return None
        else:
            return AtomicDataset.atomic_dict_to_ase_atoms(res[0])

    @staticmethod
    def atomic_dict_to_ase_atoms(a0_dict):
        cell = a0_dict["cell"]
        numbers = a0_dict["numbers"]
        positions = a0_dict["positions"]
        atoms = Atoms(numbers=numbers, positions=positions, cell=cell, pbc=True)
        atoms.info["uid"] = a0_dict["uid"]

        for key in ["forces"]:
            val = a0_dict[key]
            if val is not None:
                atoms.arrays[key] = val

        val = a0_dict["stress"]
        if val is not None:
            atoms.info["stress"] = val

        energy = a0_dict["energy"]
        if energy is not None:
            atoms.info["energy"] = energy
        return atoms

    @staticmethod
    def ase_atoms_to_atomic_dict(atoms_list):
        atoms_frame = []
        for atoms in atoms_list:
            cell = atoms.cell.array
            positions = atoms.positions
            numbers = atoms.numbers
            try:
                energy = atoms.get_potential_energy()
            except:
                energy = atoms.info.get("energy", None)
            try:
                forces = atoms.get_forces()
            except:
                forces = atoms.arrays.get("forces", None)
            try:
                stress = -atoms.get_stress() * atoms.get_volume()
            except:
                stress = atoms.info.get("stress", None)
            uid = atoms.info.get("uid", str(uuid.uuid4()))

            atoms_frame.append(
                {
                    "cell": cell,
                    "positions": positions,
                    "forces": forces,
                    "uid": uid,
                    "stress": stress,
                    "numbers": numbers,
                    "energy": energy,
                }
            )
        return atoms_frame

    def add_to_database(self, list_of_atoms):
        # check unique atoms in unique_numbers
        for item in list_of_atoms:
            numbers = set(item["numbers"])
            if not set(self.unique_numbers).issuperset(numbers):
                raise ValueError(f"unique numbers not equal, me:{self.unique_numbers}, other:{numbers}")

        if list_of_atoms is None or not list_of_atoms:
            self.logger.info(f"add 0 configurations, list_of_atoms={list_of_atoms}")
            return 0
        p0 = pd.DataFrame(list_of_atoms)
        if self._data is None:
            self._data = p0
        else:
            self._data = pd.concat([self._data, p0], ignore_index=True)
        self.logger.info(f"add {len(list_of_atoms)} configurations, total {len(self)} confs")

    def __add__(self, other):
        if isinstance(other, str):
            other = self.load(other)
        if other.unique_numbers != self.unique_numbers:
            raise ValueError(f"unique numbers not equal, me:{self.unique_numbers}, other:{other.unique_numbers}")

        data_other = other.data
        if data_other is None:
            return self.copy()
        else:
            new_dataset = self.__class__(unique_numbers=self.unique_numbers)
            new_dataset.data = pd.concat([self._data, data_other], ignore_index=True)
        return new_dataset

    def drop_duplicates(self, keep="first", ignore_index=True, inplace=True):
        if inplace:
            self._data.drop_duplicates(subset="uid", keep=keep, ignore_index=ignore_index, inplace=True)
            return self.copy()
        else:
            new_dataset = self.__class__(unique_numbers=self.unique_numbers)
            new_dataset.data = self._data.drop_duplicates(subset="uid", keep=keep, ignore_index=ignore_index)
        return new_dataset

    def copy(self):
        new_dataset = self.__class__(self.unique_numbers)
        if self._data is None:
            return new_dataset
        new_dataset.data = self.data.copy()
        return new_dataset

    @staticmethod
    def copy_ase_atoms(atoms, keep_info=False):
        new_atoms = atoms.copy()
        if keep_info:
            return new_atoms
        new_arrays = {}
        for key in new_atoms.arrays.keys():
            if key in ["numbers", "positions"]:
                new_arrays[key] = new_atoms.arrays[key]
        new_atoms.arrays = new_arrays
        new_atoms.info = {}
        return new_atoms

    def clear(self):
        self.logger.info("clear dataset")
        self.data = None

    def read_ase_file(self, fn, index=":"):
        self.logger.info(f"Reading from {fn}")
        fn = Path(fn)
        try:
            if index is None:
                atoms_list = [read(fn)]
            else:
                atoms_list = read(fn, index=index)
                if not isinstance(atoms_list, list):
                    atoms_list = [atoms_list]
            # read tag file
            if len(atoms_list) == 1:
                if atoms_list[0].info.get("uid", None) is None:
                    atoms_list[0].info["uid"] = _read_uid_tag_file(fn)
            for atoms in atoms_list:
                atoms.info["read_from"] = fn
            return atoms_list
        except Exception as e:
            # raise e
            self.logger.warning(f"Warning: read {fn} error: {e}")
            return None

    def parse_mtp_conf(self, lines):
        def index(name):
            for i, ii in enumerate(lines):
                ii = ii.lower()
                if ii.startswith(name.lower()):
                    return i
            return None

        size_i = index("Size")
        supercell_i = index("Supercell")
        energy_i = index("Energy")

        stress_i = index("PlusStress")
        atomdata_i = index("AtomData")

        natoms = int(lines[size_i + 1])
        cell = np.array([line.split() for line in lines[supercell_i + 1 : supercell_i + 4]])

        # read atoms data
        atoms_data_key = (lines[atomdata_i].split(":")[1]).split()
        atoms_data = {key: [] for key in atoms_data_key}

        for line in lines[atomdata_i + 1 : atomdata_i + natoms + 1]:
            for key, val in zip(atoms_data_key, line.split()):
                atoms_data[key].append(val)
        atomic_type = [int(i) for i in atoms_data["type"]]
        numbers = np.array([self.unique_numbers[i] for i in atomic_type], dtype=np.int)
        positions = np.zeros((natoms, 3))
        for i, key in enumerate(["cartes_x", "cartes_y", "cartes_z"]):
            positions[:, i] = [float(ii) for ii in atoms_data[key]]

        atoms = Atoms(cell=cell, positions=positions, numbers=numbers)
        forces = np.zeros((natoms, 3))
        for i, key in enumerate(["fx", "fy", "fz"]):
            if key in atoms_data_key:
                forces[:, i] = [float(ii) for ii in atoms_data[key]]
        atoms.arrays["forces"] = forces

        # read stress
        if stress_i is not None:
            stress_key = (lines[stress_i].split(":")[1]).split()
            stress_data = lines[stress_i + 1]
            # stress_data = np.array([float(i) for i in lines[stress_i+1].strip().split()])

            stress_data = np.array([float(i) for i in (lines[stress_i + 1].strip()).split()])
            stress_dict = {key: val for key, val in zip(stress_key, stress_data)}
            stress_array = np.array([float(stress_dict[i]) for i in "xx yy zz yz xz xy".split()])
            atoms.info["stress"] = stress_array
        if energy_i:
            energy = float(lines[energy_i + 1])
            atoms.info["energy"] = energy

        # read feature uid
        uid = str(uuid.uuid4())
        for line in lines:
            line = line.lower()
            if line.startswith("feature"):
                line_s = line.split()
                if line_s[1] == "uid":
                    uid = line_s[2]
        atoms.info["uid"] = uid
        return atoms

    def read_mtp_file(self, fn):
        self.logger.info(f"start read from {fn}")
        atoms_list = []
        with open(fn) as f0:
            for line in f0:
                line = line.strip()
                if line.startswith("BEGIN_CFG"):
                    single_conf = []
                    # start_conf = True
                elif line.startswith("END_CFG"):
                    # print(single_conf)
                    atoms_list.append(self.parse_mtp_conf(single_conf))
                elif line:
                    single_conf.append(line)
        # print(atoms_list)
        return atoms_list

    def read_number_of_atoms(self, fn):
        with open(fn) as f0:
            read_number_flag = False
            for line in f0:
                if "NUMBER OF ATOMS" in line:
                    read_number_flag = True
                elif read_number_flag:
                    n_atoms = int(line)
                    break
        return n_atoms

    def paser_lmp_dump(self, lines):
        lines = list(lines)
        # only read cubic cell
        bounds = np.array([[float(jj) for jj in ii.split()[:2]] for ii in lines[5:8]])
        # need to read tilt factor
        # https://lammps.sandia.gov/doc/Howto_triclinic.html
        if "xy" in lines[4] or "xz" in lines[4] or "yz" in lines[4]:
            xy, xz, yz = [float(ii.split()[2]) for ii in lines[5:8]]
        else:
            xy, xz, yz = 0, 0, 0
        xlo_bound = bounds[0, 0]
        ylo_bound = bounds[1, 0]
        zlo_bound = bounds[2, 0]
        xhi_bound = bounds[0, 1]
        yhi_bound = bounds[1, 1]
        zhi_bound = bounds[2, 1]
        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        cell = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])

        atomic_types = []
        positions = []
        forces = []
        engs = []

        for line in lines[9:]:
            ls = line.strip().split()
            atomic_types.append(int(ls[1]))
            positions.append([float(i) for i in ls[2:5]])
            forces.append([float(i) for i in ls[5:8]])
            engs.append(float(ls[-1]))
        atomic_types = np.array(atomic_types)
        positions = np.array(positions)
        forces = np.array(forces)
        engs = np.array(engs)

        numbers = np.array([self.unique_numbers[int(i) - 1] for i in atomic_types], dtype=np.int)

        atoms = Atoms(cell=cell, positions=positions, numbers=numbers)
        atoms.info["energy"] = engs.sum()
        atoms.arrays["forces"] = forces
        atoms.arrays["atomic_energy"] = engs
        return atoms

    def read_dump_file(self, fn):
        self.logger.info(f"read from {fn}")
        atoms_list = []
        n_atoms = self.read_number_of_atoms(fn)
        n = 0
        with open(fn) as f0:
            while True:
                f0.seek(0)
                n_sagment = 8 + 1 + n_atoms
                start = n * n_sagment
                end = (n + 1) * n_sagment
                lines = list(islice(f0, start, end))
                if not lines:
                    break
                atoms_list.append(self.paser_lmp_dump(lines))
                self.logger.debug(f"{len(atoms_list)} confs read.")
                # print(f"{len(atoms_list)} confs read.")
                n += 1
        return atoms_list

    def read_lmp_nc(self, fn):
        from netCDF4 import Dataset

        atoms_list = []
        rootgrp = Dataset(fn)
        types_set = np.array(rootgrp["type"][...].data, dtype=np.int)
        positions_set = np.array(rootgrp["coordinates"][...].data, dtype=np.float)
        cell_angles_set = np.array(rootgrp["cell_angles"][...].data, dtype=np.float)
        cell_angles_set = np.deg2rad(cell_angles_set)
        cell_lengths_set = np.array(rootgrp["cell_lengths"][...].data, dtype=np.float)
        # filter data
        ok_index = np.where(np.sum((types_set < 1) | (types_set > len(self.unique_numbers)), axis=1) == 0)[0]
        types_set = types_set[ok_index]
        positions_set = positions_set[ok_index]
        cell_angles_set = cell_angles_set[ok_index]
        cell_lengths_set = cell_lengths_set[ok_index]
        n_frames = types_set.shape[0]

        forces_flag = eng_flag = False
        if "forces" in rootgrp.variables:
            forces_set = np.array(rootgrp["forces"][...].data, dtype=np.float)
            forces_set = forces_set[ok_index]
            forces_flag = True
        if "c_pe" in rootgrp.variables:
            atomic_engs_set = np.array(rootgrp["c_pe"][...].data, dtype=np.float)
            atomic_engs_set = atomic_engs_set[ok_index]
            eng_flag = True

        for ii in range(n_frames):
            cell = length_angle_to_box(*cell_lengths_set[ii], *cell_angles_set[ii])
            positions = positions_set[ii]
            types = types_set[ii]
            numbers = [self.unique_numbers[t - 1] for t in types]
            atoms = Atoms(cell=cell, positions=positions, numbers=numbers)
            if eng_flag:
                atoms.arrays["atomic_energy"] = atomic_engs_set[ii]
                atoms.info["energy"] = np.sum(atoms.arrays["atomic_energy"])
            if forces_flag:
                atoms.arrays["forces"] = forces_set[ii]
            atoms_list.append(atoms)
        return atoms_list

    def format_mtp(self, index):
        row_data = self.data.iloc[index]
        # basic info
        numbers = row_data["numbers"]
        positions = row_data["positions"]
        cell = row_data["cell"]

        # extro info
        uid = row_data["uid"]
        types = np.array([self.unique_numbers.index(i) for i in numbers], dtype=np.int)

        struct_str = "BEGIN_CFG\n"
        struct_str += f"Size\n {len(numbers)}\n"
        struct_str += f"SuperCell\n"
        struct_str += f"{cell[0, 0]:20.15f} {cell[0, 1]:20.15f} {cell[0, 2]:20.15f}\n"
        struct_str += f"{cell[1, 0]:20.15f} {cell[1, 1]:20.15f} {cell[1, 2]:20.15f}\n"
        struct_str += f"{cell[2, 0]:20.15f} {cell[2, 1]:20.15f} {cell[2, 2]:20.15f}\n"

        # atom data
        forces = row_data["forces"]
        if forces is not None:
            struct_str += f"AtomData: id type cartes_x cartes_y cartes_z fx fy fz\n"
            for ii in range(len(numbers)):
                struct_str += (
                    f"{ii + 1:10d} {types[ii]:3d} "
                    f"{positions[ii][0]:20.15f} {positions[ii][1]:20.15f} {positions[ii][2]:20.15f} "
                    f"{forces[ii][0]:20.15f} {forces[ii][1]:20.15f} {forces[ii][2]:20.15f}\n"
                )
        else:
            struct_str += f"AtomData: id type cartes_x cartes_y cartes_z\n"
            for ii in range(len(numbers)):
                struct_str += (
                    f"{ii + 1:10d} {types[ii]:3d} "
                    f"{positions[ii][0]:20.15f} {positions[ii][1]:20.15f} {positions[ii][2]:20.15f} \n"
                )

        energy = row_data["energy"]
        if (energy is not None) and (not np.isnan(energy)):
            struct_str += f"Energy\n{energy:.15f}\n"
        stress = row_data["stress"]
        if stress is not None:
            struct_str += f"PlusStress: xx  yy  zz  yz  xz xy\n"
            struct_str += f"{stress[0]:.15f} {stress[1]:.15f} {stress[2]:.15f} "
            struct_str += f"{stress[3]:.15f} {stress[4]:.15f} {stress[5]:.15f}\n"
        struct_str += f"Feature uid {uid}\n"
        struct_str += f"END_CFG\n\n"
        return struct_str

    def write_mtp_data(self, fn, index_to_write=None, processes=None, step=1):
        t0 = time.time()
        if self.data is None:
            with open(fn, "w") as f:
                pass
            return None
        self.logger.info("start writing confs.")
        if index_to_write is None:
            index_set = range(0, len(self.data), step)

        if isinstance(index_to_write, int):
            if index_to_write < len(self):
                index_set = np.random.choice(range(len(self)), size=index_to_write, replace=False)
            else:
                index_set = range(len(self.data))
        elif isinstance(index_to_write, list):
            index_set = index_to_write
        # with multiprocessing.Pool(processes=processes) as pool:
        #     format_str = pool.map(self.format_mtp, index_set)
        format_str = (self.format_mtp(i) for i in index_set)
        with open(fn, "w") as f:
            for item in format_str:
                f.write(item)
        t1 = time.time()
        self.logger.info(f"Writing {len(index_set)} of {len(self)} confs.")
        self.logger.info(f"{t1 - t0:.4f} seconds to write.")

    @staticmethod
    def unique_elements_to_numbers(unique_elements):
        return sorted([atomic_numbers[i] for i in unique_elements])
