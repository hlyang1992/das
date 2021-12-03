import itertools
import logging
import os
import random
import shutil
import sys
from copy import deepcopy

import numpy as np
from path import Path
from pymtp.core import MTPCalactor, PyConfiguration

from das import AtomicDataset, Structure
from das.core import Worker
from das.dataset import atomic_masses
from das.utils import AutoRepr, natural_key, get_template


def get_force_ambiguity(force_set):
    force_set_sub_mean = force_set - np.mean(force_set, axis=0)
    force_set_sub_mean_norm = np.sum(force_set_sub_mean ** 2, axis=2)
    ee_list = np.sqrt(np.mean(force_set_sub_mean_norm, axis=0))
    # print(f"max/min/mean: {np.max(ee_list):.3f} {np.min(ee_list):.3f} {np.mean(ee_list):.3f}")
    # return np.max(ee_list)
    return ee_list


def calc_ensemble_ambiguity(atoms_list, model_fns, unique_numbers):
    calc_set = [MTPCalactor(i) for i in model_fns]
    ee_list = []
    for i, atoms in enumerate(atoms_list):
        forces_set = []
        for calc in calc_set:
            c0 = PyConfiguration.from_ase_atoms(atoms)
            calc.calc(c0)
            forces_set.append(c0.force)
        forces_set = np.array(forces_set, dtype=np.float)
        ee = get_force_ambiguity(forces_set)

        atomic_numbers = atoms.numbers
        ee_list_single = np.zeros(len(unique_numbers))
        for ii, uu in enumerate(unique_numbers):
            filter_ele_index = atomic_numbers == uu
            # if len(filter_ele_index) ==0, element uu not in this conf
            ee_list_single[ii] = np.max(ee[filter_ele_index]) if len(filter_ele_index) > 0 else 0.0
        ee_list.append(ee_list_single)

    #  ee_list NxM, N=# of confs in training set, M=# of unique elements
    ee_list = np.array(ee_list)
    return ee_list


def lmp_vars_to_params_list(lmp_vars):
    if "press" not in lmp_vars.keys():
        lmp_vars["press"] = [None]
    params_list = []
    keys = list(lmp_vars.keys())
    values = list(lmp_vars.values())
    for item in itertools.product(*values):
        param = {k: v for k, v in zip(keys, item)}
        params_list.append(param)
    return params_list


@AutoRepr
class LmpAmbiguitySampler(Worker):
    def __init__(
        self,
        repeat=1,
        af_default=0.010,
        af_limit=0.200,
        af_failed=0.500,
        over_fitting_factor=1.1,
        max_number_confs=10,
        min_number_confs=1,
        number_ignore_confs=1,
        threshold_limit_factor=1.1,
        lmp_file="lmp_af_sampler.jinja2",
        lmp_rerun_file="lmp_af_sampler_rerun.jinja2",
        run_lmp_file="run_lmp_af_sampler.jinja2",
        kind="lmp_model_sampler",
        lmp_vars=None,
        bond_hierarchy=False,
        **kwargs,
    ):
        super().__init__(kind=kind, **kwargs)
        self.logger = logging.getLogger("das.sampler.lmp_ambiguity_sampler.LmpAmbiguitySampler")
        self.repeat = repeat

        # prepare lmp template
        lmp_file = self.find_fn(lmp_file, __file__)
        self.lmp_temp = get_template(lmp_file)

        lmp_rerun_file = self.find_fn(lmp_rerun_file, __file__)
        self.lmp_rerun_temp = get_template(lmp_rerun_file)

        run_lmp_file = self.find_fn(run_lmp_file, __file__)
        self.run_template = get_template(run_lmp_file)

        self.lmp_vars = lmp_vars
        self.low_factor = over_fitting_factor
        self.max_number_confs = max_number_confs
        self.min_number_confs = min_number_confs
        self.bond_hierachy = bond_hierarchy
        self.number_ignore_confs = number_ignore_confs
        self.threshold_limit_factor = threshold_limit_factor
        if self.bond_hierachy:
            if isinstance(af_default, list):
                if len(af_default) != self.unique_numbers:
                    self.logger.error("The setting of af_default is invalid")
                    sys.exit()
                self.af_default = af_default
            else:
                self.af_default = [af_default] * len(self.unique_numbers)

            if isinstance(af_limit, list):
                if len(af_limit) != self.unique_numbers:
                    self.logger.error("The setting of af_limit is invalid")
                self.af_limit = af_limit
            else:
                self.af_limit = [af_limit] * len(self.unique_numbers)

            if isinstance(af_failed, list):
                if len(af_failed) != self.unique_numbers:
                    self.logger.error("The setting of af_failed is invalid")
                self.af_failed = af_failed
            else:
                self.af_failed = [af_failed] * len(self.unique_numbers)
        else:
            self.af_default = af_default
            self.af_limit = af_limit
            self.af_failed = af_failed

    def prepare_run(self):
        prev_model_fns = self.driver.get_prev_model_fns()

        mass_dict = {i + 1: atomic_masses[ii] for i, ii in enumerate(self.unique_numbers)}
        params_list = lmp_vars_to_params_list(self.lmp_vars)

        atoms_lists = self.dataset.ase_atoms
        self.logger.info(f"{len(atoms_lists)} need lmp sampler")
        params_list *= int(self.repeat)
        self.logger.info(f"params list: {params_list}")
        self.logger.info(f"{len(params_list)} items in params list.")

        dump_fns = [f"force.{i}.nc" for i in range(len(prev_model_fns))]
        mlip_fns = [f"mlip_{i}.ini" for i in range(len(prev_model_fns))]
        dump_mlip_fns = [(i, j) for i, j in zip(dump_fns, mlip_fns)]

        lmp_cmd = self.machine.get_cmd("lmp_cmd")
        python_cmd = self.machine.get_cmd("python_cmd")
        if not lmp_cmd or not python_cmd:
            sys.exit()
        run_str = self.run_template.render(dump_mlip_fns=dump_mlip_fns, lmp_cmd=lmp_cmd, python_cmd=python_cmd)
        with open(self.run_input_dir / "sub_run.sh", "w") as f0:
            f0.write(run_str)

        for ii, atoms in enumerate(atoms_lists):
            for jj, params in enumerate(params_list):
                tmp_params = deepcopy(params)
                # if tmp_params.get("random", None) is None:
                tmp_params["random"] = np.random.randint(0, 10000000)
                tmp_path = self.run_input_dir / f"{ii}_{jj}"
                tmp_path.mkdir_p()
                lmp_data_fn = tmp_path / "data.lmp"
                lmp_input_fn = tmp_path / "in.lmp"
                lmp_rerun_input_fn = tmp_path / "in_rerun.lmp"
                c0 = Structure.from_atoms(atoms)
                c0.write_lmp_data(lmp_data_fn, unique_numbers=self.unique_numbers)
                lmp_str = self.lmp_temp.render(mass_dict=mass_dict, **tmp_params)
                lmp_rerun_str = self.lmp_rerun_temp.render(mass_dict=mass_dict, **tmp_params)
                with open(lmp_input_fn, "w") as f0:
                    f0.write(lmp_str)
                with open(lmp_rerun_input_fn, "w") as f0:
                    f0.write(lmp_rerun_str)

                shutil.copy2(Path(__file__).abspath().dirname() / "calc_ensemble_ambiguity.py", tmp_path)

                for kk, model in enumerate(prev_model_fns):
                    mlip_fn = tmp_path / mlip_fns[kk]
                    model_fn = model.name
                    with open(mlip_fn, "w") as f0:
                        f0.write(f"mtp-filename   {model_fn}\n")
                        f0.write("select      FALSE\n")
                    shutil.copy2(model, tmp_path)

    def do_run(self):
        before_dir = Path(os.getcwd()).abspath()
        self.machine.submit_and_wait(self.run_input_dir, self.driver)
        self.machine.download(
            self.run_finished_dir, files=["force.0.nc", "*.out"],
        )
        os.chdir(before_dir)

    def post_run(self):
        prev_select_conf = self.driver.get_prev_select_conf()
        prev_model_fns = self.driver.get_prev_model_fns()
        if prev_select_conf is None:
            af_adaptive = self.af_default
        else:
            prev_select_dataset = AtomicDataset.load(prev_select_conf)
            if len(prev_select_dataset) < 1:
                af_adaptive = self.af_default
            elif self.driver.status_record.iter_i == 0:
                af_adaptive = self.af_default
            else:
                pre_select_atoms = prev_select_dataset.ase_atoms
                # train_ee: (# of confs x # of elements)
                train_ee = calc_ensemble_ambiguity(pre_select_atoms, prev_model_fns, self.unique_numbers)
                af_adaptive = np.max(train_ee, axis=0) * self.low_factor
                if self.bond_hierachy:
                    # auto increase af_limit if need
                    for ii in range(len(self.unique_numbers)):
                        if af_adaptive[ii] > self.af_limit[ii]:
                            self.af_limit[ii] *= self.threshold_limit_factor
                            self.logger.info(
                                f"{af_adaptive[ii]} > {self.af_limit[ii]/self.threshold_limit_factor}, "
                                f"increase af_limit[{ii}] to {self.af_limit[ii]}"
                            )
                    for ii in range(len(self.unique_numbers)):
                        if af_adaptive[ii] > self.af_limit[ii]:
                            af_adaptive[ii] = self.af_limit[ii]
                else:
                    af_max = np.max(af_adaptive)
                    af_min = np.min(af_adaptive)
                    if af_max > af_min * 2:
                        self.logger.warning(
                            "The ambiguity of different elements is too large,"
                            "please consider setting bond_hierachy=True."
                        )
                    af_adaptive = np.max(af_adaptive)
                    if af_adaptive > self.af_limit:
                        self.af_limit *= self.threshold_limit_factor
                        self.logger.info(
                            f"{af_adaptive} > {self.af_limit/self.threshold_limit_factor}, increase "
                            f"af_limit to {self.af_limit}"
                        )
                    if af_adaptive > self.af_limit:
                        af_adaptive = self.af_limit
        self.logger.info(f"THE AMBIGUITY THRESHOLD: {af_adaptive}")

        selected_atoms = []
        statistical_numbers = np.zeros(2, dtype=np.int)
        sub_dirs = sorted(list(self.run_finished_dir.dirs()), key=natural_key)
        for path in sub_dirs:
            tmp = AtomicDataset(self.unique_numbers)
            tmp.read_from_file(path / "force.0.nc", filetype="nc")
            tmp_atoms = tmp.ase_atoms

            tmp_numbers = np.zeros(2, dtype=np.int)
            min_index = []
            # check if lost atoms
            is_lost_atoms = self.check_lost_atoms(path / "jobs_mlip_0.ini.out")
            if not self.bond_hierachy:
                af_fn = path / "af.out"
                ee = np.loadtxt(af_fn)
                ee = ee[: len(tmp_atoms)]
                selected_index = np.where(((ee[:, 0] > af_adaptive) & (ee[:, 0] < self.af_failed)))[0]
                accurate_index = np.where(ee[:, 0] <= af_adaptive)[0]
                tmp_numbers += [len(selected_index), len(accurate_index)]
                tmp_selected_atoms = [tmp_atoms[i] for i in selected_index]
                self.logger.info(f"STATISTICS select/accurate: {tmp_numbers} in {path.name}")
                no_failed = len(accurate_index) + len(selected_index) == len(ee)
                statistical_numbers += tmp_numbers
                if (len(tmp_selected_atoms) <= self.number_ignore_confs) and (not is_lost_atoms) and no_failed:
                    self.logger.info(
                        f"the number of selected confs {len(tmp_selected_atoms)} < {self.number_ignore_confs}, skipping."
                    )
                else:
                    min_index_ii = np.argmin(np.abs(ee[:, 0] - np.mean(ee[selected_index, 0])))
                    self.logger.info(f"selected index: {min_index_ii} af:{ee[min_index_ii, 0]}")
                    min_index.append(min_index_ii)
            else:
                tmp_selected_atoms = []
                if not (len(list(path.walk("af_*.out"))) > 0):
                    self.logger.error(f"I CAN'T FIND af_*.out IN PATH: {path}")
                    continue
                for ii in range(len(self.unique_numbers)):
                    af_fn = path / f"af_{ii + 1}.out"
                    # this element not in current conf
                    if not af_fn.exists():
                        continue
                    ee = np.loadtxt(af_fn)
                    ee = ee[: len(tmp_atoms)]
                    selected_index = np.where(((ee[:, 0] > af_adaptive[ii]) & (ee[:, 0] < self.af_failed[ii])))[0]
                    accurate_index = np.where(ee[:, 0] <= af_adaptive[ii])[0]
                    tmp_numbers = [len(selected_index), len(accurate_index)]
                    tmp_selected_atoms = [tmp_atoms[i] for i in selected_index]
                    self.logger.info(f"STATISTICS select/accurate(type={ii}): {tmp_numbers} in {path.name}")
                    no_failed = len(accurate_index) + len(selected_index) == len(ee)
                    statistical_numbers += tmp_numbers
                    if (len(tmp_selected_atoms) <= self.number_ignore_confs) and (not is_lost_atoms) and no_failed:
                        self.logger.info(
                            f"the number of selected confs {len(tmp_selected_atoms)} < {self.number_ignore_confs}, skipping."
                        )
                    else:
                        min_index_ii = np.argmin(np.abs(ee[:, 0] - np.mean(ee[selected_index, 0])))
                        self.logger.info(f"selected index: {min_index_ii} af: {ee[min_index_ii, 0]}")
                        min_index.append(min_index_ii)
            min_index = list(set(min_index))
            for ii in min_index:
                selected_atoms.append(tmp_atoms[ii])
                self.logger.info(f"add {ii} to selected atoms.")
        self.logger.info(f"STATISTICS select/accurate: {statistical_numbers}")
        if len(selected_atoms) > self.max_number_confs:
            selected_atoms = random.sample(selected_atoms, k=self.max_number_confs)
        elif len(selected_atoms) <= self.min_number_confs:
            self.logger.info(f"number of seleced atoms {len(selected_atoms)} < {self.min_number_confs}, CONVERVED!")
            selected_atoms = []
        self.logger.info(f"select {len(selected_atoms)} configurations.")
        self.new_atoms_lists = selected_atoms

    @staticmethod
    def check_lost_atoms(log_fn):
        with open(log_fn) as f0:
            for line in f0:
                if "Lost atoms" in line:
                    return True
        return False
