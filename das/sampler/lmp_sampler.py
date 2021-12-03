import logging
import os
import shutil
import sys

import numpy as np
from path import Path

from das import AtomicDataset
from das import Structure
from das.core import Worker
from das.dataset import atomic_masses
from das.sampler.lmp_af_sampler import lmp_vars_to_params_list
from das.utils import AutoRepr, get_template


@AutoRepr
class LmpSampler(Worker):
    def __init__(
        self,
        repeat=1,
        kind="lmp_sampler",
        lmp_vars=None,
        run_lmp_file="run_lmp_sampler.jinja2",
        lmp_file="lmp_sampler.jinja2",
        **kwargs,
    ):
        super().__init__(kind=kind, **kwargs)
        self.logger = logging.getLogger("das.sampler.lmp_sampler.LmpSampler")
        self.repeat = repeat
        self.lmp_vars = lmp_vars

        lmp_file = self.find_fn(lmp_file, __file__)
        self.lmp_temp = get_template(lmp_file)

        run_lmp_file = self.find_fn(run_lmp_file, __file__)
        self.run_template = get_template(run_lmp_file)

    def prepare_run(self):
        prev_model_fn = self.driver.get_prev_model_fn()
        mass_dict = {i + 1: atomic_masses[ii] for i, ii in enumerate(self.unique_numbers)}

        params_list = lmp_vars_to_params_list(self.lmp_vars)

        atoms_lists = self.dataset.ase_atoms
        self.logger.info(f"{len(atoms_lists)} need lmp sampler")
        params_list *= int(self.repeat)
        self.logger.info(f"params list: {params_list}")
        self.logger.info(f"{len(params_list)} items in params list.")

        lmp_cmd = self.machine.get_cmd("lmp_cmd")
        if not lmp_cmd:
            sys.exit()
        run_str = self.run_template.render(lmp_cmd=lmp_cmd)
        with open(self.run_input_dir / "sub_run.sh", "w") as f0:
            f0.write(run_str)

        for ii, atoms in enumerate(atoms_lists):
            for jj, params in enumerate(params_list):
                # if params.get("random", None) is None:
                params["random"] = np.random.randint(0, 10000000)
                tmp_path = self.run_input_dir / f"{ii}_{jj}"
                tmp_path.mkdir_p()
                lmp_data_fn = tmp_path / "data.lmp"
                lmp_input_fn = tmp_path / "in.lmp"
                c0 = Structure.from_atoms(atoms)
                c0.write_lmp_data(lmp_data_fn, unique_numbers=self.unique_numbers)
                lmp_input_str = self.lmp_temp.render(mass_dict=mass_dict, **params)
                with open(lmp_input_fn, "w") as f0:
                    f0.write(lmp_input_str)

                with open(tmp_path / "mlip.ini", "w") as f0:
                    f0.write(f"mtp-filename   fitted.mtp\n")
                    f0.write("select      FALSE\n")
                shutil.copy2(prev_model_fn, tmp_path / "fitted.mtp")

    def do_run(self):
        before_dir = Path(os.getcwd()).abspath()
        self.machine.submit_and_wait(self.run_input_dir, self.driver)
        self.machine.download(
            self.run_finished_dir, files=["dump.nc", "*.out"],
        )

        os.chdir(before_dir)

    def post_run(self):
        finished_dataset = AtomicDataset(self.unique_numbers)
        finished_dataset.read_from_directories(self.run_finished_dir, "dump.nc", filetype="nc")
        self.new_atoms_lists = finished_dataset.ase_atoms
