import logging
import os
import sys

from ase.calculators.vasp.create_input import GenerateVaspInput
from path import Path

from das import AtomicDataset
from das.core import Worker
from das.utils import AutoRepr
from das.utils import glob_filename
from .parse_vasprun import parse_incar, count_sc_steps


def check_success(fn):
    params = parse_incar(fn)
    nelm = int(params.get("NELM", 60))
    n_sc_steps = count_sc_steps(fn)
    if n_sc_steps is None:
        return False
    else:
        return max(n_sc_steps) < nelm


def delete_failed_calculations(p0):
    p0 = Path(p0)
    files = list(glob_filename(p0, "vasprun.xml"))
    failed_files = [i for i in files if not check_success(i)]
    if len(failed_files) > 0:
        logging.warning(f"remove {len(failed_files)} failed files: {failed_files}")
    for fn in failed_files:
        fn.remove()


module_logger = logging.getLogger("das.labeler.vasp_labeler")


@AutoRepr
class VaspLabeler(Worker):
    def __init__(
        self,
        kind="vasp_labeler",
        merge_prev_conf=True,
        filter_dft_database=True,
        update_dft_database=True,
        vasp_parms=None,
        **kwargs,
    ):
        super().__init__(kind=kind, **kwargs)
        self.logger = logging.getLogger("das.labeler.vasp_labeler.VaspLabeler")
        self.merge_prev_conf = merge_prev_conf
        self.dft_database_path: Path = self.driver.root_directory / "dft_labeled.pkl"
        self.dft_database: AtomicDataset = AtomicDataset(self.unique_numbers)
        self.filter_dft_database = filter_dft_database
        self.update_dft_database = update_dft_database
        self.vasp_parms = vasp_parms
        self.gen_vasp_inputs = GenerateVaspInput()
        self.setup_vasp()

    def setup_vasp(self):
        if self.vasp_parms.get("xc", None) is None:
            self.vasp_parms["xc"] = "pbe"
        self.gen_vasp_inputs.set(**self.vasp_parms)

    def prepare_run(self):
        super().prepare_run()
        # load dft labeled dataset
        if self.dft_database_path.exists():
            self.dft_database = AtomicDataset.load(self.dft_database_path)

        if self.dataset.data is None:
            self.driver.convergence = True
            self.logger.info(f"Converaged")

        if self.driver.convergence:
            return True
        # find these need dft calculation
        todo_dft = AtomicDataset(self.unique_numbers)
        if self.dft_database.data is not None and self.filter_dft_database:
            todo_dft.data = self.dataset.data.loc[~self.dataset.data.uid.isin(self.dft_database.data.uid)]
        else:
            todo_dft.data = self.dataset.data.copy()
        ase_atoms_lists = todo_dft.ase_atoms
        self.logger.info(f"{len(ase_atoms_lists)} CONFIGURATIONS NEED LABEL!")
        if ase_atoms_lists is None or len(ase_atoms_lists) == 0:
            if self.driver is not None:
                self.driver.convergence = True
                self.logger.info("no confs need dft calculation, converaged")
        vasp_cmd = self.machine.get_cmd("vasp_cmd")
        if not vasp_cmd:
            sys.exit()
        vasp_sub_run_str = f'COMMAND="{vasp_cmd}"\n$COMMAND > jobs.out 2>&1\n\n'
        if self.driver is not None:
            if not self.driver.convergence:
                with open(self.run_input_dir / "sub_run.sh", "w") as f0:
                    f0.write(vasp_sub_run_str)
                for count, atoms in enumerate(ase_atoms_lists):
                    uid = atoms.info.get("uid", "ase")
                    tmp_path = self.run_input_dir / str(count)
                    if not tmp_path.exists():
                        tmp_path.mkdir_p()
                    with open(tmp_path / "tag", "w") as f0:
                        f0.write(uid)

                    atoms.write(
                        tmp_path / "POSCAR", format="vasp", direct=False, vasp5=True, label=uid,
                    )
                    self.gen_vasp_inputs.initialize(atoms)
                    self.gen_vasp_inputs.write_incar(atoms, directory=tmp_path)
                    self.gen_vasp_inputs.write_potcar(directory=tmp_path)
                    self.gen_vasp_inputs.write_kpoints(directory=tmp_path)

    def do_run(self):
        if self.driver is not None:
            if not self.driver.convergence:
                before_dir = Path(os.getcwd()).abspath()
                # TODO: when calculation failed, try redo
                self.machine.submit_and_wait(self.run_input_dir, self.driver)
                self.machine.download(
                    self.run_finished_dir, files=["vasprun.xml", "tag", "*.out", "OUTCAR"],
                )
                os.chdir(before_dir)

    def post_run(self):
        if self.driver.convergence:
            return True

        cur_dft_dataset = AtomicDataset(self.unique_numbers)
        # delete failed calculations
        delete_failed_calculations(self.run_finished_dir)
        cur_dft_dataset.read_from_directories(self.run_finished_dir, filename="vasprun.xml", filetype="ase")

        if self.update_dft_database:
            self.dft_database += cur_dft_dataset
            self.dft_database.save(self.dft_database_path)

        # add previous train data
        prev_train_conf = self.driver.get_prev_train_conf()
        if prev_train_conf is not None and self.merge_prev_conf:
            self.logger.info(f"add prev train conf")
            cur_dft_dataset += prev_train_conf

        self.new_atoms_lists = cur_dft_dataset.ase_atoms
