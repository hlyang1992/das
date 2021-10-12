import abc
import logging
import shutil

from das.machine import Machine
from das import AtomicDataset
from typing import Optional
from path import Path


class Worker(abc.ABC):
    dataset: Optional[AtomicDataset] = None
    machine: Optional[Machine] = None
    keep_old: bool

    def __init__(
        self,
        unique_elements=None,
        machine=None,
        driver=None,
        setup_dir=None,
        keep_old=False,
        name=None,
        kind="worker",
        debug_params=None,
    ):
        self.logger = logging.getLogger("das.core.worker.Worker")
        self.driver = driver
        self.keep_old = keep_old
        self.kind = kind
        self.machine = machine
        self.setup_dir = setup_dir
        self.name = self.kind if name is None else name
        self.debug_params = debug_params
        self.unique_numbers = AtomicDataset.unique_elements_to_numbers(unique_elements)

        self.dataset = None
        self.new_atoms_lists = None
        self.work_directory = None
        self.run_directory = None
        # directory for run input files
        self.run_input_dir = None
        # directory for run results
        self.run_finished_dir = None
        self.worker_index = None
        self.init_read_fn = None

    def run(self, worker_index=1):
        if self.driver is None:
            self.work_directory = self.debug_params["work_directory"]
        else:
            self.work_directory = self.driver.get_stage_work_directory()
        self.worker_index = worker_index

        # prepare_run
        self.init_directory_and_data()
        self.prepare_run()
        # run
        self.do_run()
        # post_run
        self.post_run()
        self.update_dataset()
        self.save_configurations()

    def init_directory_and_data(self, run_directory_name="run"):
        self.init_read_fn = self.work_directory / "conf.pkl"
        self.run_directory = self.work_directory / f"{run_directory_name}_{self.worker_index}"
        self.run_input_dir = self.run_directory / f"input_files"
        self.run_finished_dir = self.run_directory / f"finished_files"

        # TODO remove old files
        for dd in [self.run_directory, self.run_finished_dir, self.run_input_dir]:
            if dd.exists():
                self.logger.info(f"{dd} exists, remove.")
                dd.rmtree()
            self.logger.debug(f"mkdir {dd}")
            dd.mkdir_p()

        # load init dataset
        self.dataset = AtomicDataset.load(self.init_read_fn)
        if self.worker_index == 1:
            shutil.copy2(self.init_read_fn, self.run_directory / "conf_0.pkl")

        self.new_atoms_lists = []

    def update_dataset(self):
        n_new_configurations = 0 if self.new_atoms_lists is None else len(self.new_atoms_lists)
        self.logger.info(
            f"{self.__class__.__qualname__}: {self.name} " f"create {n_new_configurations} configurations."
        )
        if not self.keep_old:
            self.dataset.data = None
        self.logger.info(f"add atoms in new_atoms_lists to dataset, keep_old: {self.keep_old}")
        if self.new_atoms_lists is not None:
            self.dataset.read_ase_atoms(self.new_atoms_lists)

    def save_configurations(self):
        to_save_fn = self.run_directory / f"conf_{self.worker_index}.pkl"
        self.dataset.save(to_save_fn)
        shutil.copy2(to_save_fn, self.init_read_fn)

    def find_fn(self, fn, model_path):
        target_fns = [self.driver.config_dir / fn, Path(model_path).abspath().dirname() / fn]
        for ii in target_fns:
            if ii.exists():
                self.logger.info(f"find file: {ii}")
                return ii
        self.logger.error(f"can't find {fn} in {target_fns}, mode_path={model_path}")
        raise ValueError(f"can't find {fn}.")

    @property
    def atoms_lists(self):
        return self.dataset.ase_atoms

    def prepare_run(self):
        """Prepare input files for run."""
        return True

    def do_run(self):
        """
        do something expensive
        """
        return True

    def post_run(self):
        """
        generate self.new_atoms_set, a list of Atoms
        """
        return True
