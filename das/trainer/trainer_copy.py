import logging
import shutil

from das import AtomicDataset
from das.core import Worker
from das.utils import AutoRepr


@AutoRepr
class MTPTrainerCopy(Worker):
    def __init__(
            self,
            kind="mtp_trainer_copy",
            keep_old=True,
            copy_conf=False,
            copy_model=True,
            **kwargs,
    ):
        super().__init__(kind=kind, keep_old=keep_old, **kwargs)
        self.copy_model = copy_model
        self.copy_conf = copy_conf
        self.logger = logging.getLogger('das.trainer.trainer_copy.MTPTrainerCopy')

    def prepare_run(self):
        super().prepare_run()

    def do_run(self):
        pass

    def post_run(self):
        model_dir = self.work_directory / "model"
        if self.driver is not None:
            # copy from previous model
            if self.copy_model:
                prev_model_dir = self.driver.get_prev_model_dir()
                self.logger.info(
                    f"converaged, copy model from {prev_model_dir} to {model_dir}"
                )
                shutil.copytree(prev_model_dir, model_dir)
            # copy from previous conf
            if self.copy_conf:
                prev_train_conf = self.driver.get_prev_train_conf()
                self.logger.info(f"converaged, load conf from {prev_train_conf}")
                prev_dataset = AtomicDataset.load(prev_train_conf)
                self.new_atoms_lists = prev_dataset.ase_atoms
