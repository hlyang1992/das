import logging
import os
import re
import shutil
import sys

from path import Path

from das import AtomicDataset
from das.core import Worker
from das.trainer.mtp_trainer import write_unfitted_model
from das.utils import AutoRepr, get_template


@AutoRepr
class MTPTrainerRedo(Worker):
    def __init__(
        self,
        min_dist=1.0,
        max_dist=6.0,
        model_index=18,
        mlip_params=None,
        kind="mtp_trainer",
        keep_old=True,
        train_from_prev_model=False,
        n_models=1,
        run_mtp_file="run_mtp.jinja2",
        **kwargs,
    ):
        super().__init__(kind=kind, keep_old=keep_old, **kwargs)
        self.n_models = n_models
        if mlip_params is None:
            mlip_params = {}
        self.model_index = model_index

        run_mtp_file = self.find_fn(run_mtp_file, __file__)
        self.run_template = get_template(run_mtp_file)

        self.mlip_params = mlip_params
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.train_from_prev_model = train_from_prev_model
        self.logger = logging.getLogger("das.trainer.trainer_redo.MTPTrainerRedo")

    def gen_unfitted_model(self, fn):
        species_count = len(self.driver.global_settings["unique_elements"])
        self.logger.debug(
            f"generate unfitted model min_dist: {self.min_dist}, max_dist:{self.max_dist}, model_index:{self.model_index}"
        )
        write_unfitted_model(self.min_dist, self.max_dist, species_count, self.model_index, fn)

    def prepare_run(self):
        super().prepare_run()

        mlip_cmd = self.machine.get_cmd("mlip_cmd")
        if not mlip_cmd:
            sys.exit()
        run_str = self.run_template.render(mlip_cmd=mlip_cmd, **self.mlip_params)
        with open(self.run_input_dir / "sub_run.sh", "w") as f0:
            f0.write(run_str)

        prev_train_conf = self.driver.get_prev_train_conf()
        self.logger.info(f"converged, load conf from {prev_train_conf}")
        prev_dataset = AtomicDataset.load(prev_train_conf)

        for ii in range(self.n_models):
            tmp_path = self.run_input_dir / f"train_{ii}"
            tmp_path.mkdir_p()
            self.logger.info(f"train_from_prev_model:{self.train_from_prev_model}")

            if self.train_from_prev_model:
                try:
                    prev_model_i = self.driver.get_prev_model_fn(ii)
                    self.logger.info(f"copy unfitted.mtp from {prev_model_i}")
                    shutil.copy2(prev_model_i, tmp_path / "unfitted.mtp")
                except ValueError:
                    self.logger.info(f"can't find models, use default model")
                    self.gen_unfitted_model(tmp_path / "unfitted.mtp")
            else:
                self.gen_unfitted_model(tmp_path / "unfitted.mtp")
            prev_dataset.write_mtp_data(tmp_path / "trainset.cfg")
        self.new_atoms_lists = prev_dataset.ase_atoms

    def do_run(self):
        before_dir = Path(os.getcwd()).abspath()
        self.machine.submit_and_wait(self.run_input_dir, self.driver)
        self.machine.download(self.run_finished_dir, files=["fitted.mtp", "*.log"])
        os.chdir(before_dir)

    def post_run(self):
        model_dir = self.work_directory / "model"
        model_dir.mkdir_p()
        # sort by final loss value
        all_model_fns = sorted(list(self.run_finished_dir.walk("fitted.mtp")), key=find_final_loss)
        self.logger.info(f"found {all_model_fns} in {self.run_finished_dir}")
        for i, fn in enumerate(all_model_fns):
            dst_fn = model_dir / f"fitted_{i}.mtp"
            self.logger.info(f"copy from {fn} to {dst_fn}")
            shutil.copy2(fn, dst_fn)


def find_final_loss(fn):
    pattern = re.compile(r"BFGS iter (\d+): f=(\d+\.\d+$)")
    log_fn = fn.parent / "fitting.log"
    loss_final = None
    with open(log_fn) as f0:
        for line in f0:
            out = pattern.match(line)
            if out:
                loss_final = float(out.groups()[1])
    if loss_final is None:
        raise ValueError(f"Could not find final loss, {fn}, {log_fn}")
    return loss_final
