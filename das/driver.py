import logging
import os

from path import Path

from das import AtomicDataset, machine, sampler, labeler, trainer, selector
from das.core import StatusRecord
from das.utils import get_ascii_text

SAMPLER_MAP = {
    "scale_box": sampler.ScaleBoxSampler,
    "lmp_sampler": sampler.LmpSampler,
    "lmp_model_sampler": sampler.LmpAmbiguitySampler,
}
SELECTOR_MAP = {"merge_prev_selector": selector.MergePrevSelector}

MACHINE_MAP = {
    "machine_lsf": machine.MachineLSF,
    "machine_shell": machine.MachineShell,
}

LABELER_MAP = {"vasp": labeler.VaspLabeler}

TRAINER_MAP = {
    "mtp_trainer": trainer.MTPTrainer,
    "mtp_trainer_copy": trainer.MTPTrainerCopy,
    "mtp_trainer_redo": trainer.MTPTrainerRedo,
}

module_logger = logging.getLogger("das.driver")


class Driver:
    def __init__(
        self, params, config_dir, job_run_dir, label=None, status_flag=None, max_loops=200,
    ):
        self.logger = logging.getLogger("das.driver.Driver")
        self.max_loops = max_loops
        # self.root_directory = self.init_root_directory(job_run_dir, create=True)
        self.root_directory = self.init_root_directory(job_run_dir, create=False)
        # self.set_logger(self.root_directory/log_fn)
        self.status_file = self.root_directory / "status"
        self.status_record = StatusRecord(n_stage=4, n_interior=2, status_flag=status_flag)
        self.block_i = 0
        if self.status_file.exists():
            self.logger.info("Recover from status file")
            self.update_status_flag(False)
        self.label = label
        self.config_dir = Path(config_dir).abspath()
        self.params = params

        self.worker_template = {}
        self.machine_template = {}
        self.sampler = None
        self.selector = None
        self.labeler = None
        self.trainer = None
        self.params_iter_i = None
        self.global_settings = None
        self.init_conf_dict = {}
        self.convergence = False
        self.setup()

    def set_stage(self, stage_settings, maps):
        stage_dict = {}
        unique_elements = self.params["global_settings"]["unique_elements"]
        for key, item in stage_settings.items():
            tmp_stage_params = item.copy()
            stage_type = tmp_stage_params["kind"]
            stage_cls = maps.get(stage_type, None)
            if stage_cls is None:
                raise Exception(f"Can't find stage type: {stage_type}")
            for ii in ["setup_dir"]:
                if tmp_stage_params.get(ii):
                    tmp_stage_params[ii] = self.config_dir / item[ii]
            # setattr(current_stage, "unique_elements", unique_elements)
            # setattr(current_stage, "driver", self)
            tmp_stage_params["unique_elements"] = unique_elements
            tmp_stage_params["driver"] = self
            current_stage = stage_cls(**tmp_stage_params)
            machine_key = item.get("machine", None)
            if machine_key:
                setattr(current_stage, "machine", self.machine_template[str(machine_key)])
            stage_dict[str(key)] = current_stage
        return stage_dict

    def setup(self):
        # set machine
        for key, item in self.params["machine_settings"].items():
            machine_type = item["machine_type"]
            machine_cls = MACHINE_MAP.get(machine_type, None)
            if machine_cls is None:
                raise ValueError(f"can't find machine_type: {machine_type}")
            for ii in ["run_dir", "env_source_file", "template_file"]:
                if item.get(ii, None) is not None:
                    item[ii] = self.config_dir / item[ii]
            self.machine_template[str(key)] = machine_cls(**item)
            self.logger.debug(f"{item}, {machine}")

        self.worker_template["sampler"] = self.set_stage(self.params["sampler_settings"], SAMPLER_MAP)
        self.worker_template["selector"] = self.set_stage(self.params["selector_settings"], SELECTOR_MAP)
        self.worker_template["labeler"] = self.set_stage(self.params["labeler_settings"], LABELER_MAP)
        self.worker_template["trainer"] = self.set_stage(self.params["trainer_settings"], TRAINER_MAP)
        self.global_settings = self.params["global_settings"]
        vasp_pp_path = self.global_settings.get("vasp_pp_path", None)
        if vasp_pp_path is not None:
            os.environ["VASP_PP_PATH"] = vasp_pp_path

        for key, val in self.params["init_conf_setting"].items():
            self.init_conf_dict[key] = [self.config_dir / Path(i) for i in val]

        self.logger.debug(f"job template: {self.worker_template}")
        self.logger.info(f"init_conf_dict: {self.init_conf_dict}")

    def update_status_flag(self, write=True):
        self.logger.info(
            f"status_record, Block: {self.block_i} {self.status_record} " f"update status flag file, write: {write}"
        )
        if not write:
            with open(self.status_file) as f0:
                line = f0.readline().strip().split()
            self.block_i = int(line[0])
            self.status_record.status_flag = [int(i) for i in line[1:]]
            self.logger.info(f"read status from file. status_record: {self.block_i} {self.status_record}")
        else:
            with open(self.status_file, "w") as f0:
                new_status = str(self.block_i) + " "
                new_status += " ".join([str(i) for i in self.status_record.status_flag])
                f0.write(new_status)

    def run(self):
        n_blocks = len(self.params["iter_params"])
        self.logger.info(f"Total {n_blocks} blocks")
        for b_i in range(n_blocks):
            if self.block_i > b_i:
                self.logger.info(f"skip block: {b_i}")
                continue
            self.logger.info(get_ascii_text(f"\nBLOCK   START"))
            self.logger.info(f"Starting block {self.block_i}")
            start_iter = self.status_record.iter_i
            self.logger.info(f"Starting iteration {start_iter}")
            self.iter_params_block_i = self.params["iter_params"][b_i]
            if str(self.iter_params_block_i[-1]).lower() == "loop":
                self.iter_params_block_i.pop()
                self.iter_params_block_i *= self.max_loops
            n_iter = len(self.iter_params_block_i)
            for iter_i in range(start_iter, n_iter):
                self.status_record.iter_i = iter_i
                self.on_iter_start()
                self.run_iter()
                self.on_iter_end()
                if self.convergence:
                    self.logger.info(f"converaged at iteration: {iter_i}, move to next block.")
                    self.logger.info(get_ascii_text(f"\nConverged"))
                    break
            self.block_i += 1
            self.status_record.status_flag = [0, 0, 0]
            self.update_status_flag()

    def on_iter_start(self):
        iter_i = self.status_record.iter_i
        iter_params_i_index = self.iter_params_block_i[iter_i]
        self.params_iter_i = self.params["iter_params_template"][str(iter_params_i_index)]

        self.logger.info(f"read iter params: {self.params_iter_i}")
        # set worker
        for worker_name in ["sampler", "selector", "labeler", "trainer"]:
            worker_index = self.params_iter_i.get(worker_name, None)
            if worker_index is None:
                worker = []
            else:
                worker_index = [str(i) for i in worker_index]
                worker = [self.worker_template[worker_name][i] for i in worker_index]
            setattr(self, worker_name, worker)

        self.logger.info(f"{'*' * 10} ITER START {'*' * 10}")
        self.log()
        self.set_iteration_init_configuration()
        self.convergence = False

    def on_iter_end(self):
        self.status_record.stage_i = 0
        self.status_record.interior_i = 0
        self.update_status_flag()
        self.logger.info(f"{'*' * 10} ITER END {'*' * 10}\n\n")

    def on_stage_start(self):
        self.logger.info(f"{'=' * 10} STAGE START {'=' * 10}")
        self.status_record.interior_i = 0
        self.log()
        self.set_stage_init_configuration()

    def on_stage_end(self):
        self.status_record.interior_i = 1
        self.log()
        self.logger.info(f"{'=' * 10} STAGE END {'=' * 10}")
        self.status_record.inc_stage()
        self.update_status_flag()

    def set_iteration_init_configuration(self):
        iter_i = self.status_record.iter_i
        stage_i = self.status_record.stage_i
        if stage_i != 0:
            return 0
        init_conf_iter_i = self.params_iter_i.get("init_conf", None)
        iter_init_database = AtomicDataset.from_unique_elements(self.global_settings["unique_elements"])
        # if not init_conf given, use configuration of previous iteration
        for ii in init_conf_iter_i:
            # TODO: if conf is directory, add support
            for conf in self.init_conf_dict[str(ii)]:
                self.logger.info(f"add {conf} to init database")
                ext = conf.ext.lower()
                if ext == ".pkl":
                    iter_init_database += conf
                elif ext == ".lmp":
                    iter_init_database.read_from_file(conf, filetype="dump")
                elif ext == ".cfg":
                    iter_init_database.read_from_file(conf, filetype="mtp")
                else:
                    iter_init_database.read_from_file(conf, filetype="ase")
        iter_init_configuration_fn = self.get_stage_configuration(iter_i, stage_i=0)
        iter_init_database.save(iter_init_configuration_fn)
        self.logger.info(
            f"create iteration init configuration: {iter_init_configuration_fn} "
            f"with {len(iter_init_database)} configurations."
        )

    def set_stage_init_configuration(self):
        iter_i = self.status_record.iter_i
        stage_i = self.status_record.stage_i

        # configurations has been initialized by iteration
        if stage_i == 0:
            return 0
        stage_init_database = AtomicDataset.from_unique_elements(self.global_settings["unique_elements"])

        if stage_i > 0:
            prev_conf = self.get_stage_configuration(iter_i, stage_i - 1)
            stage_init_database += prev_conf

        stage_init_configuration_fn = self.get_stage_configuration(iter_i, stage_i)
        stage_init_database.save(stage_init_configuration_fn)
        self.logger.info(
            f"create stage init configuration: {stage_init_configuration_fn} "
            f"with {len(stage_init_database)} configurations."
        )

    def run_iter(self):
        stage_i = self.status_record.stage_i

        if stage_i < 1:
            self.on_stage_start()
            for i, h in enumerate(self.sampler, 1):
                h.run(i)
            self.on_stage_end()

        if stage_i < 2:
            self.on_stage_start()
            for i, h in enumerate(self.selector, 1):
                h.run(i)
            self.on_stage_end()

        if stage_i < 3:
            self.on_stage_start()
            for i, h in enumerate(self.labeler, 1):
                h.run(i)
            self.on_stage_end()

        if stage_i < 4:
            self.on_stage_start()
            for i, h in enumerate(self.trainer, 1):
                h.run(i)
            self.on_stage_end()

    def log(self):
        self.logger.info(
            f"Block: {self.block_i} status_record: {self.status_record} "
            f"work_directory:{self.get_stage_work_directory()}"
        )

    def get_stage_configuration(self, iter_i, stage_i, conf_name="conf.pkl"):
        return self.get_stage_work_directory(iter_i, stage_i) / conf_name

    def get_stage_work_directory(self, iter_i=None, stage_i=None, block_i=None):
        if block_i is None:
            block_i = self.block_i
        if iter_i is None:
            iter_i = self.status_record.iter_i
        if stage_i is None:
            stage_i = self.status_record.stage_i
        work_directory = Path(f"{self.root_directory}/" + f"{block_i}_{iter_i}_{stage_i}")
        if not work_directory.exists():
            self.logger.info(f"create work directory: {work_directory}")
            work_directory.mkdir_p()
        return work_directory

    def get_prev_stage_dir(self, stage):
        block_i = self.block_i
        iter_i = self.status_record.iter_i
        if block_i == iter_i == 0:
            return None
            # raise ValueError(f"can't find previous directory when block_i==iter_i==0")
        if iter_i == 0:
            prev_n_iter = max(
                [
                    int(i.name.split("_")[1])
                    for i in list(self.root_directory.dirs())
                    if len(i.name.split("_")) == 3 and int(i.name.split("_")[0]) == self.block_i - 1
                ]
            )

            return self.get_stage_work_directory(block_i=block_i - 1, iter_i=prev_n_iter, stage_i=stage)
        else:
            return self.get_stage_work_directory(iter_i=iter_i - 1, stage_i=stage)

    def get_prev_sample_dir(self):
        return self.get_prev_stage_dir(stage=0)

    def get_prev_select_dir(self):
        return self.get_prev_stage_dir(stage=1)

    def get_prev_label_dir(self):
        return self.get_prev_stage_dir(stage=2)

    def get_prev_train_dir(self):
        return self.get_prev_stage_dir(stage=3)

    def get_prev_train_conf(self, conf_name="conf.pkl"):
        prev_train_dir = self.get_prev_train_dir()
        if prev_train_dir is None:
            return None
        return prev_train_dir / conf_name

    def get_prev_select_conf(self, conf_name="conf.pkl"):
        prev_select_dir = self.get_prev_select_dir()
        if prev_select_dir is None:
            return None
        return prev_select_dir / conf_name

    def get_prev_model_dir(self):
        prev_train_dir = self.get_prev_train_dir()
        if prev_train_dir is None:
            return None
        prev_model_dir = prev_train_dir / "model"
        return prev_model_dir

    def get_prev_model_fn(self, ind=0):
        prev_model_fns = self.get_prev_model_fns()
        return prev_model_fns[ind]

    def get_prev_model_fns(self):
        prev_model_dir = self.get_prev_model_dir()
        if not prev_model_dir or not prev_model_dir.exists():
            raise ValueError(f"can't find model directory {prev_model_dir}")
        prev_model_fns = sorted(list(prev_model_dir.walk("fitted*.mtp")))
        if len(prev_model_fns) < 1:
            raise ValueError(f"can't find model: {prev_model_fns}")
        return prev_model_fns

    @staticmethod
    def init_root_directory(job_run_dir, create=False):
        directory = Path(f"{job_run_dir}")

        if directory.exists() and create:
            parent = directory.parent
            name = directory.name
            count = len(parent.glob(f"{name}_*"))
            directory_old = parent / f"{name}_{count}"
            module_logger.info(f"copy from {directory} to {directory_old}")
            directory.copytree(directory_old)
        else:
            directory.mkdir_p()
            module_logger.info(f"create {directory}")

        directory = directory.abspath()
        return directory


def read_and_check_config(config_name):
    from yaml import load

    try:
        from yaml import CLoader as Loader, CDumper as Dumper
    except ImportError:
        from yaml import Loader, Dumper

    with open(config_name, "r") as f:
        params = load(f, Loader=Loader)

    for key, value in params.items():
        module_logger.debug(f"{key}:{value}")

    # TODO: add check for params
    return params
