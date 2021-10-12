import logging
import os
import time
import uuid
from typing import List

from fabric import Connection
from jinja2 import Template
from path import Path

from das.utils import get_ascii_text


def update_bad_nodes(bad_nodes: List, root_directory: Path):
    fn: Path = root_directory / "bad_nodes"
    if fn.exists():
        with open(fn) as f0:
            for line in f0:
                node_name = line.strip()
                if node_name not in bad_nodes:
                    bad_nodes.append(node_name)
    # save to file
    with open(fn, "w") as f0:
        for item in bad_nodes:
            f0.write(f"{item}\n")
    return bad_nodes


class Machine:
    def __init__(self, run_dir=None, template_file=None, host=None, port=22, user=None, password=None,
                 machine_type="machine", n_tasks=20, extra_params=None, bad_nodes=None, attempt_times=10):
        if bad_nodes is None:
            bad_nodes = []
        if extra_params is None:
            extra_params = {}

        self.bad_nodes = bad_nodes
        self.host = host
        self.port = port
        self.user = user
        self.machine_type = machine_type
        self.password = password
        self.template_file = template_file
        self.n_tasks = n_tasks
        self.extra_params = extra_params
        self.logger = logging.getLogger('das.machine.machine.Machine')
        self.logger.debug(f"host: {self.host}, user:{self.user}, port:{self.port}, run_dir:{run_dir}")

        self._run_dir = run_dir
        self.conn = None
        self.machine_job_run_dir = None
        self.attempt_times = attempt_times

        self.setup()

    def setup(self):
        self.conn = Connection(self.host,
                               user=self.user,
                               port=self.port,
                               connect_kwargs={"password": self.password}
                               )

    @property
    def run_dir(self):
        return self._run_dir

    @run_dir.setter
    def run_dir(self, _run_dir):
        self._run_dir = _run_dir

    def test(self):
        try:
            res = self.conn.run(f"ls -l {self.run_dir}", hide=True)
            if res.stderr:
                self.logger.error(res.stderr)
                return False
            return True
        except Exception as e:
            print(e)
            return False

    def mkdir(self, directory):
        """
        create remote_dir. if remote_dir exists, remove it.
        """
        self.logger.info(f"remove and mkdir {directory}")
        self.conn.run(f"rm -rf {directory}", hide=True)
        self.conn.run(f"mkdir -p {directory}", hide=True)
        return directory

    def upload(self, local_job_dir, machine_target_dir):
        current_dir = Path(os.getcwd()).abspath()
        parent_local_job_dir = local_job_dir.parent
        os.chdir(parent_local_job_dir)

        job_dir_name = local_job_dir.name
        p_tar = job_dir_name + '.tgz'
        if p_tar.exists():
            p_tar.remove()
        self.logger.info(f"create {p_tar} at {os.getcwd()}")
        os.system(f"tar -cf {p_tar} {job_dir_name}")

        # upload
        self.logger.info(f"upload {p_tar} to {machine_target_dir}")
        res = self.conn.put(local=p_tar, remote=machine_target_dir)
        with self.conn.cd(machine_target_dir):
            # remove remote job_dir
            self.conn.run(f"tar xf {res.remote}", hide=True)
        remote_job_work_dir = machine_target_dir / job_dir_name

        os.chdir(current_dir)
        return remote_job_work_dir

    def submit(self):
        self.logger.debug(f"submit job at machine dir: {self.machine_job_run_dir}")
        with self.conn.cd(self.machine_job_run_dir):
            self.conn.run("bash run.sh")

    def submit_and_wait(self, local_job_dir, driver=None):
        self.logger.debug(f"submit job at local dir: {local_job_dir}")

        if driver:
            self.bad_nodes = update_bad_nodes(self.bad_nodes, driver.root_directory)
        self.generate_task_run_file(local_job_dir)
        tmp_path = str(uuid.uuid4())
        machine_target_dir = self.run_dir / tmp_path

        # ConnectionResetError: [Errno 104] Connection reset by peer
        for attempt in range(self.attempt_times):
            try:
                self.mkdir(machine_target_dir)
                # job_name = local_job_dir.name
                # machine_job_run_dir = machine_parent_dir / job_name
                self.machine_job_run_dir = self.upload(local_job_dir, machine_target_dir)
                self.submit()
                self.wait_job_finished()
            except ConnectionResetError as e:
                self.logger.warning(f"{e}")
                self.logger.info(f"{attempt}th attempt failed.")
                time.sleep(60)
                self.setup()
            else:
                break
        else:
            # we failed all the attempts - deal with the consequences.
            self.logger.error(f"we failed all {self.attempt_times} attempts, existing!")
            self.logger.error(get_ascii_text(f"\nERROR"))

    def wait_job_finished(self, wait_seconds=30):
        job_dir = self.machine_job_run_dir
        n_remote_jobs = len(self.get_job_sub_dirs(job_dir))

        while True:
            try:
                res_ok = self.check_file_exists_in_sub_dir(job_dir, filename="__ok__")
                res_start = self.check_file_exists_in_sub_dir(job_dir, filename="__start__")
                n_ok_jobs = 0 if res_ok is None else len(res_ok)
                n_start_jobs = 0 if res_start is None else len(res_start)
                self.logger.info(f"start/finished/total: {n_start_jobs}/{n_ok_jobs}/{n_remote_jobs} in {job_dir}")
                if n_ok_jobs == n_remote_jobs:
                    break
            except Exception as err:
                print(err)
            time.sleep(wait_seconds)

    def get_job_sub_dirs(self, job_dir):
        res = self.conn.run(f"ls -d {job_dir}/*/", hide=True)
        job_sub_dirs = [Path(i) for i in res.stdout.strip().split("\n")]
        return job_sub_dirs

    def check_file_exists_in_sub_dir(self, job_dir, filename="__ok__"):
        sub_dirs = []
        # print(f"find {remote_job_dir} -name {filename} -type f")
        res = self.conn.run(f"find {job_dir} -name {filename} -type f", hide=True)
        # print(res.stdout, res.stderr)
        if res.stdout:
            files = res.stdout.strip().split("\n")
            for fn in files:
                fn = Path(fn)
                sub_dirs.append(fn.parent)
            return sub_dirs
        return None

    def download(self, local_finished_job_dir, files, tar_fn="results.tgz"):
        work_dir = os.getcwd()
        with self.conn.cd(self.machine_job_run_dir):
            for fn in files:
                self.conn.run(f'find . -name "{fn}" -exec tar -uvf {tar_fn} {{}} +')
            remote_path = self.machine_job_run_dir / tar_fn
            local_path = local_finished_job_dir / tar_fn
            self.logger.debug(f"download {remote_path} to {local_path}")
            self.conn.get(remote=remote_path, local=local_path)
        os.chdir(local_finished_job_dir)
        os.system(f"tar xf {tar_fn}")
        os.chdir(work_dir)

    def get_cmd(self, name):
        cmd = self.extra_params.get(name, None)
        if not name:
            self.logger.error(f"{name} not found in machine setting, please check!")
            return None
        return cmd

    def generate_task_run_file(self, job_dir):
        pass

    @staticmethod
    def get_template(template_file=None, default_template=None):
        """Return default path to datafile."""
        if template_file is None:
            template_file = os.path.join(os.path.dirname(__file__), default_template)
        with open(template_file) as f0:
            job_template = Template(f0.read())
        return job_template
