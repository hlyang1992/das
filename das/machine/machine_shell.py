from path import Path

from das.machine import Machine
from das.utils import AutoRepr, split_jobs


@AutoRepr
class MachineShell(Machine):
    def __init__(self, template_file=None, env_source_file=None,
                 machine_type='machine_shell', **kwargs):
        super().__init__(machine_type=machine_type, **kwargs)
        self.env_source_file = env_source_file
        default_template = "shell_template.jinja2"
        self.job_template = self.get_template(template_file, default_template)

    def generate_task_run_file(self, job_dir):
        job_dir = Path(job_dir)
        all_sub_job_dir = sorted([p.abspath().name for p in job_dir.glob("*/")])
        n_jobs = len(all_sub_job_dir)

        job_splits_res = split_jobs(n_jobs, self.n_tasks)

        with open(job_dir / "sub_run.sh") as f0:
            sub_run_str = f0.read()
        for i, ii in enumerate(job_splits_res):
            job_paths = all_sub_job_dir[ii[0]:ii[1]]
            job_str = self.job_template.render({
                "job_paths": job_paths,
                "sub_run_str": sub_run_str,
            })
            with open(job_dir / f"jobs_{i}.job", "w") as writer:
                writer.write(job_str)
        # generate run.sh
        with open(job_dir / "run.sh", "w") as f0:
            if self.env_source_file is not None:
                with open(self.env_source_file) as f1:
                    f0.write(f1.read())
            f0.write("\nfor i in `ls *.job`;do\n bash $i &  \ndone\n")