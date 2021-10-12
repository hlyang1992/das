import re

from art import text2art
from jinja2 import Template
from path import Path


def get_ascii_text(text, **kwargs):
    art_text = text2art(text, **kwargs)
    return art_text


def glob_filename(data_paths, filename):
    if isinstance(data_paths, str):
        data_paths = [data_paths]

    all_dat_files = (jj for ii in data_paths for jj in Path(ii).walk(filename))
    return all_dat_files


def natural_key(string_):
    """See https://stackoverflow.com/questions/2545532/python-analog-of-phps-natsort-function-sort-a-list-using-a-natural-order-alg"""
    return [int(s) if s.isdigit() else s for s in re.split(r"(\d+)", string_)]


def generate_mtp_untrained_model(min_dist, max_dist, species_count, index, radial_basis_size=8):
    head_template_file = Path(__file__).abspath().dirname() / "untrained_mtps/unfitted.jinja2"
    alpha_info_file = Path(__file__).abspath().dirname() / f"untrained_mtps/{str(index).zfill(2)}.mtp"
    with open(head_template_file) as f0:
        head_data = f0.read()
    temp = Template(head_data)
    out_result = temp.render(
        {
            "min_dist": min_dist,
            "max_dist": max_dist,
            "species_count": species_count,
            "index": index,
            "radial_basis_size": radial_basis_size,
        }
    )
    with open(alpha_info_file) as f1:
        out_result += "\n" + f1.read()
    return out_result


def split_jobs(n_jobs, n_tasks):
    # TODO: make this function better
    job_splits_res = []
    freq_job = int(n_jobs / n_tasks) + 1 if n_jobs > n_tasks else 1
    for i in range(0, n_jobs, freq_job):
        if i + freq_job - 1 < n_jobs:
            n1, n2 = i, i + freq_job
        elif i < n_jobs:
            n1, n2 = i, n_jobs
        else:
            continue
        job_splits_res.append([n1, n2])
    return job_splits_res


class AutoRepr:
    """Auto generate __repr__ function for class.
    https://stackoverflow.com/questions/750908/auto-repr-method
    """

    @staticmethod
    def repr(obj):
        items = []
        for prop, value in obj.__dict__.items():
            try:
                item = "%s = %r" % (prop, value)
                assert len(item) < 500
            except AssertionError:
                item = "%s: <%s>" % (prop, value.__class__.__name__)
            items.append(item)

        return "%s(%s)" % (obj.__class__.__name__, ", ".join(items))

    def __init__(self, cls):
        cls.__repr__ = AutoRepr.repr
        self.cls = cls

    def __call__(self, *args, **kwargs):
        return self.cls(*args, **kwargs)


def get_template(template_file):
    with open(template_file) as f0:
        job_template = Template(f0.read())
    return job_template


if __name__ == "__main__":
    n_jobs, n_tasks = 30, 20
    res = split_jobs(n_jobs, n_tasks)
    print(len(res), res)
