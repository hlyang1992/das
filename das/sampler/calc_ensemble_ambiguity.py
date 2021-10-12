#!/usr/bin/env python
import sys

import numpy as np
from netCDF4 import Dataset
from path import Path


def read_force_from_nc(fn):
    rootgrp = Dataset(fn)
    atomic_types = rootgrp["type"][...].data[0]
    forces = np.array(rootgrp["/forces"][...].data, dtype=np.float64)
    return forces, atomic_types


def get_force_ambiguity(work_path):
    fns = list(work_path.walk("force.*.nc"))
    atomic_types = read_force_from_nc(fns[0])[1]
    unique_types = sorted(list(set(atomic_types.tolist())))
    unique_types.append(-1)

    force_set = np.array([read_force_from_nc(fn)[0] for fn in fns], dtype=np.float64)
    force_mean = np.mean(force_set, axis=0)
    force_diff = force_set - force_mean
    ee = np.sqrt(np.mean(np.sum(force_diff ** 2, axis=3), axis=0))
    # print(ee.shape)

    for ii in unique_types:
        filter_index = np.arange(len(atomic_types)) if ii < 0 else np.where(atomic_types == ii)[0]
        # print(ii, filter_index, atomic_types)
        res_ee = np.zeros((ee.shape[0], 3))
        res_ee[:, 0] = np.max(ee[:, filter_index], axis=1)
        res_ee[:, 1] = np.min(ee[:, filter_index], axis=1)
        res_ee[:, 2] = np.mean(ee[:, filter_index], axis=1)
        fout = work_path / f"af.out" if ii < 0 else f"af_{ii}.out"
        np.savetxt(fout, res_ee, header="af: max min mean")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        p0 = Path(sys.argv[1])
    else:
        p0 = Path(".")
    get_force_ambiguity(p0)
