import numpy as np

from das import AtomicDataset
from das.core import Worker
from das.utils import AutoRepr


def scale_atoms(atoms, factors):
    new_atoms = []
    if factors.ndim == 1:
        for scale in factors:
            tmp_atoms = AtomicDataset.copy_ase_atoms(atoms)
            tmp_atoms.set_cell(np.array(atoms.get_cell()) * scale,
                               scale_atoms=True)
            new_atoms.append(tmp_atoms)
    else:
        for i in range(factors.shape[1]):
            factor_i = np.reshape(factors[:, i], (-1, 1))
            tmp_atoms = AtomicDataset.copy_ase_atoms(atoms)
            new_cell = np.array(atoms.get_cell()) * factor_i
            tmp_atoms.set_cell(new_cell, scale_atoms=True)
            new_atoms.append(tmp_atoms)

    return new_atoms


@AutoRepr
class ScaleBoxSampler(Worker):
    def __init__(self, scale_factors=None, kind="scale_box", **kwargs):
        super().__init__(kind=kind, **kwargs)
        self.factors = np.array(scale_factors)

    def post_run(self):
        self.new_atoms_lists = []
        for atoms in self.atoms_lists:
            self.new_atoms_lists += scale_atoms(atoms, self.factors)
