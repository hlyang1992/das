from das import AtomicDataset
from das.core import Worker
from das.utils import AutoRepr


@AutoRepr
class MergePrevSelector(Worker):
    def __init__(
        self, kind="mrege_prev_selector", **kwargs,
    ):
        super().__init__(kind=kind, **kwargs)
        self.prev_train_dataset = None

    def prepare_run(self):
        prev_train_conf = self.driver.get_prev_train_conf()
        self.prev_train_dataset = AtomicDataset.load(prev_train_conf)

    def do_run(self):
        pass

    def post_run(self):
        new_dataset = AtomicDataset(self.unique_numbers)
        new_dataset += self.dataset
        new_dataset += self.prev_train_dataset
        self.new_atoms_lists = new_dataset.ase_atoms
