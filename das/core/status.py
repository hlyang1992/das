class CycleNumber:
    def __init__(self, N):
        self._val = 0
        self.N = N

    def inc(self):
        self._val += 1
        if self._val >= self.N:
            self._val = 0
            return 1
        return 0

    @property
    def val(self):
        return self._val

    @val.setter
    def val(self, val):
        if val < self.N:
            self._val = val
        else:
            self._val = val - self.N


class StatusRecord:
    def __init__(self, n_stage, n_interior, status_flag=None):
        self._status_flag = None
        self.n_stage = n_stage
        self.n_interior = n_interior
        self.status_flag = status_flag

    @property
    def status_flag(self):
        return [self._status_flag[0], self._status_flag[1].val, self._status_flag[2].val]

    @status_flag.setter
    def status_flag(self, status_flag):
        if status_flag is None:
            iter_i, stage_i, interior_i = 0, 0, 0
        else:
            iter_i, stage_i, interior_i = status_flag
        self._status_flag = [iter_i, CycleNumber(self.n_stage), CycleNumber(self.n_interior)]
        self._status_flag[1].val = stage_i
        self._status_flag[2].val = interior_i

    @property
    def iter_i(self):
        return self.status_flag[0]

    @iter_i.setter
    def iter_i(self, value):
        self._status_flag[0] = value

    @property
    def stage_i(self):
        return self.status_flag[1]

    @stage_i.setter
    def stage_i(self, value):
        self._status_flag[1].val = value

    @property
    def interior_i(self):
        return self.status_flag[2]

    @interior_i.setter
    def interior_i(self, value):
        self._status_flag[2].val = value

    def inc_iter(self):
        self._status_flag[0] += 1

    def inc_stage(self):
        flag = self._status_flag[1].inc()
        if flag:
            self.inc_iter()

    def inc_interior(self):
        flag = self._status_flag[2].inc()
        if flag:
            self.inc_stage()

    def __str__(self):
        return f"Iteration: {self.iter_i} Stage: {self.stage_i} Interior: {self.interior_i}"


if __name__ == "__main__":
    # c0 = CycleNumber(3)
    # for i in range(10):
    #     print(i, c0.inc(), c0.val)

    s0 = StatusRecord(n_stage=5, n_interior=4)
    for i in range(20):
        s0.inc_interior()
        print(s0)
    print("=" * 10)
    for i in range(10):
        s0.inc_stage()
        print(s0, s0.status_flag)

    s0 = StatusRecord(n_stage=5, n_interior=4)
    s0.inc_interior()
    print(s0)
    s0.interior_i = 0
    print(s0)
