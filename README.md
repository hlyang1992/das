# DAS

Dual Adaptive Sampling for Machine Learning Interatomic potential.

## How to cite 

If you use this code in your research, please cite this using: `Hongliang Yang, Yifan Zhu, Erting Dong, Yabei Wu, 
Jiong Yang, and Wenqing Zhang. Dual adaptive sampling and machine learning interatomic potentials for modeling materials with chemical bond hierarchy. Phys. Rev. B 104, 094310 (2021). `

## Install

### Install `pymtp`

You should first install the python interface for mtp: https://github.com/hlyang1992/pymtp

###  Install `das`

You can download the code by
 
```bash
git clone https://github.com/hlyang1992/das
cd das
cp -r <path-to-mlip-2>/untrained_mtps/*.mtp das/utils/untrained_mtps
```

Then remove the redundant settings from each mtp file. Only the following settings can be retained for each mtp file:

```yaml
radial_funcs_count = 
alpha_moments_count = 
alpha_index_basic_count = 
alpha_index_basic = 
alpha_index_times_count = 
alpha_index_times = 
alpha_scalar_moments = 
alpha_moment_mapping =
```

Install das by

```bash
cd <path-to-das>
pip install -r requirements.txt
pip install .
```

## Usage

```bash
das  config_dir  job_name
``` 

## Configuration

The configuration directory `config_dir` must contain the configuration file `conf.yaml`, which controls all sampling processes. The `conf.yaml` file should look like the following:

```yaml
"global_settings":

"machine_settings":

"selector_settings": {} 

"labeler_settings":

"trainer_settings":

"sampler_settings":

"init_conf_setting":

"iter_params_template":

"iter_params":
```

* `global_settings`: 

```yaml
"global_settings":
  # The elements in the system, the order of the elements does not matter, the program automatically numbers the 
  # atomic types according to their atomic number from smallest to largest.
  "unique_elements": [ "Co", "Sb" ]
  # path to VASP Pseudopotential Database, see detail at https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#vasp
  "vasp_pp_path": "path_to_directory" 
```

* `machine_settings`: 

All time-consuming computational tasks such as sampling, labeling, and training can be dispatched to designated machines via ssh. Currently only LSF is supported and migration to other job management systems is very easy.

```yaml
"machine_settings":
  "machine_1":
    # The supported machine types are now: `machine_lsf`, `machine_shell`
    "machine_type": "machine_lsf"
    "host": "ip address"
    "user": "username"
    "password": "password"
    # Exclude these nodes when submitting tasks.
    "bad_nodes": [ ] # #BSUB -R "hname!={{node}}"
    "port": 22
    # number of cores for each task
    "n_cores": 40 # #BSUB -n {{ncores}}
    "n_tasks": 40 # The maximum number of tasks to run simultaneously.
    "q_name": "short" # #BSUB -q {{q_name}}
    "env_source_file": "env.sh" # env.sh is in the config_dir
    "run_dir": "path-to-run-directory-in-target"
    "extra_params":
      "vasp_cmd": "mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp"
      "lmp_cmd": "mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP lmp_mlp"
      "mlip_cmd": "mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP mlp train"
      "python_cmd": "absolute path to python path"
  "machine_2":
    # setting for machchine_2
    "machine_type": "machine_lsf"
    # ...
```

You should prepare a file to set the environment variables. The program will source this file to set the environment 
variables after connecting to the machine via ssh. For technical reasons please see:
[The remote shell environment doesn’t match interactive shells](http://www.fabfile.org/faq.html#the-remote-shell-environment-doesn-t-match-interactive-shells)

* `sampler_settings`： 

```yaml
"scale_1":
  "kind": "scale_box"
  "scale_factors": [0.998, 0.9985, 0.999]
"scale_2":
  "kind": "scale_box"
  "scale_factors": [[0.998, 0.9985, 0.999, 0.997], # a
                    [1.002, 1.003, 1.004, 1.005],  # b
                    [0.997, 0.995, 0.999, 0.996]] # c
"nvt_0": 
  "kind": "lmp_model_sampler"
  "max_number_confs": 5
  "min_number_confs": 0
  "machine": "machine_1"
  "lmp_vars":
    "temp": [ 100, 150 ]
    "steps": [ 10000 ]
    "nevery": [ 20 ]
    "prev_steps": [ 0 ]
 
"npt_0": 
  "kind": "lmp_model_sampler"
  "max_number_confs": 5
  "min_number_confs": 0
  "machine": "machine_2"
  "lmp_vars":
    "temp": [ 100, 150 ]
    "steps": [ 10000 ]
    "nevery": [ 20 ]
    "press": [100, 200] # bar
    "prev_steps": [ 0 ]
```


* "labeler_settings"

We use ase to generate input files (INCAR, POTCAR, KPOINTS) for VASP calculation. Please see detail at
[Ase vasp calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#vasp)
  
```yaml
"labeler_settings":
  "vasp":
    "kind": "vasp"
    "machine": "ty_label"
    "vasp_parms":
      "xc": "pbe"
      "prec": "A"
      # other setting for vasp calculations
```

* "trainer_settings"

```yaml
"trainer_settings":
  "train_5_model":
    "kind": "mtp_trainer"
    "machine": "ty_train" 
    "model_index": 18 
    "min_dist": 1.39 
    "max_dist": 5.0
    "n_models": 5 
    "train_from_prev_model": true 
```

* `init_conf_setting`: 

```yaml
"init_conf_setting":
  "-1": [ "init_MD.cfg" ]
  "-2": [ "init_1.vasp" ]
  "-3": [ "init_2.vasp" ]
```
  
* `iter_params_template`: 

```yaml
"iter_params_template":
  "0":
    "init_conf": [ -1 ]
    "sampler": [ ]
    "selector": [ ]
    "labeler": [ ]
    "trainer": [ "train_5_model" ]
  "10":
    "init_conf": [ -2 ]
    "sampler": [ "scale_0", "nvt_0" ]
    "selector": [ ]
    "labeler": [ "vasp" ]
    "trainer": [ "train_5_model" ]
  "20":
    "init_conf": [ -3 ]
    "sampler": [ "npt_0"]
    "selector": [ ]
    "labeler": [ "vasp" ]
    "trainer": [ "train_5_model" ]
  "30":
    "init_conf": [ -2,-3 ]
    "sampler": [ "npt_0"]
    "selector": [ ]
    "labeler": [ "vasp" ]
    "trainer": [ "train_5_model" ]
```


* `iter_params`: 

```yaml
"iter_params":
  [
    [ "0" ],
    # If the last one is LOOP, repeat all the previous ones until convergence.
    ["10", "LOOP"], 
    ["30", "LOOP"],
    ["10", "10"]  
    ["20"],
  ]
```



