
import os
import numpy as np
import pandas as pd
import json
import notion_df as ndf
from typing import Iterable
from tqdm import tqdm
from itertools import product
from datetime import datetime
from shutil import copy, rmtree
from subprocess import run, PIPE


ndf.pandas()


def to_iter(x):
    return x if isinstance(x, Iterable) and not isinstance(x, str) else [x]


def to_val(value: str):
    if value.lower() == 'false':
        return False
    if value.lower() == 'true':
        return True
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def list_sum(*arglist):
    start = []
    for l in arglist:
        start += l
    return start


"""
Slurm Job Array creation script (Please pardon our dust, work in progress)
"""


# Parameters for your sweep, you only need to enter values that
#   differ from the base
config_list = [
#        
    {
    "sweeps": "20000",
    "active": "0",
    "active_epochs": "1",
    "seed": "0",
    "fixed": "False",
    "block": [False]
#         'cycle': [True]

}
]


# Specify the graph class and graph names that you want to run,
#   and how many iterations of each
_iter = np.arange(20)
nvals = np.arange(260, 501, 10)
sweep_jobs = {
    'graph_class': ['erdos_renyi'],
#     'graph': list_sum(
#         [f"erdos_renyi_bimodal_n_{n}_a_4_{i}.gset" for n, i in product(nvals, _iter)],
#         [f"barabasi_albert_bimodal_n_{n}_m_3_{i}.gset" for n, i in product(nvals, _iter)]
#     ),
#     'graph': [f"erdos_renyi_bimodal_n_{n}_a_4_{i}.gset" for n, i in product(nvals, _iter)],
    'graph': ['erdos_renyi_bimodal_n_420_a_4_15.gset', 'erdos_renyi_bimodal_n_440_a_4_5.gset'],
    'iter': range(20)
}

# Define the simulation name and short justification
# (folder will be sim{sim_no:03d}_{name})
project = 'BlockSA'
sim_no = 0
name = "er_scaling_dropoff"
simname = f'sim{sim_no:03d}_{name}'
justif = """
    Pls be more efficient
"""
timeval = datetime.now().isoformat(timespec='seconds')



# Simulation system parameters
time = '5-00:00:00'  # runtime in days-hours:minutes:seconds
partition = 'ising'
memory = 8  # req memory in GB
concurrent = 200  # max number of jobs running concurrently
sparse = True  # eigen linalg type
makejobs = 8

# Get path information from system environment variables
#  (should only need to configure once)
# SET THIS TO YOUR DATA PATH
queue_path = './queue_SA.py'
directory_path = '/scratch/mburns13/data/ISING_RESULTS'
# automatically get SA version
solver_name = os.path.abspath(__file__).split(os.sep)[-2]
bluehive = os.path.exists('/scratch/mhuang_lab')
thispath = os.path.dirname(os.path.abspath(__file__))
basepath = './base_config.json'  # Change as desired
base_config = json.loads(open(basepath).read())
base_config = dict([(u, to_val(v))
                   for u, v in base_config.items()])
cmdlineparams = set(list(base_config.keys()))
# Replace with your path to GSET folder (or keep GSET)
graph_dir = os.environ['GSET']
# Replace with your path to random inits (or don't)
random_dir = os.environ['RND']

config_list = [dict([(u, to_iter(v))
                     for u, v in i.items()]) for i in config_list]
sweep_jobs = dict([(u, to_iter(v))
                   for u, v in sweep_jobs.items()])
# __________SHOULD NOT NEED TO RECONFIGURE ANYTHING BELOW THIS LINE__________
config_list = [{**i, **sweep_jobs} for i in config_list]
njobs = 0
for i in config_list:
    njobs += np.prod([len(j) for j in i.values()])
print(njobs)
if njobs > 20000:
    print(f'Estimated {njobs} jobs. Continue? (y/n)')
    token = input()
    while (token not in {"y", "n"}):
            print('Please enter either y or n')
            token = input().lower()
    if token == 'n':
        quit()

#  Very ugly way of getting the best values in, TODO make more pythonic
cutdata = []
for gclass in sweep_jobs['graph_class']:
    bksfile = f'{graph_dir}/{gclass}/summary.txt'
    if os.path.exists(bksfile):
        delim = ' ' if gclass in {'set', 'tiny', 'small', 'K_graphs'} else ','
        df = pd.read_csv(bksfile,
                         delimiter=delim,
                         on_bad_lines='skip')
        if gclass in {'set', 'tiny', 'small', 'K_graphs'}:
            df = df.rename(
            columns={'Gset': 'graph', 'Cuts': "BKS", 'GSET': 'graph', 'CUT': 'BKS'})
        cutdata.append((gclass, df))
bksdict = {}
if len(cutdata) > 0:
    cutvals = dict(cutdata)
    bksdict = dict(
        [(k, dict([(g, c)
                   for g, c in df[['graph', 'BKS']].to_records(index=False)]))
         for k, df in cutvals.items()])


colnames = list_sum(*[list(i.keys()) for i in config_list])
paramset = set(colnames)


def main():
    print(random_dir)
    result_dir = f'{directory_path}/{solver_name}/{simname}'
    if os.path.exists(result_dir):
        print('Directory already exists, ' +
              'do you want to replace (r) or exit (q)')
        token = input().lower()
        while (token not in {"r", "q"}):
            print('Please enter either r or q')
            token = input().lower()
        if token == 'r':
            print('Removing previous directory...')
            rmtree(result_dir)
        else:
            print('Exiting...')
            exit(0)
    print('Making/building directory...')
    os.makedirs(result_dir)
    os.mkdir(result_dir+'/build')
    copy_and_compile(result_dir+'/build')
    parampath = result_dir+'/parameters'
    os.mkdir(parampath)
    # send a copy of this script and the base config to the parameter subdir
    copy(__file__, parampath)
    copy(basepath, parampath)
    jobpath = result_dir+'/jobs'

    os.mkdir(jobpath)
    jobrecords = []
    execpath = result_dir+'/build/block_sa'
    
    print('Enumerating jobs and writing launch files...')
    for ind, c in tqdm(enumerate(configurations(config_list))):
        jobno = ind + 1
        jobrecords.append(tuple([jobno] + [c[i] for i in colnames]))
        make_job_file(execpath=execpath,
                      job_no=jobno,
                      config=c,
                      target_directory=result_dir)
    numjobs = len(jobrecords)
    record_df = pd.DataFrame.from_records(jobrecords, 
                                          columns=['job']+colnames)
    record_df.to_csv(parampath+'/records.csv', index=False)
    print(f'Wrote record of {numjobs} jobs to {parampath+"/records.csv"}')
    print('Writing the sbatch launch script...')
    os.mkdir(result_dir+'/data')
    os.mkdir(result_dir+'/err')
    make_launch_file(result_dir)
    print(f'To run, enter: sbatch {result_dir}/parameters/launch_queue.sh')

    
    newrow = pd.DataFrame(data=[('In Progress',
                                 name,
                                 timeval,
                                 solver_name,
                                 project,
                                 justif,
                                 result_dir,
                                 parampath)],
                      columns=[
                          'Status',
                          'Name',
                          'Date',
                          'Solver',
                          'Project',
                          'Justification',
                          'Result Path',
                          'File Path'
                      ])
    newrow.to_notion('https://www.notion.so/bf10637ad2e74eac8d20c279fd8d6ae2?v=830e524c2432416eae1420298de0bc66',
           title="Experiment Log")
    
    
# define configurations as values within an enumerated variable space
def configurations(config_list: list):
    for parameters in config_list:
        param_names = list(parameters.keys())
        param_space = product(*list(parameters.values()))
        for config in param_space:
            newdict = dict([(name, val)
                            for name, val in
                            zip(param_names, config)])
            if not os.path.exists(
                f'{graph_dir}/{newdict["graph_class"]}/{newdict["graph"]}'
            ):
                print(os.path.join(
                    graph_dir,
                    newdict["graph_class"],
                    newdict["graph"]
                ))
                continue
            yield {**base_config, **newdict}


def copy_and_compile(target_directory: str):
    """Copy the source files to `target_directory` and compile the binary

    Args:
        target_directory (str): Directory containing the compiled/source 
            code for the experiment
    """
    if bluehive:
        os.system('module unload gcc')
        os.system('module load gcc/11.2.0')
        os.system('module load cmake')
    output = run(['cmake', f'-B{target_directory}', f'-S{thispath}'], stderr=PIPE)
    try:
        output.check_returncode()
    except:
        print("CMake error with code:\n{}".format(output.stderr.decode()))
        exit(0)
    output = run(['make', f'--directory={target_directory}', f'-j{makejobs}'], stderr=PIPE)
    try:
        output.check_returncode()
    except:
        print("Make error with code:\n{}".format(output.stderr.decode()))
        exit(0)

def make_launch_file(target_directory: str):
    script_contents = f"""#!/bin/bash
#SBATCH --job-name QUEUE_{simname}
#SBATCH --output={target_directory}/parameters/queue.run
#SBATCH --error={target_directory}/err/queue.err
#SBATCH -p {partition}
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --time={time}

python3 {queue_path} {simname}
"""
    with open(target_directory+'/parameters/launch_queue.sh', mode='w') as launch:
        launch.write(script_contents)
        
def make_job_file(execpath: str, 
                  target_directory: str, 
                  config: dict, 
                  job_no: int):
    """Generate launch files for each SA experiment

    Args:
        target_directory (str): Path to the directory containing the job    
            scripts
        config (dict)
    """
    header = f"""#!/bin/bash
#SBATCH --job-name {job_no}_BSA
#SBATCH --output={target_directory}/data/bsa_output_{job_no}.run
#SBATCH --error={target_directory}/err/bsa_output_{job_no}.err
#SBATCH -p {partition}
#SBATCH --mem={memory}G
#SBATCH --ntasks=1
#SBATCH --time={time}
    """
    filepath = f'{target_directory}/jobs/bsa_job_{job_no}.sh'
    with open(filepath, mode='w') as jobfile:
        
        jobfile.write(f'{header}\n')
        jobfile.write(f'{execpath}\\\n')
        for parameter, value in config.items():
            if parameter not in cmdlineparams:
                continue
            if isinstance(value, bool):
                if value:
                    jobfile.write(f'\t--{parameter}\\\n')
            else:
                jobfile.write(f'\t--{parameter} {value}\\\n')
        jobfile.write(
            f'\t--graph {graph_dir}/' +
            f'{config["graph_class"]}/{config["graph"]}\\\n')
        jobfile.write(f'\t--seed {config["iter"]}\\\n')
        graph_class = config["graph_class"]
        graph = config["graph"]

        if graph_class in bksdict and graph in bksdict[graph_class]:
            jobfile.write(f'\t--best {bksdict[graph_class][graph]}\\\n')



if __name__ == '__main__':
    main()