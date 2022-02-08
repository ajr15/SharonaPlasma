# script to easily run massively parallel simulations of ZDPlasKin
import argparse
import os
from itertools import product
import numpy as np
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(parent_dir)
from src.simulator import Simulator

def dict_cartesian_prod(args_dict):
    """method to take d={'a: [1, 2], 'b': [3, 4]} and return [{'a': 1, 'b': 3}, {'a': 1, 'b': 4}..."""
    res = []
    for v_tuple in product(*args_dict.values()):
        d = {k: v for k, v in zip(args_dict.keys(), v_tuple)}
        res.append(d)
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parser for command line interface to the simulator")
    parser.add_argument("parent_res_dir", type=str, help="Parent results directory, contains all raw simulation data")
    parent_res_dir = parser.parse_args("parent_res_dir").parent_res_dir
    # making results directory 
    if not os.path.isdir(parent_res_dir):
        os.makedirs(parent_res_dir)
    # initializing parameters with default values
    parameters = Simulator(1e-8, 1e-6, 5, os.path.join(parent_dir)).parameters_dict()
    # reading desired parameters from command line arguments
    for parameter, value in parameters.items():
        if not parameter == "output_directory":
            parser.add_argument("--{}_value".format(parameter), default=value, type=float, help="provide single value for {}".format(parameter))
            parser.add_argument("--{}_range".format(parameter), type=float, nargs=3, help="provide range (min, max, nsteps) for {}".format(parameter))
            args = parser.parse_args("{}_value".format(parameter), "{}_range".format(parameter))
            if not len(getattr(args, "{}_range".format(parameter))) == 0:
                parameters[parameter] = np.linspace(*getattr(args, "{}_range".format(parameter)))
            else:
                parameters[parameter] = getattr(args, "{}_value".format(parameter))
    # making arguments file for slurm array
    args_dicts = dict_cartesian_prod(parameters)
    args_str = ""
    for i, d in enumerate(args_dicts):
        # defining and creating output directory
        calc_out_dir = os.path.join(parent_res_dir, str(i))
        if not os.path.isdir(calc_out_dir):
            os.mkdir(calc_out_dir)
        # adding output_directory to parameter list
        parameters["output_directory"] = calc_out_dir
        # making Simulator object with parameters
        simulator = Simulator(**parameters)
        args_str += simulator.make_run_command() + "\n"
    # writing arguments to arguments.txt file in the slurm directory
    with open(os.path.join(parent_res_dir, "arguments.txt"), "w") as f:
        f.write(args_str)
    # making slurm submit file
    # reading slurm script tamplate
    with open(os.path.join(parent_dir, "scripts", "sbatch_template.src"), "r") as f:
        text = f.read()
    # filling template with run details
    text.replace("{njobs}", len(args_dicts))
    text.replace("{outfile}", os.path.join(parent_res_dir, "%a", "stdout.txt"))
    text.replace("{args_file}", os.path.join(parent_res_dir, "arguments.txt"))
    # writing the desired run file
    with open(os.path.join(parent_res_dir, "submit.src"), "w") as f:
        f.write(text)
    # submitting script with sbatch 
    os.system("sbatch {}".format(os.path.join(parent_res_dir, "submit.src")))