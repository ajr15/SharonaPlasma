# script to parse a simulation results to a csv
# TODO: implement the following
#       - filter results by specie
#       - filter results by location range / time range
#       - filter results by other properties (rf_power range for example) 
# possibly push all results into a large SQLite database that will be searchable for all these properties
from multiprocessing.sharedctypes import Value
import pandas as pd
import os
import json
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(parent_dir)
from src.simulator import Simulator


def res_to_df(res_dir: str) -> pd.DataFrame:
    df = pd.read_csv(os.path.join(res_dir, "zd_out.csv"))
    with open(os.path.join(res_dir, "params.json"), "r") as f:
        params = json.load(f)    
    for k, v in params.items():
        df[k] = [v for _ in range(len(df))]
    return df

def read_simulation_data(parent_res_dir: str) -> pd.DataFrame:
    all_data = pd.DataFrame()
    for res_dir in os.listdir(parent_res_dir):
        if os.path.isdir(os.path.join(parent_res_dir, res_dir)):
            df = res_to_df(os.path.join(parent_res_dir, res_dir))
            all_data = pd.concat([all_data, df])
    return all_data



if __name__ == "__main__":
    # parsing first command line argument using sys 
    # this is mandatory as one needs the output directory for the run
    import sys
    if len(sys.argv) >= 3:    
        parent_dir = sys.argv[2] # second argument is the results path
    else:
        parent_dir = "=.="
        data = pd.DataFrame()
    if os.path.isdir(parent_dir):
        data = read_simulation_data(parent_dir)
    elif not parent_dir == "=.=":
        raise ValueError("Supplied parent results directory doesn't exist !")
    # using argparse for all the other argument parsing
    import argparse
    # making command line parser
    parser = argparse.ArgumentParser(description="Command line interface to filter simulation data")
    # first reading data from a given directory
    parser.add_argument("fitered_results_path", type=str, help="path to the results csv")
    parser.add_argument("simulation_output_directory", type=str, help="path to the simulation output direrctory")
    parser.add_argument("--take_columns", type=str, default=None, nargs="+", help="(optional) specify columns to take into the results")
    # parsing filter arguments: can filter by columns and variable ranges
    for k in data.columns:
        parser.add_argument("--{}_range".format(k.strip()), type=float, default=None, nargs=2, help="minimun and maximum values for {}".format(k))
    args = parser.parse_args()
    # filtering by parameter ranges
    for col in data.columns:
        req_range = getattr(args, "{}_range".format(col.strip()).replace("-", "_"))
        if not req_range is None:
            data = data[req_range[0] <= data[col]]
            data = data[data[col] <= req_range[-1]]
    # filtering data - columns
    if not args.take_columns is None:
        data = data.loc[:, args.take_columns]
    # writing results to csv
    data.to_csv(args.fitered_results_path)
