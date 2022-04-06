from dataclasses import dataclass
import os
import json
import pandas as pd

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

@dataclass
class Simulator:

    """Class to wrap the use of zdplaskin executable into a python object.
    The object handles simulation run and I/O of the results
    ARGS:
        ..."""
    total_simulation_time: float=1e-7
    simulation_timestep: float=1e-9
    position: float=5
    output_directory: str="."
    initial_electron_density: float=1e9
    rf_power: float=100
    rf_impedance: float=50
    rf_frequency: float=1.365e7
    gas_pressure: float=0.2
    gap_length: float=10
    radius: float=13
    gas_temperature: float=300
    water_fraction: float=0.05

    def make_run_command(self):
        zd_out_file = os.path.join(self.output_directory, "zd_out.csv")
        zd_exe = os.path.join(parent_dir, "zdplaskin", "zdplaskin.exe")
        run_str = f"""{zd_exe} {self.total_simulation_time} {self.simulation_timestep} 
        {self.position} {self.initial_electron_density} {self.rf_power} {self.rf_impedance} {self.rf_frequency} {self.gas_pressure} 
        {self.gap_length} {self.radius} {self.gas_temperature} {self.water_fraction} {zd_out_file}""".replace("\n", "").replace("\t", "")
        return run_str

    def run(self):
        # running ZDPlasKin code
        run_str = self.make_run_command()
        os.system(run_str)
        # saving json with run parameters
        self.write_parameters_to_json(os.path.join(self.output_directory, "params.json"))
        print("COMPUTATION DONE")

    def parameters_dict(self):
        res = {}
        for arg in self.__dir__():
            attr = getattr(self, arg)
            if not callable(attr) and not arg.startswith("_"):
                res[arg] = attr
        return res
        
    def write_parameters_to_json(self, json_path):
        with open(json_path, "w") as f:
            json.dump(self.parameters_dict(), f)
        