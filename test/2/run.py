import os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path.append(parent_dir)
from src.simulator import Simulator

if __name__ == "__main__":
    simulator = Simulator(rf_power=100, total_simulation_time=1e-6, rf_frequency=1.1e7)
    simulator.run()