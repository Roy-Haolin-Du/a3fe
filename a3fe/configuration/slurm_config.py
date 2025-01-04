#slurm.py
from pydantic import BaseModel as _BaseModel
from .engine_config import EngineConfig
from .engine_config import GromacsConfig

class SlurmConfig(_BaseModel):
    """
    Pydantic model for holding SLURM configuration.
    """

    partition: str = "main"
    #time: str = "24:00:00" 
    nodes: int = 1
    #ntasks_per_node: int = 1
    output: str = "somd-array-gpu-%A.%a.out"
    #error: str = "slurm.err"
    gres: str = "gpu:1"

    def to_slurm_header(self) -> str:
        """
        Generates the SLURM header string based on the configuration.
        """
        return (
            f"#!/bin/bash\n"
            f"#SBATCH -o {self.output}\n"
            f"#SBATCH -p {self.partition}\n"
            f"#SBATCH -n {self.nodes}\n"
            f"#SBATCH --gres={self.gres}\n"
        )

    def generate_script(self, engine_config: EngineConfig, lambda_value: float) -> str:
        """
        Generate a SLURM script for running a simulation.

        Parameters
        ----------
        engine_config : EngineConfig
            The engine configuration specifying how to run the simulation.
        lambda_value : float
            The lambda value to use for the simulation.

        Returns
        -------
        str
            The complete SLURM script.
        """
        header = self.to_slurm_header()
        run_command = engine_config.make_run_command(lambda_value)
        return f"{header}\n\nsrun {run_command}\n"