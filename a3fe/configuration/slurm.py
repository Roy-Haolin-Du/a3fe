from pydantic import BaseModel as _BaseModel


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

    def generate_somd_script(self, lambda_value: float) -> str:
        """
        Generate the run_somd.sh script for running SOMD simulations.

        Parameters
        ----------
        lambda_value : float
            The lambda value to be used in the simulation.

        Returns
        -------
        str
            The content of the run_somd.sh script.
        """
        script_content = self.to_slurm_header()
        script_content += "\necho \"CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES\"\n\n"
        script_content += f"echo \"lambda is: {lambda_value}\"\n\n"
        script_content += f"srun somd-freenrg -C somd.cfg -l {lambda_value} -p CUDA\n"
        return script_content

    def generate_gromacs_submission_script(self) -> str:
        """
        Generates a script for submitting GROMACS simulations.
        """
        script_content = self.to_slurm_header()
        script_content += "\nmodule load gromacs\n"
        script_content += f"cd {self.output}\n"
        script_content += "echo 'Starting GROMACS simulation...'\n"
        script_content += "srun gmx mdrun -v -deffnm somd -ntomp 1\n"
        script_content += "echo 'GROMACS simulation completed.'\n"
        return script_content