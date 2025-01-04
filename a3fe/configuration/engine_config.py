from pydantic import BaseModel, Field
from typing import Dict

class EngineConfig(BaseModel):
    """
    Base class for engine configuration.

    This abstract class defines the common interface for molecular dynamics (MD) engines.
    Subclasses must implement the make_run_command method to define how simulations
    are executed for specific engines like SOMD or GROMACS.
    """

    engine_name: str = Field(..., description="Name of the MD engine (e.g., 'somd', 'gromacs').")

    def make_run_command(self, lambda_value: float) -> str:
        """
        Generate the command to run a simulation for the given lambda value.

        Parameters
        ----------
        lambda_value : float
            The lambda value to use in the simulation.

        Returns
        -------
        str
            The shell command to run the simulation.

        Notes
        -----
        This method must be implemented in subclasses for specific engines.
        """
        raise NotImplementedError("Please implement this method in a subclass for a specific engine.")


class SomdConfig(EngineConfig):
    """
    Configuration for the SOMD engine.

    This class extends the EngineConfig to include fields specific to the SOMD engine,
    such as the path to the SOMD configuration file. It also implements the make_run_command
    method to define how SOMD simulations are executed.
    """
    somd_config_file: str = Field("somd.cfg", description="Path to the SOMD configuration file.")

    def make_run_command(self, lambda_value: float) -> str:
        """
        Generate the command to run a SOMD simulation.

        Parameters
        ----------
        lambda_value : float
            The lambda value to use in the simulation.

        Returns
        -------
        str
            The shell command to run the SOMD simulation.
        """
        return f"somd-freenrg -C {self.somd_config_file} -l {lambda_value} -p CUDA"


class GromacsConfig(EngineConfig):
    """Configuration for GROMACS engine"""
    engine_name: str = "gromacs"
    mdp_file: str = Field("gromacs.mdp", description="Path to GROMACS MDP file")
    
    def generate_mdp(self, lambda_values: Dict[str, list], output_dir: str) -> None:
        """
        Generate GROMACS MDP file
        
        Parameters
        ----------
        lambda_values : Dict[str, list]
            Dictionary containing lambda values for bonded, coulomb and vdw
        output_dir : str
            Directory to write the MDP file
        """
        #temperially here, later will be put into yaml/json file
        mdp_content = [
            # Output Control
            "nstlog = 200",
            "nstenergy = 200",
            "nstxout-compressed = 1000",
            
            # Run Parameters
            "dt = 0.0020",
            "nsteps = 3000000",
            "constraints = h-bonds",
            "constraint-algorithm = LINCS",
            
            # Periodic Boundary Conditions and Cutoffs
            "pbc = xyz",
            "cutoff-scheme = Verlet",
            "ns-type = grid",
            "nstlist = 20",
            "rlist = 0.8",
            "rvdw = 0.8",
            "rcoulomb = 0.8",
            
            # Long-range Interactions
            "coulombtype = PME",
            "DispCorr = EnerPres",
            "vdwtype = Cut-off",
            
            # Pressure Coupling
            "pcoupl = C-rescale",
            "nstpcouple = 100",
            "tau-p = 4",
            "ref-p = 1.01325",
            "compressibility = 4.5e-5",
            
            # Temperature Coupling
            "integrator = sd",
            "tc-grps = system",
            "tau-t = 1.00000",
            "ref-t = 300.00",
            
            # Initial Velocities
            "gen-vel = yes",
            "gen-temp = 300.00",
            
            # Free Energy Calculation
            "free-energy = yes",
            "couple-moltype = LIG",
            "couple-lambda0 = vdw-q",
            "couple-lambda1 = none",
            "couple-intramol = yes",
            "sc-alpha = 0.5",
            "sc-power = 1",
            "sc-sigma = 0.3",
            "calc-lambda-neighbors = -1",
            "init-lambda-state = 1",
            
            # Energy Calculations
            "nstcalcenergy = 200",
            "nstdhdl = 200",
            "ld-seed = -1",
            
            # Lambda Values
            f"bonded-lambdas = {' '.join(map(str, lambda_values['bonded']))}",
            f"coul-lambdas = {' '.join(map(str, lambda_values['coul']))}",
            f"vdw-lambdas = {' '.join(map(str, lambda_values['vdw']))}"
        ]
        
        # Write MDP file
        with open(f"{output_dir}/{self.mdp_file}", "w") as f:
            f.write("\n".join(mdp_content))

    def make_run_command(self, lambda_value: float) -> str:
        """
        Generate GROMACS run command
        
        Parameters
        ----------
        lambda_value : float
            Lambda value for the simulation
            
        Returns
        -------
        str
            GROMACS run command
        """
        return (
            f"gmx grompp -f {self.mdp_file} -c gromacs.gro -p gromacs.top -o gromacs.tpr && "
            f"gmx mdrun -deffnm gromacs -dhdl dhdl.xvg"
        )


