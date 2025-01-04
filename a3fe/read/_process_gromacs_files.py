"""Functionality to manipulate GROMACS files."""

from logging import Logger as _Logger
from typing import Optional as _Optional
from typing import Tuple as _Tuple

import numpy as _np

from .exceptions import ReadError as _ReadError


def read_mdp_option(mdp_file: str, option: str) -> str:
    """read mdp option from file

    Parameters
    ----------
    mdp_file : str
        mdp file path
    option : str
        option name

    Returns
    -------
    value : str
        option value
    """
    with open(mdp_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        # skip comment and empty line
        if line.strip().startswith(';') or not line.strip():
            continue
        if line.split('=')[0].strip() == option:
            return line.split('=')[1].strip()
    raise ValueError(f"Option {option} not found in MDP file {mdp_file}")


def write_mdp_option(
    mdp_file: str, 
    option: str, 
    value: str,
    logger: _Optional[_Logger] = None
) -> None:
    """write option to gromacs mdp file

    Parameters
    ----------
    mdp_file : str
        mdp file path
    option : str
        option name
    value : str
        option value
    logger : Optional[Logger]
        logger object for logging
    """
    with open(mdp_file, "r") as f:
        lines = f.readlines()
    
    option_found = False
    for i, line in enumerate(lines):
        # skip comment and empty line
        if line.strip().startswith(';') or not line.strip():
            continue
        if line.split('=')[0].strip() == option:
            lines[i] = f"{option} = {value}\n"
            option_found = True
            break
    
    if not option_found:
        if logger:
            logger.warning(
                f"Option {option} not found in MDP file {mdp_file}. Appending to end."
            )
        # ensure last line is newline
        if lines[-1][-1] != "\n":
            lines[-1] += "\n"
        lines.append(f"{option} = {value}\n")
    
    with open(mdp_file, "w") as f:
        f.writelines(lines)


def read_xvg_data(xvg_file: str) -> _Tuple[_np.ndarray, _np.ndarray]:
    """read data from gromacs xvg file

    Parameters
    ----------
    xvg_file : str
        xvg file path

    Returns
    -------
    time : np.ndarray
        time data
    data : np.ndarray
        corresponding data value
    """
    time = []
    data = []
    
    with open(xvg_file, "r") as f:
        for line in f:
            if line.startswith(("#", "@")):  # skip comment and label line
                continue
            values = line.strip().split()
            if len(values) >= 2:
                time.append(float(values[0]))
                data.append(float(values[1]))
    
    return _np.array(time), _np.array(data)


def read_dhdl_file(dhdl_file: str) -> _Tuple[_np.ndarray, _np.ndarray, _np.ndarray]:
    """read dhdl file from gromacs

    Parameters
    ----------
    dhdl_file : str
        dhdl.xvg file path

    Returns
    -------
    time : np.ndarray
        time data
    dh : np.ndarray
        dH/dλ data
    dh_err : np.ndarray
        dH/dλ error
    """
    time, dh = read_xvg_data(dhdl_file)
    dh_err = _np.zeros_like(dh)  # GROMACS不直接提供误差估计
    
    return time, dh, dh_err
