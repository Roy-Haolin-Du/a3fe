"""Simulation engine backends."""

from ..configuration import EngineType
from ._engine import EngineBackend
from .gromacs import GromacsBackend
from .somd import SomdBackend

engine_backend_registry = {
    EngineType.SOMD: SomdBackend(),
    EngineType.GROMACS: GromacsBackend(),
}

CONFIG_FILE_SUFFIXES = tuple(
    suffix
    for backend in engine_backend_registry.values()
    for suffix in backend.config_file_suffixes
)
