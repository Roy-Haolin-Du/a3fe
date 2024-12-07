"""Utilities for the Ensemble, Window, and Simulation Classes"""

import getpass as _getpass
import glob as _glob
import logging as _logging
import os as _os
import subprocess as _subprocess
from dataclasses import dataclass as _dataclass
from time import sleep as _sleep
from typing import List as _List
from typing import Optional as _Optional

from ._logging_formatters import _A3feFileFormatter, _A3feStreamFormatter
from ._utils import retry as _retry
from .enums import JobStatus as _JobStatus


@_dataclass
class Job:
    """Class to hold information about a job"""

    virtual_job_id: int
    command_list: _List[str]
    slurm_job_id: _Optional[int] = None
    status: _JobStatus = _JobStatus.NONE  # type: ignore
    slurm_file_base: _Optional[str] = None

    def __str__(self) -> str:
        # Avoid printing the command, which may be long
        return f"Job (virtual_job_id = {self.virtual_job_id}, slurm_job_id= {self.slurm_job_id}), status = {self.status}"

    @property
    def slurm_outfile(self) -> str:
        if self.slurm_file_base is None:
            raise AttributeError(f"{self} has no slurm_outfile")
        matching_files = _glob.glob(f"{self.slurm_file_base}*")
        if len(matching_files) == 0:
            raise FileNotFoundError(f"No files matching {self.slurm_file_base}*")
        # Take the most recent file
        newest_file = max(matching_files, key=_os.path.getctime)
        return newest_file

    def has_failed(self) -> bool:
        """Check whether the job has failed"""
        with open(self.slurm_outfile, "r") as f:
            for line in f.readlines():
                error_statements = [
                    "NaN or Inf has been generated along the simulation",
                    "Particle coordinate is NaN",
                ]
                for error in error_statements:
                    if error in line:
                        return True

        return False


class VirtualQueue:
    """A virtual slurm queue which has no limit on the number
    of queued jobs, which submits slurm jobs to the real queue
    when there are few enough jobs queued. This gets round slurm
    queue limits."""

    def __init__(
        self,
        queue_len_lim: int = 2000,
        log_dir: str = "./output",
        stream_log_level: int = 10,
    ) -> None:
        """
        Initialise the virtual queue.

        Parameters
        ----------
        queue_len_lim : int, Optional, default: 2000
            The maximum number of jobs to queue in the real queue.
        log_dir : str, Optional, default: "./output"
            The directory to write the log to.
        stream_log_level : int, Optional, default: 10
            The log level for the stream handler.

        Returns
        -------
        None
        """
        self._slurm_queue: _List[Job] = []
        self._pre_queue: _List[Job] = []
        self._available_virt_job_id = 0
        self.queue_len_lim = queue_len_lim
        self.log_dir = log_dir
        self._stream_log_level = stream_log_level

        # Set up logging
        self._set_up_logging()

        # Write out initial settings
        self._update_log()

    def _set_up_logging(self) -> None:
        """Set up logging for the virtual queue"""
        # Virtual queues didn't use to have the ._stream_log_level attribute. This
        # code ensures backwards compatibility.
        if not hasattr(self, "_stream_log_level"):
            self._stream_log_level = _logging.INFO
        # If logging has already been set up, remove it
        if hasattr(self, "_logger"):
            handlers = self._logger.handlers[:]
            for handler in handlers:
                self._logger.removeHandler(handler)
                handler.close()
            del self._logger
        self._logger = _logging.getLogger(str(self))
        self._logger.setLevel(_logging.DEBUG)
        self._logger.propagate = False
        # For the file handler, we want to log everything
        file_handler = _logging.FileHandler(f"{self.log_dir}/virtual_queue.log")
        file_handler.setFormatter(_A3feFileFormatter())
        file_handler.setLevel(_logging.DEBUG)
        # For the stream handler, we want to log at the user-specified level
        stream_handler = _logging.StreamHandler()
        stream_handler.setFormatter(_A3feStreamFormatter())
        stream_handler.setLevel(self._stream_log_level)
        # Add the handlers to the logger
        self._logger.addHandler(file_handler)
        self._logger.addHandler(stream_handler)

    @property
    def stream_log_level(self) -> int:
        """The log level for the stream handler."""
        return self._stream_log_level

    @stream_log_level.setter
    def stream_log_level(self, value: int) -> None:
        """Set the log level for the stream handler."""
        self._stream_log_level = value
        self._set_up_logging()

    @property
    def queue(self) -> _List[Job]:
        """The queue of jobs, both real and virtual."""
        return self._slurm_queue + self._pre_queue

    def __str__(self) -> str:
        return self.__class__.__name__

    def submit(self, command_list: _List[str], slurm_file_base: str) -> Job:
        """
        Submit a job to the virtual queue.

        Parameters
        ----------
        command_list : _List[str]
            The command to be run by sbatch.
        slurm_file_base : str
            The base name of the slurm file to be written. This allows
            the slurm file to be checked for errors.

        Returns
        -------
        job : Job
            The Job submitted to the virtual queue."""
        virtual_job_id = self._available_virt_job_id
        # Increment this so that it is never used again for this queue
        self._available_virt_job_id += 1
        job = Job(virtual_job_id, command_list, slurm_file_base=slurm_file_base)
        job.status = _JobStatus.QUEUED  # type: ignore
        self._pre_queue.append(job)
        self._logger.info(f"{job} submitted")
        # Now update - the job will be moved to the real queue if there is space
        self.update()
        return job

    def kill(self, job: Job) -> None:
        """Kill and remove a job from the real and virtual queues."""
        # All jobs are in the joint queue, so remove it from there
        self.queue.remove(job)
        # If the job is in the real queue, kill it
        if job in self._slurm_queue:
            self._slurm_queue.remove(job)
            _subprocess.run(["scancel", str(job.slurm_job_id)])
        else:  # Job is in the pre-queue
            self._pre_queue.remove(job)
        job.status = _JobStatus.KILLED  # type: ignore

    def _read_slurm_queue(self) -> _List[int]:
        """
        Extract all running slurm job IDs from the SLURM
        queue.

        Returns
        -------
        running_slurm_job_ids: _List[int]
            List of running SLURM job IDs for the user
        """

        # Get job ids of currently running jobs. This occasionally fails when SLURM is
        # busy (e.g. 'slurm_load_jobs error: Socket timed out on send/recv operation'),
        # so retry a few times, waiting a while in between
        @_retry(times=5, exceptions=(ValueError), wait_time=120, logger=self._logger)
        def _read_slurm_queue_inner() -> _List[int]:
            """This inner function is defined so that we can pass self._logger
            to the decorator"""
            # Get job ids of currently running jobs. This assumes no array jobs.
            # Only read running, pending, suspended, and completing jobs (R, PD, S, CG).
            commands = [
                ["squeue", "-h", "-u", _getpass.getuser(), "-t", "R,PD,S,CG"],
                ["awk", "{print $1}"],
                ["grep", "-v", "-E", "'\\[|_'"],
                ["paste", "-s", "-d,", "-"],
            ]

            # Create a pipe for the first command
            process = _subprocess.Popen(commands[0], stdout=_subprocess.PIPE)

            # Chain the commands
            for cmd in commands[1:]:
                process = _subprocess.Popen(
                    cmd, stdin=process.stdout, stdout=_subprocess.PIPE
                )

            # Get the output
            output, _ = process.communicate()

            running_slurm_job_ids = [
                int(job_id)
                for job_id in output.decode("utf-8").strip().split(",")
                if job_id != ""
            ]
            return running_slurm_job_ids

        return _read_slurm_queue_inner()

    def _submit_job(self, job_command_list: _List[str]) -> int:
        """
        Submit the supplied command to slurm using sbatch.

        Parameters
        ----------
        job_command_list : _List[str]
            The comands to be run by sbatch.

        Returns
        -------
        slurm_job_id : int
            The job id for the submitted command
        """

        # Define inner loop to allow use of retry decorator with self.logger
        @_retry(
            times=15,
            exceptions=(ValueError, RuntimeError),
            wait_time=5,
            logger=self._logger,
        )
        def _submit_job_inner(job_command_list: _List[str]) -> int:
            cmds = ["sbatch"] + job_command_list
            process = _subprocess.Popen(
                cmds,
                stdin=_subprocess.PIPE,
                stdout=_subprocess.PIPE,
                stderr=_subprocess.STDOUT,
                close_fds=True,
            )
            if process.stdout is None:
                raise ValueError("Could not get stdout from process.")
            process_output = process.stdout.read()
            process_output = process_output.decode("utf-8").strip()
            try:
                slurm_job_id = int((process_output.split()[-1]))
                return slurm_job_id
            except Exception as e:
                raise RuntimeError(f"Error submitting job: {process_output}") from e

        return _submit_job_inner(job_command_list)

    def update(self) -> None:
        """Remove jobs from the queue if they have finished, then move jobs from
        the pre-queue to the real queue if there is space."""
        # First, remove jobs from the queue if they have finished
        running_slurm_job_ids = self._read_slurm_queue()
        n_running_slurm_jobs = len(running_slurm_job_ids)
        # Remove completed jobs from the queues and update their status
        for job in self._slurm_queue:
            if job.slurm_job_id not in running_slurm_job_ids:
                # Check if it has failed
                if job.has_failed():
                    job.status = _JobStatus.FAILED  # type: ignore
                else:
                    job.status = _JobStatus.FINISHED  # type: ignore
        # Update the slurm queue
        self._slurm_queue = [
            job
            for job in self._slurm_queue
            if job.slurm_job_id in running_slurm_job_ids
        ]

        # Submit jobs if possible
        if n_running_slurm_jobs < self.queue_len_lim:
            # Move jobs from the pre-queue to the real queue
            n_jobs_to_move = self.queue_len_lim - n_running_slurm_jobs
            jobs_to_move = self._pre_queue[:n_jobs_to_move]
            self._pre_queue = self._pre_queue[n_jobs_to_move:]
            self._slurm_queue += jobs_to_move
            # Submit the jobs
            for job in jobs_to_move:
                job.slurm_job_id = self._submit_job(job.command_list)

        # self._logger.info(f"Queue updated")
        # self._logger.info(f"Slurm queue slurm job ids: {[job.slurm_job_id for job in self._slurm_queue]}")
        # self._logger.info(f"Slurm queue virtual job ids: {[job.virtual_job_id for job in self._slurm_queue]}")
        # self._logger.info(f"Pre-queue virtual job ids: {[job.virtual_job_id for job in self._pre_queue]}")

    def _update_log(self) -> None:
        """Update the log file with the current status of the queue."""
        self._logger.debug("##############################################")
        for var in vars(self):
            self._logger.debug(f"{var}: {getattr(self, var)} ")
        self._logger.debug("##############################################")

    def wait(self) -> None:
        """Wait for all jobs to finish."""
        while len(self.queue) > 0:
            self.update()
            _sleep(30)

    def _flush(self) -> None:
        """Remove all the jobs from the queu, regardless of status."""
        self._slurm_queue = []
        self._pre_queue = []
        self._available_virt_job_id = 0
        self._update_log()
        self._logger.info("Queue flushed")
