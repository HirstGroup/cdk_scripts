import os
import subprocess
from io import StringIO


def run(cmd):
    """
    Run a command using subprocess and make code crash if there's an error
    and print error on console

    Parameters
    ----------
    cmd: command to run

    Returns
    -------
    output: str (and runs command)
        std output from process run 
    """

    result = subprocess.run(cmd, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print(result.stderr.decode('utf-8'))
    output = result.stdout.decode('utf-8')
    print(output)

    if result.returncode != 0:
        raise Exception(f"Command {cmd} failed!")

    return output
