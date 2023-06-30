import os
import subprocess
from io import StringIO


def run1(cmd):
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


def execute(cmd):
    """
    Execute a command and output lines one by one

    Parameters
    ----------
    cmd: str
        command to run

    Yields
    ------
    stdout_line: str
        line of output command
    """

    popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, bufsize=1)

    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 

    popen.stdout.close()

    return_code = popen.wait()

    if return_code:
        print(popen.stderr.read()) #.decode("utf-8"))
        raise subprocess.CalledProcessError(return_code, cmd)


def run(cmd):
    """
    Run a command, print output lines one by one and return all output

    Parameters
    ----------
    cmd: str
        command to run

    Returns
    -------
    output:
        standard output of command
    """

    lines = []

    for line in execute(cmd):
        print(line, end="")
        lines.append(line)
    output = ''.join(lines)

    return output