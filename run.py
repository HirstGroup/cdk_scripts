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

    popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        print(popen.stderr.read())
        raise subprocess.CalledProcessError(return_code, cmd)


def run(cmd):

    lines = []

    for line in execute(cmd):
        print(line, end="")
        lines.append(line)

    output = ''.join(lines)

    return output