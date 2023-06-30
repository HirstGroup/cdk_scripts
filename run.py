import subprocess

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

    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True)
        output = result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as err:
        print(err.output.decode("utf-8"))
        raise Exception(f'Command {cmd} failed')

    print(output)

    return output