import subprocess
import urllib2
import pyPING_supporting as support

def bowtie2():
    """Checks if bowtie2 is accessible.
    """
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')
    print('bowtie2 has been found and can be accessed.\n')
