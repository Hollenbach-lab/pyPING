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

def reference_update(force=False):
    """Checks the IPD-KIR database version, if a new version is found, or if
    'force' = True, the KIR sequence files are downloaded, and the local version
    updated. Returns boolean True if an update was made.
    """

    ## Check to see if the IPD-KIR url can be accessed
    url = urllib2.Request('http://www.ebi.ac.uk/ipd/kir/docs/version.html')
    if not internet_check(url):
        print('IPD-KIR database can not be reached, skipping update check.\n')
        return(False)

    external_database_version = support.ipd_kir_db_version()
    external_database_version_split = external_database_version.split('.')

    internal_database_version = support.kir_db_version_read()
    internal_database_version_split = internal_database_version.split('.')

    flag = False
    updated = False
    for i in range(0,len(external_database_version_split)):
        if int(external_database_version_split[i]) > int(internal_database_version_split[i]):
            flag = True
            break

    if(force):
        print("Forcing update of KIR alignments...\n")
        support.reference_download("Reference/")
        support.kir_db_version_write(external_database_version)
        updated=True
    elif(not flag):
        print('KIR alignments are up-to-date (v.' + internal_database_version + ')\n')
    else:
        print("Updating KIR alignments...\n")
        support.reference_download("Reference/")
        support.kir_db_version_write(external_database_version)
        updated=True

    return(updated)

def internet_check(url_to_check):
    try:
        urllib2.urlopen(url_to_check, timeout=5)
        return True
    except urllib2.URLError as err:
        return False
