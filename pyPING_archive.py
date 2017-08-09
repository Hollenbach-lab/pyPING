# ---------------------------------------
# pyPING_run.py archive
# ---------------------------------------

# 8/8/17 ################################
force_reference_update = True

# Check that the reference files are up to date
reference_updated = check.reference_update(force_reference_update)
if(reference_updated):
    support.build_locus_combinations('2DL1', '2DS1')
    support.build_locus_combinations('2DL2', '2DL3')
    support.build_locus_combinations('3DL1', '3DS1')
    support.build_kir_indices("Reference/")


# ---------------------------------------
# pyPING_checks.py archive
# ---------------------------------------

# 8/8/17 ################################
import urllib

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


# ---------------------------------------
# pyPING_supporting.py archive
# ---------------------------------------
import urllib

# 8/8/17 ################################

def kir_db_version_read():
    """Return the version of IPD-KIR database for the local alignment files.
    This value is stored locally in ping_info.pkl.
    """
    pkl_file = open('ping_info.pkl', 'rb') ## DEVEL: want to turn ping_info.pkl into a variable
    data = pickle.load(pkl_file)
    kir_version = data['local_kir_db_version']
    pkl_file.close()

    return(kir_version)

def kir_db_version_write(version):
    """Write the version of IPD-KIR database to a local ping_info.pkl file.
    This happens whenever the local alignments are updated from the database.
    """
    pkl_file = open('ping_info.pkl', 'rb') ## DEVEL: want to turn ping_info.pkl into a variable
    data = pickle.load(pkl_file)

    data['local_kir_db_version'] = version
    pkl_file.close()

    pkl_file = open('ping_info.pkl', 'wb') ## DEVEL: want to turn ping_info.pkl into a variable
    pickle.dump(data, pkl_file, 2)
    pkl_file.close()

    print("Local KIR alignments updated to IPD-KIR database version " + version + "\n")

def ipd_kir_db_version():
    """Return the latest IPD-KIR database version found online.
    """
    url = urllib2.Request('http://www.ebi.ac.uk/ipd/kir/docs/version.html')
    page_info = urllib2.urlopen(url)
    data = page_info.read()

    """Version finding is a hard-coded table search. There is only one table on the page,
    and the first element of that table is the latest database version.
    """
    version_find_beg = data.find('<tr align="left" valign="middle">')
    data = data[version_find_beg:]
    version_find_beg = data.find('>')
    data = data[version_find_beg+1:]
    version_find_beg = data.find('>')
    data = data[version_find_beg+1:]
    version_find_end = data.find('<')
    kir_db_version = data[:version_find_end]
    return(kir_db_version)

def reference_download(reference_directory):
    """Download the latest KIR sequences from IPD-KIR, place them into folders
    by locus name in 'reference_directory'.
    """
    database_url = 'ftp://ftp.ebi.ac.uk/pub/databases/ipd/kir/fasta/' ## This link is always to the latest version.
    url = urllib2.Request(database_url)
    page_info = urllib2.urlopen(url)
    data = page_info.read()
    data_lines = data.split("\r\n") ## Splitting the data up into line arrays.

    data_strings = [i.split(" ") for i in data_lines] ## Splitting the lines up into string arrays.
    data_strings = sum(data_strings, []) ## Removing the top level list

    data_kir = [i for i in data_strings if '_gen.fasta' in i] ## Finding all of the *_nuc.fasta files
    url_kir = [database_url + i for i in data_kir] ## Building urls for the *_nuc.fasta files

    base_reference_directory = reference_directory
    make_directory(base_reference_directory)

    """Download each KIR[LOCUS]_nuc.fasta file individually to
    'reference_directory/KIR[LOCUS]/'
    """
    for i in url_kir:
        kir_reference_directory = base_reference_directory + os.path.basename(i).split("_")[0] + '/'
        kir_referencePath = os.path.join(kir_reference_directory, os.path.basename(i))
        make_directory(kir_reference_directory)

        print("Downloading " + i + "\n")
        url = urllib2.Request(i)
        page_info = urllib2.urlopen(url)

        with open(kir_referencePath, 'wb') as local_file:
            local_file.write((page_info.read()))

def build_kir_indices(reference_directory):
    """Builds bowtie2 index files based on fasta files found in
    'reference_directory'.
    """
    kir_reference_directories = [i[0]+'/' for i in os.walk(reference_directory)]
    assert kir_reference_directories, reference_directory + ' directory not found.'

    kir_reference_directories = kir_reference_directories[1:]
    assert kir_reference_directories, 'No locus directories found in ' + reference_directory + "."

    kir_reference_files = sum(sum([[[(x+j, x+j.split("_gen.fasta")[0]) for j in i[2] if 'gen.fasta' in j] for i in os.walk(x)] for x in kir_reference_directories], []), [])

    for x in kir_reference_files:
        print('\nBuilding index for ' + x[0] + ':')
        print('\tbowtie2-build ' + x[0] + ' ' + x[1])
        subprocess.check_output(['bowtie2-build', x[0], x[1]])
        print('\n')
    print('Finished building indices.\n')

def build_locus_combinations(locus_1, locus_2):
    combination_directory = make_directory(os.path.join('Reference/', 'KIR' + locus_1 + '_' + locus_2))
    combination_file = 'KIR' + locus_1 + '_' + locus_2 + '_gen.fasta'

    ## Writing the combination file
    all_kir_file = os.path.join('Reference/', 'KIR/', 'KIR_gen.fasta')
    assert os.path.exists(all_kir_file), all_kir_file + ' not found.'

    ## Writing the combination file with both loci
    all_locus_path = os.path.join(combination_directory, 'all' + combination_file)
    with open(all_locus_path, 'w') as outfile:
        record = False
        with open(all_kir_file) as in_file:
            for line in in_file:
                if locus_1 in line or locus_2 in line:
                    record = True
                elif '>KIR:' in line:
                    record = False

                if record:
                    outfile.write(line)


    ## Writing the combination file with neither loci
    neither_locus_path = os.path.join(combination_directory, 'not' + combination_file)
    with open(neither_locus_path, 'w') as outfile:
        record = False
        with open(all_kir_file) as in_file:
            for line in in_file:
                if locus_1 in line or locus_2 in line:
                    record = False
                elif '>KIR:' in line:
                    record = True

                if record:
                    outfile.write(line)

    return(all_locus_path, neither_locus_path)
    
