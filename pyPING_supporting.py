import pickle
import os
import urllib2
import subprocess
import time
import pyPING_objects as ping



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

def make_directory(directory):
    """Make a directory if it does not already exist.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    return(directory)

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


def sample_finder(sample_directory, recursive=True):
    """Return a tuple of a list of file names containing '.fastq' within 'sample_directory',
    and a list of cooresponding file paths in relation to 'sample_directory'.
    """
    assert os.path.isdir(sample_directory), "%r is not a directory." % sample_directory

    directories = [x[0] for x in os.walk(sample_directory)]
    file_names = [x[2] for x in os.walk(sample_directory)]

    sample_path_list = []
    sample_name_list = []

    if recursive:
        for i in range(0, len(directories)):
            [sample_path_list.append(os.path.join(directories[i], j)) for j in file_names[i] if ".fastq" in j]
            [sample_name_list.append(j) for j in file_names[i] if ".fastq" in j]
        else:
            [sample_path_list.append(os.path.join(directories[0], j)) for j in file_names[0] if ".fastq" in j]
            [sample_name_list.append(j) for j in file_names[0] if ".fastq" in j]

            assert len(sample_path_list) is not 0, "No samples found in %r" % sample_directory

    return((sample_name_list, sample_path_list))

def match_finder(sample_name_list, sample_path_list, print_pairs=False):
    """Return a dictionary of matched sample names with their respective paths.
    This is an attempt to remove any fastq name input parameter. This does not check
    if the paired files are in the same directory.
    """
    print("Automatically finding paired sequences.\n")
    start = time.time()
    found_matches = {}

    while len(sample_name_list) > 0:
        sample = sample_name_list.pop(0) ## Remove first element of Name list to match
        path = sample_path_list.pop(0) ## Remove first element of Path list
        sub_str = ""
        last_str = ""

        for char in sample:
            last_str = sub_str ## Save the previous sub_string
            sub_str += char ## Build up a sub_string to match against
            match_num = sum([i.count(sub_str) for i in sample_name_list]) ## Count the number of sub_string matches in the Name list

            if match_num == 0: ## When the number of matches changes to 0 (remember the element being matched was removed)
                match_indices = [i for i in range(0, len(sample_name_list)) if last_str in sample_name_list[i]] ## Find the indices of the matches of the previous sub_string

                if len(match_indices) is not 1: ## If there was more than 1 match, a pair was not found
                    print("%s matches found for %s." % (len(match_indices), last_str))
                    print("Removing %s from match consideration.\n" % sample)
                else: ## Otherwise build up lists of matched sample names and paths
                    match_names = [sample_name_list[i] for i in match_indices]
                    match_names.append(sample)
                    match_paths = [sample_path_list[i] for i in match_indices]
                    match_paths.append(path)

                    if len(match_names[0]) is not len(match_names[1]): # If the length of the matched names is not the same, assume a pair was not found
                        print("Confused over matching %s with %s, ignoring these samples.\n" % (match_names[0], match_names[1]))
                    else: ## Otherwise add this pair to the found_matches dictionary
                        found_matches[last_str] = {match_names[i]: match_paths[i] for i in range(0,len(match_names))}

                    del sample_name_list[match_indices[0]] ## Remove the matched item from both lists
                    del sample_path_list[match_indices[0]]
                break ## Move to the next sample
    end = time.time()
    run_time = end - start
    print("It took %ss to match %s pairs.\n" % (run_time, len(found_matches.keys())))

    if print_pairs:
        for name, info in found_matches.items():
            print("Sample name: %s" % name)
            print("Fastq files: %s" % [key for key, value in info.items()])
            print("Fastq paths: %s" % [value for key, value in info.items()])
            print("")

    return(found_matches)

def sample_builder(sample_dictionary, sequence_group, results_directory):
    sample_list = []
    for sample_name, sample_info in sample_dictionary.items():
        file_name_list = [file_name for file_name, file_path in sample_info.items()]
        file_path_list = [file_path for file_name, file_path in sample_info.items()]
        if file_name_list[0] > file_name_list[1]:
            file_name_list = list(reversed(file_name_list))
            file_path_list = list(reversed(file_path_list))

        sample_list.append(ping.Sample(sample_name, sequence_group, file_name_list[0], file_name_list[1], file_path_list[0], file_path_list[1],
        results_directory=results_directory))

    return(sample_list)

def build_results_directory(results_directory, sequence_group):
    if results_directory == '':
        results_directory = os.path.normpath('Results/' + sequence_group + '_results/')
    else:
        results_directory = os.path.normpath(results_directory)
    return(results_directory)

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
