import pickle
import os
import subprocess
import time
import pyPING_objects as ping

def make_directory(directory):
    """Make a directory if it does not already exist.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    return(directory)

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
