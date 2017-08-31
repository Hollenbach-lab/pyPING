import pyPING_objects as ping
import pyPING_checks as check
import pyPING_supporting as support
import os

sample_location = '/home/common_arse/INDIGO/1_INDIGO_plates1-4_KIR_runs/INDIGO_1/'
bowtie_threads = 8
sequence_group = os.path.basename(os.path.normpath(sample_location))
results_directory = ''


# Check that bowtie2 can be accessed
check.bowtie2()

results_directory = support.build_results_directory(results_directory, sequence_group)

sample_name_list, sample_path_list = support.sample_finder(sample_location)
sample_dictionary = support.match_finder(sample_name_list, sample_path_list)
sample_list = support.sample_builder(sample_dictionary, sequence_group, results_directory)


for sample in sample_list:
    sample.kir_extract(os.path.normpath('KIR_sequences/'), bowtie_threads)
    print(sample)

#for sample in sample_list:
#    sample.locus_align('2DL1', bowtie_threads)
