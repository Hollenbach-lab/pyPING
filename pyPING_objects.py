import pyPING_checks as check
import pyPING_supporting as support
import pyPING_locus as locus
import subprocess
import os

class Sample:

    def __init__(self, name, group, name_fastq_1, name_fastq_2, path_fastq_1, path_fastq_2,
    results_directory, kir_extracted=False, kir_path_fastq_1='', kir_path_fastq_2=''):
        self.name = name
        self.group = group
        self.name_fastq_1 = name_fastq_1
        self.name_fastq_2 = name_fastq_2
        self.path_fastq_1 = path_fastq_1
        self.path_fastq_2 = path_fastq_2
        self.kir_extracted = kir_extracted
        self.kir_path_fastq_1 = kir_path_fastq_1
        self.kir_path_fastq_2 = kir_path_fastq_2
        self.results_directory = results_directory
        self.KIR_read_count = 0
        self.KIR_2DL1_read_count = 0
        self.temp_file_directory = os.path.join(results_directory, 'TempRunFiles/')

        support.make_directory(self.temp_file_directory)

    def __repr__(self):
        return '<Sample instance Name:%s Group:%s>' % (self.name, self.group)

    def __str__(self):
        return 'Name: %s\nGroup: %s\nResults Directory: %s\nFile Location: %s\nKIR extracted: %s\n' % (self.name, self.group, self.results_directory, os.path.dirname(self.path_fastq_1), self.kir_extracted)

    def kir_extract(self, kir_extracted_path, bowtie_threads):
        kir_extracted_path = support.make_directory(os.path.join(self.results_directory, kir_extracted_path))

        bt2_p = str(bowtie_threads)
        bt2_5 = '3'
        bt2_3 = '7'
        bt2_L = '20'
        bt2_i = 'S,1,0.5'
        bt2_min_score = 'L,0,-0.187'
        bt2_I = '75'
        bt2_X = '1000'
        bt2_x = 'Reference/KIR/output'
        bt2_1 = self.path_fastq_1
        bt2_2 = self.path_fastq_2
        bt2_stream = os.path.join(self.temp_file_directory, self.name + '.delete')
        bt2_al_conc = os.path.join(kir_extracted_path, self.name + '_KIR_%.fastq.gz')
        bt2_un = os.path.join(self.temp_file_directory, 'dump.delete')

        if not os.path.exists(os.path.join(kir_extracted_path, self.name + '_KIR_1.fastq.gz')):
            InputArgs = ['bowtie2', '-p', bt2_p, '--trim5', bt2_5, '--trim3', bt2_3, '-L', bt2_L, '-i', bt2_i, '--score-min', bt2_min_score,
            '-I', bt2_I, '-X', bt2_X, '-x', bt2_x, '-1', bt2_1, '-2', bt2_2, '-S', bt2_stream, '--al-conc-gz', bt2_al_conc, '--un-gz', bt2_un]

            print('\n\nExtracting KIR reads from: ' + self.name + '\n')
            InputString = ' '.join(InputArgs)
            print(InputString)
            subprocess.call(InputArgs)
            os.remove(bt2_stream)
            os.remove(bt2_un)
        else:
            print('\n\nKIR reads found for: ' + self.name + '\n')

        self.kir_extracted=True
        self.kir_path_fastq_1=os.path.join(kir_extracted_path, self.name + '_KIR_1.fastq.gz')
        self.kir_path_fastq_2=os.path.join(kir_extracted_path, self.name + '_KIR_2.fastq.gz')

        input_string = 'gunzip -c ' + self.kir_path_fastq_1 + ' | wc -l'
        line_count = int(subprocess.check_output(input_string, shell=True))

        self.KIR_read_count = line_count

    def locus_align(self, locus_to_align, bowtie_threads):
        if locus_to_align == "2DL1":
            locus.KIR_2DL1(self, bowtie_threads)
