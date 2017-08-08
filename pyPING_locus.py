import os
import subprocess
import pyPING_supporting as support

def KIR_2DL1(sample, bowtie_threads):
    results_directory = support.make_directory(os.path.join(sample.results_directory, '2DL1/'))

    # 2DL1S1 positive filter
    bt2_p = str(bowtie_threads)
    bt2_5 = '3'
    bt2_3 = '7'
    bt2_L = '20'
    bt2_i = 'S,1,0.5'
    bt2_min_score = 'L,0,-0.187'
    bt2_I = '75'
    bt2_X = '1000'
    bt2_x = 'Reference/KIR2DL1_2DS1/allKIR2DL1_2DS1'
    bt2_1 = sample.kir_path_fastq_1
    bt2_2 = sample.kir_path_fastq_2
    bt2_stream = os.path.join(sample.temp_file_directory, sample.name + '.delete')
    bt2_al_conc = os.path.join(sample.temp_file_directory, 'pos_2DL1S1_' + sample.name + '_%.fastq.gz')

    InputArgs = ['bowtie2', '-p', bt2_p, '--trim5', bt2_5, '--trim3', bt2_3, '-L', bt2_L, '-i', bt2_i, '--score-min', bt2_min_score,
    '-I', bt2_I, '-X', bt2_X, '-x', bt2_x, '-1', bt2_1, '-2', bt2_2, '-S', bt2_stream, '--al-conc-gz', bt2_al_conc]

    print('\n\nExtracting 2DL1S1 positive reads from: ' + sample.name + '\n')
    InputString = ' '.join(InputArgs)
    print(InputString)
    subprocess.call(InputArgs)
    os.remove(bt2_stream)

    # 2DL1S1 negative filter
    bt2_min_score = 'L,0,-0.155'
    bt2_x = 'Reference/KIR2DL1_2DS1/notKIR2DL1_2DS1'
    bt2_1 = os.path.join(sample.temp_file_directory, 'pos_2DL1S1_' + sample.name + '_1.fastq.gz')
    bt2_2 = os.path.join(sample.temp_file_directory, 'pos_2DL1S1_' + sample.name + '_2.fastq.gz')
    bt2_un_conc = os.path.join(sample.temp_file_directory, 'neg_2DL1S1_' + sample.name + '_%.fastq.gz')

    InputArgs = ['bowtie2', '-p', bt2_p, '--trim5', bt2_5, '--trim3', bt2_3, '-L', bt2_L, '-i', bt2_i, '--score-min', bt2_min_score,
    '-I', bt2_I, '-X', bt2_X, '-x', bt2_x, '-1', bt2_1, '-2', bt2_2, '-S', bt2_stream, '--un-conc-gz', bt2_un_conc]

    print('\n\nFiltering 2DL1S1 negative reads from: ' + sample.name + '\n')
    InputString = ' '.join(InputArgs)
    print(InputString)
    subprocess.call(InputArgs)
    os.remove(bt2_stream)
    os.remove(bt2_1)
    os.remove(bt2_2)

    # 2DL1 final alignments
    bt2_min_score = 'L,0,-0.55'
    bt2_x = 'Reference/KIR2DL1/KIR2DL1'
    bt2_1 = os.path.join(sample.temp_file_directory, 'neg_2DL1S1_' + sample.name + '_1.fastq.gz')
    bt2_2 = os.path.join(sample.temp_file_directory, 'neg_2DL1S1_' + sample.name + '_2.fastq.gz')
    bt2_al_conc = os.path.join(results_directory, sample.name + '_2DL1_%.fastq.gz')

    InputArgs = ['bowtie2', '-p', bt2_p, '--trim5', bt2_5, '--trim3', bt2_3, '-L', bt2_L, '-i', bt2_i, '--score-min', bt2_min_score,
    '-I', bt2_I, '-X', bt2_X, '-x', bt2_x, '-1', bt2_1, '-2', bt2_2, '-S', bt2_stream, '--al-conc-gz', bt2_al_conc]

    print('\n\nFinal 2DL1 alignmetn for: ' + sample.name + '\n')
    InputString = ' '.join(InputArgs)
    print(InputString)
    subprocess.call(InputArgs)
    os.remove(bt2_stream)
    os.remove(bt2_1)
    os.remove(bt2_2)
