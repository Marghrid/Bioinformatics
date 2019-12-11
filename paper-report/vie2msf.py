#!/usr/bin/python3
import glob
import os


def dehyphenator(vie_file_path, new_file_path):
    with open(vie_file_path, 'r') as file:
        lines = file.readlines()
    with open(new_file_path, 'w+') as file:
        for line in lines:
            if not line.startswith('>'):
                line = line.replace('-', '')
            file.write(line)


# files = glob.glob('**/mattbench/**/*.tfa', recursive=True)
# print(files)
#
# for file in files:
#     filename, file_extension = os.path.splitext(file)
#
#     subprocess.Popen(['seqret', file, '-osf', 'fasta', '-outseq', filename + '.fasta'])

for file in glob.glob("mattbench/**/*.fasta", recursive=True):
    filename, file_extension = os.path.splitext(file)
    dehyphenator(file, filename + ".tfa")

# files = glob.glob('**/bralibase/**/structural/*.fa', recursive=True)
# print(files)
#
# for file in files:
#     filename, file_extension = os.path.splitext(file)
#
#     subprocess.Popen(['seqret', file, '-osf', 'msf', '-outseq', filename + '.msf'])

# files = glob.glob('mattbench/mattbench/**/*.tfa', recursive=True)
# print(files)
#
# for file in files:
#     filename, file_extension = os.path.splitext(file)
#
#     subprocess.Popen(['mv', file, os.path.split(os.path.split(filename)[0])[0] + '/' + os.path.split(filename)[1] + '.tfa'])
