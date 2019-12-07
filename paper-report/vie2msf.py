#!/usr/bin/python3
import glob
import os
import subprocess

# files = glob.glob('**/bralibase/**/structural/*.fasta', recursive=True)
# print(files)
#
# for file in files:
#     filename, file_extension = os.path.splitext(file)
#
#     subprocess.Popen(['seqret', file, '-osf', 'msf', '-outseq', filename + '.msf'])

# files = glob.glob('**/bralibase/**/structural/*.fa', recursive=True)
# print(files)
#
# for file in files:
#     filename, file_extension = os.path.splitext(file)
#
#     subprocess.Popen(['seqret', file, '-osf', 'msf', '-outseq', filename + '.msf'])

# files = glob.glob('**/bralibase/**/unaligned/*.fasta', recursive=True)
# print(files)
#
# for file in files:
#     filename, file_extension = os.path.splitext(file)
#
#     subprocess.Popen(['mv', file, os.path.split(os.path.split(filename)[0])[0] + '/' + os.path.split(filename)[1] + '.tfa'])
