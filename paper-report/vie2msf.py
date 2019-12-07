#!/usr/bin/python3
import glob
import os
import subprocess

files = glob.glob('**/*.vie', recursive=True)
print(files)

for file in files:
    filename, file_extension = os.path.splitext(file)

    subprocess.Popen(['seqret', file, '-osf', 'msf', '-outseq', filename + '.msf'])
