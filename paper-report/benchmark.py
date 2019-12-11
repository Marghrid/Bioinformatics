#!/bin/python3
import os
import glob
import subprocess
import re
import csv
import threading
import sys


def output(file, aligner):
    filename, file_extension = os.path.splitext(file)
    return filename + '.out' + '_' + aligner


def solution(file):
    filename, file_extension = os.path.splitext(file)
    return filename + '.msf'


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run_aligner(aligner, inflag, outflag, extraflags):
    eprint('Running', aligner, '...')

    counter = 0
    for name, benchmark in instances.items():
        for instance in benchmark:
            counter += 1
            p = subprocess.Popen(['timeout', '600', 'militime', aligner, inflag, instance, outflag, output(instance, aligner)]
                                 + extraflags, stdout=subprocess.PIPE)
            p.wait()

            if p.returncode != 124:
                out = p.stdout.read().split()
                try:
                    time = float(out[0])
                    ram = float(out[1])
                except:
                    eprint("Error while parsing time and ram.")
                    eprint(out)
                    time = -1
                    ram = -1

                p = subprocess.Popen(['./bali_score', solution(instance), output(instance, aligner)], stdout=subprocess.PIPE)
                p.wait()
                out = p.stdout.read().decode()
                
                try:
                    SP = SP_regex.search(out).group(1)
                except:
                    SP = -1

                try:
                    TC = TC_regex.search(out).group(1)
                except:
                    TC = -1

                eprint('\t', aligner, '(' + str(counter) + '/' + str(total) + ')')
                print(aligner, instance, time, ram, SP, TC, sep=',')

            else:
                eprint('\tInstance', instance, 'on', aligner, 'timed out.')
                print(aligner, instance, -1, -1, -1, -1, sep=',')


SP_regex = re.compile("SP score= (.*)")
TC_regex = re.compile("TC score= (.*)")

eprint("Looking for benchmarks...")
benchmarks = os.listdir('benchmarks')

eprint("Found", list(benchmarks))
eprint()

instances = {}
for bench in benchmarks:
    eprint("Looking for instances in benchmarks/", bench, sep="")
    if bench[0] == 'm':
        files = glob.glob('benchmarks/' + bench + '/**/*.tfa', recursive=True)
        eprint("\tFound", len(files), "instances.")
        if len(files) != 0:
            instances[bench] = files

eprint()
total = sum(map(len, instances.values()))
eprint('Found', total, 'instances in total.')
eprint()

print('aligner', 'instance', 'time', 'ram', 'SP', 'TC', sep=',')

threads = []
for aligner in [
    ('kalign', '-i', '-o', ['--format', 'msf']),
    ('kalign2', '-i', '-o', ['--format', 'msf', '--quiet']),
    ('clustalo', '-i', '-o', ['--outfmt', 'msf', '--threads', '8', '--MAC-RAM', '48000', '--iterations', '2', '--force']),
    ('muscle', '-in', '-out', ['-msf', '-maxiters', '2', '-quiet'])
]:
    t = threading.Thread(target=run_aligner, args=aligner)
    t.start()
    threads.append(t)

for t in threads:
    t.join()

