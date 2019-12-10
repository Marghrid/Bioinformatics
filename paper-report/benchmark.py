#!/bin/python3
import os
import glob
import subprocess
import re
import csv
import threading


def output(file, aligner):
    filename, file_extension = os.path.splitext(file)
    return filename + '.out' + '_' + aligner


def solution(file):
    filename, file_extension = os.path.splitext(file)
    return filename + '.msf'


results = []


def run_aligner(aligner, inflag, outflag, extraflags):
    print('Running', aligner, '...')

    counter = 0
    for name, benchmark in instances.items():
        for instance in benchmark:
            p = subprocess.Popen(['timeout', '600', '/usr/bin/time', '-f', '%e %M', aligner, inflag, instance, outflag, output(instance, aligner)]
                                 + extraflags, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            p.wait()

            if p.returncode != 124:
                out = p.stderr.read().split()
                try:
                    time = float(out[0])
                    ram = float(out[0])
                except:
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

                counter += 1
                results.append((aligner, instance, time, ram, SP, TC))
                print('\t', aligner, '(' + str(counter) + '/' + str(total) + ')')

            else:
                results.append((aligner, instance, -1, -1, -1, -1))
                print('\tInstance', instance, 'on', aligner, 'timed out.')


SP_regex = re.compile("SP score= (.*)")
TC_regex = re.compile("TC score= (.*)")

print("Looking for benchmarks...")
benchmarks = os.listdir('benchmarks')

print("Found", list(benchmarks))
print()

instances = {}
for bench in benchmarks:
    print("Looking for instances in benchmarks/", bench, sep="")
    files = glob.glob('benchmarks/' + bench + '/**/*.tfa', recursive=True)
    print("\tFound", len(files), "instances.")
    if len(files) != 0:
        instances[bench] = files

print()
total = sum(map(len, instances.values()))
print('Found', total, 'instances in total.')
print()

threads = []
for aligner in [
    ('kalign', '-i', '-o', ['--format', 'msf']),
    ('clustalo', '-i', '-o', ['--outfmt', 'msf', '--threads', '8', '--MAC-RAM', '48000', '--iterations', '2', '--force']),
    ('muscle', '-in', '-out', ['-msf', '-maxiters', '2', '-quiet'])
]:
    t = threading.Thread(target=run_aligner, args=aligner)
    t.start()
    threads.append(t)

for t in threads:
    t.join()

with open('results.csv', 'w') as out_file:
    csv_out = csv.writer(out_file)
    csv_out.writerow(['aligner', 'instance', 'time', 'ram', 'SP', 'TC'])

    for row in results:
        csv_out.writerow(row)
