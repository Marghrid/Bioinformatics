#!/usr/bin/python3
import os
import glob
import subprocess
import re
import csv


def output(file):
    filename, file_extension = os.path.splitext(file)
    return filename + '.out'


def solution(file):
    filename, file_extension = os.path.splitext(file)
    return filename + '.msf'


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

counter = 0
with open('results.csv', 'w') as out_file:
    csv_out = csv.writer(out_file)
    csv_out.writerow(['aligner', 'instance', 'time', 'ram', 'SP', 'TC'])

    print('Running kalign...')
    for name, benchmark in instances.items():
        print("\tRunning instances of", name, '...')

        for instance in benchmark:
            p = subprocess.Popen(['/usr/bin/time', '-f', '%e %M', 'kalign', instance, '-o', output(instance), '--format', 'msf'],
                                 stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            p.wait()
            out = p.stderr.read().split()
            time = float(out[0])
            ram = float(out[0])

            p = subprocess.Popen(['./bali_score', solution(instance), output(instance)], stdout=subprocess.PIPE)
            p.wait()
            out = p.stdout.read().decode()

            regex_out = SP_regex.search(out)
            if regex_out:
                SP = regex_out.group(1)
            else:
                SP = -1

            regex_out = TC_regex.search(out)
            if regex_out:
                TC = regex_out.group(1)
            else:
                TC = -1

            csv_out.writerow(('kalign', instance, time, ram, SP, TC))
            counter += 1
            print('\t(' + str(counter) + '/' + str(total) + ')')

        out_file.flush()

