#!/usr/bin/python3
import os
import glob

tools = ['kalign', 'muscle', 'clustalo']

print("Looking for benchmarks...")
benchmarks = os.listdir('benchmarks')

print("Found", list(benchmarks))
print()

instances = {}
for bench in benchmarks:
    print("Looking for instances in benchmarks/", bench, sep="")
    files = glob.glob('benchmarks/' + bench + '/**/*.fa', recursive=True)
    files += glob.glob('benchmarks/' + bench + '/**/*.tfa', recursive=True)
    print("\tFound", len(files), "instances.")
    instances[bench] = files

