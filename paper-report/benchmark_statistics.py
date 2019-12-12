import os
import glob
from statistics import mean, median

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

for name, benchmark in instances.items():
	seq_counts = []
	for instance in benchmark:
		seq_count = 0
		for line in open(instance):
			if line.startswith('>'):
				seq_count += 1

		seq_counts.append(seq_count)
	print(name, mean(seq_counts), median(seq_counts))