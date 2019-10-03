import bl50
import random
import time
import subprocess
import os

instances_per_size = 10
min_size = 10
max_size = 1000
size_int = 10
file_name = "results.txt"

open("results.csv", "w+")

for i in range(instances_per_size):
	for size in range(min_size, max_size + 1, size_int):
		S1 = ''
		S2 = ''
		for c in range(size):
			pair = random.choice(list(bl50.scores.keys()))
			S1 += pair[0]
			S2 += pair[1]
		
		with open(os.devnull, "w") as fnull:
			start_time = time.time()
			subprocess.call(["python3", "smith_waterman.py", S1, S2, "8"], stdout = fnull)
			end_time = time.time() - start_time
		with open(file_name, "a") as fres:
			fres.write("{} {}\n".format(size, end_time))

	print(i+1, "done")
