#!/bin/python

import bl50
import random
import subprocess

executable = './smith_waterman.py'

instances_per_size = 10
min_size = 10
max_size = 1000
step = 10

samples = 18
jobs = 6

result = open("results.data", "w+")


for size in range(min_size, max_size + 1, step):
	print("Running for size", size)

	for i in range(0, samples, jobs):
		processes = []
		for j in range(0, jobs):
			S1 = ''
			S2 = ''
			for c in range(size):
				pair = random.choice(list(bl50.scores.keys()))
				S1 += pair[0]
				S2 += pair[1]

			processes.append(
				subprocess.Popen(['/usr/bin/time -f "%e %M" ' + executable + ' ' + S1 + ' ' + S2 + ' 8 > /dev/null'],
								 stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True,
								 shell=True))

		for j in range(0, jobs):
			processes[j].wait()
			time_output = processes[j].stderr.read().split()

			result.write(
				str(size) + " "
				+ "{0:.5f}".format(float(time_output[0])) + " "
				+ "{0:.5f}".format(float(time_output[1])) + "\n")
			result.flush()

result.close()
