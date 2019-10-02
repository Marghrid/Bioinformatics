# Smith-Waterman algorithm

import bl50  # file with BLOSUM 50 scoring matrix
import sys   # to read main arguments

def traceback(F, F_sources, F_max):
	# From each optimal matrix value in F_max, backtrack to find
	#  an optimal matching.
	paths = []

	for F_max_i in F_max:
		
		# Mark all cells as not visited
		F_visited = []
		for i in range(len(F)):
			F_visited.append([False] * len(F[0]))

		path = []
		# Call the recursive helper function to get all paths
		traceback_rec(F_max_i, F, F_sources, F_visited, path)
		paths.append(path)

	return paths


# A recursive helper function to print all optimal paths from source. 
# F_visited[] keeps track of vertices in current path. 
# path[] stores actual vertices.
def traceback_rec(source, F, F_sources, F_visited, path):
	# Mark the current node and store it in path[] 
	F_visited[source[0]][source[1]] = True
	#path[path_index] = source # append source to path
	path.append(source)
  
	# If current vertex is same as destination, then print 
	# current path[] 
	if F[source[0]][source[1]] == 0: # end of matching
		return path

	else:
		# Recur for all the vertices adjacent to current vertex 
		for prev in F_sources[source[0]][source[1]]:
			if not F_visited[prev[0]][prev[1]]:
				traceback_rec(prev, F, F_sources, F_visited, path)
				return
  
	# Remove current vertex from path[] and mark it as unvisited 
	path.pop()
	F_visited[source[0]][source[1]] = False

def print_table(S1, S2, F):
	# print table:

	# print S1 in first row
	print(' '*8, end = '')
	for char in S1:
		print('{:3s}'.format(char), end=' ')
	print('\n')


	for i in range(len(S2) + 1):
		# print S2 in first column
		if i > 0:
			print('{:2s}'.format(S2[i-1]), end='')
		else:
			print(' '*2, end = '')

		# table values
		for j in range(len(S1) + 1):
			print('{:3d}'.format(F[i][j]), end=' ')

		print('\n')

# Converts the paths found into a string showing the optimal local alignments.
def get_optimal_alignments(S1, S2, paths):
	optimal_alignments = ''
	for path in paths:
		M1 = ''
		M2 = ''
		match = ''
		lastM1 = -1
		lastM2 = -1
		for cell in reversed(path[:-1]):
			if cell[1] == lastM1:
				M1 += '-'
				M2 += S2[cell[0]-1]
				match += ' '
			elif cell[0] == lastM2:
				M2 += '-'
				M1 += S1[cell[1]-1]
				match += ' '
			else:
				M1 += S1[cell[1]-1]
				M2 += S2[cell[0]-1]
				match += '|'
				lastM1 = cell[1]
				lastM2 = cell[0]

		optimal_alignments += 'S1: ' + M1 + '\n'
		optimal_alignments += '    ' + match + '\n'
		optimal_alignments += 'S2: ' + M2 + '\n'

	return optimal_alignments

def smith_waterman(S1, S2, gap_pen):
	F = []          # matrix values indexed i, j, i in S2, j in S1.
	F_sources = []  # "arrows" in the matrix. Indexed like F.
	F_max = []      # indices of matrix values that are maximum
	F_max_val = 0   # maximum of matrix
	optimal_alignments = []

	# initialize all matrix F to 0
	for i in range(len(S2)+1):
		F.append([0] * (len(S1)+1))

	# initialize all matrix F_sources to empty lists
	for i in range(len(S2)+1):
		F_sources.append([])
		for j in range(len(S1)+1):
			F_sources[i].append([])


	for i in range(1, len(S2) + 1):
		for j in range(1, len(S1) + 1):
			# get score from BLOSSUM 50 scoring matrix
			score = bl50.scores.get((S1[j-1], S2[i-1]))

			# BLOSSUM 50 scoring matrix is not mirrored.
			# Try the inverse combination
			if score is None:
				score = bl50.scores.get((S2[i-1], S1[j-1]))

			# score between 2 aminoacids not in BLOSSUM 50 scoring matrix
			if score is None:
				raise Exception('Score between {} and {} not available'\
					.format(S2[i-1], S1[j-1]))

			F[i][j] = max(F[i-1][j-1] + score, F[i-1][j] - gap_pen,\
						  F[i][j-1] - gap_pen, 0)

			# save all optimal sources
			if F[i][j] == F[i-1][j-1] + score:
				F_sources[i][j].append((i-1, j-1))
			if F[i][j] == F[i-1][j] - gap_pen:
				F_sources[i][j].append((i-1, j))
			if F[i][j] == F[i][j-1] - gap_pen:
				F_sources[i][j].append((i, j-1))

			# update maximum values
			if F[i][j] == F_max_val:
				F_max.append((i, j))
			elif F[i][j] > F_max_val:
				F_max = [(i, j)]
				F_max_val = F[i][j]

	print_table(S1, S2, F)

	paths = traceback(F, F_sources, F_max)
	print(paths)

	optimal_alignments = get_optimal_alignments(S1, S2, paths)

	return F_max_val, optimal_alignments




if __name__ == "__main__":
	S1 = 'HEAGAWGHEE'
	S2 = 'PAWHEAE'
	gap_pen = 8

	if len(sys.argv) > 2:
		S1 = sys.argv[1]
		S2 = sys.argv[2]

	if len(sys.argv) > 3:
		gap_pen = int(sys.argv[3])

	print("Smith-Waterman Algorithms -- find all optimal local alignments")
	print("S1:", S1)
	print("S2:", S2)
	print("Gap penalty:", gap_pen)

	sol = smith_waterman(S1, S2, gap_pen)
	
	print('Score of optimal alignment(s):', sol[0])
	print('Optimal alignment(s):')
	print(sol[1])

