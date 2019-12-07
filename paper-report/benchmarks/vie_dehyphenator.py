def dehyphenator(vie_file_path, new_file_path):
	with open(vie_file_path, 'r') as file:
		lines = file.readlines()
	with open(new_file_path, 'w+') as file:
		for line in lines:
			if not line.startswith('>'):
				line = line.replace('-', '')
			file.write(line)

