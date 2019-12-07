from vie_dehyphenator import dehyphenator

import glob, os
os.chdir("HomFam")
for file in glob.glob("*.vie"):
	file_no_term = file.replace('.vie', '')
	dehyphenator(file, file_no_term + ".tfa")