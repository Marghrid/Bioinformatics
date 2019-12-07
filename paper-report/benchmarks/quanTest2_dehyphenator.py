from vie_dehyphenator import dehyphenator

import glob, os
os.chdir("QuanTest2/Test")
for file in glob.glob("*.vie"):
	file_no_term = file.replace('.vie', '')
	dehyphenator(file, file_no_term + ".tfa")