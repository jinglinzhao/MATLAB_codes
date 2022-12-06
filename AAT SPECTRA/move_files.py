import glob
import os
import shutil

#get file object
file1 = open("rename_targets.txt", "r")

#read all lines
lines = file1.readlines()

#traverse through lines one by one
for line in lines:
	print(line.strip())
	target = line.strip()
	FILE = glob.glob('./stars/'+ target + 'Extracted/' + '*txt')
	for i in range(len(FILE)):
		shutil.move(FILE[i], './stars/'+ target)
	shutil.rmtree('./stars/'+ target + 'Extracted', ignore_errors=True)
	shutil.rmtree('./stars/'+ target + '__MACOSX', ignore_errors=True)

file1.close