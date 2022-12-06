#get file object
file1 = open("raw_targets.txt", "r")

#read all lines
lines = file1.readlines()

#traverse through lines one by one
with open('rename_targets.txt', 'a') as file2:
	for line in lines:
		print(line.strip())
		target = line.strip()
		rename = target[:2] + target[3:] + '/'
		file2.write(rename + '\n')

file1.close
file2.close