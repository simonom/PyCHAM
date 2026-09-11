file_name = str('/Users/user/Documents/GitHub/PyCHAM/PyCHAM/input/Ng2007/' +
	'benzene_MCM_autoAPRAMfw.eqn')
fo = open(file_name, 'r+')

# get each line as a subsequent item in a list
lines = fo.readlines()
fo.close()


ic = -1 # line count
for i in lines: # loop through lines

	ic += 1 # line count

	if '<' not in i[0]: # if not an equation line
		continue
	if ':' not in i:
		
		ns = (str(': ') + i.split(' ')[-1])

		i = i.replace(i.split(' ')[-1], ns)
		
	if ';' not in i:
		
		ns = (i.split(' ')[-1]).replace('\n', ' ;\n')

		i = i.replace(i.split(' ')[-1], ns)
	
	lines[ic] = i

fo = open(file_name, 'w+')	
# write out updated text
for line in lines:
	fo.write(line)

# save and close file
fo.close()
