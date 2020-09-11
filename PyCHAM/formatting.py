'''module is a collection of functions that are used to format chemical schemes, most of them (eg, MCM) are designed for fortran, which have a different format of function names of scientific notations. They must be converted to the format that python can understand'''

import re
import collections

# Remove the comments from the total eqn list and store them in a separate list. Effective 
# comments should start with '//'. 
# This function returns the list of the equations without the 'space' and '\n' markers, along with
# source and sink spcification. In another word, a clean list
def remove_comments(total_list):
    # remove all empty lines and all the '\n' markers at the end of each lines
    naked_list = [line.strip() for line in total_list if line != '\n']
    
    # remove all the comments and use a separate list to store them
    comments = [] # list to store the comments
    comment_marker = r"//.*"
    line_num = 0
    for line_num in range(len(naked_list)):
        if (re.match(comment_marker, naked_list[line_num]) != None):
            comments.append(naked_list[line_num] + '\t' + '(Ln:' + str(line_num) + ')' + '\n')
            naked_list[line_num] = ''
    naked_list = [line for line in naked_list if line != '']

    return naked_list

# convert scientific notation from D-type to E-type
def SN_conversion(RateExp):
	SN_d_regex = r"\d+(D|d)(-?|\+?)\d+"
	RegexList = re.findall(SN_d_regex, RateExp)
	if (RegexList != []): # in case more than one 'D' type instances are in the rate expression
		for i in range(len(RegexList)): # loop the checking until all 'D/d' instances are corrected
			m = re.search(SN_d_regex, RateExp)
			TempSlice = m.group(0) #  get the slice that RE engine has extracted from the rate coef. expression
			SliceStart = m.start() # the position and length are used to replace a part of the expression
			SliceEnd = m.end()
			# replace the 'D/d' with 'e' in the SLICE
			TempSlice = TempSlice.replace('d', 'e')
			TempSlice = TempSlice.replace('D', 'e')
			# reconstruct the rate conef string
			RateExp = RateExp[:SliceStart] + TempSlice + RateExp[SliceEnd:]
			
			# convert '@' marker to '**' (power) marker
			RateExp = RateExp.replace('@', '**')
	if re.search(r'@', RateExp)!=None: # ensure using @ to raise to a power converted
		RateExp = RateExp.replace('@', '**')
	return RateExp

# This function ensures that the function names can be recognized by python
# future development.  The output rate_coef has length equal to the number of
# equations and values giving the reaction rate coefficient.
def convert_rate_mcm(rate_coef):

    # create a list that contains commonly used math functions
    # ('illegal_func_name', func_regex,'legal_func_name')
	math_func_list = [
		('exp', 'numpy.exp'),
        ('EXP', 'numpy.exp'),
        ('dsqrt', 'numpy.sqrt'),
        ('dlog', 'numpy.log'),
        ('LOG', 'numpy.log'),
        ('dabs', 'numpy.abs'),
        ('LOG10', 'numpy.log10')
    ]

    
	# if no change is made, pass the rate expression to a new variable
	# because this new variable is used in the process of photolysis rate J(n)
	new_rate = rate_coef

	# 1. replace the illegal function names
	for func_step in range(len(math_func_list)):

		if (new_rate.find(math_func_list[func_step][0]) != -1):
			new_rate = new_rate.replace(math_func_list[func_step][0], math_func_list[func_step][1]) 
		else:
			continue

	# see if the rate expression contains photolysis rate; convert 'J(n)'/'J<n>' to 'J[n]'
	photolysis_rate_str = re.findall(r"J\(\d+\)|J\<\d+\>", new_rate)
	if (photolysis_rate_str != []): 
		for j_rate in photolysis_rate_str:
			new_j_rate = j_rate.replace('(', '[')
			new_j_rate = new_j_rate.replace('<', '[')
			new_j_rate = new_j_rate.replace(')', ']')
			new_j_rate = new_j_rate.replace('>', ']')
			new_rate = new_rate.replace(j_rate, new_j_rate)


    # The output rate_coef has length equal to the number of
    # equations and values giving the reaction rate coefficient
	return new_rate
