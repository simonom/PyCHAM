'''code to remove the out string from the end of component names in the first column of a spreadsheet'''

def str_out_rem(): # define function

	# import dependencies------------
	import openpyxl
	import os
	# -------------------------------

	# state name of spreadsheet
	fname = 'Outdoor concentrations.xlsx'
	# path name
	pname = str(os.getcwd() + '/' + fname)
	
	wb = openpyxl.load_workbook(filename = pname)
	
	sheet = wb['Sheet1']
	# component names are in first column, times are in headers of first row		
	ic = -1 # count on row iteration
	
	for i in sheet.iter_rows(values_only=True): # loop through rows

		ic += 1 # count on row iteration

		if (ic == 0): # skip header
			continue
				
		# get names of components (matching chemical scheme names with out appended) 
		else:
			# component name
			con_infl_nam = (i[0])[0:-3]
				
			# replace without the out extension
			cell = str('A' + str(ic+1))
			sheet[cell] = con_infl_nam
				
	wb.save(filename = fname)
	wb.close() # close excel file


	return # end funtion

str_out_rem() # call function