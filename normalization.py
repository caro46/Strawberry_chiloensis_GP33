"""
Python script to apply normalization to data

Possible modes of execution (specify as 1st arg on command line):

means:
	Calculate the means from a set of files and store them internally
	Stores in means.db by default
	
	ex:
	python3 normalization.py means file0 file1 file2

ratio:
	Use previously calculated means to determine ratios
	methods prints results as tab-sep data to the command line

	ex:
	python3 normalization.py ratio name_file0 name_file1 file2 file3

print:
	Helper method to print the stored means

	ex:
	python3 normalization.py print

merge:
	Helper method to merge ratio results from multiple files
	Similar to cat but only one header line is saved

	ex:
	python3 normalization.py merge file0 file1 file2

window:
	Calculates window averages from ratio files
	Note: non-overlapping windows are used.

	ex:
	python3 normalization.py window window_size file0 file1 file2

"""


import sys
import shelve


MEANS_DB = "means.db"


def shelveMeans(means):
	""" Function to save calculated means to a shelve database: MEANS_DB """
	# Open MEANS_DB as a shelve for writing
	with shelve.open(MEANS_DB) as d:
		# Iterate through each name, mean in means dictionary
		for name, u in means.items():
			# Set mean into the database
			d[name] = u


def returnMeans():
	""" Function to return calculated means from the shelve database """
	# Open MEANS_DB as a shelve for reading data
	with shelve.open(MEANS_DB) as d:
		# Iterate through means dictionary and create a copy to return
		return {k: v for k, v in d.items()}


def calculateColumnMeans(filename):
	""" Function to calculate column means from a file """
	# Store calculated values
	omit = ["#CHROM", "POS"]
	header = None
	sums = None
	count = 0
	# Open file for reading an iterate through each line
	for line in open(filename):
		# Read in header or main body
		if header is None:
			# Read in the header (assumes tab-seperation)
			header = line.strip().split('\t')
			# Set column sums dictionary (initially zero for each column)
			sums = {h: 0 for h in header if h not in omit} 
		else:
			# Reading in main body fields (assumes tab-seperation)
			fields = line.strip().split('\t')
			# Iterate through each field to update sums
			for k, v in zip(header, fields):
				# Update sum for this column if column is not omitted 
				if k not in omit:
					sums[k] += int(v)
			# Update total line count
			count += 1
	# Convert column sums to averages and return dictionary of results
	return {k: v/count for k, v in sums.items()}


def extractNames(filename):
	""" Function to read in names from a file """
	# Store a set of names found (silences duplicate entries)
	names = set()
	# Open file for reading and iterate through each line
	for line in open(filename):
		# Strip new line characters and spaces from end
		line = line.strip()
		# Add line to names if not blank
		if len(line):
			names.add(line)
	# Return the set of names found in this file
	return names


def readFile(filename):
	""" Function to read through a file and return values line by line """
	# Variable to store header of file
	header = None
	# Open file for reading and iterate through each line of the file
	for line in open(filename):
		# Read in header or body
		if header is None:
			# Read in header and store as list (assumes tab-sep)
			header = line.strip().split('\t')
		else:
			# Read in row and store values in a list
			fields = []
			for value in line.strip().split('\t'):
				# Try to convert field values to an integer, 
				# If that fails, just add the value
				try:
					fields.append(int(value))
				except ValueError:
					fields.append(value)
			# Yield row data as a dictionary
			yield {k: v for k, v in zip(header, fields)}


def mergeRead(file0, file1):
	""" Function to merge through 2 files in unison """
	# Iterate through rows of each file
	for d0, d1 in zip(readFile(file0), readFile(file1)):
		# Check valid (chrom and position must match)
		if d0["#CHROM"] != d1["#CHROM"]:
			raise ValueError("Mismatch in chromosome names")
		elif d0["POS"] != d1["POS"]:
			raise ValueError("Mismatch in position number")

		# Merge dictionaries and yield
		d0.update(d1)
		yield d0


def calculateRatio(males, females, file0, file1):
	""" Calculate male and female ratios from two files """
	# Get means from database (Must be known prior to calling function)
	means = returnMeans()

	# Iterate through merged files
	mapped_males = []
	mapped_females = []
	for data in mergeRead(file0, file1):
		# If mapped lists are empty: Find mapped males and mapped females
		if not mapped_males:
			# Iterate through each male in the males set
			for m in males:
				# Iterate through each name entry in the data dictionary
				for name in data:
					# Check if male name is found and add to mapped males
					if (m in name):
						mapped_males.append(name)
			assert(len(mapped_males) == len(males))
		if not mapped_females:
			# Iterate through each female in the females set
			for f in females:
				# Iterate through each name entry in the data dictionary
				for name in data:
					# Check if female name is found and add to mapped females
					if (f in name):
						mapped_females.append(name)
			assert(len(mapped_females) == len(females))	

		# Calculate ratios or data values relative to means
		male_values = [data[m]/means[m] for m in mapped_males]
		male_avg = sum(male_values) / len(male_values)
		female_values = [data[f]/means[f] for f in mapped_females]
		female_avg = sum(female_values) / len(female_values)
		yield (data["#CHROM"], data["POS"], female_avg, male_avg)


def windowAverages(filename, window):
	""" Function to calculate a window average for a file """
	# Store sums
	male_sum = 0
	female_sum = 0
	count = 0
	current_index = 0
	# Iterate through each row of the file
	for data in readFile(filename):
		# Get window index (i.e. the floor division of the pos and window) 
		index = (data["POS"]-1) // window
		# Split on index value
		if index == current_index:
			# If calculated index is the same as the current index (same window)
			# Update running sums and count
			male_sum += float(data["MALE_AVG"])
			female_sum += float(data["FEMALE_AVG"])
			count += 1
		else:
			# Calculated index is different (Start of next window)
			# Yield prev interval
			start = window * current_index + 1
			end = window * (current_index + 1)
			yield data["#CHROM"], start, end, female_sum/count, male_sum/count
	
			# Reset values
			male_sum = float(data["MALE_AVG"])
			female_sum = float(data["FEMALE_AVG"])
			count = 1
			current_index += 1

	# Yield last
	start = window * current_index + 1
	end = window * (current_index + 1)
	yield data["#CHROM"], start, end, female_sum/count, male_sum/count	


if __name__ == "__main__":
	# Import input mode from the command line (argument 1)
	import sys
	mode = sys.argv[1]

	# Split on mode type
	if mode == "means":
		# Calculate means from input files
		files = sys.argv[2:]
		for filename in files:
			# Get means
			means = calculateColumnMeans(filename)
			# Save means in a shelve database
			shelveMeans(means)
	elif mode == "ratio":
		# Get male, female sets
		males = extractNames(sys.argv[2])
		females = extractNames(sys.argv[3])
		chroms = sys.argv[4:6]

		# Calculate ratio and print results to the terminal
		print("#CHROM", "POS", "FEMALE_AVG", "MALE_AVG", sep='\t')
		for c, p, f, m in calculateRatio(males, females, *chroms):
			print(c, p, f, m, sep='\t')
	elif mode == "print":
		# Just print the previously calculated means
		print(returnMeans())
	elif mode == "merge":
		# Merge results into a single file
		ratio_files = sys.argv[2:]
		print("#CHROM", "POS", "FEMALE_AVG", "MALE_AVG", sep='\t')
		# Iterate through each file
		for filename in ratio_files:
			# Iterate through each line in the file
			for line in open(filename):
				# Parse the line
				line = line.strip()
				if not line.startswith("#"):
					# Print the line and format precision
					line = line.split()
					print(line[0], line[1], "{:.4f}".format(float(line[2])), "{:.4f}".format(float(line[3])))
	elif mode == "window":
		# Calculating window averages
		# Get window size
		window_size = int(sys.argv[2])
		# Get ratio files
		ratio_files = sys.argv[3:]
		print("#CHROM", "POS_START", "POS_END", "FEMALE_AVG", "MALE_AVG", sep='\t')
		# Iterate throguh each file
		for filename in ratio_files:
			# Iterate throguh window averages and print to terminal and error
			for values in windowAverages(filename, window_size):
				print(*values, sep='\t')
				print(*values, sep='\t', file=sys.stderr)

