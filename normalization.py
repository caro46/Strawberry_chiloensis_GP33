import sys
import shelve


MEANS_DB = "means.db"


def shelveMeans(means):
	with shelve.open(MEANS_DB) as d:
		for name, u in means.items():
			d[name] = u


def returnMeans():
	with shelve.open(MEANS_DB) as d:
		return {k: v for k, v in d.items()}


def calculateColumnMeans(filename):
	# Get columns from file
	omit = ["#CHROM", "POS"]
	header = None
	sums = None
	count = 0
	for line in open(filename):
		if header is None:
			header = line.strip().split('\t')
			sums = {h: 0 for h in header if h not in omit} 
		else:
			fields = line.strip().split('\t')
			for k, v in zip(header, fields):
				if k not in omit:
					sums[k] += int(v)
			count += 1
	return {k: v/count for k, v in sums.items()}


def extractNames(filename):
	names = set()
	for line in open(filename):
		line = line.strip()
		if len(line):
			names.add(line)
	return names


def readFile(filename):
	header = None
	for line in open(filename):
		if header is None:
			header = line.strip().split('\t')
		else:
			fields = []
			for value in line.strip().split('\t'):
				try:
					fields.append(int(value))
				except ValueError:
					fields.append(value)
			yield {k: v for k, v in zip(header, fields)}


def mergeRead(file0, file1):
	for d0, d1 in zip(readFile(file0), readFile(file1)):
		# Check valid
		if d0["#CHROM"] != d1["#CHROM"]:
			raise ValueError("Mismatch in chromosome names")
		elif d0["POS"] != d1["POS"]:
			raise ValueError("Mismatch in position number")

		# Merge dictionaries and yield
		d0.update(d1)
		yield d0


def calculateRatio(males, females, file0, file1):
	# Get means
	means = returnMeans()

	# Iterate through merged files
	mapped_males = []
	mapped_females = []
	for data in mergeRead(file0, file1):
		# Find mapping if None
		if not mapped_males:
			for m in males:
				for name in data:
					if (m in name):
						mapped_males.append(name)
			assert(len(mapped_males) == len(males))
		if not mapped_females:
			for f in females:
				for name in data:
					if (f in name):
						mapped_females.append(name)
			assert(len(mapped_females) == len(females))	

		# Calculate ratio
		male_values = [data[m]/means[m] for m in mapped_males]
		male_avg = sum(male_values) / len(male_values)
		female_values = [data[f]/means[f] for f in mapped_females]
		female_avg = sum(female_values) / len(female_values)
		yield (data["#CHROM"], data["POS"], female_avg, male_avg)


def windowAverages(filename, window):
	male_sum = 0
	female_sum = 0
	count = 0
	current_index = 0
	for data in readFile(filename):
		# Get index 
		index = (data["POS"]-1) // window
		if index == current_index:
			male_sum += float(data["MALE_AVG"])
			female_sum += float(data["FEMALE_AVG"])
			count += 1
		else:
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
	import sys
	mode = sys.argv[1]

	# Split on mode type
	if mode == "means":
		# Calculate means from input files
		files = sys.argv[2:]
		for filename in files:
			# Get means
			means = calculateColumnMeans(filename)
			# Save in database
			shelveMeans(means)
	elif mode == "ratio":
		# Get male, female sets
		males = extractNames(sys.argv[2])
		females = extractNames(sys.argv[3])
		chroms = sys.argv[4:6]

		# Calculate ratio
		print("#CHROM", "POS", "FEMALE_AVG", "MALE_AVG", sep='\t')
		for c, p, f, m in calculateRatio(males, females, *chroms):
			print(c, p, f, m, sep='\t')
	elif mode == "print":
		print(returnMeans())
	elif mode == "merge":
		ratio_files = sys.argv[2:]
		print("#CHROM", "POS", "FEMALE_AVG", "MALE_AVG", sep='\t')
		for filename in ratio_files:
			for line in open(filename):
				line = line.strip()
				if not line.startswith("#"):
					line = line.split()
					print(line[0], line[1], "{:.4f}".format(float(line[2])), "{:.4f}".format(float(line[3])))
	elif mode == "window":
		# Get window size
		window_size = int(sys.argv[2])
		# Get ratio files
		ratio_files = sys.argv[3:]
		print("#CHROM", "POS_START", "POS_END", "FEMALE_AVG", "MALE_AVG", sep='\t')
		for filename in ratio_files:
			for values in windowAverages(filename, window_size):
				print(*values, sep='\t')
				print(*values, sep='\t', file=sys.stderr)

