def compare_files(file1, file2, ignore_first):

	with open(file1) as f1, open(file2) as f2:

		lines1 = f1.readlines()[ignore_first:]
		lines2 = f2.readlines()[ignore_first:]

		return lines1 == lines2