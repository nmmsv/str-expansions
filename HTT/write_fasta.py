with open('target.fa', 'r') as f:
	row = f.readline()
	while row != '':	
		if row.strip() == '>HTT_PRE':
			pre = f.readline().strip()
		if row.strip() == '>HTT_POST':
			post = f.readline().strip()
		row = f.readline()

#str_len_list = [0, 1, 2, 4, 8, 12, 20, 25, 31, 35, 40, 45, 50]
str_len_list = [4, 8, 16]

for l in str_len_list:
	with open('fasta/HTT_str_' + str(l) + '.fa', 'w') as f:
		f.write('>HTT_str_' + str(l)+'\n')
		f.write(pre + 'CAG' * l + post)

