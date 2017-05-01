from load_info import load_profile

def extract_pre_flank (exp_dir, length):
	arg_dict = load_profile(exp_dir)
	exp_name = arg_dict['exp_name']
	temp_fa_dir = exp_dir + '/temp/' + exp_name
	with open(temp_fa_dir + '_prefix.fa', 'r') as f:
		f.readline()
		pref = f.readline()[-length - 1:].strip()
	return pref

def extract_post_flank (exp_dir, length):
	arg_dict = load_profile(exp_dir)
	exp_name = arg_dict['exp_name']
	temp_fa_dir = exp_dir + '/temp/' + exp_name
	with open(temp_fa_dir + '_suffix.fa', 'r') as f:
		f.readline()
		post = f.readline()[:length].strip()
	return post



def extract_pre_post_flank(exp_dir, length):
	return extract_pre_flank(exp_dir, length), \
			extract_post_flank(exp_dir, length)