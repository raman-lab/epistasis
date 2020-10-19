#Nishikawa, Kyle K.
#Edited 20200826 for Epistasis manuscript

#This script uses subnetwork fold induction values to calculate types of epistasis. It is written to use *_mod_input.tsv files from the Bahadur expansion script.

import pandas as pd
import numpy as np
import argparse


def read_input(input_file):
	if input_file.endswith('.tsv'):
		binary_values = pd.read_csv(input_file,delimiter='\t',index_col=0,dtype='str')
		test = pd.read_csv(input_file,delimiter='\t',index_col=0)
		if 'binary' in test.columns.values:
			binary_loc, = np.where(test.columns.values == 'binary')
			test['binary']=binary_values['binary']
		else:
			pass
		if 'substrings' in test.columns.values:
			substring_loc, = np.where(test.columns.values == 'substrings')
			test['substrings']=binary_values['substrings']
		else:
			pass
	else:
		binary_values = pd.read_csv(input_file,index_col=0,dtype='str')
		test = pd.read_csv(input_file,index_col=0)
		if 'binary' in test.columns.values:
			binary_loc, = np.where(test.columns.values == 'binary')
			test['binary']=binary_values['binary']
		else:
			pass
		if 'substrings' in test.columns.values:
			substring_loc, = np.where(test.columns.values == 'substrings')
			test['substrings']=binary_values['substrings']
		else:
			pass
	return test

def epistasis_types(input_df,rsq_df,use_col):
	no_epi = [1,0,0,0]
	magnitude = [0,1,0,0]
	sign_epi = [0,0,1,0]
	reciprocal_sign = [0,0,0,1]
	rsq_val = rsq_df.iat[0,0]
	base_mut_ndx, = np.where(input_df['substrings']=='00')
	first_mut_ndx, = np.where(input_df['substrings']=='10')
	second_mut_ndx, = np.where(input_df['substrings']=='01')
	double_mut_ndx, = np.where(input_df['substrings']=='11')
	base_mut_val = input_df.iat[base_mut_ndx[0],use_col]
	first_mut_val = input_df.iat[first_mut_ndx[0],use_col]
	second_mut_val = input_df.iat[second_mut_ndx[0],use_col]
	double_mut_val = input_df.iat[double_mut_ndx[0],use_col]
	#set up two cases of epistasis
	if rsq_val < 0.9:
		if (first_mut_val > base_mut_val and second_mut_val < base_mut_val) or (first_mut_val < base_mut_val and second_mut_val > base_mut_val):
			if (double_mut_val < first_mut_val and double_mut_val > second_mut_val) or (double_mut_val > first_mut_val and double_mut_val < second_mut_val):
				return magnitude
			elif (double_mut_val > first_mut_val and double_mut_val > second_mut_val) or (double_mut_val < first_mut_val and double_mut_val < second_mut_val):
				return sign_epi
			else:
				print "shouldn't happen"
		elif (first_mut_val > base_mut_val and second_mut_val > base_mut_val):
			if (double_mut_val > first_mut_val and first_mut_val > second_mut_val) or (double_mut_val > second_mut_val and second_mut_val > first_mut_val):
				return magnitude
			elif (double_mut_val > first_mut_val and double_mut_val < second_mut_val) or (double_mut_val < first_mut_val and double_mut_val > second_mut_val):
				return sign_epi
			elif double_mut_val < first_mut_val and double_mut_val < second_mut_val:
				return reciprocal_sign
			else:
				print "shouldn't happen"
		elif (first_mut_val < base_mut_val and second_mut_val < base_mut_val):
			if (double_mut_val < first_mut_val and first_mut_val < second_mut_val) or (double_mut_val < second_mut_val and second_mut_val < first_mut_val):
				return magnitude
			elif (double_mut_val < first_mut_val and double_mut_val > second_mut_val) or (double_mut_val > first_mut_val and double_mut_val < second_mut_val):
				return sign_epi
			elif double_mut_val > first_mut_val and double_mut_val > second_mut_val:
				return reciprocal_sign
			else:
				print "shouldn't happen"
		else:
			print "shouldn't happen"
	else:
		return no_epi

def count_epistasis(input_list,rsq_list,use_col):
	out_df = pd.DataFrame()
	for num in range(0,len(input_list)):
		file_list = input_list[num].split('_mod_')[0].split('_')
		filename = '{0}_{1}_{2}_{3}'.format(file_list[-4],
				file_list[-3],file_list[-2],file_list[-1])
		print filename
		rsq_name = [s for s in rsq_list if filename in s]
		input_df = read_input(input_list[num])
		rsq_df = read_input(rsq_name[0])
		out_df[filename] = epistasis_types(input_df,rsq_df,use_col)
	out_df.index = ['No_epistasis','Magnitude_epistasis','Sign_epistasis','Reciprocal_sign_epistasis']
	out_df = out_df.transpose()
	counts = pd.DataFrame(out_df.sum(axis=0))
	counts.columns = ['Count']
	print counts
	return out_df,counts
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to calculate the number of epistatic subnetworks based on output from bahadur_expansion3.py with the -k subnetwork flag')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,nargs='*',
		help='One or more input files of the *mod_input.tsv from bahadur_expansion.py output')
	required.add_argument('--rsq','-r',required=True,nargs='*',
		help='An *rsq_avg.tsv file for each input file')
	parser.add_argument('--title','-t',required=False,default='bahadur_epistasis_counter',
		help='Title for output files')
	parser.add_argument('--use_col','-u',required=False,default=1,
		help='Column to use, varies depending on input (1 for no transform, 4 for transform')
	parsed = parser.parse_args()
	out_df,counts_df = count_epistasis(parsed.input,parsed.rsq,int(parsed.use_col))
	out_df.to_csv('{0}_counts_by_rsq.csv'.format(parsed.title))
	counts_df.to_csv('{0}_epistasis_types.csv'.format(parsed.title))
	