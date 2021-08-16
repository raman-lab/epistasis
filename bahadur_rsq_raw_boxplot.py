
import matplotlib
matplotlib.use('agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from collections import OrderedDict
import dask.dataframe as dd

def read_input(input):
	test = pd.read_csv(input,header=0,index_col=0,sep='\t')
	test = test.dropna(axis=1,how='all')
	test = test.dropna(axis=0,how='any')
	return test

def read_bootstrap_files(input):
	test = dd.read_csv(input,header=0)
	test = test.drop(test.columns[0],axis=1)
	test_avg = pd.DataFrame(np.mean(test,axis=0).compute())
	return test_avg

def parse_input(input_list,control):
	out_dict = OrderedDict()
	for file in input_list:
		file_list = file.split('_rsq_')[0].split('_')
		mut = '{0}_{1}_{2}_{3}'.format(file_list[-4],file_list[-3],file_list[-2],
			file_list[-1])
		out_dict[mut]=read_input(file)
	if control is not None:
		out_dict['control']=read_input(control)
	return out_dict

def parse_combinations(input_list,control,bootstrap):
	binary_str = OrderedDict()
	key_dict = OrderedDict()
	for file in input_list:
#		print file
		file_list = file.split('_rsq_')[0].split('_')
		mut_count = file_list[-4].count('1')
		filename = '{0}_{1}_{2}_{3}'.format(file_list[-4],
				file_list[-3],file_list[-2],file_list[-1])
		'''
		if ci == False:
			mut_count = file_list[-6].count('1')
			filename = '{0}_{1}_{2}_{3}'.format(file_list[-6],
				file_list[-5],file_list[-4],file_list[-3])
		else:
			mut_count = file_list[-8].count('1')
			filename = '{0}_{1}_{2}_{3}'.format(file_list[-8],
				file_list[-7],file_list[-6],file_list[-5])
		'''
		if mut_count in key_dict.keys():
			key_dict[mut_count].append(filename)
		else:
			key_dict[mut_count]=[filename]
		if bootstrap == False:
			binary_str[filename] = read_input(file)
		else:
			binary_str[filename] = read_bootstrap_files(file)
	if control is not None:
		if bootstrap == False:
			binary_str['control']=read_input(control)
		else:
			binary_str['control']=read_bootstrap_files(control)
		finalk = key_dict.keys()[-1]
		key_dict[finalk].append('control')
	else:
		pass
	return binary_str, key_dict

def concatenate_dfs_separate(input_dict,mut_dict):
	out_df = pd.DataFrame()
	for y in mut_dict:
		for x in mut_dict[y]:
			out_df[x]=input_dict[x].iloc[:,0]
	print out_df.columns
	return out_df

def box_plot(out_dict,basename,notext,single,use_range):
	hfont = {'fontname': 'Helvetica'}
	if single == False:
		to_plot = [np.array(out_dict[x].iloc[:,0]) for x in out_dict]
	else:
		only_df = out_dict[out_dict.keys()[0]]
		to_plot = [np.array(only_df.iloc[:,x]) for x in range(0,only_df.shape[1])]
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
#	print plot_vals
	if use_range == False:
		ax.boxplot(to_plot,sym='.')
	else:
		ax.boxplot(to_plot,whis='range')
	if single == False:
		if notext == False:
			ax.set_xticklabels(out_dict.keys(),fontsize=8,rotation=45,ha='right',**hfont)
		else:
			ax.tick_params(labelbottom=False)	
	else:
		if notext == False:
			ax.set_xticklabels(only_df.columns.values)
		else:
			ax.tick_params(labelbottom=False)	
	plt.tight_layout()
	if notext==False:
		plt.savefig('{0}_boxplot_labeled.pdf'.format(basename),transparent=True)
	else:
		plt.savefig('{0}_boxplot.pdf'.format(basename),transparent=True)
	plt.close()

def box_plot_separate(out_dict,mut_dict,basename,notext,use_range):
	hfont = {'fontname': 'Helvetica'}
	to_plot = []
	leg_lst = []
	for k in mut_dict:
		for x in range(0,len(mut_dict[k])):
			leg_lst.append(mut_dict[k][x])
			to_plot.append(out_dict[mut_dict[k][x]].iloc[:,0])
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	if use_range == False:
		ax.boxplot(to_plot,sym='.')
	else:
		ax.boxplot(to_plot,whis='range')
	if notext == False:
		ax.set_xticklabels(leg_lst,fontsize=8,rotation=45,ha='right',**hfont)
	else:
		ax.tick_params(labelbottom=False)
	plt.tight_layout()
	if notext==False:
		plt.savefig('{0}_reorder_boxplot_labeled.pdf'.format(basename),transparent=True)
	else:
		plt.savefig('{0}_reorder_boxplot.pdf'.format(basename),transparent=True)
	plt.close()
		

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to make *rsq_raw.tsv files into a boxplot')
	parser.add_argument('--input','-i',nargs='*',required=True,
		help='Input *rsq_raw files. Can take multiple')
	parser.add_argument('--title','-t',default='rsq_raw',
		help='Name for the output boxplot')
	parser.add_argument('--control','-c',required=False,default=None,
		help='Control rsq_raw file for comparison')
	parser.add_argument('--notext','-n',required=False,default=False,action='store_true',
		help='Call this flag to remove text')
	parser.add_argument('--single','-s',required=False,default=False,action='store_true',
		help='Call this flag to plot all columns from a single file')
	parser.add_argument('--separate','-e',required=False,default=False,action='store_true',
		help='Call this flag to reorder box plot based on number of mutants')
	parser.add_argument('--bootstrap','-b',required=False,default=False,action='store_true',
		help='Call this flag if you want to make a boxplot of bootstrap averages')
	parser.add_argument('--range','-r',required=False,default=False,action='store_true',
		help='Call this flag to force whiskers to the range of the data.')
	parsed = parser.parse_args()
	if parsed.separate==False:
		if parsed.single == False:
			in_dict = parse_input(parsed.input,parsed.control)
			box_plot(in_dict,parsed.title,parsed.notext,parsed.single,parsed.range)
		else:
			in_dict = parse_input(parsed.input,parsed.control)
			box_plot(in_dict,parsed.title,parsed.notext,parsed.single,parsed.range)
	else:
		binary_str,key_dict = parse_combinations(parsed.input,parsed.control,parsed.bootstrap)
		print key_dict
		print binary_str.keys()
		box_plot_separate(binary_str,key_dict,parsed.title,parsed.notext,parsed.range)
		
	