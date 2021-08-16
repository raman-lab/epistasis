#Python script to plot first-order rsq values from subnetworks analysis.
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import collections

def read_tsv(input_file):
	if input_file.lower().endswith('.tsv'):
		return pd.read_csv(input_file,sep='\t',header=0,index_col=0)
	elif input_file.lower().endswith('.csv'):
		return pd.read_csv(input_file,header=0,index_col=0)

#binary_str is the dict containing dataframes while key_dict matches the keys to the filename
def parse_combinations(input_list,control,ci):
	binary_str = collections.OrderedDict()
	key_dict = collections.OrderedDict()
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
		binary_str[filename] = read_tsv(file)
	if control is not None:
		binary_str['control']=read_tsv(control)
		for k in key_dict:
			key_dict[k].append('control')
	else:
		pass
	return binary_str, key_dict

def add_ci(binary_dict,key_dict,ci_dict):
	for k in binary_dict:
		binary_dict[k]['ci_low']=binary_dict[k].iloc[:,0]-ci_dict[k].iat[0,2]
		binary_dict[k]['ci_high']=ci_dict[k].iat[1,2]-binary_dict[k].iloc[:,0]
		print binary_dict[k]
	return binary_dict

#function for plotting whole bahadur expansion rsq with bias corrected CI
def read_full_rsq(input_list,ci_list):
	input_df = read_tsv(input_list[0])
	low = []
	high = []
	for ci in range(0,len(ci_list)):
		ci_df = read_tsv(ci_list[ci])
		low.append(input_df.iat[ci,0]-ci_df.iat[0,2])
		high.append(ci_df.iat[1,2]-input_df.iat[ci,0])
	input_df['ci_low']=low
	input_df['ci_high']=high
	return input_df
	
def concatenate_dfs(input_list,binary_str):
	out_df = pd.DataFrame()
	for x in range(0,len(input_list)):
		out_df[binary_str[x]] = input_list[x].iloc[0,:]
	print out_df
	return out_df

def concatenate_dfs_separate(input_dict,mut_list,rearrange,ci):
	out_df = pd.DataFrame()
	if rearrange == False:
		for x in mut_list:
			out_df[x]=input_dict[x].iloc[0,:]
	else:
		if 'control' in mut_list:
			for x in mut_list[:-1]:
				out_df[x]=input_dict[x].iloc[0,:]
			out_df.sort_values(axis=1,by='avg_0',ascending=True,inplace=True)
			out_df[mut_list[-1]]=input_dict[mut_list[-1]].iloc[0,:]
		else:
			for x in mut_list:
				out_df[x]=input_dict[x].iloc[0,:]
			out_df.sort_values(axis=1,by='avg_0',ascending=True,inplace=True)
	if ci == True:
		out_df.index = ['avg','stdev','num','ci_low','ci_high']
	else:
		out_df.index = ['avg','stdev','num']
	print out_df
	return out_df


def plot_rsq(out_df,xaxis_label,yaxis_label,separate,title,notext,use_ci):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	rownames = out_df.columns.values
	xvals = np.arange(0,out_df.shape[1])
	if use_ci == False:
		ax.bar(xvals+0.5,out_df.iloc[0,:],yerr=out_df.iloc[1,:],
			width=1,edgecolor='black',color='#77DD77',ecolor='black',
			linewidth=1,capsize=5,align='center',
			error_kw={'elinewidth': 1,'ecapthick': 2}
			)
	else:
		ax.bar(xvals+0.5,out_df.iloc[0,:],yerr=np.array(out_df.iloc[3:5,:]),
			width=1,edgecolor='black',color='#77DD77',ecolor='black',
			linewidth=1,capsize=5,align='center',
			error_kw={'elinewidth': 1,'ecapthick': 2}
			)
	ax.set_xlim(xmin=0)
	ax.set_xticks(xvals+0.5)
	if notext == False:
		ax.set_xticklabels(rownames,rotation=45,ha='right')
	else:
		ax.tick_params(labelbottom=False)
	ax.set_xlabel(xaxis_label)
	ax.set_ylabel(yaxis_label)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.suptitle('Rsq_values')
	plt.tight_layout()
	if notext == True:
		if separate == False:
			plt.savefig(title+'_selected_rsq_values.pdf')
		else:
			plt.savefig(title+'_rsq_bar.pdf')
	else:
		if separate == False:
			plt.savefig(title+'_selected_rsq_values_labeled.pdf')
		else:
			plt.savefig(title+'_rsq_bar_labeled.pdf')

def plot_same_axis(out_dict,xaxis_label,yaxis_label,title,notext,use_ci):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	rownames = []
	start = [0]
	end = []
	total_length = 0
	xtick_loc = []
	for x in range(0,len(out_dict.keys())):
		total_length += out_dict[x].shape[1]
		end.append(total_length)
		total_length +=1
		start.append(total_length)
		for row in out_dict[x].columns.values:
			rownames.append(row)
		xvals = np.arange(start[x],end[x])
		for test in xvals:
			xtick_loc.append(test+0.5)
		if use_ci == False:
			ax.bar(xvals+0.5,out_dict[x].iloc[0,:],yerr=out_dict[x].iloc[1,:],
				width=1,edgecolor='black',color='#77DD77',ecolor='black',
				linewidth=1,capsize=5,align='center',
				error_kw={'elinewidth': 1,'ecapthick': 2}
				)
		else:
			ax.bar(xvals+0.5,out_dict[x].iloc[0,:],yerr=np.array(out_dict[x].iloc[3:5,:]),
				width=1,edgecolor='black',color='#77DD77',ecolor='black',
				linewidth=1,capsize=5,align='center',
				error_kw={'elinewidth': 1,'ecapthick': 2}
				)
	ax.set_xticks(xtick_loc)
	if notext == False:
		ax.set_xticklabels(rownames,rotation=45,ha='right')
	else:
		ax.tick_params(labelbottom=False)		
	ax.set_xlim(xmin=0)
	ax.set_xlabel(xaxis_label)
	ax.set_ylabel(yaxis_label)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.suptitle('Rsq_values')
	plt.tight_layout()
	if notext == True:
		plt.savefig(title+'_combined_rsq_bar.pdf')
	else:
		plt.savefig(title+'_combined_rsq_bar_labeled.pdf')

def plot_full_rsq(in_df,xaxis_label,yaxis_label,figname,notext):
	fig,ax = plt.subplots()
	fig.set_size_inches(8,6)
	xvals = np.arange(in_df.shape[0])
	ax.bar(xvals+0.5,in_df.iloc[:,0],yerr=np.array(in_df.iloc[:,3:5].transpose()),width=1,
		edgecolor='black',color='#77DD77',ecolor='black',
		linewidth=1,capsize=5,align='center',
		error_kw={'elinewidth': 1,'ecapthick': 2}
		)
	ax.set_xticks(xvals+0.5)
	if notext == False:
		ax.set_xticklabels(in_df.index.values)
	else:
		ax.tick_params(labelbottom=False)
	ax.set_xlim(xmin=0)
	ax.set_xlabel(xaxis_label)
	ax.set_ylabel(yaxis_label)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.suptitle('Rsq_values')
	plt.tight_layout()
	if notext == True:
		plt.savefig(figname+'_full_expansion_rsq.pdf',transparent=True)
	else:
		plt.savefig(figname+'_full_expansion_rsq_labeled.pdf')
	
def combine_funcs(input_list,control,separate,rearrange,title,same_axis,notext,conf_int,
	ci_control,full):
	if full == False:
		if separate == True:
			if same_axis == False:
				binary_str,key_dict = parse_combinations(input_list,control,False)
				if conf_int == None:
					for k in key_dict:
						print key_dict[k]
						concat_df = concatenate_dfs_separate(binary_str
							,key_dict[k],rearrange,False)
						plot_rsq(concat_df,'Mutants','Rsquared',separate,
							'{0}_{1}_set'.format(title,k),notext,False)
				else:
					conf_int_dict,ci_key_dict = parse_combinations(conf_int,
						ci_control,True)
					try:
						key_dict.keys() == ci_key_dict.keys()
					except:
						print "Fix the ci dict matching the input dict"
						quit()
					binary_str = add_ci(binary_str,key_dict,conf_int_dict)
					for k in key_dict:
						print key_dict[k]
						concat_df = concatenate_dfs_separate(binary_str,
							key_dict[k],rearrange,True)
						plot_rsq(concat_df,'Mutants','Rsquared',separate,
							'{0}_{1}_set'.format(title,k),notext,True)
				
			else:
				binary_str,key_dict = parse_combinations(input_list,None,False)
				out_dict = collections.OrderedDict()
				if conf_int is not None:
					conf_int_dict,ci_key_dict = parse_combinations(conf_int,
						ci_control,True)
					binary_str = add_ci(binary_str,key_dict,conf_int_dict)
					try:
						key_dict.keys() == ci_key_dict.keys()
					except:
						print "Fix the ci dict matching the input dict"
						quit()
					print key_dict
					for k in key_dict:
						out_dict[k] = concatenate_dfs_separate(binary_str,
							key_dict[k],rearrange,True)
				else:
					for k in key_dict:
						out_dict[k] = concatenate_dfs_separate(binary_str,
							key_dict[k],rearrange,False)
				if control is not None:
					control_df = pd.DataFrame(read_tsv(control).iloc[0,:])
					if ci_control is not None:
						cic_df = pd.DataFrame(read_tsv(ci_control))
						control_df.loc['ci_low']=control_df.iat[0,0]-cic_df.iat[0,2]
						control_df.loc['ci_high']=cic_df.iat[1,2]-control_df.iat[0,0]
					else:
						pass
					out_dict[len(key_dict.keys())]=control_df
					out_dict[len(key_dict.keys())].transpose()
					out_dict[len(key_dict.keys())].columns=['control']
					if ci_control is not None:
						out_dict[len(key_dict.keys())].index = ['avg','stdev','num',
							'ci_low','ci_high']
					else:
						out_dict[len(key_dict.keys())].index = ['avg','stdev','num']
					print out_dict[len(key_dict.keys())]
				else:
					pass
	#			print out_dict
				if conf_int is not None:
					plot_same_axis(out_dict,'Mutants','Rsquared',title,notext,True)
				else:
					plot_same_axis(out_dict,'Mutants','Rsquared',title,notext,False)
		else:
			binary_str = []
			for file in input_list:
				file_list = file.split('_')
				filename = '{0}-{1}-{2}-{3}'.format(file_list[-6],
				file_list[-5],file_list[-4],file_list[-3])
				binary_str.append(filename)
			input_dfs = [read_tsv(x) for x in parsed.input]
			if parsed.control is not None:
				control_df = read_tsv(parsed.control)
				binary_str.append('Control')
				input_dfs.append(control_df)
			else:
				pass
			out_df = concatenate_dfs(input_dfs,binary_str)
			plot_rsq(out_df,'Mutants','Rsquared',separate,title,notext)
	else:
		if conf_int is not None and len(input_list)==1:
			pass
		else:
			print "See --full description for details on use"
			quit()
		in_df = read_full_rsq(input_list,conf_int)
		print in_df
		plot_full_rsq(in_df,'Order','Rsquared',title,notext)
		
		
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description='Script to make bar graph of selected Rsquared values')
	parser.add_argument('--input','-i',nargs='*',required=True,
		help='Any number of rsq_avg plots from bahadur_analysis.py.  Must have binary code in name of file')
	parser.add_argument('--control','-c',required=False,default=None,
		help='Control pdb to show Rsq values of a no-epistasis case.')
	parser.add_argument('--separate','-s',required=False,default=False,
		action='store_true',
		help='Separate subnetworks by number of mutations in starting node')
	parser.add_argument('--rearrange','-r',required=False,default=False,
		action='store_true',
		help='Sort bar-graph by value low to high. The control will always be last.')
	parser.add_argument('--title','-t',required=False,default=None,
		help='Title for bar graph output file')
	parser.add_argument('--oneaxis','-o',required=False,default=False,action='store_true',
		help='Call this flag to force a single axis')
	parser.add_argument('--notext','-n',required=False,default=False,action='store_true',
		help='Call this flag to turn off x-axis label')
	parser.add_argument('--conf_int','-ci',required=False,nargs='*',default=None,
		help='Rsq confidence interval outfiles from bootstrap_se.py')
	parser.add_argument('--conf_int_control','-cic',required=False,default=None,
		help='Control csv to parse confidence intervals')
	parser.add_argument('--full','-f',required=False,default=False,action='store_true',
		help='Call this flag if there is a full bahadur expansion rsq to plot. Requires that a single input is passed and multiple ci files are passed')
	parsed = parser.parse_args()
	binary_str = []
	combine_funcs(parsed.input,parsed.control,parsed.separate,
		parsed.rearrange,parsed.title,parsed.oneaxis,parsed.notext,parsed.conf_int,
		parsed.conf_int_control,parsed.full)
	