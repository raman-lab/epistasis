#Nishikawa, Kyle K.
#Edited 20200826 for epistasis manuscript
"""
Python script to be used in conjunction with flowjo_parse2.py to normalize for fluorescence differences arising from ligand concentration in GFP controls.

This script requires the *_fold_ind.tsv file from flowjo_parse2.py for a GFP set and the *_baseline_avg.tsv file from the experimental set.
"""

import numpy as np
import pandas as pd
import argparse
from collections import OrderedDict

def read_file(input_file):
	return pd.read_csv(input_file,sep='\t',index_col=0)

def error_calc(avg1,avg2,stdev1,stdev2):
	error = np.float(np.sqrt(((avg1/avg2)**2)*(((stdev1/avg1)**2)+((stdev2/avg2)**2))))
	return error

def parse_avg_std(input_df):
	stdev_dict = OrderedDict()
	avg_dict = OrderedDict()
	rownames = input_df.index.values
	
	for row in range(0,input_df.shape[0]):
		avgs = input_df.iloc[row,:].tolist()[0::2]
		stdevs = input_df.iloc[row,:].tolist()[1::2]
		avg_dict[rownames[row]] = avgs
		stdev_dict[rownames[row]] = stdevs
	avg_df = pd.DataFrame.from_dict(avg_dict,orient='index')
	stdev_df = pd.DataFrame.from_dict(stdev_dict,orient='index')
	return avg_df,stdev_df
	
#need a function to add a 1 to the start of the row
def add_one(data_frame):
	temp_list =[1,0]
	for num in data_frame.iloc[0,:].tolist():
		temp_list.append(num)
	temp_dict = {data_frame.index.values[0]: temp_list}
	return pd.DataFrame.from_dict(data=temp_dict,orient='index')

def calc_fi(base_avg,base_std,ctrl_avg,ctrl_std):
	averages = OrderedDict()
	stdevs = OrderedDict()
	rownames = base_avg.index.values
	for row in range(0,base_avg.shape[0]):
		stdev_list = []
		average_row = base_avg.iloc[row,:]/ctrl_avg.iloc[0,:]
		for avg in range(0,len(average_row)):
			stdev_list.append(error_calc(base_avg.iat[row,avg],ctrl_avg.iat[0,avg],
				base_std.iat[row,avg],ctrl_std.iat[0,avg]))
		averages[rownames[row]]=average_row
		stdevs[rownames[row]]=stdev_list
	avg_df = pd.DataFrame.from_dict(averages,orient='index')
	stdev_df = pd.DataFrame.from_dict(stdevs,orient='index')
	return avg_df,stdev_df

def stitch_dfs(df1,df2):
	rownames = df1.index.values
	new_dict = OrderedDict()
	colnames = []
	for row in range(0,len(rownames)):
		temp_list = []
		for col in range(0,df1.shape[1]):
			temp_list.append(df1.iat[row,col])
			temp_list.append(df2.iat[row,col])
		new_dict[rownames[row]] = temp_list
		temp_list = []
	for col in range(0,df1.shape[1]):
		colnames.append('avg_{0}'.format(col))
		colnames.append('stdev_{0}'.format(col))
	out_df = pd.DataFrame.from_dict(data=new_dict,orient='index')
	out_df.columns = colnames
	return out_df
	
			
	
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to scale baseline_avg values by fold induction of control')
	required = parser.add_argument_group('Required')
	required.add_argument('--baseline','-b',required=True,help='Baseline_avg file')
	required.add_argument('--fold_ind','-f',required=True,help='fold_ind file')
	parser.add_argument('--no_zero','-n',required=False,default=False,action='store_true',
		help='Call this flag to NOT ADD a zero to the first column of the normalization. this is required for fold induction normalizations but not baseline normalization')
	parsed = parser.parse_args()
	
	basename = parsed.baseline.split('/')[-1].split('.')[0]
	
	baseline_df = read_file(parsed.baseline)
	print baseline_df
	fold_ind_df = read_file(parsed.fold_ind)
	print fold_ind_df
	if parsed.no_zero == True:
		pass
	else:
		fold_ind_df = add_one(fold_ind_df)
	print fold_ind_df
	baseline_avg,baseline_stdev = parse_avg_std(baseline_df)
	fold_ind_avg,fold_ind_stdev = parse_avg_std(fold_ind_df)
	print 'baseline_avg'
	print baseline_avg
	print 'baseline_stdev'
	print baseline_stdev
	print fold_ind_avg
	print fold_ind_stdev
	corrected_baseline_avg,corrected_baseline_stdev = calc_fi(baseline_avg,baseline_stdev,
		fold_ind_avg,fold_ind_stdev)
	print corrected_baseline_avg
	print corrected_baseline_stdev
	final_out = stitch_dfs(corrected_baseline_avg,corrected_baseline_stdev)
	final_out.to_csv(basename+'_corrected.tsv',sep='\t')
	
	
	