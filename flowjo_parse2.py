#Nishikawa,Kyle K.
#Edit 20200825 for epistasis manuscript

#Script is designed to work on data tables written by FlowJo consisting of median RFU values for wells in a 96-well format.

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import collections
import re
from collections import OrderedDict

#--------argparse checks--------
def xcel_file(file):
	if file.endswith('.xls') or file.endswith('.xlsx'):
		return file
	else:
		msg = "Please use a valid Excel file"
		raise argparse.ArgumentTypeError(msg)

def is_int(interger):
	try:
		lel = int(interger)
		return lel
	except:
		msg = "Please use an interger"
		raise argparse.ArgumentTypeError(msg)
def is_well(string):
	well_coords = re.split('(\d+)',string)
	if well_coords[0].upper() in list('ABCDEFGH'):
#		print "yes"
		if int(well_coords[1]) in list(range(1,13)):
			return string
		else:
			msg = "Please use a valid well ID"
			raise argparse.ArgumentTypeError(msg)
	else:
		msg = "Please use a valid well ID"
		raise argparse.ArgumentTypeError(msg)




#Function to read the data in pandas
def grab_df(input_file):
	fluorescence_df = pd.read_excel(input_file,header=None,usecols=[0,1])
	fluorescence_df.columns=['sample','median']
	return fluorescence_df

def trunchate_df(dataframe,take_log):
	new_dict = {}
	dataframe = dataframe.drop([0,dataframe.shape[0]-1,dataframe.shape[0]-2])
	if parsed.log == True:
		dataframe['median']=np.log10(np.float64(dataframe['median']))
	sample = dataframe.iloc[:,0]
	labels = []
	for item in sample:
#		print item
		if "Specimen" in item:
			name = item.split('_')[2]
		elif "_" in item:
			name = item.split('_')[0]
		else:
			name = item.split('.')[0]
		labels.append(name)
	dataframe.index=labels
	sample_list = dataframe.index.values
	for sam in sample_list:
		sam_index, = np.where(sample_list==sam)
		sam.encode('ascii','ignore')
		sample_as_string = dataframe.iloc[sam_index,0].to_string()
		if sample_as_string.find(sam)==-1:
			print 'do not match: ', sam, sample_as_string
			quit()
		else:
			pass
	return dataframe

def total_cols_rows(total,vert,xbyy_list,start_well_coords):
	if xbyy_list != None:
		xbyy_rows = int(xbyy_list[0])
		xbyy_cols = int(xbyy_list[1])
	total_samples = num_var*num_sam*num_rep
	if vert == False:
		if xbyy_list == None:
			rows_filled = (total_samples+start_well_coords[1]) / 12
			leftovers = total_samples % 12
			rows_filled += start_well_coords[0]
			if leftovers != 0:
				rows_filled += 1
			else:
				pass
			if (12-start_well_coords[1])<leftovers:
				rows_filled +=1
			else:
				pass
			cols_filled = None
		else:
			rows_filled = xbyy_rows + start_well_coords[0]
			cols_filled = xbyy_cols + start_well_coords[1]
	else:
		if xbyy_list == None:
			cols_filled = (total_samples+start_well_coords[0])/8
			leftovers = total_samples % 8
			cols_filled += start_well_coords[1]
			if leftovers != 0:
				cols_filled += 1
			else:
				pass
			if (8-start_well_coords[0])<leftovers:
				cols_filled += 1 
			else:
				pass
			rows_filled = None
		else:
			cols_filled = xbyy_cols+start_well_coords[1]
			rows_filled = xbyy_rows+start_well_coords[0]
	return rows_filled,cols_filled

def get_samples(row,col,dataframe,samples):
	well_row = rows_plate[row]
	well_col = cols_plate[col]
	well_name = str(well_row)+str(well_col)
	label_index, = np.where(samples == well_name)
	if len(label_index) != 0:
		return float(dataframe.iloc[label_index,1]),label_index
	else:
		print "help"
		return None

def attune_grab_baselines(df):
	samples = df.index.values
	baseline_dict = collections.OrderedDict()
	total_samples = num_sam * num_var * num_rep
	rep_values = []
	rows_filled,cols_filled = total_cols_rows(total_samples,
		vertical,xbyy_list,start_well_coords)
	if vertical == False:
		if xbyy == None:
			for row in range(start_well_coords[0],rows_filled):
				if row != start_well_coords[0]:
					#grabs nonstarting rows
					for col in range(0,12):
						base_val,label_index = get_samples(row,col,df,samples)
						rep_values.append(base_val)
						if len(rep_values)==num_sam:
							baseline_dict[samples[int(label_index)-(num_sam-1)]]=rep_values
							rep_values = []
						else:
							pass
				else:
					for col in range(start_well_coords[1],12):
						base_val,label_index = get_samples(row,col,df,samples)
						rep_values.append(base_val)
						if len(rep_values)==num_sam:
							baseline_dict[samples[int(label_index)-(num_sam-1)]] = rep_values
							rep_values = []
						else:
							pass
		else:
			for row in range(start_well_coords[0],rows_filled):
				for col in range(start_well_coords[1],cols_filled):
					base_val,label_index = get_samples(row,col,df,samples)
					rep_values.append(base_val)
					if len(rep_values)==num_sam:
						baseline_dict[samples[int(label_index)-(num_sam-1)]] = rep_values
						rep_values = []
					else:
						pass
	else:
		if xbyy == None:
			for col in range(start_well_coords[1],cols_filled):
				if col != start_well_coords[1]:
					for row in range(0,8):
						base_val,label_index = get_samples(row,col,df,samples)
						rep_values.append(base_val)
						if len(rep_values)==num_sam:
							baseline_dict[samples[int(label_index)-(num_sam-1)]] = rep_values
							rep_values = []
						else:
							pass
				else:
					for row in range(start_well_coords[0],8):
						base_val,label_index = get_samples(row,col,df,samples)
						rep_values.append(base_val)
						if len(rep_values)==num_sam:
							baseline_dict[samples[int(label_index)-(num_sam-1)]] = rep_values
							rep_values = []
						else:
							pass
		else:
			for col in range(start_well_coords[1],cols_filled):
				for row in range(start_well_coords[0],rows_filled):
					base_val,label_index = get_samples(row,col,df,samples)
					rep_values.append(base_val)
					if len(rep_values)==num_sam:
						baseline_dict[samples[int(label_index)-(num_sam-1)]] = rep_values
						rep_values = []
					else:
						pass
	baseline_df = pd.DataFrame.from_dict(data=baseline_dict,orient='index')
	return baseline_df

#the baseline_df should be a data frame of variants and their respective fluorescence values for each sample containing that variant.

def avg_stdev(df,delcol_list):
	samples = df.index.values
	averages = OrderedDict()
	stdevs = OrderedDict()
	total_reps = [x for x in range(0,num_var*num_rep)]
	total_rep_sep = total_reps[::num_rep]
	#rep_sep defines the subgroup on the y-axis
	start = 0
	for rep_sep in total_rep_sep:
		temp_avg = []
		temp_stdev = []
		#num_sam iterates through the different averages for a particular group
		if num_rep == 1:
			for x in range(0,num_sam):
				test = df.iat[rep_sep,x]
				temp_avg.append(np.mean(test))
				temp_stdev.append(np.std(test))
			if len(temp_avg)==len(temp_stdev)==num_sam:
				averages[samples[rep_sep]]=temp_avg
				stdevs[samples[rep_sep]]=temp_stdev
				temp_avg = []
				temp_stdev = []
			else:
				pass
		else:
			for x in range(0,num_sam):
#				print x
				test = df.iloc[rep_sep:rep_sep+num_rep,x]
				print test
				temp_avg.append(np.mean(test))
				temp_stdev.append(np.std(test))
			if len(temp_avg)==len(temp_stdev)==num_sam:
				averages[samples[rep_sep]]=temp_avg
				stdevs[samples[rep_sep]]=temp_stdev
				temp_avg = []
				temp_stdev = []
			else:
				pass
	averages_df = pd.DataFrame.from_dict(data=averages,orient='index')
	stdevs_df = pd.DataFrame.from_dict(data=stdevs,orient='index')
	if delcol_list is not None:
		dropcol_list = [int(x) for x in delcol_list.split(',')]
		averages_df.drop(averages_df.columns[dropcol_list],axis=1,inplace=True)
		stdevs_df.drop(stdevs_df.columns[dropcol_list],axis=1,inplace=True)
	else:
		pass
	return averages_df,stdevs_df

def error_calc(avg1,avg2,stdev1,stdev2):
	error = np.float(np.sqrt(((avg1/avg2)**2)*(((stdev1/avg1)**2)+((stdev2/avg2)**2))))
	return error

def fold_ind_calcs(averages_df,stdevs_df):
	fold_ind_dict = OrderedDict()
	fi_stdev_dict = OrderedDict()
	samples = averages_df.index.values
	sample_indices = [x for x in range(0,num_sam)][::2]
#	print sample_indices
	for var in range(0,averages_df.shape[0]):
		if doser == False:
			temp_avgs = []
			temp_stds = []
			for base in sample_indices:
				temp_avgs.append(np.float(averages_df.iat[var,base+1]/
					averages_df.iat[var,base]))
				temp_stds.append(error_calc(averages_df.iat[var,base+1],
					averages_df.iat[var,base],stdevs_df.iat[var,base+1],
					stdevs_df.iat[var,base]))
			fold_ind_dict[samples[var]]=temp_avgs
			fi_stdev_dict[samples[var]]=temp_stds
		else:
			vals = averages_df.iloc[var,:].tolist()
			stds = stdevs_df.iloc[var,:].tolist()
			baseline = vals[0]
			baseline_std = stds[0]
			temp_fi = []
			temp_std = []
			for val in range(1,len(vals)):
				temp_fi.append(np.float(vals[val]/baseline))
				temp_std.append(error_calc(vals[val],baseline,stds[val],baseline_std))
			fold_ind_dict[samples[var]] = temp_fi
			fi_stdev_dict[samples[var]] = temp_std
	fold_ind_df = pd.DataFrame.from_dict(data=fold_ind_dict,orient='index')
	fi_stdev_df = pd.DataFrame.from_dict(data=fi_stdev_dict,orient='index')		
	return fold_ind_df,fi_stdev_df

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
	parser = argparse.ArgumentParser(description='Script flowjo data tables in the format of two columns: name \ fluorescence')
	required = parser.add_argument_group('Required')

	required.add_argument('--input','-i',type=xcel_file,action='store',required=True,
		dest='input_file',help='Input excel file')
	required.add_argument('--variants','-v',type=is_int,action='store',required=True,
		dest='num_var',help='The number of variants tested')
	required.add_argument('--replicates','-r',type=is_int,action='store',required=True,
		dest='num_rep',help='The number of replicates per variant')
	required.add_argument('--samples','-s',type=is_int,action='store',required=False,
		default=4,dest='num_sam',help='The number of samples per replicate')
	parser.add_argument('--vertical','-l',action='store_true',required=False,
		default=False,dest='orientation',
		help='Call this flag if the samples in a replicate were arranged vertically on the plate.  Default is false (samples oriented horizontally for a replicate)')
	parser.add_argument('--start','-p',type=is_well,action='store',required=False,
		default='A1',dest='start_well',help='Specify the well to start')
	parser.add_argument('--doseresponse','-d',action='store_true',required=False,
		default=False,dest='doseresponse',
		help='Call this flag if there are multiple induced samples and one control sample')
	parser.add_argument('--xbyy','-x',action='store',required=False,dest='xbyy',
		help='Call this flag if there is no offset and and the selected wells are in a grid of x rows and y columns.  Give in the form of rows,cols')
	parser.add_argument('--baselines','-b',action='store_true',required=False,
		default=False,dest='blines',
		help='Call this flag if you do not want any fold induction output or baseline average output')
	parser.add_argument('--suffix','-sx',action='store',required=False,default=None,
		dest='suffix',help='Define a suffix to be added to the end of all output files, default is None')
	parser.add_argument('--delcol','-c',action='store',required=False,default=None,
		help='Interger or comma-separated list for columns to remove from the final dataframes')
	parser.add_argument('-log',action='store_true',required=False,default=False,
		help='Call this flag to take the log of all fluorescence values')
	parsed = parser.parse_args()

	input_file = parsed.input_file
	num_var = parsed.num_var
	num_rep = parsed.num_rep
	num_sam = parsed.num_sam
	vertical = parsed.orientation
	start_well = parsed.start_well.upper()
	doser = parsed.doseresponse
	xbyy = parsed.xbyy
	blines = parsed.blines
	suffix = parsed.suffix

	filename = input_file.split('/')[-1].split('.')[0]
	rows_plate = list('ABCDEFGH')
	cols_plate = list(range(1,13))
	plate_dict = {}
	for item in rows_plate:
		plate_dict[item]=cols_plate
	start_well_list = re.split('(\d+)',start_well)
	start_well_coords = (rows_plate.index(start_well_list[0]),cols_plate.index(int(start_well_list[1])))
	print start_well_coords

	if xbyy:
		xbyy_list = xbyy.split(',')
		xbyy_rows = int(xbyy_list[0])
		xbyy_cols = int(xbyy_list[1])
		if 12 < start_well_coords[1] + xbyy_cols:
			print "Use a valid number of columns"
			quit()
		elif 8 < start_well_coords[0] + xbyy_rows:
			print "Use a valid number of rows"
			quit()
		else:
			pass
	else:
		xbyy_list = None

	input_df = grab_df(input_file)
	trunch_df = trunchate_df(input_df,parsed.log)
	baseline_df = attune_grab_baselines(trunch_df)
	averages_df,stdevs_df = avg_stdev(baseline_df,parsed.delcol)
	print averages_df
	fold_ind_df,fi_stdev_df = fold_ind_calcs(averages_df,stdevs_df)

	baseline_avg_df = stitch_dfs(averages_df,stdevs_df)

	fold_ind_avg_df = stitch_dfs(fold_ind_df,fi_stdev_df)

	if parsed.delcol == None:
		pass
	else:
		dropcol_list = [int(x) for x in parsed.delcol.split(',')]
		baseline_df.drop(baseline_df.columns[dropcol_list],axis=1,inplace=True)
	
	if suffix == None:
		if parsed.log == False:
			baseline_df.to_csv(filename+'_baseline.tsv',sep='\t')
			baseline_avg_df.to_csv(filename+'_baseline_avg.tsv',sep='\t')
			fold_ind_avg_df.to_csv(filename+'_fold_ind_avg.tsv',sep='\t')
		else:
			baseline_df.to_csv(filename+'_baseline_log.tsv',sep='\t')
			baseline_avg_df.to_csv(filename+'_baseline_avg_log.tsv',sep='\t')
			fold_ind_avg_df.to_csv(filename+'_fold_ind_avg_log.tsv',sep='\t')
	else:
		if parsed.log == False:
			baseline_df.to_csv(filename+'_baseline_'+suffix+'.tsv',sep='\t')
			baseline_avg_df.to_csv(filename+'_baseline_avg_'+suffix+'.tsv',sep='\t')
			fold_ind_avg_df.to_csv(filename+'_fold_ind_avg_'+suffix+'.tsv',sep='\t')
		else:
			baseline_df.to_csv(filename+'_baseline_log_'+suffix+'.tsv',sep='\t')
			baseline_avg_df.to_csv(filename+'_baseline_avg_log_'+suffix+'.tsv',sep='\t')
			fold_ind_avg_df.to_csv(filename+'_fold_ind_avg_log_'+suffix+'.tsv',sep='\t')




	
			
	

	