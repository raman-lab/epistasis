#Nishikawa, Kyle K.
#For epistasis manuscript
#This script performs the Bahadur expansion on a data set (provided as averages and standard deviations).

import argparse
import numpy as np
import pandas as pd
import itertools
from collections import OrderedDict
import collections
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def grab_file(input_file,binary):
	if binary == False:
		if input_file.endswith('.tsv'):
			return pd.read_csv(input_file,delimiter='\t')
		elif input_file.endswith('.csv'):
			return pd.read_csv(input_file)
		elif input_file.endswith('.xlsx'):
			return pd.read_excel(input_file)
		else:
			print "Please use a .xlsx, .csv, or .tsv file"
			quit()
	else:
		if input_file.endswith('.tsv'):
			temp = pd.read_csv(input_file,delimiter='\t',header=0,dtype='str')
		elif input_file.endswith('.csv'):
			temp = pd.read_csv(input_file,header=0,dtype='str')
		elif input_file.endswith('.xlsx'):
			temp = pd.read_excel(input_file,header=0,dtype='str')
		print temp
		for col in temp.columns:
			print col
			if "binary" != col:
				temp[col] = pd.to_numeric(temp[col])
			else:
				pass
		return temp

#Generates binary 0/1 combinations based on an interger length
def generate_binary(length):
	return list(itertools.product([0,1],repeat=length))

#For generating subnetworks based on binary list
def four_node_sets(binary_list,length,verbose):
	by_count = collections.OrderedDict()
	accepted_list = []
	for binary in binary_list:
		if binary.count(1) <= length-2:
			accepted_list.append(binary)
		else:
			pass
	print "Starting node_ids:\n",accepted_list
	for num in range(0,length+1):
		for binary in binary_list:
			if num == binary.count(1):
				if num in by_count.keys():
					by_count[num].append(binary)
				else:
					by_count[num] = [binary]
			else:
				pass
	print "Binaries separated by number of mutations:\n",by_count
	out_dict = collections.OrderedDict()
	for num in range(0,length-1):
		single_change_list = list(itertools.combinations(by_count[num+1],2))
		double_change_list = list(itertools.combinations(by_count[num+2],1))
		for binary in by_count[num]:
			one_index, = np.where(np.array(binary)==1)
			for single_pairs in single_change_list:
				single_set = []
				if set(np.where(np.asarray(single_pairs[0])==1)[0]) > set(one_index) and set(np.where(np.asarray(single_pairs[1])==1)[0]) > set(one_index):
					single_set.append(single_pairs)
					one_counter, = np.where(np.asarray(single_pairs[0])==1)
					two_counter, = np.where(np.asarray(single_pairs[1])==1)
					union = collections.Counter(one_counter) | collections.Counter(two_counter)
					union_list = list(union.elements())
					for bin in double_change_list:
						bin_ones, = np.where(np.asarray(bin[0])==1)
						if set(bin_ones) == set(union_list):
							single_set.append(bin[0])
						else:
							pass
					if binary in out_dict.keys():
						out_dict[binary].append(single_set)
					else:
						out_dict[binary] = [single_set]
				else:
					pass
	if verbose == True:
		for k in out_dict:
			print '------'
			print k
			for v in out_dict[k]:
				print v
	else:
		pass
	return out_dict

#returns a list of dataframes with binary / data / std for every group in the four_node_dict
def associate_data(four_node_dict,input_df,data_col,stdev_col):
	binary_col = input_df.columns.values.tolist().index('binary')
	out_dict = OrderedDict()
	for k in four_node_dict:
		binary_str = ''.join(str(x) for x in k)
		binary_loc, = np.where(input_df.iloc[:,binary_col]==binary_str)
		for binary_set in four_node_dict[k]:
			temp_dict = OrderedDict()
			binary_list = [binary_str]
			data_list = [input_df.iat[binary_loc[0],data_col]]
			std_list = [input_df.iat[binary_loc[0],stdev_col]]
			mut_list = [''.join(str(y) for y in x) for x in binary_set[0]]
			mut_list.append(''.join(str(y) for y in binary_set[1]))
			for mut in mut_list:
				mut_ndx, = np.where(input_df.iloc[:,binary_col]==mut)
				binary_list.append(mut)
				data_list.append(input_df.iat[mut_ndx[0],data_col])
				std_list.append(input_df.iat[mut_ndx[0],stdev_col])
			temp_dict['binary']=binary_list
			temp_dict['averages'] = data_list
			temp_dict['stdevs'] = std_list
			temp_df = pd.DataFrame.from_dict(temp_dict)
			if k not in out_dict.keys():
				out_dict[k] = [temp_df]
			else:
				out_dict[k].append(temp_df)
	return out_dict

#this function matches a larger string to a 2-mutation subnetwork via constant values in the binary strings
def determine_equivalency(binary_array):
	binary_list = [list(x) for x in binary_array]
	binary_array = np.array(binary_list)
	print binary_array
	true_vals = np.all(binary_array == binary_array[0,:],axis=0)
	actual_index = []
	for x in range(0,true_vals.shape[0]):
		if true_vals[x] == False:
			actual_index.append(x)
		else:
			pass
#	print actual_index
	out_binary = OrderedDict()
	for item in binary_list:
		bruh = ''.join(item[x] for x in actual_index)
		out_binary[''.join(x for x in item)] = bruh
	out_df = pd.DataFrame.from_dict(out_binary, orient='index')
	return out_df
		

#Converts quadruple mutant data to binary sequence. Keys are based on synthesis approach and can be altered.
def map_to_quad(input_df):
	b_dict = {
		0: '1010',
		1: '1001',
		2: '1011',
		3: '1000',
		4: '0110',
		5: '0101',
		6: '0111',
		7: '0100',
		8: '1110',
		9: '1101',
		10: '1111',
		11: '1100',
		12: '0010',
		13: '0001',
		14: '0011',
		15: '0000',
		}
	temp_df = pd.DataFrame.from_dict(b_dict,orient='index',dtype='str')
	input_df['binary'] = temp_df.iloc[:,0]
	return input_df

#Converts binary to z values	
def binary_to_z(binary_list):
	z_dict = OrderedDict()
	for binary_tuple in binary_list:
		binary_str = ''.join(str(x) for x in binary_tuple)
		temp_list = []
		for num in binary_tuple:
			if num == 0:
				temp_list.append(-1)
			else:
				temp_list.append(1)
		z_dict[binary_str]=temp_list		
	return z_dict

#gives dictionary of psi values for every binary possible. A single binary term gives 2^n psi terms
def orthonormal_basis(z_dict,length):
	final_dict = OrderedDict()
	name_dict = OrderedDict()
	interaction_groups = [x for x in range(1,length+1)]
	for binary in z_dict:
		temp_list = []
		interaction_list = []
		for L in range(0,length+1):
			if L == 0:
				temp_list.append(1)
				interaction_list.append(0)
			else:
				test = list(itertools.combinations(z_dict[binary],L))
				temp_int = list(itertools.combinations(interaction_groups,L))
				if L == 1:
					#for some reason, placing this on a single line gives a generator obj
					for ndx in range(0,len(test)):
						temp_list.append(test[ndx][0])
						interaction_list.append(temp_int[ndx][0])
				else:
					for ndx in range(0,len(test)):
						temp_list.append(np.prod(test[ndx]))
						int_string = ''.join(str(x) for x in temp_int[ndx])
						interaction_list.append(int_string)
		final_dict[binary] = temp_list
		name_dict[binary] = interaction_list
	
	return final_dict,name_dict

#Verifies that the psi dataframe is orthonormal
def check_orthonormality(orthonormal_df,length):
	doubles = list(itertools.combinations(range(0,2**length),2))
	col_failures = []
	row_failures = []
	for double in doubles:
		col_out = np.dot(orthonormal_df.iloc[:,double[0]],
			orthonormal_df.iloc[:,double[1]])
		if col_out != 0:
			col_failures.append(double)
		else:
			pass
		row_out = np.dot(orthonormal_df.iloc[double[0],:],
			orthonormal_df.iloc[double[1],:])
		if row_out != 0:
			row_failures.append(double)
		else:
			pass
	if len(col_failures) == len(row_failures) == 0:
		return None,None
	else:
		return col_failures,row_failures

#Calculates Bahadur coefficients
def calc_w(psi_df,fluor_df,binary_row,fluor_row,verbose):
	final_w_vals_dict = OrderedDict()
	for row in range(0,psi_df.shape[0]):
		w_list = []
		for col in range(0,psi_df.shape[1]):
			colname = psi_df.columns.values[col]
			col_list = fluor_df.iloc[:,binary_row].tolist()
			if colname in col_list:
				col_ndx = col_list.index(colname)
				temp_w = np.float(fluor_df.iat[col_ndx,fluor_row] * psi_df.iat[row,col])
				w_list.append(temp_w)
			else:
				pass
		final_w_val = np.float((np.sum(w_list))/(2**length))
		final_w_vals_dict['w_{0}'.format(psi_df.index.values[row])]=final_w_val
		if verbose == True:
			print w_list
			sys.stdout.write('w_{0}\t{1}\n'.format(psi_df.index.values[row],final_w_val))
	final_wvals_df = pd.DataFrame.from_dict(data=final_w_vals_dict,orient='index')
	return final_wvals_df

#Calculates R^2 based on Bahadur coefficients
def calc_rsq(psi_df,w_vals_df,fluor_df,length,psi_groups_df,binary_row,fluor_row,phrase,
	need_graphs,verbose):
	terms_dict = OrderedDict()
	terms_list = psi_groups_df.iloc[:,0].tolist()
	for x in range(0,length+1):
		temp_list = []
		if x == 0:
			terms_dict[0] = [0]
		elif x != 0:
			for item in range(1,len(terms_list)):
				if x == len(str(terms_list[item])):
					temp_list.append(item)
				terms_dict[x] = temp_list
		else:
			pass
	psi_cols = psi_df.columns.values
	w_cols = w_vals_df.columns.values
	out_dict = OrderedDict()
	fluor_names = fluor_df.iloc[:,binary_row].tolist()
	rsq_dict = OrderedDict()
	for k in terms_dict:
		if verbose == True:
			print '------'
			print '{0}_order'.format(k)
		fluor_calc_vals = OrderedDict()
		end_point = terms_dict[k][-1]
		fluor_vals = fluor_df.iloc[:,fluor_row]
		w_vals = w_vals_df.iloc[:,0][:end_point+1]
		reordered_input = OrderedDict()
		#THis for loop is made under the assumption that the psi and w-vals are in the same order for different interaction terms (which they should be)
		for x in range(0,psi_df.shape[1]):
			binary = psi_cols[x]
			fluor_ndx = fluor_names.index(binary)
			psi_vals = psi_df.iloc[:,x][:end_point+1]
			#Uses the psi vals for a particular binary string and the wvals to calculate
			fluor_calc_vals[binary]=np.float(np.dot(psi_vals,w_vals))
			reordered_input[binary]=fluor_df.iat[fluor_ndx,fluor_row]
		out_dict['{0}_order'.format(str(k))] = fluor_calc_vals
		temp_df = pd.DataFrame.from_dict(data=fluor_calc_vals,orient='index')
		reordered_df = pd.DataFrame.from_dict(data=reordered_input,orient='index')
		if verbose == True:
			print "Reordered binary input:\n",reordered_df
			print "Calculated fluorescence:\n",temp_df
		p_coef = np.corrcoef(reordered_df.iloc[:,0],temp_df.iloc[:,0])[1,0]
		r_sq = p_coef**2
		rsq_dict[k]=r_sq
		if need_graphs == True:
			fig,ax = plt.subplots(1,1)
			ax.scatter(reordered_df.iloc[:,0],temp_df.iloc[:,0],
				c='#545454',edgecolor='black',
				label='Rsq_value: {0}'.format(r_sq))
			lims = [
				np.min([ax.get_xlim(),ax.get_ylim()]),
				np.max([ax.get_xlim(), ax.get_ylim()])
				]
			ax.plot(lims,lims,'k-',alpha=0.5,zorder=0)
			ax.set_ylabel('predicted_{0}'.format(phrase))
			ax.set_xlabel('actual_{0}'.format(phrase))
			ax.legend()
			plt.savefig(basename+'_'+str(phrase)+'_{0}_order_rsq.png'.format(k))
			plt.clf()
	return rsq_dict,out_dict

def graph_dfs(df,title,xaxis_label,yaxis_label,end,psi_groups_df,color,stdevs):
	rownames = df.index.values
	ind = np.arange(df.shape[0])
	breaks = range(2,20,2)
	if color == None:
		barcolor='#d3d3d3'
	else:
		barcolor = color
	plot_y = df.iloc[:,0]
	if np.amax(plot_y) > 1000*np.amin(plot_y) and np.amax(plot_y) > 4:
		fig,(ax1,ax) = plt.subplots(2,1,sharex=True)
		if stdevs is None:
			graph1 = ax1.bar(ind+0.5,df.iloc[:,0],width=1,linewidth=1,align='center',
				color=barcolor)
			graph2 = ax.bar(ind+0.5,df.iloc[:,0],width=1,linewidth=1,align='center',
				color=barcolor)
		else:
			graph1 = ax1.bar(ind+0.5,df.iloc[:,0],yerr=stdevs.iloc[:,0],width=1,
				linewidth=1,align='center',color=barcolor,
				ecolor='black',capsize=10,
				error_kw={'elinewidth': 1,'ecapthick': 2}
				)
			graph2 = ax.bar(ind+0.5,df.iloc[:,0],yerr=stdevs.iloc[:,0],width=1,
				linewidth=1,align='center',color=barcolor,
				ecolor='black',capsize=10,
				error_kw={'elinewidth': 1,'ecapthick': 2}
				)
		for item in range(0,len(breaks)):
			if np.amax(plot_y) < breaks[item]:
				ymin = (np.amax(plot_y))*.9
				ax1.set_ylim(ymin,breaks[item])
				ax1.set_xlim(0,ind.size)
				break
			else:
				pass
		for item in range(0,len(breaks)):
			if np.amin(plot_y) < breaks[item]:
				ax.set_ylim(ymax = breaks[item])
				ax.set_xlim(0,ind.size)
				break
			else:
				pass
		ax1.spines['bottom'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax1.xaxis.tick_top()
		ax1.tick_params(labeltop=False,length=0)
		ax.xaxis.tick_bottom()
		ax.set_xticks(ind+0.5)
		ax.tick_params(length=0)
	else:
		fig,ax = plt.subplots(1,1)
		fig.set_size_inches(8,6)
		if stdevs is None:
			ax.bar(ind+0.5,df.iloc[:,0],
				width=1,color=barcolor,edgecolor='black',
				linewidth=1,capsize=10,
				)
		else:
			ax.bar(ind+0.5,df.iloc[:,0],yerr=stdevs.iloc[:,0],width=1,
				linewidth=1,align='center',color=barcolor,
				ecolor='black',capsize=10,
				error_kw={'elinewidth': 1,'ecapthick': 2}
				)

	if psi_groups_df is not None:
		groups_list = [list(str(x)) for x in psi_groups_df.iloc[:,0].tolist()]
		name_strings = ['-'.join(x) for x in groups_list]
		ax.set_xticklabels(name_strings,rotation=45,ha='right')
	else:
		ax.set_xticklabels(rownames)
	ax.set_xlim(xmin=0)
	ax.set_xticks(ind+0.5)
	ax.set_xlabel(xaxis_label)
	ax.set_ylabel(yaxis_label)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.suptitle('{0}_{1}_{2}_{3}'.format(title,xaxis_label,yaxis_label,end),fontsize=8)
	plt.tight_layout()
	plt.savefig('{0}_{1}_{2}_{3}.pdf'.format(title,xaxis_label,yaxis_label,end))
	plt.close(fig)




'''
This function is called bootstrap but it's actually a Monte Carlo sampling of a normal distribution around the mean of the data using the standard deviation column.
'''
def bootstrap_data(input_df,bootstrap,data_col,stdev_col,ln,psi_df,psi_groups_df,
	length,mut7b,barcolor,binary_row,label,verbose,take_log,iterate,threshold):
	rsq_dict_fin = OrderedDict()
	wvals_df_fin = pd.DataFrame()
	wvals_df_avgs = pd.DataFrame()
	wvals_df_stdevs = pd.DataFrame()
	if iterate==True:
		test_range = 0
		diff_df_avgs = pd.DataFrame()
		diff_df_stdev = pd.DataFrame()
		while True:
			for x in range(test_range,test_range+bootstrap):
				temp_sam = np.random.normal(input_df.iloc[:,data_col],input_df.iloc[:,stdev_col])
				if ln == True:
					sampled = np.log(temp_sam)
					while True:
						sampled = np.log(temp_sam)
						if np.isnan(sampled).any() == True:
							temp_sam = np.random.normal(input_df.iloc[:,data_col],input_df.iloc[:,stdev_col])
						else:
							break
					act_sampled = np.log(input_df.iloc[:,data_col])
					untransformed = temp_sam
				elif take_log == True:
					while True:
						sampled = np.log10(temp_sam)
						if np.isnan(sampled).any() == True:
							temp_sam = np.random.normal(input_df.iloc[:,data_col],input_df.iloc[:,stdev_col])
						else:
							break
					act_sampled = np.log10(input_df.iloc[:,data_col])
					untransformed = temp_sam
				else:
					sampled = temp_sam
					act_sampled = input_df.iloc[:,data_col].copy()
					untransformed = temp_sam
				temp_df = input_df.copy()
				temp_df['rand_sam']=sampled
				temp_df['act_sam']=act_sampled
				temp_df['untransformed']=untransformed
				if mut7b == True:
					temp_df = map_to_quad(temp_df)
				else:
					pass
				if verbose == True:
					print temp_df
				else:
					pass
				if binary_row is None:
					binary_row = temp_df.columns.values.tolist().index('binary')
				else:
					pass
				data_row=temp_df.columns.values.tolist().index('rand_sam')
				act_row = temp_df.columns.values.tolist().index('act_sam')
				final_wvals_df = calc_w(psi_df,temp_df,binary_row,data_row,verbose)
				wvals_df_fin[x] = final_wvals_df.iloc[:,0]
				rsq_dict,predicted_fluors = calc_rsq(psi_df,final_wvals_df,
					temp_df,length,psi_groups_df,binary_row,data_row,'pass',False,verbose)
				for k in rsq_dict:
					if k not in rsq_dict_fin:
						rsq_dict_fin[k]=[rsq_dict[k]]
					else:
						rsq_dict_fin[k].append(rsq_dict[k])
			wvals_df_avgs[test_range] = wvals_df_fin.apply(np.nanmean,axis=1)
			wvals_df_stdevs[test_range] = wvals_df_fin.apply(np.nanstd,axis=1)
			if wvals_df_avgs.shape[1] > 1:
				diff_col_avgs = np.absolute(wvals_df_avgs[test_range-bootstrap]-wvals_df_avgs[test_range])
				diff_col_stdev = np.absolute(wvals_df_stdevs[test_range-bootstrap]-wvals_df_stdevs[test_range])
				if (np.absolute(diff_col_avgs / wvals_df_avgs[test_range])<threshold).all() == True and (np.absolute(diff_col_stdev / wvals_df_stdevs[test_range])<threshold).all() == True:
					print 'Threshold Met'
					print diff_col_avgs / wvals_df_avgs[test_range]
					print (np.absolute(diff_col_avgs / wvals_df_avgs[test_range])<threshold)
					print diff_col_stdev / wvals_df_stdevs[test_range]
					print (np.absolute(diff_col_stdev / wvals_df_stdevs[test_range])<threshold)
					diff_df_avgs[test_range]=diff_col_avgs / wvals_df_avgs[test_range]
					diff_df_stdev[test_range]=diff_col_stdev / wvals_df_stdevs[test_range]
					if label == None:
						diff_df_avgs.to_csv('bahadur_expansion3_iterate_diff_avgs.csv')
						diff_df_stdev.to_csv('bahadur_expansion3_iterate_diff_stdevs.csv')
					else:
						diff_df_avgs.to_csv(
							'bahadur_expansion3_iterate_{0}_diff_avgs.csv'.format(label))
						diff_df_stdev.to_csv(
							'bahadur_expansion3_iterate_{0}_diff_stdevs.csv'.format(label))
					break
				else:
					diff_df_avgs[test_range]=diff_col_avgs / wvals_df_avgs[test_range]
					diff_df_stdev[test_range]=diff_col_stdev / wvals_df_stdevs[test_range]
					if label == None:
						diff_df_avgs.to_csv('bahadur_expansion3_iterate_diff_avgs.csv')
						diff_df_stdev.to_csv('bahadur_expansion3_iterate_diff_stdevs.csv')
					else:
						diff_df_avgs.to_csv(
							'bahadur_expansion3_iterate_{0}_diff_avgs.csv'.format(label))
						diff_df_stdev.to_csv(
							'bahadur_expansion3_iterate_{0}_diff_stdevs.csv'.format(label))
					test_range = test_range+bootstrap
			else:
				test_range = test_range+bootstrap
	else:
		bootstrap = int(bootstrap)
		for x in range(0,bootstrap):
			temp_sam = np.random.normal(input_df.iloc[:,data_col],input_df.iloc[:,stdev_col])
			if ln == True:
				sampled = np.log(temp_sam)
				while True:
					sampled = np.log(temp_sam)
					if np.isnan(sampled).any() == True:
						temp_sam = np.random.normal(input_df.iloc[:,data_col],input_df.iloc[:,stdev_col])
					else:
						break
				act_sampled = np.log(input_df.iloc[:,data_col])
				untransformed = temp_sam
			elif take_log == True:
				while True:
					sampled = np.log10(temp_sam)
					if np.isnan(sampled).any() == True:
						temp_sam = np.random.normal(input_df.iloc[:,data_col],input_df.iloc[:,stdev_col])
					else:
						break
				act_sampled = np.log10(input_df.iloc[:,data_col])
				untransformed = temp_sam
			else:
				sampled = temp_sam
				act_sampled = input_df.iloc[:,data_col].copy()
				untransformed = temp_sam
			temp_df = input_df.copy()
			temp_df['rand_sam']=sampled
			temp_df['act_sam']=act_sampled
			temp_df['untransformed']=untransformed
			if mut7b == True:
				temp_df = map_to_quad(temp_df)
			else:
				pass
			print temp_df
			if binary_row is None:
				binary_row = temp_df.columns.values.tolist().index('binary')
			else:
				pass
			data_row=temp_df.columns.values.tolist().index('rand_sam')
			act_row = temp_df.columns.values.tolist().index('act_sam')
			final_wvals_df = calc_w(psi_df,temp_df,binary_row,data_row,verbose)
			wvals_df_fin[x] = final_wvals_df.iloc[:,0]
			rsq_dict,predicted_fluors = calc_rsq(psi_df,final_wvals_df,
				temp_df,length,psi_groups_df,binary_row,data_row,'pass',False,verbose)
			for k in rsq_dict:
				if k not in rsq_dict_fin:
					rsq_dict_fin[k]=[rsq_dict[k]]
				else:
					rsq_dict_fin[k].append(rsq_dict[k])
		wvals_df_avgs[bootstrap] = wvals_df_fin.apply(np.mean,axis=0)
		wvals_df_stdevs[bootstrap] = wvals_df_fin.apply(np.std,axis=0)
	rsq_df = pd.DataFrame.from_dict(rsq_dict_fin)
	rsq_avgs = OrderedDict()
	rsq_stdevs = OrderedDict()
	for col in range(1,rsq_df.shape[1]):
		rsq_avgs[col]=np.nanmean(rsq_df.iloc[:,col])
		rsq_stdevs[col] = np.nanstd(rsq_df.iloc[:,col])
	rsq_avgs_df = pd.DataFrame.from_dict(rsq_avgs,orient='index')
	rsq_stdevs_df = pd.DataFrame.from_dict(rsq_stdevs,orient='index')
	if label is None:
		label = input_df.columns.values.tolist()[data_col]
	else:
		pass
	graph_dfs(rsq_avgs_df,basename,'Orders','Rsquared','all_muts_{0}'.format(label),None,
		barcolor,rsq_stdevs_df)
	if ln == True:
		input_df['ln'] = np.log(input_df.iloc[:,data_col].copy())
	elif take_log == True:
		input_df['log10'] = np.log10(input_df.iloc[:,data_col].copy())
	else:
		pass
	if mut7b == True:
		map_to_quad(input_df)
	else:
		pass
#	input_df = map_to_quad(input_df)
	return rsq_avgs_df,rsq_stdevs_df,rsq_df,input_df,wvals_df_avgs,wvals_df_stdevs,wvals_df_fin


def subnetwork_bootstrap(input_df,bootstrap,data_col,stdev_col,ln,
	length,mut7b,barcolor,basename,write,verbose,take_log,iterate,threshold):
	print '----'
	temp_df = input_df.copy()
	if mut7b == True:
		temp_df = map_to_quad(temp_df)
	else:
		pass
	print temp_df
	rsq_dict_fin = OrderedDict()
	new_bin_list = generate_binary(2)
	sub_z_dict = binary_to_z(new_bin_list)
	sub_psi_dict,sub_psi_groups = orthonormal_basis(sub_z_dict,length)
	sub_psi_df = pd.DataFrame.from_dict(sub_psi_dict)
	sub_psi_groups_df = pd.DataFrame.from_dict(sub_psi_groups)
	four_node_dict = four_node_sets(binary_dict,int(parsed.length),verbose)

	mini_network = associate_data(four_node_dict,temp_df,
		data_col,stdev_col)
	for k in mini_network:
		for df in mini_network[k]:
			rsq_dict = OrderedDict()
			binary_strings = df.loc[:,'binary']
			print "binary strings:\n",binary_strings
			binary_label = '_'.join(x for x in binary_strings)
			substrings = determine_equivalency(binary_strings)
			df['substrings']=substrings.values
			print df
			subcol = df.columns.values.tolist().index('substrings')
			rsq_avgs_df,rsq_stdevs_df,rsq_df,mod_input_df,wvals_df_avgs,wvals_df_stdevs,wvals_df_fin = bootstrap_data(df,
					bootstrap,1,2,ln,sub_psi_df,sub_psi_groups_df,2,False,
					barcolor,subcol,binary_label,verbose,take_log,iterate,threshold
					)
			rsq_final_df = stitch_dfs(rsq_avgs_df,rsq_stdevs_df)
			rsq_final_df['num']=np.repeat(rsq_df.shape[0],rsq_avgs_df.shape[0])
			wvals_df = stitch_dfs(wvals_df_avgs,wvals_df_stdevs)
			if write == True:
				mod_input_df.to_csv('{0}_{1}_mod_input.tsv'.format(basename,binary_label),
					sep='\t')
				rsq_final_df.to_csv('{0}_{1}_rsq_avgs.tsv'.format(basename,binary_label),
					sep='\t')
				rsq_df.to_csv('{0}_{1}_rsq_raw.tsv'.format(basename,binary_label),
					sep='\t')
				wvals_df.to_csv('{0}_{1}_wvals_avgs.csv'.format(basename,binary_label))
				wvals_df_fin.to_csv(
					'{0}_{1}_wvals_raw.csv'.format(basename,binary_label))
			else:
				print mod_input_df
				print rsq_final_df
				print rsq_df
	if write == True:
		sub_psi_df.to_csv('{0}_subgroups_psi.tsv'.format(basename),sep='\t')
		sub_psi_groups_df.to_csv('{0}_subgroups_psi_groups.tsv'.format(basename),sep='\t')
	else:
		print sub_psi_df
		print sub_psi_groups_df
				
			
				
	
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
		
		

def coefficient_calc(num,fluorescence,psi_vals,length):
	return (1/(2**length))*np.sum(fluorescence*psi_vals)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='script for bahadur analysis')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,
		help='CSV, TSV, or excel file of fluorescences')
	parser.add_argument('--length','-l',required=False,default=4,
		help='Number of variants (length of binary string)')
	parser.add_argument('-7b',dest='mut7b',required=False,default=False,
		action='store_true',help='Call this flag if mutants are from 7B')
	parser.add_argument('-ln',required=False,default=False,action='store_true',
		help='Call this flag to apply a natural log to all fluorescence values')
	parser.add_argument('--write','-w',default=False,action='store_true',
		help='Call this flag to force output to different files')
	parser.add_argument('--binary','-b',default=False,action='store_true',
		help='Call this flag if the input file lists variants in a binary form')
	parser.add_argument('--use_col','-u',default=1,action='store',
		help='Specify a column of data to use for analysis, 0-indexed')
	parser.add_argument('--std','-s',default=4,action='store',
		help='Specify column of data standard deviation, 0-indexed')
	parser.add_argument('--bootstrap','-x',default=500,action='store',
		help='Number of times to bootstrap values')
	parser.add_argument('--subnetwork','-k',default=False,action='store_true',
		help='Call flag to force subnetwork analysis')
	parser.add_argument('--color','-y',dest='barcolor',required=False,default='#77FF77',
		help='Specify bar color, default = #77FF77')
	parser.add_argument('--verbose','-v',required=False,default=False,
		action='store_true',
		help='Call this flag for verbose output to terminal')
	parser.add_argument('-log',default=False,action='store_true',required=False,
		help='Call this flag to take the log10 of all data')
	parser.add_argument('--iterate','-e',default=False,action='store_true',required=False,
		help='Call this flag to iterate in increments of --bootstrap until convergence (mean and standard deviation of w-values change by less than threshold percent')
	parser.add_argument('--threshold','-d',default=0.001,required=False,
		help='Call this flag to change threshold to which the iterate process goes to. Given as a percent difference (provide decimal) between current and previous iterations. Default is 0.001 (0.1 percent).')
	
	parsed = parser.parse_args()
	length = int(parsed.length)
	basename = parsed.input.split('/')[-1].split('.')[0]
	input_df = grab_file(parsed.input,parsed.binary)
	print input_df
	binary_dict = generate_binary(length)
	z_dict = binary_to_z(binary_dict)
	psi_dict,psi_groups = orthonormal_basis(z_dict,length)
	psi_df = pd.DataFrame.from_dict(data=psi_dict)
	psi_groups_df = pd.DataFrame.from_dict(data=psi_groups)
	col_fails,row_fails = check_orthonormality(psi_df,length)
	np.set_printoptions(precision=5)
	if col_fails is None and row_fails is None:
		pass
	else:
		binary_list = psi_df.columns.values
		print binary_list
		print col_fails
		quit()
	print "Data col:\n",input_df.iloc[:,int(parsed.use_col)]
	print "StDev col:\n",input_df.iloc[:,int(parsed.std)]

	if parsed.subnetwork == True:
		subnetwork_bootstrap(input_df,int(parsed.bootstrap),int(parsed.use_col),
			int(parsed.std),parsed.ln,int(parsed.length),parsed.mut7b,parsed.barcolor,
			basename,parsed.write,parsed.verbose,parsed.log,parsed.iterate,
			np.float64(parsed.threshold))

	else:
		rsq_avgs,rsq_stdevs,rsq_raw_df,out_df,wvals_df_avgs,wvals_df_stdevs,wvals_df_fin = bootstrap_data(input_df,
			int(parsed.bootstrap),
			int(parsed.use_col),int(parsed.std),parsed.ln,psi_df,psi_groups_df,
			int(parsed.length),parsed.mut7b,parsed.barcolor,None,None,parsed.verbose,
			parsed.log,parsed.iterate,np.float64(parsed.threshold)
			)
		rsq_df = stitch_dfs(rsq_avgs,rsq_stdevs)
		rsq_df['num'] = np.repeat(rsq_raw_df.shape[0],rsq_df.shape[0])
		print wvals_df_avgs
		print wvals_df_stdevs
		wvals_df = stitch_dfs(wvals_df_avgs,wvals_df_stdevs)
		if parsed.write == True:
			out_df.to_csv(basename+'_input.tsv',sep='\t')
			psi_df.to_csv(basename+'_psi.tsv',sep='\t')
			psi_groups_df.to_csv(basename+'_psi_groups.tsv',sep='\t')
			rsq_df.to_csv(basename+'_rsq_avgs.tsv',sep='\t')
			rsq_raw_df.to_csv(basename+'_rsq_raw.tsv',sep='\t')
			wvals_df.to_csv(basename+'_wvals_avgs.csv')
			wvals_df_fin.to_csv(basename+'_wvals_raw.csv')
		else:
			print input_df
			print psi_df
			print psi_groups_df
			print rsq_df 
			print wvals_df



	
		
		
		