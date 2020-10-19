#Nishikawa, Kyle K.
#Edited 20200826 for Epistasis manuscript

#This script creates dose response curves for a set 
import scipy.optimize
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from collections import OrderedDict
import re
import traceback



def read_tsv(input_file):
	return pd.read_csv(input_file,index_col=0,delimiter='\t')

def func(conc,n,kd,offset):
	return (conc**n)/(kd+(conc**n)) + offset

def func2(conc,n,kd,fmax):
	return baseline + ((fmax - baseline) * ((conc**n)/((kd**n)+(conc**n))))

def func3(conc,n,kd):
	return baseline + ((fmax - baseline) * ((conc**n)/((kd**n)+(conc**n))))
	
def ratio_data(input_df):
	rownames = input_df.index.names
	ratio_dict = OrderedDict()
	for row in range(0,input_df.shape[0]):
		fold_ind_list_max = max(input_df.iloc[row,:].tolist())
		ratio_list = [x/fold_ind_list_max for x in input_df.iloc[row,:].tolist()]
		ratio_dict[rownames[row]] = ratio_list
	return pd.DataFrame.from_dict(data=ratio_dict,orient='index')

def trunchate_data(input_df,conc_list):
	temp_dict = OrderedDict()
	rownames = input_df.index.names
	for row in range(0,len(rownames)):
		temp_dict[rownames[row]]=input_df.iloc[row,len(conc_list)]
	return pd.DataFrame.from_dict(data=temp_dict,orient='index')

def rmscalc(yact,yfit):
	return np.sqrt(np.sum((yfit-yact)**2)/yfit.shape[0])
	
def fit_data(input_df,conc_list,yaxis_label,xunits,force_fmax):
	rownames = input_df.index.values
	minconc = float(min(conc_list))
	if '-' and '_' in str(rownames[0]):
		rownames = [re.sub('\_|\-','',x) for x in rownames]
	else:
		pass
	output_dict = OrderedDict()
	fit_xvals = np.arange(0,max(conc_list)*10)
	for row in range(0,input_df.shape[0]):
		print row
		fig,ax=plt.subplots(1,1)
		fig.set_size_inches(8,6)
		fold_ind = input_df.iloc[row,:].tolist()[0:][::2]
		fmax_approx = max(fold_ind)
		global baseline
		baseline = fold_ind[conc_list.index(minconc)]
		if force_fmax == True:
			global fmax
			fmax = max(fold_ind)
			print fmax
		else:
			pass
		stdevs = input_df.iloc[row,:].tolist()[1:][::2]
		baseline_stds = stdevs[conc_list.index(minconc)]
		conc = conc_list
		y_tups = [(conc_list[x],fold_ind[x],stdevs[x]) for x in range(0,len(conc_list))]
		y_tups = sorted(y_tups,key=lambda x: x[0])
		y_zipped = zip(*y_tups)
		print y_zipped[1]
		print y_zipped[2]
		failed_fits = []
		try:
			if force_fmax == False:
				popt,pcov = scipy.optimize.curve_fit(func2,y_zipped[0],y_zipped[1],
					p0=[1,100,fmax_approx],
					method='lm',sigma=y_zipped[2],absolute_sigma=True)
			else:
				popt,pcov = scipy.optimize.curve_fit(func3,y_zipped[0],y_zipped[1],
					p0=[1,100],
					method='lm',sigma=y_zipped[2],absolute_sigma=True)
		except:
			failed_fits.append(rownames[row])
			traceback.print_exc()
			temp_list = [np.nan for x in range(0,7)]
			temp_list.append(baseline)
			temp_list.append(baseline_stds)
			output_dict[rownames[row]]=temp_list
		else:
			if force_fmax == False:
				temp_list = [x for x in popt]
				for num in np.sqrt(np.diag(pcov)):
					temp_list.append(num)
			else:
				temp_list = [x for x in popt]
				temp_list.append(fmax)
				for num in np.sqrt(np.diag(pcov)):
					temp_list.append(num)
				temp_list.append(stdevs[y_zipped[1].index(fmax)])
			temp_list.append(baseline)
			temp_list.append(baseline_stds)
			if force_fmax == False:
				rmserror = rmscalc(y_zipped[1],func2(y_zipped[0],*popt))
				ax.errorbar(fit_xvals,func2(fit_xvals,*popt),color='red',label='fit: n = %5.3f\n kd = %5.3f\n fmax=%5.3f' % tuple(popt))
			else:
				rmserror = rmscalc(y_zipped[1],func3(y_zipped[0],*popt))
				ax.errorbar(fit_xvals,func3(fit_xvals,*popt),color='red',label='fit: n = %5.3f\n kd = %5.3f\n' % tuple(popt))
			temp_list.append(rmserror)
			output_dict[rownames[row]]=temp_list
		ax.errorbar(y_zipped[0],y_zipped[1],yerr=y_zipped[2],fmt='o',
			color='blue',capsize=5, elinewidth=1, markeredgewidth=1,
			label='data\n {0}'.format(rmserror)
			)
		ax.set_xlabel('Concentration ('+xunits+')')
		ax.set_ylabel(yaxis_label)
		ax.set_xscale('log')
		ax.set_xlim(right=max(y_zipped[0])*5)
		if ttgr7b == True:
			rowname = rownames[row]
			plt.title(str(rownames[row]) + ' - ' + node_id_dict[rowname])
		else:
			plt.title(rownames[row])
		plt.legend(loc=2)
		if force_fmax == True:
			plt.savefig(base_name+'_'+str(rownames[row])+'_forcefmax_doseresponse_fit.pdf')
		else:
			plt.savefig(base_name+'_'+str(rownames[row])+'_'+'doseresponse_fit.pdf')
	print failed_fits
#	print output_dict
	output_df = pd.DataFrame.from_dict(data=output_dict,orient='index')
	colnames = ['n','kd','fmax','n_std','kd_std','fmax_std',
		'baseline','baseline_std','rms']
	output_df.columns = colnames
	return output_df,failed_fits

def round_df(input_df,round):
	test_df = input_df.astype('str')
	rownames = test_df.index.values
	new_dict = OrderedDict()
	colnames = input_df.columns.values
	for row in range(0,test_df.shape[0]):
		temp_list = []
		for item in test_df.iloc[row,:]:
			test = '{0}'.format(float("{0:.{1}g}".format(np.float64(item),round)))
			temp_list.append(test)
		new_dict[rownames[row]]=temp_list
	rounded_df = pd.DataFrame.from_dict(data=new_dict,orient='index',dtype=np.float64)
	rounded_df.columns = colnames
	return rounded_df

def graph_output(output_df,node_id_dict):
	barnames = [node_id_dict[x] for x in output_df.index.values]
	ind = np.arange(output_df.shape[0])
	for x in range(0,3):
		fig,ax = plt.subplots(1,1)
		fig.set_size_inches(10,7.5)
		coords = range(0,6)[x::3]
		induction_val = output_df.iloc[:,coords[0]]
		induction_val[induction_val>500000]=0
		error_vals = output_df.iloc[:,coords[1]]
		error_vals[error_vals>500000]=0
		ax.bar(ind+0.5,induction_val,
			width=1,color='#d3d3d3',edgecolor='black',yerr=output_df.iloc[:,coords[1]],
			linewidth=1,capsize=10,error_kw={'elinewidth': 1,'ecapthick': 2}
			)
		ax.set_xlim(xmin=0)
		ax.set_xlabel('Mutants')
		if np.amax(output_df.iloc[:,coords[0]])*5 < np.amax(output_df.iloc[:,coords[1]]):
			ax.set_ylim(0,np.amax(output_df.iloc[:,coords[0]])*5)
		else:
			pass
		ax.set_ylabel(output_df.columns.values[x])
		ax.set_xticks(ind+0.5)
		ax.set_xticklabels(barnames,rotation=60,ha='right',fontdict={'fontsize': 12})
		
		plt.title(output_df.columns.values[x])
		plt.tight_layout()
		plt.savefig(base_name+'_{0}.pdf'.format(output_df.columns.values[x]))
		
		plt.clf()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to fit dose response data to the hill equation')
	required = parser.add_argument_group('required')
	required.add_argument('--input','-i',dest='input',required=True,action='store',
		help='Input tsv file containing dose response data and standard deviations in alternating columns.')
	required.add_argument('--concentrations','-c',dest='concentrations',required=True,
		action='store',help='Comma-separated list of concentrations')
	parser.add_argument('--round','-r',default=None,action='store',required=False,
		help='Specify number of significant figures to round to.  Default is none')
	parser.add_argument('--7b',default=False,action='store_true',
		required=False,dest='ttgr7b',
		help='Call this flag to specify that the output naming scheme should follow TtgR 7B naming conventions')
	parser.add_argument('--yaxis','-y',default='fold_induction',action='store',
		required=False,dest='yaxis_label',help='This flag specifies the yaxis label for any graph output')
	parser.add_argument('--xunits','-x',default='mM',action='store',required=False,
		dest='units',help='Labels used for graphs')
	parser.add_argument('--force_fmax','-f',action='store_true',required=False,
		default=False,
		help='Call this flag to force fmax to be the maximum fluorescence during the fit process')

	parsed = parser.parse_args()
	input_file = parsed.input
	concentrations = parsed.concentrations
	round = parsed.round
	ttgr7b = parsed.ttgr7b
	yaxis_label = parsed.yaxis_label
	xunits = parsed.units

	base_name = input_file.split('/')[-1].split('.')[0]
	doser_list = concentrations.split(',')
	conc_list = [float(x) for x in doser_list]

	node_id_dict = {
		1: 'C137I M167L',
		2: 'C137I F168Y',
		3: 'C137I\n M167L F168Y',
		4: 'C137I',
		5: 'I141W M167L',
		6: 'I141W F168Y',
		7: 'I141W\n M167L F168Y',
		8: 'I141W',
		9: 'C137I\n I141W M167L',
		10: 'C137I\n I141W F168Y',
		11: 'C137I I141W\n M167L F168Y',
		12: 'C137I I141W',
		13: 'M167L',
		14: 'F168Y',
		15: 'M167L F168Y',
		16: 'WT',
		}
	input_df = read_tsv(input_file)
	if len(conc_list) != input_df.shape[1]/2:
		input_df = trunchate_data(input_df,conc_list)
	else:
		pass
	output_df,failed_fits = fit_data(input_df,conc_list,yaxis_label,xunits,parsed.force_fmax)

	if round:
		rounded_df = round_df(output_df,round)
		rounded_df.to_csv(base_name+'_fit_out_rounded_{0}.tsv'.format(round),sep='\t',
			na_rep='NaN')
		output_df.to_csv(base_name+'_fit_out.tsv',sep='\t',na_rep='NaN')
	else:
		output_df.to_csv(base_name+'_fit_out.tsv',sep='\t',na_rep='NaN')
	print failed_fits
	graph_output(output_df,node_id_dict)
