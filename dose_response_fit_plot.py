#Nishikawa, Kyle K.
#Edit 20200826 for Epistasis manuscript
#Python script to re-plot dose_response_fit.py output as pdfs in a different format for figures

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_parameters(input_file):
	return pd.read_csv(input_file,index_col=0,delimiter='\t')

def func2(conc,n,kd,fmax):
	return baseline + ((fmax - baseline) * ((conc**n)/((kd**n)+(conc**n))))

def func3(conc,n,kd):
	return baseline + ((fmax - baseline) * ((conc**n)/((kd**n)+(conc**n))))
	
def rmscalc(yact,yfit):
	return np.sqrt(np.sum((yfit-yact)**2)/yfit.shape[0])
	
def plot_input(parameters_df,input_vals_df,use_row,conc_list,xunits,ylabel,base_name,
	data_colors,fit_colors,marker_list,no_stds):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	minconc = float(min(conc_list))
	rownames = input_df.index.values
	temp_title = '_'.join(str(x) for x in use_row)
	fit_xvals = np.arange(0,max(conc_list)*10)
	for row in range(0,len(use_row)):
		params = parameters_df.iloc[use_row[row],0:3]
		print params
		fold_ind = input_df.iloc[use_row[row],:].tolist()[0:][::2]
		fmax_approx = max(fold_ind)	
		global baseline
		baseline = fold_ind[conc_list.index(minconc)]
		stdevs = input_df.iloc[use_row[row],:].tolist()[1:][::2]
		baseline_stds = stdevs[conc_list.index(minconc)]
		print params
		ax.errorbar(fit_xvals,func2(fit_xvals,*params),color=fit_colors[row],
			label='{0}_Fit'.format(use_row[row]),zorder=2*row+1
			)
		if '-' and '_' in str(rownames[0]):
			rownames = [re.sub('\_|\-','',x) for x in rownames]
		else:
			pass
		y_tups = [(conc_list[x],fold_ind[x],stdevs[x]) for x in range(0,len(conc_list))]
		y_tups = sorted(y_tups,key=lambda x: x[0])
		y_zipped = zip(*y_tups)
		print y_zipped[1]
		print y_zipped[2]
		rmserror = rmscalc(y_zipped[1],func2(np.asarray(y_zipped[0]),*params))
		if no_stds == False:
			ax.errorbar(y_zipped[0],y_zipped[1],yerr=y_zipped[2],fmt=marker_list[row],
				color=data_colors[row],capsize=5, elinewidth=1, markeredgewidth=1,
				label='{0}_Experimental'.format(use_row[row]),zorder=2*row+2
				)
		else:
			ax.errorbar(y_zipped[0],y_zipped[1],fmt=marker_list[row],
				color=data_colors[row],capsize=5, elinewidth=1, markeredgewidth=1,
				label='{0}_Experimental'.format(use_row[row]),zorder=2*row+2
				)
	ax.set_xlabel('Concentration ({0})'.format(xunits))
	ax.set_ylabel(ylabel)
	ax.set_xscale('log')
	ax.set_xlim(right=max(conc_list)*5)
#	ax.set_ylim(ymin=0)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.legend(loc=2)
	plt.title(temp_title)
	plt.tight_layout()
	if no_stds == False:
		plt.savefig(base_name+'_'+temp_title+'_'+'doseresponse_fit.pdf')
	else:
		plt.savefig(base_name+'_'+temp_title+'_'+'nostd_doseresponse_fit.pdf')
	plt.close(fig)
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''
		Python script to re-plot dose_response_fit.py output 
		as pdfs in a different format for figures
		'''
		)
	required=parser.add_argument_group('Required')
	required.add_argument('--input','-i',required=True,
		help='Input file with fluorescences and standard deviations from flowjo_parse2.py'
		)
	required.add_argument('--concs','-c',required=True,
		help='''
			Comma-separated list of concentrations in the order that they were listed for 
			the input file
			'''
		)
	required.add_argument('--params','-p',required=True,
		help='TSV file with fit parameters')
	parser.add_argument('--use_row','-u',default=0,required=False,
		help='''
			Row (0-indexed) to use for both the parameters and the input file.
			Can also be a comma-separated list of rows to plot on the same graph
			''')
	parser.add_argument('--xunits','-x',default='uM',required=False,
		help='Specify units for x-axis, default uM')
	parser.add_argument('--ylabel','-y',default='Fluorescence',required=False,
		help='Y-axis label')
	parser.add_argument('--nostds','-s',default=False,required=False,action='store_true',
		help='Call this flag to Not use the standard deviations in the plot')
	parsed = parser.parse_args()
	input_df = read_parameters(parsed.input)
	basename=parsed.input.split('/')[-1].split('.')[0]
	params_df = read_parameters(parsed.params)
	data_colors = ['#000000','#478ebf','#fa5340','#28f93d','#5340fa','#f9a628']
	fit_colors = ['#858585','#80b1d3','#fb8072','#72fb80','#8072fb','#fbc572']
	marker_shapes = ['o','v','s','P','x','D']
	if ',' in parsed.use_row:
		use_row = [int(x) for x in parsed.use_row.split(',')]
	else:
		use_row = [int(parsed.use_row)]
	conc_list = [np.float64(x) for x in parsed.concs.split(',')]
	print conc_list
	print input_df
	plot_input(params_df,input_df,use_row,conc_list,parsed.xunits,parsed.ylabel,basename,
		data_colors,fit_colors,marker_shapes,parsed.nostds)