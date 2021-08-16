#Function to plot bahadur rsq_raw values as a histogram.  Uses rsq_raw.tsv files from bahadur_expansion2.py


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse


def read_input(input):
	return pd.read_csv(input,header=0,sep='\t',index_col=0)

def plot_hist(random_rsq_df,col,name):
	bruh = np.array(random_rsq_df.iloc[:,col].copy())
	bruh = bruh[~np.isnan(bruh)]
	print bruh
	cycles=random_rsq_df.shape[0]
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	test = ax.hist(bruh)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('$R^2$')
	ax.set_ylabel('Count')
	plt.title('rsq_hist_cycles_{0}'.format(cycles))
	plt.savefig(name+'_{0}_hist_cycles_{1}.pdf'.format(col,cycles),transparent=True)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Bahadur script to make histograms from rsq_raw values based on bahadur_expansion2.py')
	parser.add_argument('--input','-i',required=True,
		help='Input tsv of raw_rsq values')
	parsed=parser.parse_args()
	base = parsed.input.split('/')[-1].split('.')[0]
	input_df = read_input(parsed.input)
#	print input_df
	for col in range(1,input_df.shape[1]):
		plot_hist(input_df,col,base)