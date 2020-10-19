#Nishikawa, Kyle K.
#Edited 20200827 for Epistasis Manuscript
#Script to graph the rawdata files from ITC experiments.  Files should be in the format time,sec

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np

#checks start
def is_csv(string):
	if string.lower().endswith('.csv'):
		return string
	else:
		msg = 'Please use a csv file'
		raise argparse.ArgumentTypeError(msg)
#---------------

def read_csv(input_file):
	if nano:
		temp_df = pd.read_csv(input_file,header=0,sep='\t')
	else:
		temp_df = pd.read_csv(input_file,header=0)
	return temp_df

def plot_heats(df):
	fig,ax = plt.subplots()
	time = df.iloc[:,0]
	heats = df.iloc[:,1]
	plt.plot(time,heats)
	if nano:
		ax.set_ylabel("Heat(uJ/s)",rotation=90)
		ax.set_xlabel("Time(s)")
	else:
		ax.set_ylabel("Heat(uCal/s)",rotation=90)
		ax.set_xlabel("Time(s)")
	plt.savefig(name+'.png')

if __name__== '__main__':
	parser = argparse.ArgumentParser(description='Script to graph raw data files from ITC experiments as heat vs sec.')

	parser.add_argument('--input','-i',required=True,action='store',dest='input_file',
		help='Input CSV file')
	parser.add_argument('--nano','-n',required=False,action='store_true',default=False,
		dest='nano',help='Use flag if raw data from nano ITC')

	parsed = parser.parse_args()
	input_file = parsed.input_file
	nano = parsed.nano
	no_slash = input_file.split('/')[-1]
	name = no_slash.split('.')[0]
	in_df = read_csv(input_file)
	#print in_df
	plot_heats(in_df)