#Nishikawa, Kyle K.
#Edited 20200827 for Epistasis manuscript
#Python script to adjust DH files based on a baseline injection heats file.

import argparse
import numpy as np
import itertools
import sys
#checks start
def is_DH(string):
	if string.lower().endswith('.dh'):
		return string
	else:
		msg = 'Please use a csv file'
		raise argparse.ArgumentTypeError(msg)




def grab_heats(baseline_file):
	heat_data = []
	with open(baseline_file,'r') as the_file:
		for line in itertools.islice(the_file,5,None):
			heat = float(line.split(',')[-1].rstrip())
			heat_data.append(heat)
	return heat_data

def average_heats(baseline_list,range_list):
	parsed_list = []
	if all_baselines == True:
		parsed_list = baseline_list
	else:
		for item in range_list:
			parsed_list.append(baseline_list[item])
	parsed_array = np.asarray(parsed_list)
	parsed_mean = np.mean(parsed_array)
	parsed_std = np.std(parsed_array)
	return parsed_mean,parsed_std

#Output needs to be in a DOS format (lines end with \r\n)
def subtract_heats(pmean,input_file):
	temp_heats = []
	with open(input_file,'r') as the_file:
		header = [next(the_file).rstrip() for x in xrange(6)]
		for line in itertools.islice(the_file,0,None):
			heat = float(line.split(',')[-1].rstrip())
			adj_heat = heat-pmean
			temp_heats.append('{0},{1}'.format(line.split(',')[0],adj_heat))
	for line in header:
		sys.stdout.write(line+'\r\n')
	for line in temp_heats:
		sys.stdout.write(line+'\r\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to graph raw data files from ITC experiments as heat vs sec.')

	parser.add_argument('--input','-i',required=True,action='store',dest='input_file',
		help='Input DH file')
	parser.add_argument('--control','-c',required=True,action='store',dest='control_file',
		help='Control DH file to subtract baseline heats from')
	parser.add_argument('--range','-r',required=True,nargs="*",action='store',
		dest='the_range',default="all",
		help='Select range of baselines to average (count from 0).  default is all.  Comma separated selects all between commas (inclusive).  Space means only written')
	parsed = parser.parse_args()
	input_file = parsed.input_file
	control_file = parsed.control_file
	start_range = parsed.the_range
	base_name = input_file.split('/')[-1].split('.')[0]

	final_range = []
	all_baselines = False
	for item in start_range:
		if item.lower() == 'all':
			all_baselines = True
		elif ',' in item:
			temp = item.split(',')
			start = int(temp[0])
			end = int(temp[1])
			test = list(range(int(temp[0]),int(temp[1])+1))
			for item in test:
				final_range.append(item)
		else:
			try:
				int(item)
			except:
				print "Please use either 'all' or an interger"
				quit()
			else:
				number = int(item)
				final_range.append(number)
	basename = input_file.split('/')[-1].split('.')[0]
	baseline_heats = grab_heats(control_file)
	pmean,pstd = average_heats(baseline_heats,final_range)
	subtract_heats(pmean,input_file)
		