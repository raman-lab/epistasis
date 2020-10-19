#Nishikawa, Kyle K.
#Edited 20200827 for Epistasis Manuscript
#Generates a network graph using Plotly.
#Note: This uses python3

import networkx as nx
import itertools
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib as mpl
import collections
import numpy as np
#import plotly.plotly as py
#import plotly.graph_objs as go
from mpl_toolkits.mplot3d import Axes3D
import mayavi.mlab as mlab
import sys


#---Argparse checks start
def tsv(file):
	if file.endswith('.tsv') or file.endswith('.tsv'):
		return file
	else:
		msg = "Please use a valid tsv file"
		raise argparse.ArgumentTypeError(msg)
#---Argparse checks end

	
#normalizing based on log2
def normalization(df,base):
	temp_dict = collections.OrderedDict()
	temp_list = df.iloc[:,base].tolist()
	temp_array = np.array(temp_list)
	temp_dict[0] = np.log2(temp_array)
	temp_df=pd.DataFrame.from_dict(temp_dict)
	return temp_df

#We will need an input .xlsx file with three columns: node name \ avg induction \ st dev
def read_tsv(input_file):
	input_df = pd.read_csv(input_file,delimiter='\t',index_col=0)
	return input_df

def tsv_check(input_df):
	rows = input_df.shape[0]
	if rows != len(node_id_dict):
		sys.stdout.write("Not the right number of nodes for the TtgR mutant.  Quitting...")
		quit()
	else:
		pass

def make_graph(g,df,color,data_col,cbar_label,cmap_norm,labeling):
	dict_keys = node_pos.keys()
	max_vals = df.max()
	max_val = max_vals[0]
#	sys.stdout.write(max_val)
	val_dict = collections.OrderedDict()
	g.add_nodes_from(node_pos.keys())
	for n,p in node_pos.items():
		g.nodes[n]['pos'] = p
	for n,p in node_id_dict.items():
		g.nodes[n]['name']=p
	for n in range(0,df.shape[0]):
		g.nodes[n]['data'] = df.iloc[n,data_col]
		val_dict[n]=df.iloc[n,data_col]
	val_df = pd.DataFrame.from_dict(val_dict,orient='index')
	tests = val_df.max()
	test = tests[0]
	mins = val_df.min()[0]
#	mins=0
	sys.stdout.write("val_df: \n")
	print(val_df)
	fig,ax1 = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	fold_induction = input_df.iloc[:,0]
	virid = plt.get_cmap(color)
	if cmap_norm == False:
		cnorm = colors.Normalize(vmin=mins,vmax=test)
	else:
		cnorm = colors.DivergingNorm(vmin=mins,vmax=test,vcenter=0)
	scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=virid)
	networkx_colors = [scalarmap.to_rgba(val_dict[val]) for val in val_dict]
	
	for val in val_dict:
		ax1.plot([0],[0],color=scalarmap.to_rgba(val_dict[val]),label=node_id_dict[val],
			linewidth=1
		)
	net = nx.draw_networkx_nodes(g,node_pos,node_shape='o',
		node_color=networkx_colors,
		with_labels=False,ax=ax1,node_size=1100)
	net.set_edgecolor('black')
	nx.draw_networkx_edges(g,node_pos,linewidth=0.5)
	plt.axis('off')
	fig.set_facecolor('w')
#	plt.legend(loc='best')
	scalarmap._A=[]
	
	cbar = plt.colorbar(scalarmap)
	if parsed.log2 == True:
		cbar.ax.set_title('{0}_{1}'.format('log2',cbar_label))
	else:
		cbar.ax.set_title(cbar_label)
	
	if txt == True:
		for num in range(0,len(node_binary_dict)):
			plt.text(node_pos[num][0],node_pos[num][1]+3.5,node_binary_dict[num],ha = 'center',fontsize = 10)
			plt.text(node_pos[num][0],node_pos[num][1]-4.5,str(val_dict[num])[0:5],ha='center',fontsize=7)
	else:
		pass
	plt.tight_layout()
	if labeling == False:
		plt.savefig(filename+'_network.pdf',transparent=True)
	else:
		plt.savefig(filename+'_network_labeled.pdf',transparent=True)



if __name__ == '__main__':
	parser=argparse.ArgumentParser(description='Script to generate a network graph keyed for the TtgR_7B mutation intermediates')
	parser.add_argument('--input','-i',type=tsv,action='store',dest='input_file',
		required=True,help='Input TSV file to generate the network graphs')
	parser.add_argument('--3d','-d',type=tsv,action='store',dest='threed',
		required=False,default=None,
		help='Optional TSV file to generate network graphs incorporating baseline fluorescence and fold induction')
	parser.add_argument('--baseline','-b',action='store',dest='base',required=False,
		default=0,help='Optional Argument to specify baseline value')
	parser.add_argument('--text','-t',dest='txt',required=False,default=True,
		action='store_false',help='Call this flag to turn labeling OFF')
	parser.add_argument('--log2','-l',required=False,default=False,action='store_true',
		help='Call this input to force log2 scaling of input')
	parser.add_argument('--colormap','-c',default='Blues',action='store',
		help='Colormap to use')
	parser.add_argument('--use_col','-u',default=0,action='store',
		help='Column to use in the input file.  Default 0, 0-indexed (first column is rownames, second column is index 0)')
	parser.add_argument('--label','-x',default='RFU',action='store',
		help='Label for colorbar')
	parser.add_argument('--normalize','-n',default=False,action='store_true',
		help='Call this flag to normalize the color bar to x=0')
	parser.add_argument('--log10','-log10',required=False,default=False,
		action='store_true',
		help='Call this argument to force log10 scaling of input')

	parsed = parser.parse_args()
	input_file = parsed.input_file
	threed = parsed.threed
	base=parsed.base
	txt = parsed.txt

	#Set up node names and locations
	node_id_dict = {
		0: 'C137I M167L',
		1: 'C137I F168Y',
		2: 'C137I\n M167L F168Y',
		3: 'C137I',
		4: 'I141W M167L',
		5: 'I141W F168Y',
		6: 'I141W\n M167L F168Y',
		7: 'I141W',
		8: 'C137I\n I141W M167L',
		9: 'C137I\n I141W F168Y',
		10: 'C137I I141W\n M167L F168Y',
		11: 'C137I I141W',
		12: 'M167L',
		13: 'F168Y',
		14: 'M167L F168Y',
		15: 'WT',
		}
	node_binary_dict = {
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
	g = nx.Graph()

	node_pos = {
		0: (10,32.5),
		1: (10,20),
		2: (15,20),
		3: (5,20),
		4: (10,7.5),
		5: (10,-7.5),
		6: (15,7.5),
		7: (5,7.5),
		8: (15,-7.5),
		9: (15,-20),
		10: (20,0),
		11: (10,-20),
		12: (5,-7.5),
		13: (5,-20),
		14: (10,-32.5),
		15: (0,0)
		}
	input_df = read_tsv(input_file)
	#input_df = input_df.iloc[:,int(parsed.use_col)]
	print(input_df)
	tsv_check(input_df)
	filenames = input_file.split('/')[-1].split('.')
	filename = filenames[0]
	norm_df = normalization(input_df,0)
	if threed:
		try:
			int(base)
		except:
			sys.stdout.write("Use an interger\n")
			quit()
		base=int(base)
		threed_df = read_tsv(threed)
	#	sys.stdout.write threed_df
	#	sys.stdout.write threed_df.iloc[:,base]
		tsv_check(threed_df)
		base_norm_df = normalization(threed_df,base)
		if int(base) in [x for x in range(0,threed_df.shape[1])]:
			pass
		else:
			sys.stdout.write('Int not in range')
			quit()
		make_3d_graph(input_df,threed_df,base)
	
	else:
		if parsed.log2 == True and parsed.log10 == False:
			test = np.log2(input_df.iloc[:,0])
			input_df['log_data']=test
			print(input_df)
			filename = filename+'_log2'
			data_col = input_df.columns.values.tolist().index('log_data')
		elif parsed.log10 == True and parsed.log2 == False:
			test = np.log10(input_df.iloc[:,0])
			input_df['log10_data'] = test
			print(input_df)
			filename = filename+'_log10'
			data_col = input_df.columns.values.tolist().index('log10_data')
		else:
			data_col = 0
		make_graph(g,input_df,parsed.colormap,data_col,parsed.label,parsed.normalize,txt)
