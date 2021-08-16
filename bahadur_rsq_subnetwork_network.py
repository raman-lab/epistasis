#fancy network graphs for rsq values from bahadur_expansion2.py.

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import argparse
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib as mpl
import collections
import itertools

def read_tsv(input_file):
	return pd.read_csv(input_file,sep='\t',header=0,index_col=0)

def concatenate_dfs(input_dict,mut_list):
	out_df = pd.DataFrame()
	for x in mut_list:
		out_df[x]=input_dict[x].iloc[0,:]
	out_df.index = ['avg','stdev','num']
	return out_df

def get_shape_coords(num,centerx,centery,radius,control):
	angle=0
	
	if control is not None:
		n=num-1
	else:
		n=num
	angle_increment = (2*np.pi)/n
	coords = collections.OrderedDict()
	for bruh in range(0,n):
		x=centerx + (radius*np.cos(angle))
		y=centery + (radius*np.sin(angle))
		angle += angle_increment
		coords[bruh]=(x,y)
	if control is not None:
		coords[num-1]=(centerx,centery)
	return coords

#this function assumes *rsq_avgs.tsv files are used with a subnetwork size of 2 mutations
def parse_combinations(input_list,control):
	binary_str = collections.OrderedDict()
	key_dict = collections.OrderedDict()
	for file in input_list:
		file_list = file.split('_')
		mut_count = file_list[-6].count('1')
		filename = '{0}_{1}_{2}_{3}'.format(file_list[-6],
			file_list[-5],file_list[-4],file_list[-3])
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

def make_scatterplot(g,df,coords,name_dict,color,cbar_label,cmap_norm,txt,filename,largefig):
	val_dict = collections.OrderedDict()
	g.add_nodes_from(coords.keys())
	for n,p in coords.items():
		g.nodes[n]['pos'] = p
	for n,p in name_dict.items():
		g.nodes[n]['name']=p
	for n in range(0,df.shape[1]):
		g.nodes[n]['data'] = df.iloc[0,n]
		val_dict[n]=df.iloc[0,n]
	fig,ax1 = plt.subplots(1,1)
	fig.set_size_inches(8,8)
	virid = plt.get_cmap(color)
	if cmap_norm == False:
#		cnorm = colors.Normalize(vmin=np.amin(df.iloc[0,:]),
#			vmax=np.amax(df.iloc[0,:]))
		cnorm = colors.Normalize(vmin=0,
			vmax=1)
	else:
		cnorm = colors.DivergingNorm(vmin=np.amin(df.iloc[0,:]),
			vmax=np.amax(df.iloc[0,:]),vcenter=0)
	scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=virid)
	networkx_colors = [scalarmap.to_rgba(df.iat[0,n]) for n in range(0,df.shape[1])]
#	print(len(networkx_colors))
	for n in range(0,df.shape[1]):
		colnames = df.columns.values.tolist()
		ax1.plot([0],[0],color=scalarmap.to_rgba(df.iat[0,n]),label=colnames[n],
			linewidth=1
		)
	if largefig == True:
		nsize = 2000
		xaxscale=1.25
		labeldist = 0.15
	else:
		nsize = 3000
		xaxscale=1.5
		labeldist=0.2
	bruh = nx.draw_networkx_nodes(g,coords,node_shape='o',
		node_color=networkx_colors,
		with_labels=False,ax=ax1,node_size=nsize)
	bruh.set_edgecolor('black')
	nx.draw_networkx_edges(g,coords,linewidth=0.5)
	ax1.set_xlim(xmin=-xaxscale,xmax=xaxscale)
	ax1.set_ylim(ymin=-xaxscale,ymax=xaxscale)
	plt.axis('off')
	fig.set_facecolor('w')
	scalarmap._A=[]
	cbar = plt.colorbar(scalarmap)
	cbar.ax.set_title(cbar_label)
	
	if txt == True:
		for num in range(0,df.shape[1]):
			plt.text(coords[num][0],coords[num][1]+labeldist,
				df.columns.values[num],ha = 'center',fontsize = 10)
			plt.text(coords[num][0],coords[num][1]-labeldist,
				str(df.iat[0,num])[0:5],ha='center',fontsize=7)
	else:
		pass
	plt.tight_layout()
	if txt == True:
		plt.savefig(filename+'_rsq_network_{0}_labeled.pdf'.format(color),transparent=True)
	else:
		plt.savefig(filename+'_rsq_network_{0}.pdf'.format(color),transparent=True)
	plt.clf()

		
		
def combine_funcs(input_file_list,control,colormap,cbar_label,normalize,txt,filename):
	binary_str,key_dict = parse_combinations(input_file_list,control)
#	print binary_str
#	print key_dict
	for k in key_dict:
		concat_df = concatenate_dfs(binary_str,key_dict[k])
		print(concat_df)
		print(concat_df.shape[1])
		
		network_coords = get_shape_coords(concat_df.shape[1],0,0,1,control)
		print(network_coords)
		print(len(network_coords.keys()))
		name_dict = collections.OrderedDict()
		for x in range(0,len(key_dict[k])):
			name_dict[x] = key_dict[k][x]
		print(name_dict)
		g=nx.Graph()
		test=np.amin(concat_df.iloc[0,:])
		if concat_df.shape[1]>7:
			make_scatterplot(g,concat_df,
				network_coords,name_dict,colormap,cbar_label,False,txt,
				'{0}_{1}_set'.format(filename,k),True)
		else:
			make_scatterplot(g,concat_df,
				network_coords,name_dict,colormap,cbar_label,False,txt,
				'{0}_{1}_set'.format(filename,k),False)
	
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='fancy network graphs for rsq values from bahadur_expansion2.py.')
	required = parser.add_argument_group('Required')
	required.add_argument('--input','-i',nargs='*',
		help='1 or more *rsq_avg.tsv files generated from bahadur_expansion2.py')
	parser.add_argument('--control','-r',default=None,
		help='Control rsq_avg.tsv file for randomized data OR non-epistatic data')
	parser.add_argument('--text','-t',required=False,default=True,
		action='store_false',help='Call this flag to turn labeling OFF')
	parser.add_argument('--colormap','-c',default='Greens',action='store',
		help='Colormap to use')
	parser.add_argument('--label','-x',default='Rsquared',action='store',
		help='Label for colorbar')
	parser.add_argument('--normalize','-n',default=False,action='store_true',
		help='Call this flag to normalize the color bar to x=0')
	parser.add_argument('--title',default='selected',
		help='Define title')
	parser.add_argument('--bar','-b',required=False,action='store_true',default=False,
		help='Call this flag to force a bar graph for each set of subnetworks.  Replaces network graph')
		
	parsed = parser.parse_args()
	combine_funcs(parsed.input,parsed.control,parsed.colormap,
		parsed.label,parsed.normalize,parsed.text,parsed.title)