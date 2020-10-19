#Nishikawa, Kyle K.
#Edited 20200826 for epistasis manuscript
#This script does bootstrap analysis of an input file

import argparse
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import decimal


def read_input(file):
	if file.endswith('.csv'):
		test = pd.read_csv(file,index_col=0)
		test = test.dropna(how='all',axis=0)
		return test
	elif file.endswith('.tsv'):
		test = pd.read_csv(file,index_col=0,sep='\t')
		test = test.dropna(how='all',axis=0)
		return test
	else:
		pass


def get_sample(arr,bootstrap):
	test = np.random.choice(arr,(arr.shape[0],bootstrap))
	out_df = pd.DataFrame(test)
	return out_df

def bootse(df):
	means = df.apply(np.mean,axis=0)
	se = np.std(means)
	return means,se

def jack_sample(arr):
	out_df = pd.DataFrame()
	for x in range(0,arr.shape[0]):
		out_df[x] = np.concatenate((arr[:x],arr[x+1:]),axis=None)
	return out_df
	
def jackse(df):
	means = np.mean(df)
	jse =  np.sqrt(np.float64((means.shape[0]-1.0)/means.shape[0])*np.sum((means-np.mean(means))**2))
	return means,jse


#Takes a 1D array as an input.
#Gives a corrected 100(1-2alpha) confidence interval
def bca_calc(arr,nboot,alpha,converge):
	if converge == False:
		theta_hat = np.nanmean(arr)
		nx = arr.shape[0]
		bsamp = get_sample(arr,nboot)
		print bsamp.shape
		bmeans,bse = bootse(bsamp)
		print bmeans
		jsamp = jack_sample(arr)
		jmeans,jse = jackse(jsamp)
	else:
		theta_hat = np.nanmean(arr)
		nx = arr.shape[0]
		bsamp = get_sample(arr,nboot)
		jsamp = jack_sample(arr)
		jmeans,jse = jackse(jsamp)
		while True:
			old_bmeans,old_bse = bootse(bsamp)
			bsamp = pd.concat([bsamp,get_sample(arr,nboot)])
			bmeans,bse = bootse(bsamp)
			break
	print bse, jse
	if bse == 0 or jse == 0:
		alpha1=alpha
		alpha2=1-alpha
		z0 = None
		ahat = None
	else:
		z0 = decimal.Decimal(
			scipy.stats.norm.ppf(np.sum(bmeans<theta_hat)/np.float64(nboot))
			)
		atop = decimal.Decimal(np.sum((np.mean(jmeans)-jmeans)**3.0))
		abot = decimal.Decimal(6.0*(((jse**2.0)*(nx/(nx-1.0)))**(3/2.0)))
		ahat = atop/abot
		print ahat,abot
		alpha1 = scipy.stats.norm.cdf(np.float64(
			z0+(z0+decimal.Decimal(scipy.stats.norm.ppf(alpha)))/(
			1-ahat*(z0+decimal.Decimal(scipy.stats.norm.ppf(alpha))))))
		alpha2 = scipy.stats.norm.cdf(np.float64(
			z0+(z0+decimal.Decimal(scipy.stats.norm.ppf(1.0-alpha)))/(
			1-ahat*(z0+decimal.Decimal(scipy.stats.norm.ppf(1-alpha))
			))))
	print alpha1,alpha2
	confpoint = np.quantile(bmeans,[alpha1,alpha2])
	percentiles = [alpha1,alpha2]
	u=jmeans-theta_hat
	return confpoint,z0,ahat,u,bmeans,bse,bsamp,percentiles

def ci_calc(arr,alpha):
	return np.quantile(arr,[alpha,1.0-alpha])

def graph_mean_rsq(bmeans,name,usecol):
	fig,ax = plt.subplots(1,1)
	fig.set_size_inches(8,6)
	ax.hist(bmeans)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.title('rsq_hist_cycles_{0}'.format(bmeans.shape[0]))
	plt.savefig(name+'_col_{0}_rsq_hist_bootstrap_{1}.pdf'.format(usecol,bmeans.shape[0]))

def format_out(perc,confpoint,z0,ahat,bse,alpha):
	out_df = pd.DataFrame()
	out_df['actual_percentile']=[alpha,1.0-alpha]
	out_df['adjusted_percentiles']=perc
	out_df['confidence_interval']=confpoint
	out_df['bias_correction_z0']=[z0,None]
	out_df['acceleration_factor_ahat']=[ahat,None]
	out_df['bootstrap_standard_error']=[bse,None]
	return out_df

#no bias correction if the bse == jse == 0
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Bootstrap functions for resampling and calculating standard error of the mean. No bias correction if bse == jse == 0')
	parser.add_argument('--input','-i',required=True,
		help='TSV or CSV input file')
	parser.add_argument('--alpha','-a',required=False,default=0.025,
		help='Float to calculate the 100(1-2*alpha) confidence interval')
	parser.add_argument('--nboot','-n',required=False,default=10000,
		help='Number of times to bootstrap the samples')
	parser.add_argument('--usecol','-u',required=False,default=0,
		help='Specify column to use. Default is 0')
	parser.add_argument('--write','-w',action='store_true',required=False,default=False,
		help='Call flag to write bootstrap output')
	parser.add_argument('--conf_int','-ci',required=False,action='store_true',
		default=False,
		help='Call this flag to calculate confidence interval without bootstrap')
	parser.add_argument('--converge','-c',required=False,action='store_true',
		default=False,
		help='Call this flag to bootstrap until convergence. Bootstrap iteration group size determined by the nboot flag')
	parser.add_argument('--threshold','-t',required=False,default=0.001,
		help='Call this flag to change threshold to which the iterate process goes to. Given as a percent difference (provide decimal) between current and previous iterations. Default is 0.001 (0.1 percent).')
	parsed = parser.parse_args()
	
	input_df = read_input(parsed.input)
	basename = parsed.input.split('/')[-1].split('.')[0]
	in_arr = np.array(input_df.iloc[:,int(parsed.usecol)])
	if parsed.conf_int == False:
		confpoint,z0,ahat,u,bmeans,bse,bsamp,perc = bca_calc(
			in_arr,np.int64(parsed.nboot),
			np.float64(parsed.alpha),parsed.converge)
		out_df = format_out(perc,confpoint,z0,ahat,bse,np.float64(parsed.alpha))
		if parsed.write==True:
			bsamp.to_csv('{0}_col_{1}_bootstrap_{2}.csv'.format(
				basename,parsed.usecol,parsed.nboot),sep=',')
			graph_mean_rsq(bmeans,basename,parsed.usecol)
			u = pd.DataFrame(u)
			u.to_csv('{0}_col_{1}_jknife_influence.csv'.format(
				basename,parsed.usecol),sep=',')
			out_df.to_csv('{0}_col_{1}_bootstrap_log.csv'.format(
				basename,parsed.usecol),sep=',')
		else:
			pass
			print 'Percentiles,{0}'.format(perc)
			print 'Confidence Interval,{0}'.format(confpoint)
			print 'Bias correction factor z0,{0}'.format(z0)
			print 'Acceleration factor ahat,{0}'.format(ahat)
			print 'Bootstrap standard error,{0}'.format(bse)
	else:
		conf_int = ci_calc(in_arr,np.float64(parsed.alpha))
		out_df = format_out(None,conf_int,None,None,None,np.float64(parsed.alpha))
		if parsed.write == True:
			out_df.to_csv('{0}_conf_int.csv'.format(basename),sep=',')
		else:
			print 'Percentiles,{0}'.format(parsed.alpha)
			print 'Confidence Interval,{0}'.format(conf_int)
	