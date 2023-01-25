import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path='/u/mcurie/Output/' #end with '/'
dir_list = os.listdir(path)
print(dir_list)

#dir_list=dir_list[:30]


mode=2 	#mode 1: merge all the files and plot
		#mode 2: calculate the average and std, then plot

df=pd.read_csv(path+'114100.csv')
keys=df.keys()
#df_merge=df

if mode==1:
	for i in dir_list:
		if i[0]=='0':
			continue

		df=pd.read_csv(path+i)
		df_merge=pd.concat([df_merge, df], axis=0)

	df_good=df_merge[(df_merge.Lref>0.2)]

	plt.clf()
	plt.scatter(df_merge['coll_ei'],df_merge['a_LTe'],alpha=0.3)
	plt.xlabel(r'$/nu_{ei}$ (kHz)')
	plt.ylabel(r'$a/L_{Te}$')

if mode==2:
	d={}
	for key in keys:
			d[key]=[]
			d[key+'_err']=[]

	for i in dir_list:
		print(i)
		#print(i[0])
		if i[0]=='1' or i[0]=='2':
			pass
		else:
			continue
		df=pd.read_csv(path+i)
		for key in keys:
			d[key].append(np.mean(df[key]))
			d[key+'_err'].append(np.std(df[key]))

	df=pd.DataFrame.from_dict(d)
	df_work=df[(df.Lref>0.2)]
	print(df)
	plt.clf()
	#plt.plot(df_work['coll_ei'],df_work['a_LTe'],,alpha=0.01)
	plt.errorbar(df_work['coll_ei'],df_work['a_LTe'],\
			xerr=df_work['coll_ei_err'],\
			yerr=df_work['a_LTe_err'],\
			marker='s',ms=10,\
			linestyle='none',alpha=0.02)
	plt.xlabel(r'$\nu_{ei}$ (kHz)')
	plt.xscale('log')
	plt.ylabel(r'$a/L_{Te}$')
	plt.savefig(path+'0nu_Lt_log_err_plot.png')

	plt.clf()
	#plt.plot(df_work['coll_ei'],df_work['a_LTe'],,alpha=0.01)
	plt.errorbar(df_work['coll_ei'],df_work['a_LTe'],\
			xerr=df_work['coll_ei_err'],\
			yerr=df_work['a_LTe_err'],\
			marker='s',ms=10,\
			linestyle='none',alpha=0.02)
	plt.xlabel(r'$\nu_{ei}$ (kHz)')
	#plt.xscale('log')
	plt.ylabel(r'$a/L_{Te}$')
	plt.savefig(path+'0nu_Lt_err_plot.png')

	plt.clf()
	plt.scatter(df_work['coll_ei'],df_work['a_LTe'],s=1,alpha=0.2)
	#plt.errorbar(df_work['coll_ei'],df_work['a_LTe'],\
	#		xerr=df_work['coll_ei_err'],\
	#		yerr=df_work['a_LTe_err'],\
	#		marker='s',ms=10,\
	#		linestyle='none',alpha=0.02)
	plt.xlabel(r'$\nu_{ei}$ (kHz)')
	plt.xscale('log')
	plt.ylabel(r'$a/L_{Te}$')
	plt.savefig(path+'0nu_Lt_log_plot.png')

	plt.clf()
	plt.scatter(df_work['coll_ei'],df_work['a_LTe'],s=1,alpha=0.2)
	#plt.errorbar(df_work['coll_ei'],df_work['a_LTe'],\
	#		xerr=df_work['coll_ei_err'],\
	#		yerr=df_work['a_LTe_err'],\
	#		marker='s',ms=10,\
	#		linestyle='none',alpha=0.02)
	plt.xlabel(r'$\nu_{ei}$ (kHz)')
	#plt.xscale('log')
	plt.ylabel(r'$a/L_{Te}$')
	plt.savefig(path+'0nu_Lt_plot.png')


	plt.clf()
	#plt.plot(df_work['coll_ei'],df_work['a_LTe'],,alpha=0.01)
	plt.errorbar(df_work['shot_num'],df_work['Lref'],\
			xerr=df_work['shot_num_err'],\
			yerr=df_work['Lref_err'],\
			marker='s',ms=10,\
			linestyle='none',alpha=0.02)
	plt.xlabel('shot_num')
	#plt.xscale('log')
	plt.ylabel('Lref')
	plt.savefig(path+'0shot_num_Lref_err_plot.png')

