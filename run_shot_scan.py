from OMFITlib_MDS_obj import MDS_obj
import pandas as pd

test=False
scan=True
debug=False

output_dir='/u/mcurie/Output/' #finish with '/'

if test:
	shot_num=132588
	MDS_obj=MDS_obj()
	df=MDS_obj.Auto_scan(shot_num=shot_num,plot=True)
	
	root['OUTPUTS']['DataFrame']=df
	df.to_csv('file.csv')
	root['OUTPUTS']['CSV'] = OMFITcsv('file.csv')
	root['OUTPUTS']['CSV'].deploy('/u/mcurie/Output/output.csv')
	
	#dic=df.to_dict()
	#root['OUTPUTS']['CSV']=OMFITcsv(filename='file.csv',dic)


if scan:
	shot_min=100000
	shot_max=250000

	#shot_min=132588
	#shot_max=132590

	shot_list=np.arange(shot_min,shot_max,dtype=int)
	work_list=np.zeros(len(shot_list))

	MDS_obj=MDS_obj()

	for i in range(len(shot_list)):
		shot_num=shot_list[i]
		print('shot_num='+str(shot_num))
		if debug:
			df=MDS_obj.Auto_scan(shot_num=shot_num,plot=False)
			
			root['OUTPUTS']['DataFrame']['df_tmp']=df
			df.to_csv('file.csv')
			root['OUTPUTS']['CSV'] = OMFITcsv('file.csv')
			root['OUTPUTS']['CSV'].deploy(output_dir+str(shot_num)+'.csv')
			work_list[i]=1
		else:
			try:
				df=MDS_obj.Auto_scan(shot_num=shot_num,plot=False)
				
				root['OUTPUTS']['DataFrame']['df_tmp']=df
				df.to_csv('file.csv')
				root['OUTPUTS']['CSV'] = OMFITcsv('file.csv')
				root['OUTPUTS']['CSV'].deploy(output_dir+str(shot_num)+'.csv')
				work_list[i]=1
			except:
				work_list[i]=0

	d={}
	d['shot_num']=shot_list
	d['work_list']=work_list

	df=pd.DataFrame.from_dict(d)
	root['OUTPUTS']['DataFrame']['summary']=df
	df.to_csv('file.csv')
	root['OUTPUTS']['CSV'] = OMFITcsv('file.csv')
	root['OUTPUTS']['CSV'].deploy(output_dir+'0summary.csv')

