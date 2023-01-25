from OMFITlib_MDS_obj import MDS_obj
import pandas as pd

scan=True
debug=False

output_dir='/u/mcurie/Output/' #finish with '/'

shot_min_list=[114100,206700]
shot_max_list=[142500,206800]
#shot_min_list=[132588]
#shot_max_list=[132595]

#clear the OUTPUTS
root['OUTPUTS']['DataFrame']=OMFITtree("")
root['OUTPUTS']['CSV']=OMFITtree("")


if scan:
	shot_list=[]
	for i in range(len(shot_min_list)):
		shot_min=shot_min_list[i]
		shot_max=shot_max_list[i]
		for shot_num in np.arange(shot_min,shot_max,dtype=int):
			shot_list.append(shot_num)
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

