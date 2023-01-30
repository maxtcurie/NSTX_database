import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

from OMFITlib_MDS_obj import MDS_obj

shot_num=132591 #132588, 132591 
device='nstx' #'nstx', 'd3d'
MDS_obj_0=MDS_obj(device=device,shot_num=shot_num)
mode='dB' #'all', 'dB'

if mode=='all':
	#get all the quantities
	quant_list=MDS_obj_0.get_all_quant(plot=False)
	MDS_obj_0.chose_time(plot=True)
	#cut the R from psi=0 to psi=0.99 ish
	MDS_obj_0.cut_R()
	#cut time to from 0.3s to 0.8s
	MDS_obj_0.cut_t()
	#MDS_obj_0.test_save_worthy()
	#interpolation all the quantities to uniform grid
	quant_list=MDS_obj_0.interp_all_quant(interp_factor=1.2,plot=False)
	#calculate the average Lref
	Lref,Lref_err=MDS_obj_0.calc_Lref_avg(plot=False)
	#calculate the dn_dR=dne/dR
	dn_dR=MDS_obj_0.calc_dn_dR(plot=False)
	#calculate the dT_dR=dTe/dR
	dT_dR=MDS_obj_0.calc_dT_dR(plot=False)
	#calculate the R, a/Lne, ne in dne/dR peak
	R_ne_list,a_Lne_list,ne_ped_list,ne_ped_index_list=MDS_obj_0.calf_Ln_peak(plot=False)
	#calculate the R, a/LTe, Te in dTe/dR peak
	R_Te_list,a_LTe_list,Te_ped_list,Te_ped_index_list=MDS_obj_0.calf_Lt_peak(plot=False)
	#calculate the collision
	coll_ei=MDS_obj_0.calc_coll_ei(plot=False)
	coll_ei_ped_value=[coll_ei[i,ne_ped_index_list[i]] for i in range(len(ne_ped_index_list))]
	#calculate the beta
	beta=MDS_obj_0.calc_beta(plot=False)
	beta_ped_value=[beta[i,ne_ped_index_list[i]] for i in range(len(ne_ped_index_list))]
	t=MDS_obj_0.quant_list['Te']['time']



	d={}
	d['shot_num']=[shot_num]*len(t)
	d['time']=t
	d['Lref']=MDS_obj_0.quant_list['Lref']['data']
	d['beta']=beta_ped_value
	d['a_Lne']=a_Lne_list
	d['a_LTe']=a_LTe_list
	d['coll_ei']=coll_ei_ped_value
	df=pd.DataFrame(d)

	print(df)

if mode=='dB':
	MDS_obj_0.get_dB1(plot=False)
	MDS_obj_0.find_dB1_peak(plot=True)