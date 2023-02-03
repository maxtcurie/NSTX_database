import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

from OMFITlib_FFT_general import spectral_density

test=True

mref=2.
m_SI = mref *1.6726*10**(-27)
me_SI = 9.11*10**(-31)
c  = 1.
qref = 1.6*10**(-19)

#Reference on NSTX MSD database: 
#https://nstx.pppl.gov/nstx/Software/FAQ/signallabels.html


#to future self: #GENERAL needed require attention for general cases
class MDS_obj:
	#profile_type=('pfile','ITERDB')	

	def __init__(self,device='nstx',shot_num=132588):
		self.shot_num=shot_num
		self.device=device
		self.quant_list={}
		self.coord={}
		if self.device=='nstx':
			self.keys_t_r=['Bt', 'Te', 'ne']
			self.keys_psi=['q0psi']
			self.keys_t=['Bref', 'Lref']
			self.R_min=0.90
			self.R_max=1.45
			self.t_min=0.3
			self.t_max=0.8

		elif self.device=='d3d':
			self.keys_t_r=['Te', 'ne']
			self.keys_t=['Bref', 'Lref']
			self.R_min=0.90
			self.R_max=1.45
			self.t_min=0.3
			self.t_max=0.8

	def quant_report(quant):
		#print(quant['info'])
		print(np.shape(quant['data']))
		print(np.shape(quant['time']))
		try:
			print(np.shape(quant['R']))
		except:
			pass
	
	def finite_differences(self,var,grid):

		def get_mat_fd_d1_o4(size,dx,plot_matrix=False):
			"""Creates matrix for centered finite difference, first derivative, 4th order.
			size: size of (number of elements in) quantity to be differentiated
			dx: grid spacing (for constant grid)."""

			prefactor=1.0/(12.0*dx)
			mat=np.zeros((size,size),dtype='float')	
			for i in range(size):
				if i-1 >= 0:
					mat[i,i-1]=-8
				if i-2 >= 0:
					mat[i,i-2]=1
				if i+1 <= size-1:
					mat[i,i+1]=8
				if i+2 <= size-1:
					mat[i,i+2]=-1
		   
			mat=prefactor*mat

			if plot_matrix:
				plt.contourf(mat,50)
				plt.colorbar()
				plt.show()

			return mat

			
		def fd_d1_o4(var,grid,mat=False):
			"""Centered finite difference, first derivative, 4th order.
			var: quantity to be differentiated.
			grid: grid for var 
			mat: matrix for the finite-differencing operator. if mat=False then it is created"""

			if not mat:
				mat=get_mat_fd_d1_o4(len(var),grid[1]-grid[0])

			dvar=-np.dot(mat,var)
			dvar[0]=0.0
			dvar[1]=0.0
			#dvar[2]=0.0
			dvar[-1]=0.0
			dvar[-2]=0.0
			#dvar[-3]=0.0
			return -dvar 

		
		dvar=fd_d1_o4(var,grid)
		return dvar


	#z2=interp_2D(z,x,y,x2,y2,plot=False)
	def interp_2D(self,z,x,y,x2,y2,plot=False): #for z(x,y)
		(nx,ny)=np.shape(z)

		z_tmp=np.zeros((nx,len(y2)))
		#print(len(x2))
		#print(len(y2))
		for i in range(nx):
			z_tmp[i,:]=np.interp(y2,y,z[i,:])

		z2=np.zeros((len(x2),len(y2)))
		for i in range(len(y2)):
			z2[:,i]=np.interp(x2,x,z_tmp[:,i])
		if plot:
			plt.clf()
			plt.plot(y,z[0,:],label='original',alpha=0.1)
			plt.plot(y2,z2[0,:],label='interp',alpha=0.1)
			for i in range(len(x)):
				plt.plot(y,z[i,:],alpha=0.1)
			for i in range(len(x2)):
				plt.plot(y2,z2[i,:],alpha=0.6)
			plt.legend()
			plt.show()
		return z2

	def smooth(self,y, box_pts):
		box = np.ones(box_pts)/box_pts
		y_smooth = np.convolve(y, box, mode='same')
		return y_smooth

	def set_shot_num(self,shot_num):
		self.shot_num=shot_num
		self.t_min=0.3
		self.t_max=0.8

	def set_device(self,device):
		self.device=device
		self.t_min=0.3
		self.t_max=0.8
		if device=='nstx':
			self.R_min=0.90
			self.R_max=1.45

	def get_Zeff(self,plot=False):
		pass
		#\ACTIVESPEC::TOP.CHERS.ANALYSIS.CT1:ZEFF
	
	#quant['ne']=MDS_obj.get_ne(plot=True) #ne(time,r)
	def get_ne(self,plot=False):
		quant={}
		quant['name']='ne'
		entry_tmp=OMFITmdsValue(server=self.device,treename='activespec',shot=self.shot_num,TDI='\\ACTIVESPEC::TOP.CHERS.ANALYSIS.CT1:DEN')
		quant['data']=entry_tmp.data() #ne(time,r)
		(nt,nr)=np.shape(quant['data'])
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(1)
		quant['R']=entry_tmp.dim_of(0)*0.01

		if plot:
			plt.plot(quant['R'],quant['data'][int(0.5*nt),:])   #?????????
			plt.xlabel('R (m)')
			plt.ylabel('ne (cm-3)')
		self.quant_list[quant['name']]=quant
		return quant

	#quant['Te']=MDS_obj.get_Te(plot=True) #Te(time,r)
	def get_Te(self,plot=False):
		quant={}
		quant['name']='Te'
		entry_tmp=OMFITmdsValue(server=self.device,treename='activespec',shot=132588,TDI='\\ACTIVESPEC::TOP.MPTS.OUTPUT_DATA.BEST:FIT_TE	')
		quant['data']=entry_tmp.data() #te(time,r)
		quant['data']=quant['data'].T
		(nt,nr)=np.shape(quant['data'])
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)
		quant['R']=entry_tmp.dim_of(1)*0.01

		if plot:
			plt.plot(quant['R'],quant['data'][int(0.5*nt),:])   
			#plt.plot(quant['time'],quant['data'][int(0.5*nr),:])   #?????????
			plt.xlabel('R (m)')
			plt.ylabel('Te (keV)')
			plt.ylim(0,1)
		self.quant_list[quant['name']]=quant
		return quant

	#quant['q0psi']=MDS_obj.get_qpsi(plot=True)  #q(time,r)
	def get_q0psi(self,plot=False):
		quant={}
		quant['name']='q0psi'
		entry_tmp=OMFITmdsValue(server=self.device,treename='efit01',shot=self.shot_num,TDI='\\EFIT01::QPSI')
		quant['data']=entry_tmp.data() #q(time,r)
		(nt,nr)=np.shape(quant['data'])
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)
		quant['r_t']=entry_tmp.dim_of(1) #r(time)=r[time,r]
		quant['R']=entry_tmp.dim_of(1)[0,:] #r(time)=r[time,r]

		if plot:
			plt.scatter(quant['r_t'][int(0.5*nt),:],\
						quant['data'][int(0.5*nt),:])
			plt.xlabel('psi')
			plt.ylabel('q0')
		self.quant_list[quant['name']]=quant
		return quant

	#quant['q0']=MDS_obj.get_q0(plot=True)  #q(time) q0 on magnetic axis
	def get_q0(self,plot=False):
		quant={}
		quant['name']='q0'
		entry_tmp=OMFITmdsValue(server=self.device,treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:Q0')
		quant['data']=entry_tmp.data() #q(time)
		#(nt,nr)=np.shape(quant['data'])
		#print(np.shape(quant['data']))
		
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)

		if plot:
			plt.scatter(quant['time'],quant['data'])
			plt.xlabel('time (s)')
			plt.ylabel('q0')
		self.quant_list[quant['name']]=quant
		
		return quant

	#quant['psin']=MDS_obj.get_psi_R(plot=True)  #psi, R
	def get_psi_R(self,plot=False):
		quant={}
		quant['name']='psin'
		entry_tmp=OMFITmdsValue(server=self.device,treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.GEQDSK:PSIN')

		quant['data']=entry_tmp.data()
		quant['psi']=entry_tmp.dim_of(0)
		quant['R_psi']=entry_tmp.dim_of(1)*0.01
	
		if plot:
			plt.scatter(quant['psi'],quant['R'])
			plt.xlabel('psi')
			plt.ylabel('R (m)')
			
		#self.coord=quant
		
		return quant

	#quant['Bt']=MDS_obj.get_Bt(plot=True) #Bt(time,r)
	def get_Bt(self,plot=False):
		quant={}
		quant['name']='Bt'
		entry_tmp=OMFITmdsValue(server=self.device,treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.DERIVED:BTZ0')
		quant['data']=entry_tmp.data() #Bt(time,r)
		(nt,nr)=np.shape(quant['data'])
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)
		quant['R']=entry_tmp.dim_of(1)[int(0.5*nt),:] #GENERAL

		if plot:
			plt.plot(quant['R'][int(0.5*nt),:],quant['data'][int(0.5*nt),:])   #?????????
			plt.xlabel('R (m)') 
			plt.ylabel('Bt (T)')
		self.quant_list[quant['name']]=quant
		return quant

	#quant['Bref']=MDS_obj.get_Bref(plot=True) #Bref(time) Bt at axis
	def get_Bref(self,plot=False):
		quant={}
		quant['name']='Bref'
		entry_tmp=OMFITmdsValue(server=self.device,treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:BT0')
		quant['data']=entry_tmp.data() #Bt0(time)
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)

		if plot:
			plt.plot(quant['time'],quant['data'])   
			plt.xlabel('time(s)') # ?
			plt.ylabel('B0t (T)')

		self.quant_list[quant['name']]=quant
		return quant

	#dB1=get_dB1(plot=True)
	def get_dB1(self,window_for_FFT='hann',plot=False):
		entry_tmp=OMFITmdsValue(server=self.device,treename='operations',shot=self.shot_num,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNLF')
		t_L=entry_tmp.dim_of(0)
		dB_L=entry_tmp.data()

		entry_tmp=OMFITmdsValue(server=self.device,treename='operations',shot=self.shot_num,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNMF')
		t_M=entry_tmp.dim_of(0)
		dB_M=entry_tmp.data()
		#print(np.shape(dB_L))
		#print(np.shape(dB_M))


		frequency,dB_L_frequency_sq=spectral_density(dB_L,t_L,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
		frequency_kHZ=frequency/1000.
		dB_L_frequency=abs(np.sqrt(dB_L_frequency_sq))

		frequency,dB_M_frequency_sq=spectral_density(dB_M,t_M,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
		frequency_kHZ=frequency/1000.
		dB_M_frequency=abs(np.sqrt(dB_M_frequency_sq))
		quant={}
		quant['name']='dB'
		quant['data_low_freq']=dB_L_frequency
		quant['data_mid_freq']=dB_M_frequency
		quant['freq']=frequency_kHZ
		self.dB=quant
		if plot:
			plt.clf()
			plt.scatter(frequency_kHZ,dB_M_frequency+dB_L_frequency)
			plt.xlabel('kHz')
			plt.ylabel('Gauss?/sqrt(Hz)')
			plt.yscale('log')
			plt.grid()
			plt.show()
		return quant

	def get_dB1_t(self,window_for_FFT='hann',plot=False):
		entry_tmp=OMFITmdsValue(server=self.device,treename='operations',shot=self.shot_num,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNLF')
		t_L=entry_tmp.dim_of(0)
		dB_L=entry_tmp.data()

		entry_tmp=OMFITmdsValue(server=self.device,treename='operations',shot=self.shot_num,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNMF')
		t_M=entry_tmp.dim_of(0)
		dB_M=entry_tmp.data()
		#print(np.shape(dB_L))
		#print(np.shape(dB_M))

		quant_list['Lref']

		frequency,dB_L_frequency_sq=spectral_density(dB_L,t_L,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
		frequency_kHZ=frequency/1000.
		dB_L_frequency=abs(np.sqrt(dB_L_frequency_sq))

		frequency,dB_M_frequency_sq=spectral_density(dB_M,t_M,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
		frequency_kHZ=frequency/1000.
		dB_M_frequency=abs(np.sqrt(dB_M_frequency_sq))
		quant={}
		quant['name']='dB'
		quant['data_low_freq']=dB_L_frequency
		quant['data_mid_freq']=dB_M_frequency
		quant['freq']=frequency_kHZ
		self.dB=quant
		if plot:
			plt.clf()
			plt.scatter(frequency_kHZ,dB_M_frequency+dB_L_frequency)
			plt.xlabel('kHz')
			plt.ylabel('Gauss?/sqrt(Hz)')
			plt.yscale('log')
			plt.grid()
			plt.show()
		return quant


	
	def get_dn1(self,window_for_FFT='hann',plot=False):
		entry_tmp=OMFITmdsValue(server=self.device,treename='operations',shot=self.shot_num,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNLF')
		t_L=entry_tmp.dim_of(0)
		dB_L=entry_tmp.data()

		entry_tmp=OMFITmdsValue(server=self.device,treename='operations',shot=self.shot_num,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNMF')
		t_M=entry_tmp.dim_of(0)
		dB_M=entry_tmp.data()
		#print(np.shape(dB_L))
		#print(np.shape(dB_M))


		frequency,dB_L_frequency_sq=spectral_density(dB_L,t_L,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
		frequency_kHZ=frequency/1000.
		dB_L_frequency=abs(np.sqrt(dB_L_frequency_sq))

		frequency,dB_M_frequency_sq=spectral_density(dB_M,t_M,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
		frequency_kHZ=frequency/1000.
		dB_M_frequency=abs(np.sqrt(dB_M_frequency_sq))
		quant={}
		quant['name']='dB'
		quant['data_low_freq']=dB_L_frequency
		quant['data_mid_freq']=dB_M_frequency
		quant['freq']=frequency_kHZ
		self.dB=quant
		if plot:
			plt.clf()
			plt.scatter(frequency_kHZ,dB_M_frequency+dB_L_frequency)
			plt.xlabel('kHz')
			plt.ylabel('Gauss?/sqrt(Hz)')
			plt.yscale('log')
			plt.grid()
			plt.show()
		return quant


	#quant['Lref']=MDS_obj.get_Lref(plot=True) #Lref(time) minor radius
	def get_Lref(self,plot=False):
		quant={}
		quant['name']='Lref'
		entry_tmp=OMFITmdsValue(server=self.device,treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:AMINOR')
		quant['data']=entry_tmp.data() #minor_radius(time)
		#quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)
		if self.device=='d3d':
			quant['time']=quant['time']*0.001

		if plot:
			plt.plot(quant['time'],quant['data'])   
			plt.xlabel('time(s)') 
			plt.ylabel('r (m)')

		self.quant_list[quant['name']]=quant
		return quant

	def calc_ome(self,plot=False):
		kymin=n0*q0*rhoref/(Lref*x0_center)
		kyGENE =kymin * (q/q0) * np.sqrt(te_u/te_mid) * (x0_center/uni_rhot) #Add the effect of the q varying
		omMTM = kyGENE*(tprime_e+nprime_e)

	#Lref,Lref_err=calc_Lref_avg()
	def calc_Lref_avg(self,plot=False):
		Lref_t=self.quant_list['Lref']['data']
		t=self.quant_list['Lref']['time']
		Lref_avg=np.mean(Lref_t)
		Lref_err=np.std(Lref_t)/Lref_avg
		if plot:
			plt.plot(t,Lref_t)
			plt.axhline(Lref_avg)
		self.Lref=Lref_avg
		self.Lref_err=Lref_err
		return self.Lref,self.Lref_err


	
	def calc_grad(self,data,t,R):
		dprime=np.zeros(np.shape(data))

		for i in range(len(t)):

			d=data[i,:]
			d = self.smooth(d, 5) 

			#print(len(d))
			#print(len(R))
			dd = np.gradient(d,R)
			#print()
			#print(len(d))
			#print(len(dd))
			dprime[i,:] = -dd/d
		
		return dprime

	#dT_dR=calc_dT_dR(plot=False)
	def calc_dT_dR(self,plot=False):
		Te=self.quant_list['Te']['data']
		t=self.quant_list['Te']['time']
		R=self.quant_list['Te']['R']
		(nt,nr)=np.shape(Te)
		dT_dR=np.zeros((nt,nr))

		for i in range(nt):
			Te_r=self.smooth(Te[i,:],5)
			dT_dR[i,:]=-self.finite_differences(Te_r,R)

		if plot:
			Te_tmp=self.smooth(Te[int(0.5*nt),:],5)
			plt.scatter(R,Te_tmp/np.max(Te_tmp),label='Te')
			plt.plot(R,dT_dR[int(0.5*nt),:]/np.max(abs(dT_dR[int(0.5*nt),:])),label='dT/dR')
			plt.axvline(R[np.argmax(dT_dR[int(0.5*nt),:])],color='red')
			plt.legend()
			plt.xlabel('r (m)',fontsize=15)
			plt.ylabel('a.u.',fontsize=15)
		self.dT_dR=dT_dR

		return self.dT_dR

	#dn_dR=calc_dn_dR(plot=False)
	def calc_dn_dR(self,plot=False):
		ne=self.quant_list['ne']['data']
		t=self.quant_list['ne']['time']
		R=self.quant_list['ne']['R']
		(nt,nr)=np.shape(ne)
		dn_dR=np.zeros((nt,nr))

		for i in range(nt):
			ne_r=self.smooth(ne[i,:],5)
			dn_dR[i,:]=-self.finite_differences(ne_r,R)

		if plot:
			ne_tmp=self.smooth(ne[int(0.5*nt),:],5)
			plt.scatter(R,ne_tmp/np.max(ne_tmp))
			plt.plot(R,dn_dR[int(0.5*nt),:]/np.max(abs(dn_dR[int(0.5*nt),:])))

		self.dn_dR=dn_dR

		return self.dn_dR

	#R_list,a_Lne_list,ne_ped_list=calc_Ln_peak(plot=False)
	def calc_Ln_peak(self,plot=False):
		R=self.quant_list['ne']['R']
		t=self.quant_list['ne']['time']
		dn_dR=self.dn_dR
		(nt,nr)=np.shape(dn_dR)
		
		R_list=[]
		a_Lne_list=[]
		ne_ped_list=[]
		ne_ped_index_list=[]

		for i in range(nt):
			dn_dR_tmp=dn_dR[i,:]
			R_location,dn_dR_max,max_index=self.find_peak(R[3:],dn_dR_tmp[3:],plot=False)
			ne_ped=(self.quant_list['ne']['data'][i,3:])[max_index]
			#a/Lne = (a/ne) * (dne/da)
			a_Lne=(dn_dR_max/ne_ped)*self.Lref

			R_list.append(R_location)
			a_Lne_list.append(a_Lne)
			ne_ped_list.append(ne_ped)
			ne_ped_index_list.append(max_index)

		if plot:
			plt.clf()
			plt.plot(t,a_Lne_list)
			plt.xlabel('time')
			plt.ylabel('a/Lne')
			plt.show()

		self.R_ne_mid_ped=R_list
		self.a_Lne_mid_ped=a_Lne_list
		self.ne_mid_ped=ne_ped_list
		self.ne_ped_index=ne_ped_index_list
		return R_list,a_Lne_list,ne_ped_list, ne_ped_index_list

	#R_list,a_LTe_list,Te_ped_list=calc_Lt_peak(plot=False)
	def calc_Lt_peak(self,plot=False):
		R=self.quant_list['Te']['R']
		t=self.quant_list['Te']['time']
		dT_dR=self.dT_dR
		(nt,nr)=np.shape(dT_dR)
		
		R_list=[]
		a_LTe_list=[]
		Te_ped_list=[]
		Te_ped_index_list=[]

		for i in range(nt):
			dT_dR_tmp=dT_dR[i,:]
			R_location,dT_dR_max,max_index=self.find_peak(R[3:],dT_dR_tmp[3:],plot=False)
			Te_ped=(self.quant_list['Te']['data'][i,3:])[max_index]
			#a/LTe = (a/Te) * (dTe/da)
			a_LTe=(dT_dR_max/Te_ped)*self.Lref

			R_list.append(R_location)
			a_LTe_list.append(a_LTe)
			Te_ped_list.append(Te_ped)
			Te_ped_index_list.append(max_index)

		if plot:
			if 1==0:
				t_index=int(0.5*nt)
				plt.clf()
				plt.plot(R,(dT_dR/Te_ped)*self.Lref)
				plt.axvline(R_location[t_index])
				plt.xlabel('time (s)')
				plt.ylabel('a/LTe')
				plt.show()

			if 1==1:
				plt.clf()
				plt.plot(t,a_LTe_list)
				plt.xlabel('time (s)')
				plt.ylabel('a/LTe')
				plt.show()

		self.R_Te_mid_ped=R_list
		self.a_LTe_list=a_LTe_list
		self.Te_ped_list=Te_ped_list

		return R_list,a_LTe_list,Te_ped_list,Te_ped_index_list

	#coll_ei=calc_coll_ei(plot=False)
	def calc_coll_ei(self,plot=False):
		Lref=self.Lref
		#ne (/cm^3) ---> ne (10^19 /m^3)
		ne=self.quant_list['ne']['data']*10.**(6.-19.)
		#Te (keV)
		te=self.quant_list['Te']['data']

		(nt,nr)=np.shape(ne)

		#coll_ei (cs/a)
		coll_c=2.3031*10**(-5)*Lref*ne/(te)**2*(24-np.log(np.sqrt(ne*10**13)/(te*1000)))
		#coll_ei 
		coll_ei=4.*coll_c*np.sqrt(te*1000.*qref/me_SI)/Lref
		#coll_ei in kHz
		coll_ei=coll_ei/(2.*np.pi*1000.)
		self.coll_ei=coll_ei
		if plot:
			t=self.quant_list['Te']['time']
			R=self.quant_list['Te']['R']
			plt.clf()
			plt.plot(R,coll_ei[int(0.5*nt),:])
			plt.xlabel('R (m)')
			plt.ylabel('coll_ei (kHz)')
			plt.show()
		return coll_ei
	
	#beta=calc_beta(plot=False)
	def calc_beta(self,plot=False):
		quant={}
		#ne (/cm^3) ---> ne (10^19 /m^3)
		ne=self.quant_list['ne']['data']*10.**(6.-19.)
		#Te (keV)
		te=self.quant_list['Te']['data']

		time=self.quant_list['Te']['time']
		R=self.quant_list['Te']['R']
		(nt,nr)=np.shape(ne)

		Bref=self.quant_list['Bref']['data']
		Bref=np.tile(Bref, (nr, 1)).T

		beta=403.*10**(-5)*np.divide(ne*te,(Bref**2.))

		quant['name']='beta'
		quant['data']=beta
		quant['R']=R
		quant['time']=time
		self.quant_list['beta']=quant
		if plot:
			plt.clf()
			plt.plot(quant_list['beta']['data'][int(0.5*nt),:])
			plt.xlabel('R (m)')
			plt.ylabel('beta')
			plt.show()
		return beta

	def calc_ome(self,plot=False):
		quant={}
		#ne (/cm^3) ---> ne (10^19 /m^3)
		ne=self.quant_list['ne']['data']*10.**(6.-19.)
		#Te (keV)
		te=self.quant_list['Te']['data']

		time=self.quant_list['Te']['time']
		R=self.quant_list['Te']['R']
		(nt,nr)=np.shape(ne)

		Bref=self.quant_list['Bref']['data']
		Bref=np.tile(Bref, (nr, 1)).T

		beta=403.*10**(-5)*np.divide(ne*te,(Bref**2.))

		quant['name']='beta'
		quant['data']=beta
		quant['R']=R
		quant['time']=time
		self.quant_list['beta']=quant
		if plot:
			plt.clf()
			plt.plot(quant_list['beta']['data'][int(0.5*nt),:])
			plt.xlabel('R (m)')
			plt.ylabel('beta')
			plt.show()
		return beta


	#x_location,y_max=find_peak(x,y,plot=False)
	def find_peak(self,x,y,plot=False):
		max_index=np.argmax(abs(y))
		if plot:
			plt.clf()
			plt.plot(x,y)
			plt.axvline(x[max_index])
			plt.show()
		return x[max_index],y[max_index],max_index

	
	def find_dB1_peak(self,plot=False):
		quant=self.dB
		
		dB_L_frequency=quant['data_low_freq']
		dB_M_frequency=quant['data_mid_freq']
		frequency_kHZ=quant['freq']

		dB_frequency=dB_M_frequency+dB_L_frequency
		
		dB_frequency=self.smooth(dB_frequency, box_pts=50)
		peaks,_=signal.find_peaks(dB_frequency,prominence=abs(np.mean(dB_frequency)-3*np.std(dB_frequency)))
		#peaks,_=signal.find_peaks(dB_frequency, prominence=1)
		#print(peaks)
		f_list=frequency_kHZ[peaks]
		dB_list=dB_frequency[peaks]
		if plot:
			plt.plot(frequency_kHZ,dB_frequency,alpha=0.3)
			plt.scatter(f_list,dB_list,s=10)
			plt.xlabel('kHz')
			plt.ylabel('Gauss?/sqrt(Hz)')
			#plt.xscale('log')
		return f_list,dB_list


	#H_mode_or_not=self.judge_H_mode(plot=True)
	def judge_H_mode(self,plot=True):
		#1: H-mode, #0: L-mode, #-1: do not use
		dn_dR=self.dn_dR
		dT_dR=self.dT_dR
		R=self.quant_list['ne']['R']
		t=self.quant_list['ne']['time']
		(nt,nr)=np.shape(dT_dR)
		confinment_mode_list=np.zeros(nt,dtype=int)
		for i in range(nt):
			confinment_mode=0
			confinment_mode_list[i]=confinment_mode
		if plot:
			plt.clf()
			#plt.plot(R,dn_dR[int(0.5*nt),:]/np.max(abs(dn_dR[int(0.5*nt),:])),label='n')
			#plt.plot(R,dT_dR[int(0.5*nt),:]/np.max(abs(dT_dR[int(0.5*nt),:])),label='T')
			#plt.plot(R,dT_dR[int(0.5*nt),:]*dn_dR[int(0.5*nt),:]/np.max(abs(dT_dR[int(0.5*nt),:]*dn_dR[int(0.5*nt),:])),label='p')
			
			plt.plot(R,dn_dR[int(0.5*nt),:],label='n')
			#plt.plot(R,dT_dR[int(0.5*nt),:],label='T')
			#plt.plot(R,dT_dR[int(0.5*nt),:]*dn_dR[int(0.5*nt),:],label='p')
			
			plt.legend()
			plt.show()
		self.confinment_mode_list=confinment_mode_list
		return confinment_mode_list

	def interp_all_quant(self,interp_factor=2,plot=False):

		nt_max=0
		for key in self.keys:
			try:
				time=self.quant_list[key]['time']
			except:
				continue

			tmp=len(time)
			if nt_max<tmp:
				nt_max=tmp
		

 		#GENERAL needed
		nR_max=0
		for key in self.keys:
			try:
				R=self.quant_list[key]['R']
			except:
				continue
			
			tmp=len(R)
			if nR_max<tmp:
				nR_max=tmp
			
		t_u=np.linspace(self.t_min,self.t_max,int(nt_max*interp_factor))
		R_u=np.linspace(self.R_min,self.R_max,int(nR_max*interp_factor))

		self.coord['time']=t_u
		self.coord['R']=R_u
		
		#interperlation
		for key in self.keys_t_r:
			data=self.quant_list[key]['data']
			t=self.quant_list[key]['time']
			R=self.quant_list[key]['R']
			#print(key)
			#print(np.shape(data))
			#print(len(t))
			#print(len(R))
			#print(len(t_u))
			#print(len(R_u))
			data_u=self.interp_2D(data,t,R,t_u,R_u,plot=plot)
			
			self.quant_list[key]['data']=data_u
			self.quant_list[key]['time']=t_u
			self.quant_list[key]['R']=R_u
			

		
		
		for key in self.keys_t:
			data=self.quant_list[key]['data']
			t=self.quant_list[key]['time']
			data_u=np.interp(t_u,t,data)
			
			self.quant_list[key]['data']=data_u
			self.quant_list[key]['time']=t_u

		return self.quant_list

	
	def chose_time(self,plot=False):
		Lref=self.quant_list['Lref']['data']
		time=self.quant_list['Lref']['time']
		t_min_index=np.argmin(abs(time-self.t_min))
		Lref=Lref[t_min_index:]
		time=time[t_min_index:]

		#def moving_avg(x, w):
		#	return np.convolve(x, np.ones(w), 'valid') / w
		def moving_avg(x, w):
			return np.array([np.mean(x[i:i+w]) for i in range(len(x)-w-1)])
		
		def moving_std(x, w):
			return np.array([np.std(x[i:i+w]) for i in range(len(x)-w-1)])
		MACD_5=moving_avg(Lref, 5)
		time_5=moving_avg(time, 5)
		M_std_5=moving_std(Lref, 5)
		quant=M_std_5/MACD_5
		crit=0.02
		#print(quant)
		try:
			index=next(x[0] for x in enumerate(quant) if x[1] > crit)
			t_max=time_5[index]
			self.t_max=t_max
		except:
			t_max=self.t_max

		if plot:
			plt.plot(time,Lref,label='Data')
			plt.plot(time_5,MACD_5,label='Avgerage')
			plt.plot(time_5,M_std_5,label='Standard Deviation')
			plt.xlabel('time (s)',fontsize=15)
			plt.ylabel('a (m)',fontsize=15)
			plt.axvline(t_max,color='red')
			plt.legend()

	def cut_R(self,plot=False):
		if plot:
			plt.plot(self.quant_list['Te']['R'],\
					self.quant_list['Te']['data'][10,:])
			plt.xlabel('r (m)',fontsize=15)
			plt.ylabel(r'$T_e (keV)$',fontsize=15)
			plt.axvline(self.R_min,color='red')
			plt.axvline(self.R_max,color='red')
			#plt.legend()
		for key in self.keys_t_r:
			data=self.quant_list[key]['data']
			R=self.quant_list[key]['R']
			
			R_min_index=np.argmin(abs(R-self.R_min))
			R_max_index=np.argmin(abs(R-self.R_max))

			#print(R_min_index)
			#print(R_max_index)

			self.quant_list[key]['data']=data[:,R_min_index:R_max_index+1]
			self.quant_list[key]['R']=R[R_min_index:R_max_index+1]
			


	def cut_t(self):
		for key in self.keys_t_r:
			data=self.quant_list[key]['data']
			t=self.quant_list[key]['time']
			
			t_min_index=np.argmin(abs(t-self.t_min))
			t_max_index=np.argmin(abs(t-self.t_max))

			self.quant_list[key]['data']=data[t_min_index:t_max_index+1,:]
			self.quant_list[key]['time']=t[t_min_index:t_max_index+1]
		
		for key in self.keys_t:
			data=self.quant_list[key]['data']
			t=self.quant_list[key]['time']
			
			t_min_index=np.argmin(abs(t-self.t_min))
			t_max_index=np.argmin(abs(t-self.t_max))

			self.quant_list[key]['data']=data[t_min_index:t_max_index+1]
			self.quant_list[key]['time']=t[t_min_index:t_max_index+1]
	

	def get_all_quant(self,plot=False):
		if self.device=='nstx':
			quant=self.get_Bt(plot)		#Bt(time,r)
			#quant=self.get_q0psi(plot)  	#q(time,r)
			quant=self.get_Te(plot) 	#Te(r,time)
			quant=self.get_ne(plot) 	#ne(time,r)
			quant=self.get_Bref(plot) 	#Bref(time) Bt at axis
			quant=self.get_Lref(plot) 	#Lref(time) minor radius
			#quant=self.get_psi_R(plot)
		elif self.device=='d3d':
			quant=self.get_Te(plot) 	#Te(r,time)
			quant=self.get_ne(plot) 	#ne(time,r)
			quant=self.get_Bref(plot) 	#Bref(time) Bt at axis
			quant=self.get_Lref(plot) 	#Lref(time) minor radius
		
		self.keys=self.quant_list.keys()
		return self.quant_list


	def Auto_scan(self,device='nstx',shot_num=132588,plot=False):
		self.set_shot_num(shot_num=shot_num)
		self.set_device(device=device)
		#get all the quantities
		quant_list=self.get_all_quant(plot=False)
		self.chose_time(plot=False)
		#cut the R from psi=0 to psi=0.99 ish
		self.cut_R(plot=False)
		#cut time to from 0.3s to 0.8s
		self.cut_t()
		#self.test_save_worthy()
		#interpolation all the quantities to uniform grid
		quant_list=self.interp_all_quant(interp_factor=1.2,plot=False)

		#calculate the average Lref
		Lref,Lref_err=self.calc_Lref_avg(plot=False)
		#calculate the dn_dR=dne/dR
		dn_dR=self.calc_dn_dR(plot=False)
		#calculate the dT_dR=dTe/dR
		dT_dR=self.calc_dT_dR(plot=True)

		#calculate the R, a/Lne, ne in dne/dR peak
		R_ne_list,a_Lne_list,ne_ped_list,ne_ped_index_list=self.calc_Ln_peak(plot=False)
		#calculate the R, a/LTe, Te in dTe/dR peak
		R_Te_list,a_LTe_list,Te_ped_list,Te_ped_index_list=self.calc_Lt_peak(plot=False)
		
		#calculate if that is H-mode
		#H_mode_or_not=self.judge_H_mode(plot=False)

		#calculate the collision
		coll_ei=self.calc_coll_ei(plot=False)
		coll_ei_ped_value=[coll_ei[i,ne_ped_index_list[i]] for i in range(len(ne_ped_index_list))]
		beta=self.calc_beta(plot=False)
		beta_ped_value=[beta[i,ne_ped_index_list[i]] for i in range(len(ne_ped_index_list))]
		t=self.quant_list['Te']['time']

		d={}
		d['shot_num']=[shot_num]*len(t)
		d['time']=t
		d['Lref']=self.quant_list['Lref']['data']
		d['beta']=beta_ped_value
		d['a_Lne']=a_Lne_list
		d['a_LTe']=a_LTe_list
		d['coll_ei']=coll_ei_ped_value
		d['R_ne_list']=R_ne_list
		d['R_Te_list']=R_Te_list

		df=pd.DataFrame(d)

		if plot:
			if 1==0:
				plt.clf()
				plt.plot(t,R_ne_list,label=r'R($ne_{mid ped}$)')
				plt.plot(t,R_Te_list,label=r'R($Te_{mid ped}$)')
				plt.legend()
				plt.xlabel('t (s)')
				plt.ylabel('R (m)')
				plt.show()
			if 1==0:
				plt.clf()
				plt.plot(t,coll_ei_ped_value)
				plt.xlabel('t (s)')
				plt.ylabel('coll_ei_ped_value (kHz)')
				plt.show()
			if 1==0:
				plt.clf()
				plt.plot(t,beta_ped_value)
				plt.xlabel('t (s)')
				plt.ylabel('beta')
				plt.show()
		return df

	
if test:
	shot_num=132588 #132588 (H-mode), 141716 (L-mode)
	device='nstx' #'nstx', 'd3d'
	MDS_obj=MDS_obj(device=device,shot_num=shot_num)
	#dB=MDS_obj.get_dB1(plot=False)
	#f_list,dB_list=MDS_obj.find_dB1_peak(plot=True)
	#print(f_list)
	#print(dB_list)
	df=MDS_obj.Auto_scan(device=device,shot_num=shot_num,plot=False)
	print(df)
	if 1==0:
		plt.plot(df.Lref)
		plt.plot(df.beta)
