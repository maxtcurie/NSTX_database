import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


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

	def __init__(self,shot_num=132588):
		self.shot_num=shot_num
		self.quant_list={}
		self.coord={}
		self.keys_t_r=['Bt', 'Te', 'ne']
		self.keys_psi=['q0psi']
		self.keys_t=['Bref', 'Lref']
		self.R_min=0.90
		self.R_max=1.45
		self.t_min=0.3
		self.t_max=0.8

	def quant_report(quant):
		print(quant['info'])
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

	def get_Zeff(self,plot=False):
		pass
		#\ACTIVESPEC::TOP.CHERS.ANALYSIS.CT1:ZEFF
	
	#quant['ne']=MDS_obj.get_ne(plot=True) #ne(time,r)
	def get_ne(self,plot=False):
		quant={}
		quant['name']='ne'
		entry_tmp=OMFITmdsValue(server='nstx',treename='activespec',shot=self.shot_num,TDI='\\ACTIVESPEC::TOP.CHERS.ANALYSIS.CT1:DEN')
		quant['data']=entry_tmp.data() #ne(time,r)
		(nt,nr)=np.shape(quant['data'])
		quant['info']=entry_tmp.xarray()
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
		entry_tmp=OMFITmdsValue(server='nstx',treename='activespec',shot=132588,TDI='\\ACTIVESPEC::TOP.MPTS.OUTPUT_DATA.BEST:FIT_TE	')
		quant['data']=entry_tmp.data() #te(time,r)
		quant['data']=quant['data'].T
		(nt,nr)=np.shape(quant['data'])
		quant['info']=entry_tmp.xarray()
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
		entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=self.shot_num,TDI='\\EFIT01::QPSI')
		quant['data']=entry_tmp.data() #q(time,r)
		(nt,nr)=np.shape(quant['data'])
		quant['info']=entry_tmp.xarray()
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
		entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:Q0')
		quant['data']=entry_tmp.data() #q(time)
		#(nt,nr)=np.shape(quant['data'])
		print(np.shape(quant['data']))
		
		quant['info']=entry_tmp.xarray()
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
		entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.GEQDSK:PSIN')

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
		entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.DERIVED:BTZ0')
		quant['data']=entry_tmp.data() #Bt(time,r)
		(nt,nr)=np.shape(quant['data'])
		quant['info']=entry_tmp.xarray()
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
		entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:BT0')
		quant['data']=entry_tmp.data() #Bt0(time)
		quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)

		if plot:
			plt.plot(quant['time'],quant['data'])   
			plt.xlabel('time(s)') # ?
			plt.ylabel('B0t (T)')

		self.quant_list[quant['name']]=quant
		return quant

	#quant['Lref']=MDS_obj.get_Lref(plot=True) #Lref(time) minor radius
	def get_Lref(self,plot=False):
		quant={}
		quant['name']='Lref'
		entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=self.shot_num,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:AMINOR')
		quant['data']=entry_tmp.data() #minor_radius(time)
		quant['info']=entry_tmp.xarray()
		quant['time']=entry_tmp.dim_of(0)

		if plot:
			plt.plot(quant['time'],quant['data'])   
			plt.xlabel('time(s)') 
			plt.ylabel('r (m)')

		self.quant_list[quant['name']]=quant
		return quant

	def get_all_quant(self,plot=False):
		quant=self.get_Bt(plot)		#Bt(time,r)
		quant=self.get_q0psi(plot)  	#q(time,r)
		quant=self.get_Te(plot) 	#Te(r,time)
		quant=self.get_ne(plot) 	#ne(time,r)
		quant=self.get_Bref(plot) 	#Bref(time) Bt at axis
		quant=self.get_Lref(plot) 	#Lref(time) minor radius
		quant=self.get_psi_R(plot)
		
		self.keys=self.quant_list.keys()
		return self.quant_list

	def interp_all_quant(self,inter_factor=2,plot=False):

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
			
		t_u=np.linspace(self.t_min,self.t_max,int(nt_max*inter_factor))
		R_u=np.linspace(self.R_min,self.R_max,int(nR_max*inter_factor))

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

	def smooth_all_quant(self,box_pts=3,plot=False):
		pass

	def chose_time(self,box_pts=5,plot=False):
		t_min=0
		t_max=1
		pass

	def cut_R(self):
		for key in self.keys_t_r:
			data=self.quant_list[key]['data']
			R=self.quant_list[key]['R']
			
			R_min_index=np.argmin(abs(R-self.R_min))
			R_max_index=np.argmin(abs(R-self.R_max))

			#print(R_min_index)
			#print(R_max_index)

			self.quant_list[key]['data']=data[:,R_min_index:R_max_index+1]
			self.quant_list[key]['R']=R[R_min_index:R_max_index+1]

		#plt.plot(self.quant_list[key]['R'],\
		#		self.quant_list[key]['data'][10,:])


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


	def calc_beta(self):
		beta=403.*10**(-5)*ne*te/Bref**2.

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
			plt.scatter(R,Te_tmp/np.max(Te_tmp))
			plt.plot(R,dT_dR[int(0.5*nt),:]/np.max(abs(dT_dR[int(0.5*nt),:])))

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

	#R_list,a_Lne_list,ne_ped_list=calf_Ln_peak(plot=False)
	def calf_Ln_peak(self,plot=False):
		R=self.quant_list['ne']['R']
		t=self.quant_list['ne']['time']
		dn_dR=self.dn_dR
		(nt,nr)=np.shape(dn_dR)
		
		R_list=[]
		a_Lne_list=[]
		ne_ped_list=[]

		for i in range(nt):
			dn_dR_tmp=dn_dR[i,:]
			R_location,dn_dR_max,max_index=self.find_peak(R[3:],dn_dR_tmp[3:],plot=False)
			ne_ped=(self.quant_list['ne']['data'][i,3:])[max_index]
			#a/Lne = (a/ne) * (dne/da)
			a_Lne=(dn_dR_max/ne_ped)*self.Lref

			R_list.append(R_location)
			a_Lne_list.append(a_Lne)
			ne_ped_list.append(ne_ped)

		if plot:
			plt.clf()
			plt.plot(t,a_Lne_list)
			plt.xlabel('time')
			plt.ylabel('a/Lne')
			plt.show()

		self.R_ne_mid_ped=R_list
		self.a_Lne_mid_ped=a_Lne_list
		self.ne_mid_ped=ne_ped_list

		return R_list,a_Lne_list,ne_ped_list

	#R_list,a_LTe_list,Te_ped_list=calf_Lt_peak(plot=False)
	def calf_Lt_peak(self,plot=False):
		R=self.quant_list['Te']['R']
		t=self.quant_list['Te']['time']
		dT_dR=self.dT_dR
		(nt,nr)=np.shape(dT_dR)
		
		R_list=[]
		a_LTe_list=[]
		Te_ped_list=[]

		for i in range(nt):
			dT_dR_tmp=dT_dR[i,:]
			R_location,dT_dR_max,max_index=self.find_peak(R[3:],dT_dR_tmp[3:],plot=False)
			Te_ped=(self.quant_list['Te']['data'][i,3:])[max_index]
			#a/LTe = (a/Te) * (dTe/da)
			a_LTe=(dT_dR_max/Te_ped)*self.Lref

			R_list.append(R_location)
			a_LTe_list.append(a_LTe)
			Te_ped_list.append(Te_ped)

		if plot:
			plt.clf()
			plt.plot(t,a_LTe_list)
			plt.xlabel('time')
			plt.ylabel('a/LTe')
			plt.show()

		self.R_Te_mid_ped=R_list
		self.a_LTe_list=a_LTe_list
		self.Te_ped_list=Te_ped_list

		return R_list,a_LTe_list,Te_ped_list


	#x_location,y_max=find_peak(x,y,plot=False)
	def find_peak(self,x,y,plot=False):
		max_index=np.argmax(abs(y))
		if plot:
			plt.clf()
			plt.plot(x,y)
			plt.axvline(x[max_index])
			plt.show()
		return x[max_index],y[max_index],max_index

	def calc_grad(self,data,t,R):
		dprime=np.zeros(np.shape(data))

		for i in range(len(t)):

			d=data[i,:]
			d = self.smooth(d, 5) 

			print(len(d))
			print(len(R))
			dd = np.gradient(d,R)
			print()
			print(len(d))
			print(len(dd))
			dprime[i,:] = -dd/d
		
		return dprime

	def calc_nu(self):
		coll_c=2.3031*10**(-5)*Lref*ne/(te)**2*(24-np.log(np.sqrt(ne*10**13)/(te*1000)))
		coll_ei=4.*coll_c*np.sqrt(te*1000.*qref/me_SI)/Lref



	def Auto_scan(self,shot_num=132588):
		self.set_shot_num(shot_num=shot_num)
		#get all the quantities
		quant_list=self.get_all_quant(plot=False)
		#cut the R from psi=0 to psi=0.99 ish
		self.cut_R()
		#cut time to from 0.3s to 0.8s
		self.cut_t()
		#interpolation all the quantities to uniform grid
		quant_list=self.interp_all_quant(inter_factor=1.2,plot=False)
		#calculate the average Lref
		Lref,Lref_err=self.calc_Lref_avg(plot=False)
		#calculate the dn_dR=dne/dR
		dn_dR=self.calc_dn_dR(plot=False)
		#calculate the dT_dR=dTe/dR
		dT_dR=self.calc_dT_dR(plot=False)
		R_ne_list,a_Lne_list,ne_ped_list=self.calf_Ln_peak(plot=False)
		R_Te_list,a_LTe_list,Te_ped_list=self.calf_Lt_peak(plot=True)

	
if test:
	shot_num=132588
	MDS_obj=MDS_obj()
	MDS_obj.Auto_scan(shot_num=shot_num)