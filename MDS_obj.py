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

	def __init__(self,shot_num):
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

	def calc_Lt(self,plot=False):
		Te=self.quant_list['Te']['data']
		t=self.quant_list['Te']['time']
		R=self.quant_list['Te']['R']
		R_2=(R[:-1]+R[1:])*0.5
		(nt,nr)=np.shape(ne)
		dTe_2=np.zeros((nt,nr-1))
		for i in range(nt):
			ne_r=self.smooth(ne[i,:],5)
			dTe_2[i,:]=(ne_r[:-1]-ne_r[1:])/(R[:-1]-R[1:])
		R_2=R_2[2:]
		dTe_2=dTe_2[:,2:]
		if plot:
			#plt.plot(R,ne)
			plt.scatter(R,ne[int(0.5*nt),:]/np.max(ne[int(0.5*nt),:]))
			plt.plot(R_2,dne_2[int(0.5*nt),:]/np.max(dne_2[int(0.5*nt),:]))
			#plt.ylim(0,2)
		self.Lne=dne_2
		self.R_2=R_2
		return self.R_2,self.Lne


	def calc_Ln(self,plot=False):
		ne=self.quant_list['ne']['data']
		t=self.quant_list['ne']['time']
		R=self.quant_list['ne']['R']
		R_2=(R[:-1]+R[1:])*0.5
		(nt,nr)=np.shape(ne)
		dne_2=np.zeros((nt,nr-1))
		for i in range(nt):
			ne_r=self.smooth(ne[i,:],5)
			dne_2[i,:]=(ne_r[:-1]-ne_r[1:])/(R[:-1]-R[1:])
		R_2=R_2[2:]
		dne_2=dne_2[:,2:]
		if plot:
			#plt.plot(R,ne)
			plt.scatter(R,ne[int(0.5*nt),:]/np.max(ne[int(0.5*nt),:]))
			plt.plot(R_2,dne_2[int(0.5*nt),:]/np.max(dne_2[int(0.5*nt),:]))
			#plt.ylim(0,2)
		self.Lne=dne_2
		self.R_2=R_2
		return self.R_2,self.Lne


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



if test:
	MDS_obj=MDS_obj(shot_num=132588)
	quant_list=MDS_obj.get_all_quant(plot=False)

	Lref,Lref_err=MDS_obj.calc_Lref_avg(plot=True)

	MDS_obj.cut_R()
	MDS_obj.cut_t()
	quant_list=MDS_obj.interp_all_quant(inter_factor=1.2,plot=False)
	
	R,Lne_t_r=MDS_obj.calc_Ln(plot=False)

	Lref,Lref_err=MDS_obj.calc_Lref_avg(plot=True)
	print(Lref,Lref_err)


	