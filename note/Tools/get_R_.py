#working

entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::TOP.RESULTS.GEQDSK:R')
#time=entry_tmp.dim_of(0)

R=entry_tmp.data() #R(time,r)
#print(entry_tmp.xarray())

time=entry_tmp.dim_of(0)
psi=entry_tmp.dim_of(1)

plt.plot(time,psi)

#print(R)
print(np.shape(R))
#print(psi)
print(np.shape(psi))

#plt.plot(np.arange(len(R[60,:])),R[60,:])
#plt.plot(psi,R[60,:])
#plt.ylim(0,2)
#print(info)

#psi=entry_tmp.dim_of(0) #??????

#print(np.shape(psi))

#plt.plot(psi,R)
'''
info=entry_tmp.xarray()
r_t=entry_tmp.dim_of(1) #r(time)=r[time,r]

#print(dim1)
#print(np.shape(q0))
#print(np.shape(dim1))

plt.clf()
plt.scatter(dim1[64,:],q0[64,:])
plt.xlabel('psi')
plt.ylabel('q0')
plt.show()

'''