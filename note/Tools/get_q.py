#working

entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::QPSI')
time=entry_tmp.dim_of(0)

q0=entry_tmp.data() #q(time,r)

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