#working

entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:AMINOR')
time=entry_tmp.dim_of(0)

r=entry_tmp.data() #q(time,r)
info=entry_tmp.xarray()

#print(r)
print(np.shape(r))

print(np.shape(time))

plt.plot(time,r)
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