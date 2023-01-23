#working
entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::TOP.RESULTS.DERIVED:BTZ0')

Bt=entry_tmp.data() #ne(time,r)
info=entry_tmp.xarray()

time=entry_tmp.dim_of(1)
r=entry_tmp.dim_of(0)

print(info)
print(np.shape(Bt))
print(np.shape(r))
print(np.shape(time))
print(time[1,:])


plt.plot(time[1,:],Bt[60,:])   #?????????
plt.xlabel('R (m)') # ?
plt.ylabel('Bt (T)')
