#working
entry_tmp=OMFITmdsValue(server='nstx',treename='activespec',shot=132588,TDI='\\ACTIVESPEC::TOP.MPTS.OUTPUT_DATA.BEST:FIT_TE	')

te=entry_tmp.data() #te(time,r)
info=entry_tmp.xarray() 

r=entry_tmp.dim_of(1)
time=entry_tmp.dim_of(0)
units=entry_tmp.units()


print(info)
print(np.shape(te))
print(np.shape(r))
print(np.shape(time))

plt.plot(r,te[:,10])
plt.xlim(90,165)
plt.xlabel('R (cm)')
plt.ylabel('Te (keV)')
