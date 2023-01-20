entry_tmp=OMFITmds(server='nstx',treename='activespec',shot=132588)

print(entry_tmp.keys())
print(entry_tmp['MPTS'].keys())
print(entry_tmp['MPTS']['OUTPUT_DATA'].keys())
print(entry_tmp['MPTS']['OUTPUT_DATA']['BEST'].keys())

entry_tmp=entry_tmp['MPTS']['OUTPUT_DATA']['TS2']['FIT_TE']

te=entry_tmp.data() #ne(time,r)



time=entry_tmp.dim_of(1)
r=entry_tmp.dim_of(0)
info=entry_tmp.xarray()
units=entry_tmp.units()
#print(info)
#print(np.shape(te))
#print(unit)

#plt.plot(te[:,30])
#plt.xlim(0,0.95)