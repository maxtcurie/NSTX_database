entry_tmp=OMFITmds(server='nstx',treename='activespec',shot=132588)

a=entry_tmp['CHERS']['ANALYSIS']['CT1']
#print(a.keys())
#print(a['DEN'])

ne=a['DEN'].data() #ne(time,r)
time=a['DEN'].dim_of(0)
r=a['DEN'].dim_of(1)

info=a['DEN'].xarray()
#print(type(ne))
print(np.shape(ne))
print(np.shape(dim0))
print(np.shape(dim1))

plt.plot(r,ne[30,:])