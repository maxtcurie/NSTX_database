entry_tmp=OMFITmds(server='nstx',treename='activespec',shot=132588)

a=entry_tmp['CHERS']['ANALYSIS']['CT1']

ne=a['DEN'].data() #ne(time,r)
r=a['DEN'].dim_of(1)
time=a['DEN'].dim_of(0)

info=a['DEN'].xarray()
print(info)
print(np.shape(ne))
print(np.shape(r))
print(np.shape(time))

plt.plot(time,ne[30,:])   #?????????