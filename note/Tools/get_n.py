#working
entry_tmp=OMFITmds(server='nstx',treename='activespec',shot=132588)

a=entry_tmp['CHERS']['ANALYSIS']['CT1']

ne=a['DEN'].data() #ne(time,r)
time=a['DEN'].dim_of(1)
r=a['DEN'].dim_of(0)

info=a['DEN'].xarray()
print(info)
print(np.shape(ne))
print(np.shape(r))
print(np.shape(time))

plt.plot(r,ne[30,:])   #?????????
plt.xlabel('R (cm)')
plt.ylabel('ne (cm-3)')