#working

#https://omfit.io/tutorial/SCRIPTS_ADAS_demo_read_adf12.html

entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:BETAN')
#time=entry_tmp.dim_of(0)

beta=entry_tmp.data() #R(time,r)
info=entry_tmp.xarray()


time=entry_tmp.dim_of(0)

print(np.shape(beta))
print(np.shape(time))
#psi=entry_tmp.dim_of(1)

plt.plot(time,beta)

