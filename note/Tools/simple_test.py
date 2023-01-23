
entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::QPSI')
x=entry_tmp.dim_of(0)
q0=entry_tmp.data()
print(np.shape(q0))

entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::TOP.RESULTS.AEQDSK:PSI0')
x=entry_tmp.dim_of(0)
psi=entry_tmp.data()


entry_tmp.data()

plt.clf()
plt.scatter(x,psi)
plt.xlabel('psi')
plt.ylabel('q0')
plt.show()

'''

entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::TOP.RESULTS.GEQDSK:QPSI')
q0_2D=entry_tmp.data()

x=np.arange(len(psi))

print(np.shape(psi))

plt.clf()
plt.scatter(x,psi)
plt.show()

plt.clf()
plt.scatter(x,q0)
plt.show()


'''