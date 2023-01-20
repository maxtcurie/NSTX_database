entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::QPSI')
time=entry_tmp.dim_of(0)

q0=entry_tmp.data()

print(np.shape(x))
print(np.shape(q0))

x=np.arange(65)

plt.clf()
plt.scatter(x,q0[64,:])
plt.xlabel('psi')
plt.ylabel('q0')
plt.show()