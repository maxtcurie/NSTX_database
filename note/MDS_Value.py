entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=132588,TDI='\\EFIT01::q0')

x=entry_tmp.dim_of(0)
y=entry_tmp.data()
unit=entry_tmp.units()
info=entry_tmp.xarray()

print(info)

