Signal= OMFITmdsValue(server=device,shot=shot,treename=treename,TDI=f'\\{treename}::{node_name}')
device=root['SETTINGS']['EXPERIMENT']['device'],
treename=root['SETTINGS']['EXPERIMENT'].get('treename_plotDerive','EFIT01'),
node_name =root['SETTINGS']['EXPERIMENT'].get('node_name_plotDerive', 'wmhd'),
shots = [root['SETTINGS']['EXPERIMENT']['shot']],
times_start = root['SETTINGS']['EXPERIMENT'].get('times_start_plotDerive',[0.38]),
times_end = root['SETTINGS']['EXPERIMENT'].get('times_end_plotDerive',[0.42])
