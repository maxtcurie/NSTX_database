import numpy as np 

working_show_num_list=[]
base_name='shot_survey'
shot_num_list=np.arange(163300,163310,dtype=int)
shot_num_list=[132588]

OMFIT[base_name]={}

#https://nstx.pppl.gov/nstx/Software/FAQ/signallabels.html

TDI_list=['\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNLF',\
            '\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNMF',\
            '\\EFIT01::TOP.RESULTS.AEQDSK:PSI0',\
            '\\EFIT01::q0',\
            '\\EFIT01::TOP.RESULTS.AEQDSK:DENSITY'\
        ]

name_list=['dB_low_f','df_mid_f','psi','q0','density']

TDI_list=TDI_list[2:]
name_list=name_list[2:]

def read_data_from_db(shot_num):
    db_dic_item={}
    for i in range(len(TDI_list)):   
        entry_tmp=OMFITmdsValue(server='nstx',treename='efit01',shot=shot_num,TDI=TDI_list[i])
        print(entry_tmp)
        db_dic_item[name_list[i]]=entry_tmp
    return db_dic_item

db_dic={}
for shot_num in shot_num_list:
    #try:
    #except:
    #    pass
    if 1==1:
        a=read_data_from_db(shot_num)
        db_dic[str(shot_num)]=a
        working_show_num_list.append(shot_num)
    

print(working_show_num_list)
print(db_dic)

OMFIT[base_name]=db_dic

name='dB_low_f'
plt.clf()
for shot_num in working_show_num_list:
    #print(OMFIT[base_name][str(shot_num)].keys())
    print(OMFIT[base_name][str(shot_num)][name])
    plt.plot(\
        #OMFIT[base_name][str(shot_num)]['psi'],\
        OMFIT[base_name][str(shot_num)]['q0'])
plt.show()


#? available_efits_from_mds