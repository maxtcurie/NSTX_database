
#doc: https://omfit.io/_modules/omfit_classes/omfit_rdb.html


shot_min=163300
shot_max=163310



#OMFITrdb['select * from summaries where shot>'\
#        +str(shot_min)+'and shot<'+str(shot_max)]


set_rdb_password('nstxrdb', username='mcurie', password='C!!urie1312467')

OMFIT["rdb"]=OMFITrdb("select * from sys.databases ", \
                        server='NSTX',db='nstxlogs')


#a=OMFIT["rdb"].get_db_structure()
#available_efits_from_rdb
#print(OMFIT["rdb"])
#print(OMFIT["rdb"].keys())
