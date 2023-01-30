import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

path='./Output/' #end with '/'
dir_list = os.listdir(path)
print(dir_list)


shot_list=[]
for i in dir_list:
	if i[0]=='1' or i[0]=='2':
		pass
	else:
		continue
	print(i[:6])
	shot_num=int(i[:6])
	shot_list.append(shot_num)

print(shot_list)
