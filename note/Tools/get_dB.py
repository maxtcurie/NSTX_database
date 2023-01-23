#working
window_for_FFT='hann'

from scipy import signal
#from FFT_general import spectral_density

def sort_x_f(x_unsort,f_unsort): 
   
    arr_unsort=[x_unsort,f_unsort]
    f_x_unsort=tuple(map(tuple, np.transpose(arr_unsort)))
      
    f_x_sort=sorted(f_x_unsort, key=lambda f_x_unsort: f_x_unsort[0])
    f_x_sort=np.array(f_x_sort)
    f_x_sort=np.transpose(f_x_sort)

    x_sort=f_x_sort[0,:]
    f_sort=f_x_sort[1,:]
    x_sort=x_sort.astype(type(x_unsort[0]))
    f_sort=f_sort.astype(type(f_unsort[0]))

    #print('2')
    #print('type(x_sort[0])'+str(type(x_sort[0])))
    #print('type(f_sort[0])'+str(type(f_sort[0])))

    return x_sort,f_sort

#About welch method: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.periodogram.html
#'boxcar' Also known as a rectangular window or Dirichlet window, this is equivalent to no window at all.
#window types: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window
#intruction video of Welch's method: https://youtu.be/YK1F0-3VvQI
def spectral_density(function,time,percent=0.5,window_for_FFT='hann',plot=False):
    time=np.array(time)
    
    time,function=sort_x_f(time,function)
    dt=time[1:]-time[:-1]
    dt_min=np.mean(dt)

    if abs(np.std(dt))>=np.min(dt)*0.01:
        print('time step is NOT uniform. interperlating')
        uni_time = np.linspace(min(time),max(time),int(abs((max(time)-min(time))/dt_min)*1.5)) #uniform time
        uni_function = np.interp(uni_time,time,function)
    else:
        uni_time=time
        uni_function=function


    fs=1./np.mean(abs(uni_time[1:]-uni_time[:-1]))
    print('avg_dt='+str(np.mean(abs(uni_time[1:]-uni_time[:-1]))))
    print('std_dt='+str(np.std(abs(uni_time[1:]-uni_time[:-1]))))
    #f, Pxx_den = signal.welch(uni_function, fs, nperseg=len(uni_function), window=window_for_FFT) #, scaling='spectrum')

    f, Pxx_den = signal.welch(uni_function, fs, nperseg=int(percent*len(uni_function)), window=window_for_FFT,return_onesided=False, scaling='density')
    #f, Pxx_den = signal.periodogram(uni_function, fs)

    #Sort frequency to monotonic increase
    f, Pxx_den=sort_x_f(f, Pxx_den)

    if plot==True:
        plt.plot(f, Pxx_den,label='Pxx_den')
        plt.plot(f, np.sqrt(Pxx_den),label='sqrt(Pxx_den)')
        #plt.semilogy(f, Pxx_den)
        #plt.ylim([1e-7, 1e2])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [V**2/Hz]')
        plt.grid()
        plt.legend()
        plt.show()
    return f, Pxx_den


entry_tmp=OMFITmdsValue(server='nstx',treename='operations',shot=132588,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNLF')
t_L=entry_tmp.dim_of(0)
dB_L=entry_tmp.data()

entry_tmp=OMFITmdsValue(server='nstx',treename='operations',shot=132588,TDI='\\OPERATIONs::TOP.MAGNETICS.MIRNOV.B.BMHDODDNMF')
t_M=entry_tmp.dim_of(0)
dB_M=entry_tmp.data()

print(np.shape(dB_L))
print(np.shape(dB_M))


frequency,dB_L_frequency_sq=spectral_density(dB_L,t_L,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
frequency_kHZ=frequency/1000.
dB_L_frequency=abs(np.sqrt(dB_L_frequency_sq))

frequency,dB_M_frequency_sq=spectral_density(dB_M,t_M,percent=0.1,window_for_FFT=window_for_FFT,plot=False)
frequency_kHZ=frequency/1000.
dB_M_frequency=abs(np.sqrt(dB_M_frequency_sq))

plt.clf()
plt.scatter(t_L,dB_L)
plt.xlabel('s')
plt.xlabel('Gauss?')
plt.show()

plt.clf()
plt.scatter(t_M,dB_M)
plt.xlabel('s')
plt.xlabel('Gauss?')
plt.show()

plt.clf()
plt.scatter(t_L,t_M)
plt.xlabel('s')
plt.xlabel('Gauss?')
plt.grid()
plt.show()


plt.clf()
plt.scatter(frequency_kHZ,dB_L_frequency)
plt.xlabel('kHz')
plt.ylabel('Gauss?/sqrt(Hz)')
plt.show()

plt.clf()
plt.scatter(frequency_kHZ,dB_M_frequency+dB_L_frequency)
plt.xlabel('kHz')
plt.ylabel('Gauss?/sqrt(Hz)')
plt.yscale('log')
plt.grid()
plt.show()


#mode number? 
#how to process? 
#standard process?