import numpy as np
import casac
import json


######### loading json parameters######
with open('params.json') as json_data:
    d = json.load(json_data)

path = d["Simulator_params"]["working_directory"].encode("ascii","ignore")

    
def add_receiver_noise(msname, sefd, spw_id=0, load=None):
        """ baseline dependent thermal noise only. Calculated from SEFDs, tint, dnu """
        
        if load:
            thermal_noise = np.load(path+'/receiver_noise.npy')
        else:
            
            tb = casac.casac.table()
            
            # get data and meta-data from MS
            tb.open(path+msname) 
            data = tb.getcol('DATA').T # all the visibilities.
            A0 = tb.getcol('ANTENNA1') # ant1 in baseline y
            A1 = tb.getcol("ANTENNA2") # ant2 in baseline y  
            ant_unique = np.unique(np.hstack((A0, A1))) 
            time_unique = np.unique(tb.getcol('TIME')) # all unique times (i.e. the number the time samples in the observation
            tint = time_unique[1]-time_unique[0] # time resolution 
            baseline_dict = dict([((x, y), np.where((A0 == x) & (A1 == y))[0]) 
                                      for x in ant_unique for y in ant_unique if y > x]) # create baseline dictionary - very useful!
            tb.close()

            tb.open(path+msname+'/SPECTRAL_WINDOW')
            chan_freq = tb.getcol('CHAN_FREQ').flatten()       # centre freqeuncy of all channels
            chan_width = tb.getcol('CHAN_WIDTH').flatten()[0]  # chan width of all channels
            tb.close()
            
            tb.open(path+msname+'/ANTENNA')
            station_names = tb.getcol('NAME') # get antenna/station names
            Nant = len(station_names)
            tb.close()
            
            
            thermal_noise = np.zeros(data.shape, dtype='complex')
            receiver_rms = np.zeros(data.shape, dtype='float')
            size = (time_unique.shape[0], chan_freq.shape[0], 4)
            for a0 in range(Nant):
                for a1 in range(Nant):
                    if a1 > a0:
                        rms = sefd/np.sqrt(abs(2 * tint * chan_width))

                        thermal_noise[baseline_dict[(a0, a1)]] =\
                            np.random.normal(0.0, rms, size=size) + 1j * np.random.normal(0.0, rms, size=size)
                        receiver_rms[baseline_dict[(a0,a1)]] = rms 

            np.save(path+'receiver_noise', thermal_noise)
        data = np.add(data, thermal_noise)
        tb.open(path+msname,nomodify=False)
        tb.putcol('DATA',data.T)
        tb.flush()
        tb.unlock()
        tb.close()
        print("Thermal noise added")
        
print("Adding Thermal noise ...............")
add_receiver_noise(d["Simulator_params"]["msfile_name"].encode("ascii","ignore"), 
                   d["Simulator_params"]["SEFD_Jy"], spw_id=0, load=None)