import pandas as pd

spectra_df = pd.read_csv("spectra.csv")

list=[]
for i in range(100):
    list.append(1205.0+5*i)

## removing wavlengths after 1200 since they have low contribution in absorption by si and pk
spectra_df = spectra_df.set_index("wavelength")
spectra_df = spectra_df.drop(list)
spectra_df = spectra_df.reset_index()
spectra_df['time'] = pd.to_datetime(spectra_df["time"], format='%Y-%m-%d %H:%M:%S')
spectra_df.time = (spectra_df.time - pd.Timestamp("2014-01-01")) // pd.Timedelta('1m')
##N is number of minutes for the interval to split the entire year into
N = 240
weather_df = pd.read_csv("weather.csv")
datetime_list=weather_df.time
weather_df['time'] = pd.to_datetime(weather_df["time"], format='%Y-%m-%d %H:%M:%S')
weather_df.time = (weather_df.time - pd.Timestamp("2014-01-01")) // pd.Timedelta('1m')
weather_df["DateTime"] = datetime_list
##files that will be there after breaking into 4 hour intervals- this was already calculated
files = 1130
datetime = []
azimuth = []
zenith = []
gni=[]
ghi=[]
temp =[]
temp_c=[]
for i in range(files):
    ## removing times with sun below horizon
    ze = weather_df.iloc[i * int(N / 5), 12]
    if ze >= 90.0:
        continue
    
    else:
        
        ## making splitted spectrum file
        temp_df = spectra_df.set_index("time")
        df = temp_df.loc[weather_df.iloc[i*int(N/5), 0]]
        ## naming the file
        filename_str = str(weather_df.iloc[i*int(N/5), 15])
        filename_str = filename_str.replace(":", "-")
        filename_str = filename_str.replace(" ", "-")
        filename_str = 'spectrum_files_real_time\\'+filename_str
        df.to_csv(filename_str+".spc", index=False, sep='\t',header=False)
        
        ## now appending hyperparameters at these times and ghi and gni values
        datetime.append(filename_str)
        az = weather_df.iloc[i * int(N / 5), 11]
        gn= weather_df.iloc[i * int(N / 5), 4]
        gh = weather_df.iloc[i * int(N / 5), 2]
        T= weather_df.iloc[i * int(N / 5), 7]
        T_C= weather_df.iloc[i * int(N / 5), 14]
        #az = az.replace(".", "p")
        #ze = ze.replace(".", "p")
        gni.append(gn)
        azimuth.append(az)
        zenith.append(ze)
        ghi.append(gh)
        temp.append(T)
        temp_c.append(T_C)
        



## saving file with hyperparameters : this will be used during optimization
df=pd.DataFrame({'DateTime': datetime, 'zenith': zenith, 'azimuth': azimuth,'Temperature':temp,'Temperature_cell':temp_c, 'GNI':gni, 'GHI':ghi})
df.to_csv("datetime_zenith_azimuth_temperature_gni_ghi.csv")



