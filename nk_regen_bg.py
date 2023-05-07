import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', -1)
file_pero="jw_pero3cat_2to1_match_interp.nk"
file_si="c-Si_interp.nk"
temp_pero="integer_template_pero.nk"

##temp is the file containing values for interpolation purpose
df_temp_pero=pd.read_csv(temp_pero,sep="\t",index_col=False)

f="datetime_zenith_azimuth_temperature.csv"
dat=pd.read_csv(f)

## create shifted bandgap nk file
def wav_at(T):
    bg_p=1.710344827586206896551724137931 + 0.00031*(T-25)
    bg_s=1.11-0.000273*(T-27)
    wl_p=1240/bg_p
    wl_s=1240/bg_s
    a=pd.read_csv(file_pero,sep="\t",index_col=False)
    # plt.plot(a['nm'], a['k'],label='o')
    a['nm']=a['nm']+(wl_p-725)
    res=np.interp(df_temp_pero['nm'],a['nm'],a['k'])
    res_n = np.interp(df_temp_pero['nm'], a['nm'], a['n'])
    # plt.plot(a['nm'],res,label='n')
    # plt.legend()
    # plt.show()
    a['nm'] = df_temp_pero['nm'].astype(int)
    a['k'] = res
    a['n'] =res_n

    a.to_csv("folder_nk_bg_shifted\\"+'pero_n_k_' + str(T) + '.nk',index=False,sep='\t')

    b = pd.read_csv(file_si, sep="\t", index_col=False)
    # plt.plot(b['nm'], b['k'], label='o')
    b['nm'] = b['nm'] + (wl_s - 1117.177)
    res = np.interp(df_temp_pero['nm'], b['nm'], b['k'])
    res_n = np.interp(df_temp_pero['nm'], a['nm'], a['n'])
    # plt.plot(b['nm'], res, label='n')
    # plt.legend()
    # plt.show()
    b['nm'] = df_temp_pero['nm'].astype(int)
    b['k'] = res
    b['n'] =res_n
    b.to_csv("folder_nk_bg_shifted\\"+'si_n_k_' + str(T) + '.nk',index=False,sep='\t')

arr=dat['Temperature_cell']

for i in range(len(arr)):
    wav_at(arr[i])





