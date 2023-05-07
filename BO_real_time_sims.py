## Copyright. All rights reserved.
## This code is the intellectual property of Shaurya Anand.
import matlab.engine
from os.path import exists
import pandas
import datetime
import os
import sys
import pandas as pd
from numpy import savetxt
import sklearn.gaussian_process as gp
import numpy as np
from numpy import genfromtxt
from sklearn import preprocessing
import matplotlib.pyplot as plt
import time
from scipy.stats import norm
from sklearn.gaussian_process.kernels import ExpSineSquared


# "N is no. of layers"
N = 15

df_n = pandas.read_csv('data_frame_tandem_real_time_sims.csv')
df_n.columns = ['start', 'end', 'total points']

##training set - first argument here
X = np.empty(N,dtype=object)
for i in range(N):
    X[i] = np.random.randint(df_n.iloc[i, 0], df_n.iloc[i, 1], 1).astype(int)

#Copying the set
X_new = X.astype(int)

#Converting the array into a data frame
df = pandas.DataFrame()
for i in range(N):
    df[i] = X[i]

##preprocessing
df = df.astype(float)

##scanning space set
for i in range(N):
    d = df[i].to_numpy()
    d = d.reshape(-1, 1)
    scaler_x = preprocessing.StandardScaler().fit(d)
    df[i] = scaler_x.transform(d)


## Function to generate the test set or scanning set.
def scan_set_gen():
    dftest = pandas.DataFrame()
    temp = np.empty(N, dtype=object)
    for i in range(N):
        temp[i] = np.random.randint(df_n.iloc[i, 0], df_n.iloc[i, 1], df_n.iloc[i, 2]).astype(int)
    for i in range(N):
        dftest[i] = temp[i].flatten()
    dftest = dftest.astype(float)
    scaler = np.empty(N, dtype=object)
    for i in range(N):
        d = dftest[i].to_numpy()
        d = d.reshape(-1, 1)
        scaler[i] = preprocessing.StandardScaler().fit(d)
        dftest[i] = scaler[i].transform(d)
    return dftest, scaler

c = 0
iter = 1000
maximum_tracker = np.zeros(iter)
start_time = time.time()

F_train=[]
dat_frame=[]
Table=[]

print("Enter run name")
rname=str(input())

dza = pd.read_csv("datetime_zenith_azimuth.csv")
files = 1092
for t in range(iter):
    P_max_year = 0.0
    ## Input files are created in the correct folder from where a Matlab script can pick them up and simulate it using CROWM in terminal. The same script transfers the file with values of current etc. back to another folder
    for k in range(files):
        st = '_'
        for i in range(N):
            st = st + str(X_new[i]) + "_"
        original = "PK_Si_tandem_DST_reference_real_time_sims.txt"  ## An example input file is put in the repository 
        f = open(r"C:\Users\Shaurya\Crowm_sims_BO_shaurya\In\PK_Si_tandem" + st + str(dza.DateTime[k])+"_" + str(dza.zenith[k])+"_" + str(dza.azimuth[k])+ ".txt", "w")
        ct = 0
        ## Reformatting original input file for our simulation
        for line in open(original):
            ct = ct+1
            li=line.strip()
            if ct == 10:
                li = "Illumination spectrum file:            "+str(dza.DateTime[k])+".spc"
                f.writelines(li + "\n")
            elif ct == 13:
                li = "Incident zenith angle [deg]:           "+str(dza.zenith[k])
                f.writelines(li + "\n")
            elif ct == 14:
                li = "Incident azimuth angle [deg]:          "+str(dza.azimuth[k])
                f.writelines(li + "\n")
            elif ct == 47:
                li = "Base name of the output files:         PK_Si_tandem" + st + str(dza.DateTime[k])+"_" + str(dza.zenith[k])+"_" + str(dza.azimuth[k])
                f.writelines(li + "\n")
                ## base name is the name of the output file
            else:
                f.writelines(li + "\n")

        ##putting thickness numbers
        for i in range(N):
            l = str(X_new[i]) + "\n"
            f.write(l)

        l = str(0)
        f.write(l)
        f.close()
    ## Waiting for file after a matlab script transfers it back to the correct folder
    for k in range(files):
        print("Got File no. ",k)
        file_exists=False
        print("Waiting... at",datetime.datetime.now())
        while(file_exists == False):
            path_to_file = r"C:\Users\Shaurya\Crowm_sims_BO_shaurya\Out\PK_Si_tandem" + st + str(dza.DateTime[k])+"_" + str(dza.zenith[k])+"_" + str(dza.azimuth[k])+"_Jsc.txt"
            file_exists = exists(path_to_file)
            if file_exists==True:
                print("File found")
                break
            else:
                time.sleep(10)
        ##reading current density values
        ct = 0
        f1 = open(r"C:\Users\Shaurya\Crowm_sims_BO_shaurya\Out\PK_Si_tandem" + st + str(dza.DateTime[k])+"_" + str(dza.zenith[k])+"_" + str(dza.azimuth[k])+"_Jsc.txt", "r")
        for line in f1:
            ct= ct + 1
            if ( ct == 10 ):
                l1 = float(line[:11]) ##Jsc_pk
            if ( ct == 14 ):
                l2 = float(line[:11]) ##Jsc_si
        ##MATLAB electrical simulation code start with 4 inputs : thickness,finger number, Jsc_pk,Jsc_si
        eng = matlab.engine.start_matlab()
        thickness = float(X_new[2])
        finger_no = float(2)
        Jsc_Pk = float(l1*0.001)
        Jsc_Si = float(l2*0.001)
        print(Jsc_Si, Jsc_Pk, thickness)
        Jsc_device, Jsc_CROWM, Jsc_total_unshadowed, Jsc_total_shadowed, Jsc_top, Jsc_bot, Voc, FF, Eff, V_mpp, J_mpp = eng.Tandem3Diode_v2(thickness, finger_no, Jsc_Pk, Jsc_Si, nargout=11)
        eng.quit()
        ##Power calculation
        P_max_year = P_max_year + (V_mpp * J_mpp * 0.01 )  #

    ##copying the thickness and corresponding current vlaues into a file for later use
    F_train = np.append(F_train, P_max_year)
    X_new_copy = X_new
    dat_frame.append(X_new_copy)
    savetxt("df_"+rname+".csv", dat_frame)
    savetxt("F_train_"+rname+".txt", F_train)

    f_max = np.max(F_train)
    ##GPR model fitting
    dftest, scaler = scan_set_gen()
    gpr = gp.GaussianProcessRegressor(kernel=None, optimizer = " fmin_l_bfgs_b", n_restarts_optimizer=10, alpha=1e-3,
                                      normalize_y=True, random_state=2).fit(df,F_train)

    mu, sigma = gpr.predict(dftest, return_std=True)

    p = (mu - f_max) / sigma
    CDF = norm.cdf(p)
    PDF = norm.pdf(p)

    EI = (mu - f_max) * CDF + sigma * PDF

    EI[np.where(sigma < 1e-4)] = 0

    newsample = dftest.iloc[np.argmax(EI)]

    X_new = np.empty(N, dtype=int)

    for i in range(N):
        X_new[i] = scaler[i].inverse_transform((newsample[i]).reshape(-1, 1))

    ##thickness values into the dataset , Ftrain already has current values

    df.loc[len(df)] = newsample

    #F_train = np.append(F_train, (F[tuple(X_new[i] for i in range(len(X_new)))]))
    print(P_max_year, "          P_max_year at   ", X_new_copy)
    
    ## convergence criteria for maxima found
    if (abs(Eff - f_max)/f_max) < 1e-3:
        tick = True
        print("Location         " + "Maxima")
        print(*X_new_copy, "         ", P_max_year)
        last_iter = t

    else:
        tick = False

    if  tick == True :
        c = c + 1
    else:
        c = 0

    if c == 15:
        break

    maximum_tracker[t] = np.max(F_train)

##plots
plt.plot(maximum_tracker[0:last_iter], marker='+')
plt.title("Z_value vs iterations")
plt.xlabel("Iterations")
plt.ylabel("Z_value")
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
