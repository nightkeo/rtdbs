# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:28:05 2023

@author: Kaine
"""

import numpy as np

def TS_read(filename):
    data=np.loadtxt(filename, dtype='str',delimiter=',')
    
    r=data[1:,0].astype(np.float64)
    
    t_1=data[0,1:]
    t=np.zeros(int(len(t_1)/2))
    m=list()
    dm=list()
    j=0
    for i in range(0,len(t_1),2):
        t[j]=t_1[i][0:t_1[i].find('s')]
        j=j+1
        m.append(data[1:,i+1].astype(np.float64))
        dm.append(data[1:,i+2].astype(np.float64))
    t=t.astype(np.float64)
    
    #delete bad data
    m_new=list()
    r_new=list()
    dm_new=list()
    
    for j in range(0,len(t)):
        num_bad=(m[j]==1)
        buf_r=list()
        buf_m=list()
        buf_dm=list()
        for i in range(0,len(r)):
            if not num_bad[i]:
                buf_r.append(r[i])
                buf_m.append(m[j][i])
                buf_dm.append(dm[j][i])
        
        r_new.append(np.array(buf_r))
        m_new.append(np.array(buf_m))
        dm_new.append(np.array(buf_dm))

    
    # Create dictionary of values to return
    result = {
              'ntime': len(t),
              'time':t,
              'r':r_new,
              'array':m_new,
              'darray':dm_new
              }
    return result