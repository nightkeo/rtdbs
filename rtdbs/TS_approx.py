# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:56:59 2023

@author: Kaine
"""

import numpy as np

from scipy.optimize import least_squares 
from functools import partial
from scipy import interpolate
from rtdbs.equ_approx import find_rho_theta_with_R_in_Z_in, polinom_fun

def TS_approx(geom_approx,r_ts,z_ts,ne_ts,dne_ts,alpha,rho_start,max_it=10_000):
    
    #copy data and "save" initial data
    r_ts_in=r_ts.copy()
    z_ts_in=z_ts.copy()
    ne_ts_in=ne_ts.copy()
    dne_ts_in=dne_ts.copy()
    
    #extract points in approx region
    r_geo_approx_alpha=geom_approx['r_geo_approx_alpha']
    z_geo_approx_alpha=geom_approx['z_geo_approx_alpha']

    n1=np.linspace(0,1,r_geo_approx_alpha.shape[0])
    n2=np.linspace(0,1,r_geo_approx_alpha.shape[1])

    n11=np.linspace(0,1,r_geo_approx_alpha.shape[0]*20)
    n22=np.linspace(0,1,r_geo_approx_alpha.shape[1]*20)
    
    r_geo_sec_approx2=interpolate.RectBivariateSpline( n1, n2, r_geo_approx_alpha  )(n11,n22)
    z_geo_sec_approx2=interpolate.RectBivariateSpline( n1, n2, z_geo_approx_alpha  )(n11,n22)
    
    tol=1e-3    

    l_ts_region=list()
    for i in range(len(r_ts)):
        l_loc=np.sqrt( (r_ts[i]-r_geo_sec_approx2)**2 + (z_ts[i]-z_geo_sec_approx2)**2 )
        l_ts_region.append( np.min(l_loc) )
    l_ts_region=np.array(l_ts_region)
    n_ts_true=np.argwhere(l_ts_region<tol)[:,0]

    r_ts=r_ts[n_ts_true]
    z_ts=z_ts[n_ts_true]
    ne_ts=ne_ts[n_ts_true]
    dne_ts=dne_ts[n_ts_true]

    rho_ts=list()
    theta_ts=list()
    for i in range(len(r_ts)):
        out=find_rho_theta_with_R_in_Z_in(geom_approx,alpha,r_ts[i],z_ts[i],max_it=max_it)
        rho_ts.append( out[0] )
        theta_ts.append( out[1] )
    rho_ts=np.array(rho_ts)
    theta_ts=np.array(theta_ts)

    #sort
    n_sort_ts=np.argsort(rho_ts)
    rho_ts=rho_ts[n_sort_ts]
    r_ts=r_ts[n_sort_ts]
    z_ts=z_ts[n_sort_ts]
    ne_ts=ne_ts[n_sort_ts]
    dne_ts=dne_ts[n_sort_ts]

    #append density at start point
    rho_ts=np.append(rho_ts,rho_start)
    r_ts=np.append(r_ts,0)
    z_ts=np.append(z_ts,0)
    ne_ts=np.append(ne_ts,np.min(ne_ts)*1e-6)
    # ne_ts=np.append(ne_ts,1e9/1e13)
    dne_ts=np.append(dne_ts,0.001)
    
    #append density at start2 point (EXP!)
    rho_start2=1.5
    rho_ts=np.append(rho_ts,rho_start2)
    ne_ts=np.append(ne_ts,np.min(ne_ts)*1e-3)
    dne_ts=np.append(dne_ts,np.min(dne_ts)*1e-3)
    # rho_start2=1.5
    # rho_ts=np.append(rho_ts,rho_start2)
    # ne_ts=np.append(ne_ts,-2.0)
    # dne_ts=np.append(dne_ts,0.001)
    
    #append density at magnetic axis
    rho_ts=np.insert(rho_ts,0,0)
    r_ts=np.insert(r_ts,0,0)
    z_ts=np.insert(z_ts,0,0)
    ne_ts=np.insert(ne_ts,0,np.max(ne_ts))
    dne_ts=np.insert(dne_ts,0,0.001)
    
    print(rho_ts)
    print(ne_ts)
    
    nrho_ts_int=100
    rho_ts_int=np.linspace(0,rho_start2,nrho_ts_int)

    ne_ts_int_up=np.interp(rho_ts_int,rho_ts,ne_ts+dne_ts)
    ne_ts_int_low=np.interp(rho_ts_int,rho_ts,ne_ts-dne_ts)
    ne_ts_int=(ne_ts_int_up+ne_ts_int_low)/2

    # def ne_fun(ne_koeff,rho):
    #     return (ne_koeff[0]-ne_koeff[1])*(1-(rho/np.max(rho))**ne_koeff[2])**ne_koeff[3]+ne_koeff[1]

    def min_fun(ne_koeff,ne,rho):
        ne_approx=polinom_fun(ne_koeff,rho)
        # ne_approx=ne_fun(ne_koeff,rho)
        out=(ne-ne_approx)**2
        # dne_approx_norm=np.gradient(ne_approx,rho)
        # # dne_approx_norm=np.nan_to_num(dne_approx_norm,nan=1)
        # dne_pol=np.where(dne_approx_norm>0,1,0)*dne_approx_norm
        # return out+np.abs(dne_pol)
        return out

    #approx
    func_to_optimize = partial( 
                                    min_fun,
                                    ne=ne_ts_int,
                                    rho=rho_ts_int
                                    )

    initial=np.zeros(5)+2
    
    output=least_squares(func_to_optimize,
                        initial,                     
                        method='lm',
                        max_nfev=max_it
                        ) 

    ne_koeff=output['x']

    ne_approx=polinom_fun(ne_koeff,rho_ts_int)
    # ne_approx=ne_fun(ne_koeff,rho_ts_int)
    
    result = {'ne_koeff':ne_koeff,

              #initial data
              'r_ts_in':r_ts_in,
              'z_ts_in':z_ts_in,
              'ne_ts_in':ne_ts_in,
              'dne_ts_in':dne_ts_in,
              
              'r_ts':r_ts,
              'z_ts':z_ts,
              'ne_ts':ne_ts,
              'dne_ts':dne_ts,
              
              'rho_ts':rho_ts,
              'theta_ts':theta_ts,
              'rho_ts_int':rho_ts_int,
              'ne_ts_int':ne_ts_int,
              'ne_ts_int_low':ne_ts_int_low,
              'ne_ts_int_up':ne_ts_int_up,
              'ne_approx':ne_approx
              }
              
    return result

def save_dens(dens,filename='input_test3.dat'):
    out=np.vstack((dens['rho_ts_int'],dens['ne_ts_int'])).T
    np.savetxt(filename,out)


























