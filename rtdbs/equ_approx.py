# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 17:24:59 2023

@author: Kaine
"""

import numpy as np
from scipy import interpolate
from functools import partial
from scipy.optimize import least_squares 
from scipy import integrate

def polinom_fun(x,rho):
    n=len(x)
    y=0
    for i in range(0,n):
        y=y+x[i]*rho**i
    return y

def approx_magn(geom_koeff,rho_a,theta,a0,r0,z0):
    delta=geom_koeff[0]
    lam=geom_koeff[1]
    gamma=geom_koeff[2]
    bxx=geom_koeff[3]
    bzz=geom_koeff[4]
    
    r=( -delta+rho_a*np.cos(theta)-gamma*np.sin(theta)**2-bxx*np.sin(2*theta)**2*np.cos(theta) )*a0+r0
    z=( lam*np.sin(theta)-bzz*np.sin(2*theta) )*a0+z0

    return r,z

def approx_magn_from_pol(geom_koeff_pol,rho_a,theta,a0,r0,z0):
    delta_koeff=geom_koeff_pol[0:4+1]
    lam_koeff=geom_koeff_pol[5:9+1]
    gamma_koeff=geom_koeff_pol[10:14+1]
    bxx_koeff=geom_koeff_pol[15:19+1]
    bzz_koeff=geom_koeff_pol[20:24+1]
    
    delta=polinom_fun(delta_koeff,rho_a)
    lam=polinom_fun(lam_koeff,rho_a)
    gamma=polinom_fun(gamma_koeff,rho_a)
    bxx=polinom_fun(bxx_koeff,rho_a)
    bzz=polinom_fun(bzz_koeff,rho_a)
    
    geom_koeff=[delta,lam,gamma,bxx,bzz]
    
    r,z=approx_magn(geom_koeff,rho_a,theta,a0,r0,z0)
    
    return r,z

def equ_approx(geom,
               max_it=10_000,
               rho_a_sec=np.linspace(0.3,1,10),
               theta_sec=np.linspace(-np.pi/3,+np.pi/3,40),
               alpha=1.217
               ):
    

    r_geo=geom['r_geo']
    z_geo=geom['z_geo']
    
    r_sep=geom['r_sep']
    z_sep=geom['z_sep']

    rho_a=geom['rho_a']
    theta=geom['theta']

    theta2d,rho_a2d=geom['theta2d'],geom['rho_a2d']
    
    #approx in predifined region
    r_geo=interpolate.RectBivariateSpline( theta, rho_a, r_geo  )(theta_sec,rho_a_sec)
    z_geo=interpolate.RectBivariateSpline( theta, rho_a, z_geo  )(theta_sec,rho_a_sec)
    rho_a=rho_a_sec
    theta=theta_sec
    
    #
    a0=(r_sep.max()-r_sep.min())/2
    r0=geom['rmaxis']
    z0=geom['zmaxis']

    def error_fun(x1,x2):
        return (x1-x2)**2

    def fit_fun(x,rho,theta,a0,r0,z0,r_in,z_in):
        r_out,z_out=approx_magn(x,rho,theta,a0,r0,z0)
        return error_fun(r_in,r_out)+error_fun(z_in,z_out)

    geom_koeff=list()

    for jj in range(0,len(rho_a)):
        rho_a_loc=rho_a[jj]
        r_geo_loc=r_geo[:,jj]
        z_geo_loc=z_geo[:,jj]
        
        #approx
        func_to_optimize = partial( 
                                        fit_fun,
                                        rho=rho_a_loc,
                                        theta=theta,
                                        a0=a0,
                                        r0=r0,
                                        z0=z0,
                                        r_in=r_geo_loc,
                                        z_in=z_geo_loc
                                  )
        
        initial=np.ones(4+1)
        
        output=least_squares(func_to_optimize,
                            initial,                     
                            method='lm',
                            max_nfev=max_it
                            ) 
        
        geom_koeff_loc=output['x']
        
        geom_koeff.append( geom_koeff_loc )

    geom_koeff=np.array(geom_koeff)

    geom_koeff_pol=list()
    for ii in range(0,geom_koeff.shape[1]):
        
        #approx
        def min_fun(geom_koeff_pol,rho_a,input_array):
            return ( polinom_fun(geom_koeff_pol,rho_a) - input_array )**2
        
        func_to_optimize = partial( 
                                        min_fun,
                                        rho_a=rho_a,
                                        input_array=geom_koeff[:,ii]
                                  )
        
        initial=[0,0,0,0,0]
        output=least_squares(func_to_optimize,
                            initial,                     
                            method='lm',
                            max_nfev=max_it
                            ) 
        
        geom_koeff_pol_loc=output['x']
        
        geom_koeff_pol.append( geom_koeff_pol_loc )
    geom_koeff_pol=np.array(geom_koeff_pol).flatten()

    r_geo_approx,z_geo_approx=approx_magn_from_pol(geom_koeff_pol,rho_a2d,theta2d,a0,r0,z0)
    r_geo_approx_alpha,z_geo_approx_alpha=approx_magn_from_pol(geom_koeff_pol,rho_a2d,theta2d,a0*alpha,r0,z0)

    #save
    result = {
    'a0':a0,
    'r0':r0,
    'z0':z0,

    'rho_a_sec':rho_a,
    'theta_sec':theta,

    'r_geo_approx':r_geo_approx,
    'z_geo_approx':z_geo_approx,

    'r_geo_approx_alpha':r_geo_approx_alpha,
    'z_geo_approx_alpha':z_geo_approx_alpha,

    'geom_koeff_pol':geom_koeff_pol
    }
    
    return result

def save_geom_koeff_pol(geom_approx,geom,filename_geom='input_test.dat',filename_Btor='Btor.txt'):
    #save data
    
    geom_koeff_pol=geom_approx['geom_koeff_pol']
    a0=geom_approx['a0']
    r0=geom_approx['r0']
    z0=geom_approx['z0']
    
    #input_test.dat
    f = open(filename_geom, 'w')
    
    f.write('           5\n')
    f.write('  '+str(a0*100)+'\n')
    f.write('  '+str(r0*100)+'\n')
    f.write('  '+str(z0*100)+'\n')
    for i in geom_koeff_pol:
        f.write('  '+str(i)+'\n')
    f.write('  '+str(0.0000)+'\n')
    f.write('  '+str(1.0000)+'\n')
    f.write('  '+str(2.0000)+'\n')
    f.write('  '+str(3.0000)+'\n')
    f.write('  '+str(4.0000)+'\n')
    
    f.close()
    
    #Btor.txt
    f = open(filename_Btor, 'w')
    
    f.write(str(geom['bcentr'])+'\n')
    f.write(str(geom['rcentr'])+'\n')

    f.write('//\n')
    f.write('*************************************************\n')
    f.write('в gфайле смотри строчки\n')
    f.write('BTOR =    3.941500000000000E-001\n')
    f.write('RCENTR =    3.600000000000000E-001\n')
    f.write('*************************************************\n')
    f.write('BTOR=Btor_0 Торроидальное магнитное поле при R=R0 в системе си! смотри g файл!\n')
    f.write('RCENTER=R0_gfile  //геом центр камеры смотри g файл\n')
    f.write(' \n')

    f.close()
    
    
    
def find_rho_theta_with_R_in_Z_in(geom_approx,alpha,R_in,Z_in,max_it=10_000):
    #koeff approx
    a0=geom_approx['a0']
    r0=geom_approx['r0']
    z0=geom_approx['z0']
    geom_koeff_pol=geom_approx['geom_koeff_pol']
    
    def min_fun(x,geom_koeff_pol,a0,alpha,r0,z0,R_in,Z_in):
        rho=x[0]
        theta=x[1]
        r,z=approx_magn_from_pol(geom_koeff_pol,rho,theta,a0*alpha,r0,z0)
        
        error_out=( (r-R_in)**2+(z-Z_in)**2 )*np.ones(2)
        
        if rho<0 or theta <-np.pi or theta>np.pi:
            error_out=error_out+1
        return error_out

    func_to_optimize = partial( 
                                    min_fun,
                                    geom_koeff_pol=geom_koeff_pol,
                                    a0=a0,
                                    alpha=alpha,
                                    r0=r0,
                                    z0=z0,
                                    R_in=R_in,
                                    Z_in=Z_in
                              )

    initial=[0.01,0.0]
    output=least_squares(func_to_optimize,
                        initial,                     
                        method='lm',
                        max_nfev=max_it
                        ) 

    rho_out,theta_out=output['x']
    
    return rho_out,theta_out

















def find_rho_sep(geom_approx,alpha,max_it=10_000):
    r_sep_approx=geom_approx['r_geo_approx'][:,-1]
    z_sep_approx=geom_approx['z_geo_approx'][:,-1]
    
    num_r_sep_approx_max=np.argmax(r_sep_approx)
    
    r_sep_approx_max=r_sep_approx[num_r_sep_approx_max]
    z_sep_approx_max=z_sep_approx[num_r_sep_approx_max]
    
    rho_sep,_=find_rho_theta_with_R_in_Z_in(geom_approx,alpha,r_sep_approx_max,z_sep_approx_max,max_it=max_it)
    
    return rho_sep

def find_dpsi_drho(g,geom_approx,alpha,rho_start,rho_min=0.15,max_it=10_000):
    #grid
    r=g['r'] #m
    z=g['z'] #m

    #lim
    r_lim=g['lim'][:,0] #m
    z_lim=g['lim'][:,1] #m
        
    #magnetic axis
    rmaxis=g['rmaxis'] #m
    zmaxis=g['zmaxis'] #m
    
    #psi
    psi=g['psirz'].T #Wb/rad
    
    #psi axis and boundary
    psi_bndry=g['ssibry'] #Wb/rad
    psi_axis=g['ssimag'] #Wb/rad
    
    #norm!
    psi = psi-psi_axis #strart from zero
    psi_bndry = psi_bndry - psi_axis
    psi_axis=0
    
    psi = np.abs(psi)
    psi_bndry=np.abs(psi_bndry)
    
    #increase resolution
    nr_new=500
    nz_new=1000
    
    r_new=np.linspace(r_lim.min()*0.95,r_lim.max()*1.05,nr_new)
    z_new=np.linspace(z_lim.min()*1.05,z_lim.max()*1.05,nz_new)
    
    psi=interpolate.RectBivariateSpline( r, z, psi  )(r_new,z_new)
    
    r=r_new
    z=z_new
    
    r2d,z2d=np.meshgrid(r,z,indexing="ij")
    
    # psi=geom['psi']
    # r2d=geom['r2d']
    # z2d=geom['z2d']

    # #find psi(rho) dependence
    # rmaxis=geom['rmaxis']
    # zmaxis=geom['zmaxis']
    
    #find point of magn. axis
    l_maxis_grid=np.sqrt( (r2d-rmaxis)**2 + (z2d-zmaxis)**2 )
    maxis_ind=np.argmin(l_maxis_grid)
    rmaxis_ind,zmaxis_ind=np.unravel_index(maxis_ind,psi.shape)

    psi_line=psi[rmaxis_ind:,zmaxis_ind]
    r_line=r2d[rmaxis_ind:,zmaxis_ind]
    z_line=z2d[rmaxis_ind:,zmaxis_ind]
    
    #for test
    psi_z0=psi[:,zmaxis_ind]
    r_z0=r2d[:,zmaxis_ind]
    z_z0=z2d[:,zmaxis_ind]
    
    rho_line=list()
    for i in range(0,len(r_line)):
        rho_l,theta_l=find_rho_theta_with_R_in_Z_in(geom_approx,alpha,r_line[i],z_line[i],max_it=max_it)
        rho_line.append( rho_l )
    rho_line=np.array(rho_line)

    rho_line_int=np.linspace(rho_min,rho_start,100)
    psi_line_int=np.interp(rho_line_int,rho_line,psi_line)    
    
    #dpsi_drho
    dpsi_line_drho=np.gradient(psi_line_int,rho_line_int)
    
    #approx dpsi_drho
    def min_fun(dpsi_drho_koeff,dpsi_line_drho,rho_line_int):
        return ( polinom_fun(dpsi_drho_koeff,rho_line_int) - dpsi_line_drho )**2

    func_to_optimize = partial( 
                                    min_fun,
                                    dpsi_line_drho=dpsi_line_drho,
                                    rho_line_int=rho_line_int
                              )

    initial=[0,0,0,0,0]
    output=least_squares(func_to_optimize,
                        initial,                     
                        method='lm',
                        max_nfev=max_it
                        ) 

    dpsi_drho_koeff_pol=output['x']

    dpsi_line_drho_approx=polinom_fun(dpsi_drho_koeff_pol,rho_line_int)
    
    #for testing and trouble shooting
    psi_line_int_approx=integrate.cumtrapz(dpsi_line_drho_approx, rho_line_int, initial=0)
    
    result = {'dpsi_drho_koeff_pol':dpsi_drho_koeff_pol,
    'rho_line':rho_line,
    'r_line':r_line,
    'z_line':z_line,
    'psi_line':psi_line,
    
    'rho_line_int':rho_line_int,
    'psi_line_int':psi_line_int,
    'psi_line_int_approx':psi_line_int_approx,
    
    'dpsi_line_drho':dpsi_line_drho,
    'dpsi_line_drho_approx':dpsi_line_drho_approx,
    
    'psi_z0':psi_z0,
    'r_z0':r_z0,
    'z_z0':z_z0
    
    }
    return result

def save_dpsi_drho_koeff_pol(dpsi_drho_data,rho_start,theta_start,n_start_vect_rot,output_filename='input_test2.dat'):
    dpsi_drho_koeff_pol=dpsi_drho_data['dpsi_drho_koeff_pol']
    f = open(output_filename, 'w')

    for i in dpsi_drho_koeff_pol:
        f.write(str(i)+'\n')

    f.write(' \n')
    f.write(str(rho_start)+' '+str(theta_start)+'\n')
    f.write(str(n_start_vect_rot[0])+' '+str(n_start_vect_rot[1])+' '+str(n_start_vect_rot[2])+'\n')
    f.write('5 coef dpsi/drho; w_start; theta_start\n')
    f.write('Nx Ny Nz start for ray tracing\n')

    f.close()
    



