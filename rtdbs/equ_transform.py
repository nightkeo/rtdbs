# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:13:28 2023

@author: Kaine
"""

import numpy as np
from scipy import interpolate
from skimage import measure

def generate_points(rr, zz, n_point=100):
    """
    Generate points along the line that may be used to check
    if the plasma is limited or not.
    """

    # Interpolate line points.
    # Make an interpolator for point location as function of normalised distance
    # along the wall
    points = np.array([rr, zz]).T
    distance = np.cumsum(np.sqrt(np.sum(np.diff(points, axis=0) ** 2, axis=1)))
    distance = np.insert(distance, 0, 0) / distance[-1]
    
    interpolator = interpolate.interp1d(distance, points, kind="linear", axis=0)
    new_distances = np.linspace(0, 1, n_point, endpoint=True)
    interpolated_points = interpolator(new_distances)

    R = np.asarray(interpolated_points[:, 0])
    Z = np.asarray(interpolated_points[:, 1])

    return R, Z

def interp_line(x,y,x0,y0,ntheta):
    
    #interpolate and sort with theta
    th=np.arctan2(y-y0,x-x0)
    
    n_th=np.argsort(th)
    th=th[n_th]
    x=x[n_th]
    y=y[n_th]
    
    th_new=np.linspace(th[0],th[-1],ntheta)
    
    x_new = interpolate.interp1d(th,x,kind='linear')(th_new)
    y_new = interpolate.interp1d(th,y,kind='linear')(th_new)
    
    #interpolate and sort with theta - TRUE
    lpol_l=( np.gradient(x_new)**2+np.gradient(y_new)**2 )**0.5
    lpol_f=np.sum(lpol_l)
    lpol=np.cumsum(lpol_l)
    
    th=lpol/lpol_f*2*np.pi-np.pi

    th_new=np.linspace(th[0],th[-1],ntheta)
    
    x_new2 = interpolate.interp1d(th,x_new,kind='linear')(th_new)
    y_new2 = interpolate.interp1d(th,y_new,kind='linear')(th_new)
    
    return x_new2,y_new2

def equ_transform(g,nrho=41,ntheta=151):
    
    nr=len(g['r'])
    nz=len(g['z'])
    
    r=g['r'] #m
    z=g['z'] #m
    
    rmin=r.min()
    rmax=r.max()
    
    zmin=z.min()
    zmax=z.max()
    
    dr=np.diff(r)[0]
    dz=np.diff(z)[0]
    
    r2d,z2d=np.meshgrid(r,z,indexing="ij")
    
    r2d=r2d #m
    z2d=z2d #m
    
    #psi
    psi=g['psirz'].T #Wb/rad
    
    #sep
    n_sep=g['nbdry']
    r_sep=g['bdry'][:,0] #m
    z_sep=g['bdry'][:,1] #m
    
    #lim
    n_lim=g['limitr']
    r_lim=g['lim'][:,0] #m
    z_lim=g['lim'][:,1] #m
        
    #vacuum toroidal fileld
    bcentr=g['bcentr'] #T
    
    #limiter centr (for bt calc!)
    rcentr=g['rcentr'] #m
        
    #plasma current
    current=g['current'] #A
        
    #magnetic axis
    rmaxis=g['rmaxis'] #m
    zmaxis=g['zmaxis'] #m
    
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

    #extend number of limiter points    
    r_sep_new,z_sep_new=generate_points(r_sep,z_sep,n_point=1000)
    r_sep=r_sep_new
    z_sep=z_sep_new
    n_sep=len(r_sep)
    
    #extend number of limiter points
    r_lim_new,z_lim_new=generate_points(r_lim,z_lim,n_point=1000)
    r_lim=r_lim_new
    z_lim=z_lim_new
    n_lim=len(r_lim)
    
    # r_new=np.linspace(r_lim.min()*0.95,r_lim.max()*1.05,nr_new)
    # z_new=np.linspace(z_lim.min()*1.05,z_lim.max()*1.05,nz_new)

    r_new=np.linspace(r_sep.min()*0.95,r_sep.max()*1.05,nr_new)
    z_new=np.linspace(z_sep.min()*1.05,z_sep.max()*1.05,nz_new)
    
    psi_fun=interpolate.RectBivariateSpline( r2d[:, 0], z2d[0, :], psi  )
    psi_new=psi_fun(r_new,z_new)
    
    r=r_new
    z=z_new
    
    nr=nr_new
    nz=nz_new
    
    psi=psi_new
    
    dr=np.diff(r)[0]
    dz=np.diff(z)[0]
    
    r2d,z2d=np.meshgrid(r,z,indexing="ij")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #grid for psi
    rho=np.linspace(0,1,nrho)
    
    #psi norm
    psi_norm = (psi - psi_axis) / (psi_bndry - psi_axis)
    psi_norm=np.clip(psi_norm, 0.0, 1.0)
    
    theta=np.linspace(-np.pi,+np.pi,ntheta)
    
    r_geo=list()
    z_geo=list()
    for i in range(1,len(rho)-1):
        cl=[p for p in measure.find_contours(psi_norm, rho[i]) if (p[0,0] == p[-1,0])][0]
    
        x=np.interp(cl[:,0], range(0,len(r)), r)
        y=np.interp(cl[:,1], range(0,len(z)), z)
        
        x_new,y_new=interp_line(x,y,rmaxis,zmaxis,ntheta)
        
        r_geo.append(x_new)
        z_geo.append(y_new)
    
    #append magnetic axis
    r_geo.insert(0,rmaxis*np.ones(len(theta)))
    z_geo.insert(0,zmaxis*np.ones(len(theta)))
    
    # #append plasma boundary
    # x_new,y_new=interp_line(r_sep_new,z_sep_new,rmaxis,zmaxis,ntheta)
    # r_geo.append(x_new)
    # z_geo.append(y_new)
    
    #ALTERNATIVE FOR PLASMA BOUNDARY
    cl=[p for p in measure.find_contours(psi_norm, 0.999) if (p[0,0] == p[-1,0])][0]
    x=np.interp(cl[:,0], range(0,len(r)), r)
    y=np.interp(cl[:,1], range(0,len(z)), z)
    x_new,y_new=interp_line(x,y,rmaxis,zmaxis,ntheta)
    r_geo.append(x_new)
    z_geo.append(y_new)
    
    r_geo=np.array(r_geo).T
    z_geo=np.array(z_geo).T
    
    # сшивка кривых при theta=-pi и +pi
    r_stich = (r_geo[0].copy()+r_geo[-1].copy())/2
    z_stich = (z_geo[0].copy()+z_geo[-1].copy())/2
    r_geo[0] = r_geo[-1] = r_stich
    z_geo[0] = z_geo[-1] = z_stich
    
    rho_a=(np.max(r_geo,axis=0)-np.min(r_geo,axis=0))/2
    rho_a=rho_a/rho_a[-1]
    
    #2d grid 
    theta2d,rho_a2d=np.meshgrid(theta,rho_a,indexing='ij')
    
    # Create dictionary of values to return
    result = {'nr': nr, 'nz':nz,
              'r':r, 'z':z,
              'r2d':r2d, 'z2d':z2d,
              'rmin':rmin, 'rmax':rmax,
              'zmin':zmin, 'zmax':zmax,
              'dr':dr, 'dz':dz,
              'psi':psi,
              'psi_norm':psi_norm,
              'psi_bndry':psi_bndry,
              'psi_axis':psi_axis,
              'rmaxis':rmaxis, 'zmaxis':zmaxis,
              'bcentr':bcentr, 'rcentr':rcentr,
              'current':current,
              'n_sep':n_sep, 'r_sep':r_sep, 'z_sep':z_sep,
              'n_lim':n_lim, 'r_lim':r_lim, 'z_lim':z_lim,
              
              'r_geo':r_geo, 'z_geo':z_geo,
              'ntheta':ntheta, 'nrho':nrho,
              
              'theta':theta,
              'theta2d':theta2d,
              
              'rho':rho,
              
              'rho_a':rho_a,              
              'rho_a2d':rho_a2d,  
              }
    
    return result


