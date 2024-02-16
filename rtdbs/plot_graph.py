# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:45:41 2023

@author: Kaine
"""

import matplotlib.pylab as plt
import numpy as np
from scipy import interpolate

def equ_plot(geom,geom_approx,dens,output_path_all,plot_type='full_geo_vs_full_approx'):
    r_lim=geom['r_lim']
    z_lim=geom['z_lim']

    r_geo=geom['r_geo']
    z_geo=geom['z_geo']

    r_sep=geom['r_sep']
    z_sep=geom['z_sep']

    r_geo_approx=geom_approx['r_geo_approx']
    z_geo_approx=geom_approx['z_geo_approx']

    r_geo_approx_alpha=geom_approx['r_geo_approx_alpha']
    z_geo_approx_alpha=geom_approx['z_geo_approx_alpha']
    
    r_ts_in=dens['r_ts_in']
    z_ts_in=dens['z_ts_in']
    
    #for good view
    rho_a_sec=np.linspace(0.0,1,10),
    theta_sec=np.linspace(-np.pi,+np.pi,60)
    rho_a=geom['rho_a']
    theta=geom['theta']
    r_geo=interpolate.RectBivariateSpline( theta, rho_a, r_geo  )(theta_sec,rho_a_sec)
    z_geo=interpolate.RectBivariateSpline( theta, rho_a, z_geo  )(theta_sec,rho_a_sec)

    r_geo_approx=interpolate.RectBivariateSpline( theta, rho_a, r_geo_approx  )(theta_sec,rho_a_sec)
    z_geo_approx=interpolate.RectBivariateSpline( theta, rho_a, z_geo_approx  )(theta_sec,rho_a_sec)

    r_geo_approx_alpha=interpolate.RectBivariateSpline( theta, rho_a, r_geo_approx_alpha  )(theta_sec,rho_a_sec)
    z_geo_approx_alpha=interpolate.RectBivariateSpline( theta, rho_a, z_geo_approx_alpha  )(theta_sec,rho_a_sec)
    
    fig, ax = plt.subplots(figsize=(8,10))
    ax.set_aspect('equal')
    
    if plot_type=='full_geo_vs_full_approx':
        ax.set_title('full_geo_vs_full_approx')
        ax.plot(r_geo,z_geo,'-k')
        ax.plot(r_geo.T,z_geo.T,'-k')
        ax.plot(r_geo_approx,z_geo_approx,'-r')
        ax.plot(r_geo_approx.T,z_geo_approx.T,'-r')
    elif plot_type=='full_geo':
        ax.set_title('full_geo')
        ax.plot(r_geo,z_geo,'-k')
        ax.plot(r_geo.T,z_geo.T,'-k')
    elif plot_type=='full_approx':
        ax.set_title('full_approx')
        ax.plot(r_geo_approx,z_geo_approx,'-r')
        ax.plot(r_geo_approx.T,z_geo_approx.T,'-r')   
    elif plot_type=='full_approx_alpha':
        ax.set_title('full_approx_alpha')
        ax.plot(r_geo_approx_alpha,z_geo_approx_alpha,'-r')
        ax.plot(r_geo_approx_alpha.T,z_geo_approx_alpha.T,'-r')  
    else:
        print('Wrong settings!!!')
    
    #lim
    ax.plot(r_lim,z_lim,'-g',label='lim')

    #sep
    ax.plot(r_sep,z_sep,'-m',label='sep')
    
    #ts
    ax.plot(r_ts_in,z_ts_in,'ob')
    
    ax.set_ylabel('Z, m')
    ax.set_xlabel('R, m')
    ax.set_xlim(xmin=r_lim.min()*0.7,xmax=r_lim.max()*1.2)
    ax.set_ylim(ymin=z_lim.min()*1.2,ymax=z_lim.max()*1.2)
    ax.legend()
    fig.tight_layout()
    plt.show()
    
    fig.savefig(output_path_all+'//'+plot_type+'.png', dpi=fig.dpi)
    
    plt.close()

def dpsi_drho_graph(dpsi_drho_data,output_path_all):       
    rho_line_int=dpsi_drho_data['rho_line_int']
    
    psi_line_int=dpsi_drho_data['psi_line_int']
    psi_line_int_approx=dpsi_drho_data['psi_line_int_approx']

    dpsi_line_drho=dpsi_drho_data['dpsi_line_drho']
    dpsi_line_drho_approx=dpsi_drho_data['dpsi_line_drho_approx']
    
    #Psi vs rho
    fig, ax = plt.subplots(figsize=(8,10))
    ax.set_title('Psi vs rho')
    ax.plot(rho_line_int,psi_line_int,'-.b',label='psi')
    ax.plot(rho_line_int,psi_line_int_approx,'--r',label='psi_approx')
    ax.set_ylabel('psi, Wb/rad')
    ax.set_xlabel('rho')
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig(output_path_all+'//'+'psi_vs_rho'+'.png', dpi=fig.dpi)
    plt.close()
    
    #dPsi_drho vs rho
    fig, ax = plt.subplots(figsize=(8,10))
    ax.set_title('dPsi_drho vs rho')
    ax.plot(rho_line_int,dpsi_line_drho,'-.b',label='dpsi_drho')
    ax.plot(rho_line_int,dpsi_line_drho_approx,'--r',label='dpsi_drho_approx')
    ax.set_ylabel('dpsi_drho, Wb/rad')
    ax.set_xlabel('rho')
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig(output_path_all+'//'+'dpsi_drho_vs_rho'+'.png', dpi=fig.dpi)
    plt.close()

def dens_plot(dens,output_path_all):
    rho_ts=dens['rho_ts']
    ne_ts=dens['ne_ts']
    dne_ts=dens['dne_ts']
    
    rho_ts_int=dens['rho_ts_int']
    ne_ts_int=dens['ne_ts_int']
    ne_ts_int_low=dens['ne_ts_int_low']
    ne_ts_int_up=dens['ne_ts_int_up']
    
    #ne vs rho
    fig, ax = plt.subplots(figsize=(7,8))
    ax.set_title('ne vs rho')
    ax.errorbar(rho_ts, ne_ts, yerr=dne_ts,fmt='ob')
    ax.plot(rho_ts_int,ne_ts_int,'-r',label='mean')
    ax.plot(rho_ts_int,ne_ts_int_low,'-r',label='low')
    ax.plot(rho_ts_int,ne_ts_int_up,'-r',label='up')
    ax.set_ylabel('$n_e$, $10^{19}$ $m^{-3}$')
    ax.set_xlabel('rho')
    ax.legend()
    fig.tight_layout()
    plt.show()
    fig.savefig(output_path_all+'//'+'ne_vs_rho'+'.png', dpi=fig.dpi)
    plt.close()
    



def ray_plot(geom,ray_data,output_path_all):
    r_lim=geom['r_lim']
    z_lim=geom['z_lim']

    r_geo=geom['r_geo']
    z_geo=geom['z_geo']

    r_sep=geom['r_sep']
    z_sep=geom['z_sep']
    
    r_ray=ray_data['R']
    z_ray=ray_data['Z']
    
    #for good view
    rho_a_sec=np.linspace(0.0,1,10),
    theta_sec=np.linspace(-np.pi,+np.pi,60)
    rho_a=geom['rho_a']
    theta=geom['theta']
    r_geo=interpolate.RectBivariateSpline( theta, rho_a, r_geo  )(theta_sec,rho_a_sec)
    z_geo=interpolate.RectBivariateSpline( theta, rho_a, z_geo  )(theta_sec,rho_a_sec)
    
    fig, ax = plt.subplots(figsize=(8,10))
    ax.set_aspect('equal')
    
    ax.set_title('ray')
    ax.plot(r_geo,z_geo,'-k')
    ax.plot(r_geo.T,z_geo.T,'-k')

    #lim
    ax.plot(r_lim,z_lim,'-g',label='lim')

    #sep
    ax.plot(r_sep,z_sep,'-m',label='sep')

    #ray
    ax.plot(r_ray,z_ray,'-b',label='ray')
    
    ax.set_ylabel('Z, m')
    ax.set_xlabel('R, m')
    ax.set_xlim(xmin=r_lim.min()*0.7,xmax=r_lim.max()*1.2)
    ax.set_ylim(ymin=z_lim.min()*1.2,ymax=z_lim.max()*1.2)
    ax.legend()
    fig.tight_layout()
    plt.show()
    
    fig.savefig(output_path_all+'//'+'ray_plot'+'.png', dpi=fig.dpi)
    
    plt.close()
    
