# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 19:57:29 2024

@author: user
"""

import numpy as np

def rotate_wave_vector(r_start_vect,n_start_vect,R_start):
    xs=r_start_vect[0]
    ys=r_start_vect[1]
    zs=r_start_vect[2]

    nxs=n_start_vect[0]
    nys=n_start_vect[1]
    nzs=n_start_vect[2]
    
    #решаем квадр ур-е
    diskr=4*(xs*nxs+ys*nys)*(xs*nxs+ys*nys)-4*(nxs*nxs+nys*nys)*(xs*xs+ys*ys-R_start*R_start)

    # h1=(-2*xs*nxs-2*ys*nys+np.sqrt(diskr))/(2*(nxs*nxs+nys*nys))
    h2=(-2*xs*nxs-2*ys*nys-np.sqrt(diskr))/(2*(nxs*nxs+nys*nys))

    #подходит второе
    X_start=xs+h2*nxs
    Y_start=ys+h2*nys
    Z_start=zs+h2*nzs

    #далее поворачиваем систему координат в торроидальном направление чтобы Y_start=0
    #при таком повороте координаты старта превращаются в X_start=R_start=0.62 Y_start=0 Z_start=zs+h2*nzs;
    Z_start=zs+h2*nzs
    phi_rot=np.arctan(Y_start/X_start)    
    R_start_check=X_start*np.cos(phi_rot)+Y_start*np.sin(phi_rot)
    #для проверки что X'_start действительно=R_start
    #R_start=X'_start поскольку дальше в программе используется именно R_start

    #a компоненты волнового вектора преобразуются как:
    nxs_rot=nxs*np.cos(phi_rot)+nys*np.sin(phi_rot)
    nys_rot=-nxs*np.sin(phi_rot)+nys*np.cos(phi_rot)
    nzs_rot=nzs
    
    #pack to vectors
    n_start_vect_rot=[nxs_rot,nys_rot,nzs_rot]
    r_start_vect_rot=[X_start,Y_start,Z_start]
    
    return n_start_vect_rot,r_start_vect_rot,R_start_check,phi_rot

def create_initial_data_file(exe_path,alpha,r_start_vect,n_start_vect,frequency,moda,InversT):
    f = open(exe_path+'/initial_data.txt', 'w')

    f.write(str(alpha)+'\n')

    for k in r_start_vect:
        f.write(str(k)+'\n')

    for k in n_start_vect:
        f.write(str(k)+'\n')

    f.write(str(frequency)+'E9\n')

    f.write(str(moda)+'\n')

    f.write(str(InversT)+'\n')

    f.write(' \n')
    f.write('//alpha (масштабный коэффицент )\n')
    f.write('//x_start(m)\n')
    f.write('//y_start(m)\n')
    f.write('//z_start(m)\n')
    f.write('//Nx_start\n')
    f.write('//Ny_start\n')
    f.write('//Nz_start\n')
    f.write('//F(hz)(линейная)\n')
    f.write('//moda  0-обыкновенная 1-необыкновенная \n')
    f.write('//InversT  в какую сторону интегрировать 0 или 1 \n')

    f.close()


























