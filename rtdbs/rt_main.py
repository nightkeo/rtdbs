import numpy as np
import time
import os
from shutil import copyfile

from rtdbs.wave_data_preparation import rotate_wave_vector, create_initial_data_file

from rtdbs.equ_read import readg
from rtdbs.nerv_geom import transform_equ
from rtdbs.equ_transform import equ_transform
from rtdbs.equ_approx import equ_approx, save_geom_koeff_pol, find_rho_sep, find_rho_theta_with_R_in_Z_in, find_dpsi_drho, save_dpsi_drho_koeff_pol

from rtdbs.TS_approx import TS_approx,save_dens

from rtdbs.plot_graph import equ_plot, dpsi_drho_graph, dens_plot, ray_plot

def rt_main(
            alpha,r_start_vect,n_start_vect,R_start,frequency,moda,InversT,
            filename_equ,output_path_all,
            r_ts,z_ts,ne_ts,dne_ts
            ):
    
    #path for exe
    exe_path='.//exe//'
    
    
    
    #prepare wave data
    ###############################################################################
    n_start_vect_rot,r_start_vect_rot,R_start_check,phi_rot=rotate_wave_vector(r_start_vect,n_start_vect,R_start)
    
    #create initial data file
    create_initial_data_file(exe_path,alpha,r_start_vect,n_start_vect,frequency,moda,InversT)
    ###############################################################################
    
    
    
    #read equilibrium data
    ###############################################################################
    g=readg(filename_equ)
    geom=equ_transform(g,nrho=41,ntheta=151)
    
    #read and transform equilibrium data with old library (new [currently] - WRONG!!!)
    geom_data_old=transform_equ(filename_equ)
    geom['r_geo']=geom_data_old['R_geo']
    geom['z_geo']=geom_data_old['Z_geo']
    
    #normalized minor radius
    rho_a=(np.max(geom['r_geo'],axis=0)-np.min(geom['r_geo'],axis=0))/2
    geom['rho_a']=rho_a/rho_a[-1]
    
    #EFIT files rcentr is WRONG!!!
    geom['rcentr']=0.36
    
    ##############################################################################
    
    
    
    #approx equilibrium
    ###############################################################################
    geom_approx=equ_approx(geom,
                    max_it=100_000,
                    rho_a_sec=np.linspace(0.1,1,10),
                    theta_sec=np.linspace(-np.pi/3,np.pi/3,50),
                    alpha=alpha
                    )
    #save data
    save_geom_koeff_pol(geom_approx,
                                    geom,
                                    filename_geom=exe_path+'/input_test.dat',
                                    filename_Btor=exe_path+'/Btor.txt'
                                    )
    ##############################################################################
    
    
    
    ###############################################################################
    #find rho_sep at max R approx
    rho_sep=find_rho_sep(geom_approx,alpha,max_it=100_000)
    print('rho_sep=',rho_sep)
    #find rho_start and theta_start
    rho_start,theta_start=find_rho_theta_with_R_in_Z_in(geom_approx,alpha,r_start_vect_rot[0],r_start_vect_rot[2],max_it=10_000)
    ###############################################################################
    
    
    
    #dpsi/drho поиск
    ###############################################################################
    dpsi_drho_data=find_dpsi_drho(g,geom_approx,alpha,rho_start,rho_min=0.1,max_it=100_000)
    
    #save data
    save_dpsi_drho_koeff_pol(dpsi_drho_data,rho_start,theta_start,n_start_vect_rot,output_filename=exe_path+'/input_test2.dat')
    ###############################################################################
    
    
    
    
    #dens approx
    ###############################################################################
    dens=TS_approx(geom_approx,r_ts,z_ts,ne_ts,dne_ts,alpha,rho_start,max_it=10_000)
    #save dens data
    save_dens(dens,filename=exe_path+'/input_test3.dat')
    ###############################################################################
    
    
    
    
    #plot
    ###############################################################################
    arg_list=['full_geo_vs_full_approx','full_geo','full_approx','full_approx_alpha']
    for arg in arg_list:
        equ_plot(geom,geom_approx,dens,output_path_all,plot_type=arg)
    dpsi_drho_graph(dpsi_drho_data,output_path_all)
    dens_plot(dens,output_path_all)
    ###############################################################################
    
    
    
    #START!!!
    ##############################################################################
    dir = os.path.abspath(os.curdir)
    os.chdir(dir+exe_path)
    os.startfile('Ray_tracing_1.exe')
    os.chdir(dir)
    ##############################################################################
    
    # copy and delete files
    #############################################################################
    
    #wait 10 sec before file copy and deleting
    time.sleep(30)
    
    filename_list=[
        'input_test.dat',
        'input_test2.dat',
        'input_test3.dat',
        'raytracing.log',
        'Btor.txt',
        'cutoff.txt',
        'initial_data.txt',
        'ne_log_kuk.txt',
        'out1.txt',
        'out3.txt',
        'ray_kuk.txt'
                    ]
    
    for filename in filename_list:
        cfn=exe_path+'//'+filename
        if os.path.exists(cfn):
            print('Delete and copy '+cfn)
            copyfile(cfn,output_path_all+'//'+filename)
            os.remove(cfn)
    
    ##############################################################################
    
    
    
    #open ray_data
    #############################################################################
    ray_filename=output_path_all+'//'+'ray_kuk.txt'
    if os.path.exists(ray_filename):
        ray_data_buf=np.loadtxt(ray_filename,skiprows=2)
        try:
            ray_data = {'rho':ray_data_buf[:,0],
                      'X':ray_data_buf[:,1],
                      'Y':ray_data_buf[:,2],
                      'Z':ray_data_buf[:,3],
                      'R':ray_data_buf[:,4],
                      'Npar':ray_data_buf[:,5],
                      'Nper':ray_data_buf[:,6],
                      'Bpol/Btor':ray_data_buf[:,7],
                      'Step':ray_data_buf[:,8]
                      }                
            #plot ray
            ray_plot(geom,ray_data,output_path_all)
        except:
            print('bad ray.txt !!!')
    ##############################################################################

    




