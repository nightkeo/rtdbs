# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 18:17:40 2022

@author: Kaine
"""

import numpy as np
import re
from scipy.interpolate import interp1d, interp2d
import scipy.integrate
import matplotlib.pylab as plt

def file_numbers(fp):
    """Generator to get numbers from a text file"""
    toklist = []
    while True:
        line = fp.readline()
        if not line: break
        # Match numbers in the line using regular expression
        pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?'
        toklist = re.findall(pattern, line)
        for tok in toklist:
            yield tok

def readg(f):
    """ Reads a G-EQDSK file

    Parameters
    ----------

    f = Input file. Can either be a file-like object,
        or a string. If a string, then treated as a file name
        and opened.

    Returns
    -------

    """
    #TODO: Make this into a class and precalculate psirz spline
    if isinstance(f, str):
        # If the input is a string, treat as file name
        with open(f) as fh: # Ensure file is closed
            return readg(fh) # Call again with file object

    # Read the first line, which should contain the mesh sizes
    desc = f.readline()
    if not desc:
        raise IOError("Cannot read from input file")

    s = desc.split() # Split by whitespace
    if len(s) < 3:
        raise IOError("First line must contain at least 3 numbers")

    idum = int(s[-3])
    nw = int(s[-2])
    nh = int(s[-1])

    # Use a generator to read numbers
    token = file_numbers(f)

    rdim   = float(next(token))
    zdim   = float(next(token))
    rcentr = float(next(token))
    rleft  = float(next(token))
    zmid   = float(next(token))

    rmaxis = float(next(token))
    zmaxis = float(next(token))
    simag  = float(next(token))
    sibry  = float(next(token))
    bcentr = float(next(token))

    current= float(next(token))
    simag  = float(next(token))
    xdum   = float(next(token))
    rmaxis = float(next(token))
    xdum   = float(next(token))

    zmaxis = float(next(token))
    xdum   = float(next(token))
    sibry  = float(next(token))
    xdum   = float(next(token))
    xdum   = float(next(token))

    # Read arrays
    def read_array(n, name="Unknown"):
        data = np.zeros([n])
        try:
            for i in np.arange(n):
                data[i] = float(next(token))
        except:
            raise IOError("Failed reading array '"+name+"' of size ", n)
        return data

    # read 2d array
    def read_2d(nw, nh, name="Unknown"):
        data = np.zeros([nw, nh])
        for j in np.arange(nh):
            for i in np.arange(nw):
                data[i,j] = float(next(token))
        return data

    fpol   = read_array(nw, "fpol")
    pres   = read_array(nw, "pres")
    ffprim = read_array(nw, "ffprim")
    pprime = read_array(nw, "pprime")
    psirz  = read_2d(nw, nh, "psirz")
    qpsi   = read_array(nw, "qpsi")

    # Read boundary and limiters, if present
    nbbbs  = int(next(token))
    limitr = int(next(token))

    if nbbbs > 0:
        rbbbs = np.zeros([nbbbs])
        zbbbs = np.zeros([nbbbs])
        for i in range(nbbbs):
            rbbbs[i] = float(next(token))
            zbbbs[i] = float(next(token))
    else:
        rbbbs = [0]
        zbbbs = [0]

    if limitr > 0:
        rlim = np.zeros([limitr])
        zlim = np.zeros([limitr])
        for i in range(limitr):
            rlim[i] = float(next(token))
            zlim[i] = float(next(token))
    else:
        rlim = [0]
        zlim = [0]

    # Construct R-Z mesh
    r = np.linspace(rleft, rleft + rdim, nw)
    z = np.linspace(zmid - 0.5*zdim, zmid + 0.5*zdim, nh)

    # Create dictionary of values to return
    result = {'nw': nw, 'nh':nh,        # Number of horizontal and vertical points
              'r':r, 'z':z,                     # Location of the grid-poinst
              'rdim':rdim, 'zdim':zdim,         # Size of the domain in meters
              'rcentr':rcentr, 'bcentr':bcentr, # Reference vacuum toroidal field (m, T)
              'rleft':rleft,                  # R of left side of domain
              'zmid':zmid,                      # Z at the middle of the domain
              'rmaxis':rmaxis, 'zmaxis':zmaxis,     # Location of magnetic axis
              'ssimag':simag, # Poloidal flux at the axis (Weber / rad)
              'ssibry':sibry, # Poloidal flux at plasma boundary (Weber / rad)
              'current':current,
              'psirz':psirz.T,    # Poloidal flux in Weber/rad on grid points
              'fpol':fpol,  # Poloidal current function on uniform flux grid
              'ffprim':ffprim, # derivative of poloidal flux
              'pres':pres,  # Plasma pressure in nt/m^2 on uniform flux grid
              'pprime':pprime, # derivative of pressure
              'qpsi':qpsi,  # q values on uniform flux grid
              'nbdry':nbbbs, 'bdry':np.vstack((rbbbs,zbbbs)).T, # Plasma boundary
              'limitr':limitr, 'lim':np.vstack((rlim,zlim)).T} # Wall boundary

    return result

def transform_equ(filename):
    #EFIT
    g=readg(filename)
    
    r_efit=g['bdry'][:,0]
    z_efit=g['bdry'][:,1]
    rlim_efit=g['lim'][:,0]
    zlim_efit=g['lim'][:,1]
    
    bcentr=g['bcentr'] #T toroidal magnetic field
    current=g['current']/1000 #kA plasma current
    a=np.abs(np.max(r_efit)-np.min(r_efit))/2 #m minor radius
    R0=np.abs(np.max(r_efit)+np.min(r_efit))/2 #m major radius
    bcentr=g['bcentr'] # T toroidal magnetic field
    rcentr=g['rcentr'] # Reference vacuum toroidal field (m, T)
    rmaxis=g['rmaxis'] # Location of magnetic axis
    zmaxis=g['zmaxis'] # Location of magnetic axis
    ssimag=g['ssimag'] # Poloidal flux at the axis (Weber / rad)
    ssibry=g['ssibry'] # Poloidal flux at plasma boundary (Weber / rad)
    
    fpol=g['fpol']
    qpsi=g['qpsi']
    
    
    #additional
    rgrid=g['r']
    zgrid=g['z']
    psirz=g['psirz']
    
    
    
    #ssimag=np.min(psirz)
    psirz = psirz-ssimag #strart from zero
    ssibry = (ssibry - ssimag)*0.999
    ssimag=0
    
    psirz = np.abs(psirz)
    ssibry=np.abs(ssibry)
    
    #interpolate
    rgrid_b1=np.linspace(0.8*np.min(rlim_efit), 1.2*np.max(rlim_efit), num=109)
    zgrid_b1=np.linspace(1.2*np.min(zlim_efit), 1.2*np.max(zlim_efit), num=109)
    psirz_b1 = interp2d(rgrid, zgrid, psirz, kind='cubic')(rgrid_b1,zgrid_b1)
    psirz_b1=psirz_b1-ssimag
    
    #interesting field
    rbdy_max=1.1*np.max(r_efit)
    rbdy_min=0.9*np.min(r_efit)
    
    zbdy_max=1.1*np.max(z_efit)
    zbdy_min=1.1*np.min(z_efit)
    
    #number of r and z of int field
    n_r_max=np.argmin(np.abs(rgrid-rbdy_max))
    n_r_min=np.argmin(np.abs(rgrid-rbdy_min))
    
    n_z_max=np.argmin(np.abs(zgrid-zbdy_max))
    n_z_min=np.argmin(np.abs(zgrid-zbdy_min))
    
    rgrid=rgrid[n_r_min:n_r_max]
    zgrid=zgrid[n_z_min:n_z_max]
    psirz=psirz[n_z_min:n_z_max,n_r_min:n_r_max]
    
    #interpolate
    rgrid_b=np.linspace(np.min(rgrid), np.max(rgrid), num=2222)
    zgrid_b=np.linspace(np.min(zgrid), np.max(zgrid), num=2222)
    psirz_b = interp2d(rgrid, zgrid, psirz, kind='cubic')(rgrid_b,zgrid_b)
    psirz_b=psirz_b-np.min(psirz_b)
    
    
    
    # psirz_b=(psirz_b-ssimag)/(ssibry-ssimag)
    # ssibry=1
    # ssimag=0
    
    
    
    psimesh=np.linspace(ssimag, ssibry, num=65)
    psimesh_b=np.linspace(ssimag, ssibry, num=500)
    
    qpsi_b=interp1d(psimesh, qpsi, kind='cubic')(psimesh_b)
    fpol_b=interp1d(psimesh, fpol, kind='cubic')(psimesh_b)
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rho=np.linspace(0, 1, num=100)
    
    rho_efit=np.zeros(len(psimesh_b))
    for i in range(1,len(psimesh_b)):
        rho_efit[i]=np.sqrt(2*np.pi*scipy.integrate.simps(qpsi_b[0:i], psimesh_b[0:i]))
    rho_efit=rho_efit/np.max(rho_efit)
    
    psipol=interp1d(rho_efit[1:len(rho_efit)], psimesh_b[1:len(rho_efit)], kind='cubic')(rho)
    qpsi=interp1d(rho_efit[1:len(rho_efit)], qpsi_b[1:len(rho_efit)], kind='cubic')(rho)
    g_eq_b=interp1d(rho_efit[1:len(rho_efit)], fpol_b[1:len(rho_efit)], kind='cubic')(rho)
    
     ##GRAPH
    
     #     plt.rcParams.update({'font.size':14})
     #     plt.rc('font',family='serif')
    
     #Poloidal magnetic flux
    fig,ax=plt.subplots(figsize=(8,10))
    ax.set_aspect('equal')
    ax.set_title('Poloidal magnetic flux, Wb/rad')
    ax.set_xlabel('R, m')
    ax.set_ylabel('Z, m')
    tcf=plt.contour(rgrid_b,zgrid_b,psirz_b,psipol,cmap='jet')
    tcf=plt.contourf(rgrid_b,zgrid_b,psirz_b,100,cmap='jet')
    
    #boundary
    qq_x=list()
    qq_y=list()
    check = list()
    cs=plt.contour(rgrid_b,zgrid_b,psirz_b,psipol)
    for i in range(0,len(psipol)):
        vin = cs.collections[i].get_paths()
        for j in range(0,len(vin)):
            w = cs.collections[i].get_paths()[j].vertices
            if abs(w[0,1]-w[-1,1])+abs(w[0,0]-w[-1,0]) < 1e-10:
                x = w[:,0]
                y = w[:,1]
                plt.plot(x,y,'.y')
                plt.plot(x[0],y[0],'ob')                
                qq_x.append(x)
                qq_y.append(y)
                break
          
                 
    plt.plot(rlim_efit,zlim_efit,'-g',label='limiter')
    plt.plot(r_efit,z_efit,'-ok',label='boundary')
    plt.plot(rmaxis,zmaxis,'or')
     
    plt.colorbar(tcf)
    fig.tight_layout()
     
    plt.show()
    
    th_new2=np.linspace(0.0, np.pi, num=2000)
    th_new=np.linspace(-np.pi, np.pi, num=1000)
    r_g=np.zeros((len(qq_x),len(th_new)))
    z_g=np.zeros((len(qq_x),len(th_new)))

    
    for j in range(0,len(qq_x)):
        x=0
        y=0
        d=0
        lpol=0
        th1=0
        th2=0
        th=0
        x_new=0
        y_new=0
        ii=0
        a=0
        b=0
        x_new2=0
        y_new2=0
        buf_x=0
        buf_y=0
        
        x=np.flipud(qq_x[j])
        y=np.flipud(qq_y[j])
        
        ii=np.argsort(np.abs(y-zmaxis))
        a = max(ii[0:2])#a,b=ii[0:2]
        if len(x)/a > 2:
            for a in ii[2:]:
                if len(x)/a < 2: break
        x_2=x[0:-1].copy()
        y_2=y[0:-1].copy()
        x_split=np.array_split(x_2,[a])
        x_2 = np.hstack([x_split[1],x_split[0],x_split[1][0]])
        y_split=np.array_split(y_2,[a])
        y_2 = np.hstack([y_split[1],y_split[0],y_split[1][0]])
        
        #interp
    
        jj = np.argsort(np.abs(y_2-zmaxis))
        b = max(jj[0:2])
        if len(x_2)/b < 1.5:
            for b in jj[2:]:
                if len(x_2)/b > 1.5: break
        x_c1, x_c2= np.array_split(x_2,[b])
        y_c1, y_c2 = np.array_split(y_2,[b])
            
            
        d=np.sqrt( np.diff(x_c1)**2+np.diff(y_c1)**2 )
        lpol=np.hstack((0,np.cumsum(d)))
        th1=(np.pi/np.max(lpol)*lpol - np.pi)
    
        d=np.sqrt( np.diff(x_c2)**2+np.diff(y_c2)**2 )
        lpol=np.hstack((0,np.cumsum(d)))+d[0]
        th2=np.pi/np.max(lpol)*lpol
        th = np.hstack([th1, th2])  
    
    
        th[-1] = np.pi
    
    
        x_new = interp1d(th,x_2,kind='linear')(th_new)
        y_new = interp1d(th,y_2,kind='linear')(th_new)
        r_g[j,:]=x_new.copy()
        z_g[j,:]=y_new.copy()
    
    r_g2=r_g.T
    z_g2=z_g.T
    
    r_g2[:,0]=np.ones(len(th_new))*rmaxis
    z_g2[:,0]=np.ones(len(th_new))*zmaxis
    
    th_eq=np.linspace(-np.pi, +np.pi, num=151)
    rho_eq=np.linspace(0.0, 1.0, num=41)
    #rho_new=np.linspace(0.0, 1.0, num=len(rho))
    rho_new=np.linspace(0.0, 1.0, num=len(r_g))
        
    r_geo = interp2d(rho_new,th_new, r_g2, kind='cubic')(rho_eq,th_eq)
    z_geo = interp2d(rho_new,th_new, z_g2, kind='cubic')(rho_eq,th_eq)
    r_stich = (r_geo[0].copy()+r_geo[-1].copy())/2
    z_stich = (z_geo[0].copy()+z_geo[-1].copy())/2
    r_geo[0] = r_geo[-1] = r_stich
    z_geo[0] = z_geo[-1] = z_stich
    
    rho_eq=np.linspace(0.0, 1.0, num=41)
    psipol_eq=interp1d(rho, psipol, kind='cubic')(rho_eq)
    psipol_eq=psipol_eq-psipol_eq[0]
    q_eq=interp1d(rho, qpsi, kind='cubic')(rho_eq)
    g_eq=interp1d(rho, g_eq_b, kind='cubic')(rho_eq)
         
    #close graph (graph is used for equilibrium transforming)
    plt.close()
        
    # Create dictionary of values to return
    result = {'rho_eq': rho_eq,
              'psipol':psipol_eq,
              'q_eq':q_eq,
              'th_eq':th_eq,
              'g_eq':g_eq,
              'R_geo':r_geo,
              'Z_geo':z_geo,
              'R_grid':rgrid_b1,
              'Z_grid':zgrid_b1,
              'PsiRZ':psirz_b1
    } 
    
    return result

def derivatives_geom(r_geo, z_geo, rho_eq, th_eq):
    #berem vse prostranstvennye proisvodnye
    dR_dTH, dR_dRHO = np.gradient(r_geo,th_eq, rho_eq)
    dZ_dTH, dZ_dRHO= np.gradient(z_geo,th_eq, rho_eq)
    dR_dRHO[:,0] = dR_dRHO[:,1]
    dR_dRHO[:,-1] = dR_dRHO[:,-2]
    dR_dTH[0,:] = dR_dTH[1,:]
    dR_dTH[-1,:] = dR_dTH[-2,:]
    dZ_dRHO[:,0] = dZ_dRHO[:,1]
    dZ_dRHO[:,-1] = dZ_dRHO[:,-2]
    dZ_dTH[0,:] = dZ_dTH[1,:]
    dZ_dTH[-1,:] = dZ_dTH[-2,:]
    d2R_dRHOdTH = np.gradient(dR_dTH,th_eq, rho_eq)[1]
    d2Z_dRHOdTH = np.gradient(dZ_dTH,th_eq, rho_eq)[1]
    d2R_dRHOdTH[:,0] = d2R_dRHOdTH[:,1]
    d2R_dRHOdTH[:,-1] = d2R_dRHOdTH[:,-2]
    d2Z_dRHOdTH[:,0] = d2Z_dRHOdTH[:,1]
    d2Z_dRHOdTH[:,-1] = d2Z_dRHOdTH[:,-2]
    result = {'dR_dTH':dR_dTH,
              'dR_dRHO':dR_dRHO,
              'dZ_dTH':dZ_dTH,
              'dZ_dRHO':dZ_dRHO,
              'd2R_dRHOdTH':d2R_dRHOdTH,
              'd2Z_dRHOdTH':d2Z_dRHOdTH
             }
    return result 

def limiter_calc():
    # limiter
    ew1 = [0.121, 0.12191139240506328, 0.12282278481012658, 0.12373417721518987, 0.12464556962025317, 0.1270886075949367,
               0.13050632911392407, 0.1339240506329114, 0.13734177215189874, 0.14141772151898735, 0.1477974683544304,
               0.15417721518987343, 0.16055696202531647, 0.1669367088607595, 0.17369620253164558, 0.18053164556962026,
               0.18736708860759496, 0.19420253164556964, 0.20154430379746838, 0.20951898734177216, 0.21749367088607596,
               0.22546835443037974, 0.23343037974683545, 0.24117721518987342, 0.24892405063291143, 0.2566708860759494,
               0.26441772151898735, 0.27216455696202535, 0.2799113924050633, 0.2876582278481013, 0.29540506329113925,
               0.3052405063291139, 0.32050632911392407, 0.3357721518987342, 0.35103797468354436, 0.3663037974683544,
               0.3742784810126582, 0.38134177215189874, 0.3884050632911393, 0.39546835443037975, 0.4023037974683545,
               0.4089113924050633, 0.4155189873417721, 0.422126582278481, 0.42870886075949366, 0.4350886075949367,
               0.44146835443037974, 0.4478481012658228, 0.4542278481012659, 0.46044303797468356, 0.4665949367088608,
               0.47274683544303797, 0.4788987341772152, 0.4848987341772152, 0.49059493670886073, 0.49629113924050633,
               0.5019873417721519, 0.5076835443037975, 0.5131645569620253, 0.5186329113924051, 0.5241012658227848,
               0.5295696202531647, 0.5346582278481014, 0.5394430379746836, 0.5442278481012659, 0.5490126582278482,
               0.5537594936708862, 0.5583164556962026, 0.562873417721519, 0.5674303797468354, 0.5719873417721518,
               0.5761898734177215, 0.5802911392405063, 0.5843924050632912, 0.588493670886076, 0.5927721518987341,
               0.5973291139240506, 0.601886075949367, 0.6064430379746835, 0.611, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62, 0.62,
               0.62, 0.62, 0.62, 0.611, 0.6054117647058823, 0.5998235294117646, 0.5942352941176471, 0.5888823529411764,
               0.5838529411764706, 0.5788235294117646, 0.5737941176470588, 0.5682941176470588, 0.5627058823529412,
               0.5571176470588236, 0.5514558823529412, 0.5455882352941177, 0.5397205882352941, 0.5338529411764706,
               0.5274117647058824, 0.5207058823529412, 0.514, 0.507264705882353, 0.5002794117647059, 0.4932941176470588,
               0.48630882352941174, 0.47902941176470587, 0.47148529411764706, 0.46394117647058825, 0.45639705882352943,
               0.4485882352941177, 0.44076470588235295, 0.4329411764705882, 0.4250147058823529, 0.41691176470588237,
               0.4088088235294118, 0.4007058823529412, 0.3921617647058824, 0.3835, 0.37483823529411764, 0.3661176470588235,
               0.3571764705882353, 0.3482352941176471, 0.33929411764705886, 0.32982352941176474, 0.32004411764705887,
               0.3102647058823529, 0.3005, 0.2909999999999999, 0.2815, 0.2719999999999999, 0.26249999999999996, 0.253, 0.2435,
               0.23400000000000004, 0.22425, 0.21447058823529416, 0.20469117647058824, 0.19535294117647062, 0.1869705882352941,
               0.17858823529411766, 0.17020588235294115, 0.16223529411764706, 0.1544117647058823, 0.14658823529411763,
               0.13933823529411762, 0.1351470588235294, 0.13095588235294114, 0.12676470588235292, 0.12435294117647058,
               0.12323529411764705, 0.12211764705882353, 0.121, 0.1215, 0.12361111111111112, 0.12572222222222224,
               0.12783333333333333, 0.12994444444444445, 0.13205555555555556, 0.13416666666666668, 0.13627777777777778,
               0.1383888888888889, 0.1405, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.1405,
               0.139, 0.1375, 0.136, 0.1345, 0.133, 0.1315, 0.13, 0.1285, 0.127, 0.1265, 0.1265, 0.1265, 0.1265, 0.1265,
               0.1265, 0.1265, 0.1265, 0.1265, 0.1265, 0.141, 0.141, 0.121]
    ew2 = [0.44, 0.4491139240506329, 0.4582278481012658, 0.4673417721518987, 0.47645569620253164, 0.4833417721518987,
           0.48881012658227846, 0.4942784810126582, 0.499746835443038, 0.5049113924050633, 0.5090126582278481, 0.5131139240506329,
           0.5172151898734177, 0.5213164556962026, 0.5238987341772152, 0.5261772151898735, 0.5284556962025316, 0.5307341772151899,
           0.5321012658227848, 0.5323291139240507, 0.5325569620253164, 0.5327848101265823, 0.5329493670886076, 0.5320379746835443,
           0.5311265822784811, 0.5302151898734178, 0.5293037974683544, 0.5279367088607595, 0.5263417721518987, 0.524746835443038,
           0.5231518987341772, 0.5205443037974684, 0.5153037974683544, 0.5100632911392405, 0.5048227848101265, 0.49958227848101266,
           0.49575949367088606, 0.4921139240506329, 0.4884683544303797, 0.4848227848101266, 0.48083544303797465,
           0.47650632911392404, 0.47217721518987343, 0.4678481012658228, 0.4634683544303798, 0.4586835443037975,
           0.4538987341772152, 0.44911392405063294, 0.4443291139240506, 0.43921518987341773, 0.43397468354430374,
           0.42873417721518986, 0.4234936708860759, 0.4180253164556962, 0.4121012658227848, 0.40617721518987343,
           0.400253164556962, 0.3943291139240507, 0.3881898734177215, 0.3820379746835443, 0.37588607594936707,
           0.3697341772151898, 0.3635822784810126, 0.3574303797468354, 0.3512784810126582, 0.34512658227848103,
           0.33886075949367084, 0.3320253164556962, 0.3251898734177215, 0.3183544303797468, 0.31151898734177214,
           0.30468354430379746, 0.2978481012658228, 0.29101265822784805, 0.28417721518987343, 0.27716455696202535,
           0.26987341772151896, 0.2625822784810127, 0.2552911392405064, 0.248, 0.23, 0.1788888888888889, 0.12777777777777777,
           0.07666666666666666, 0.025555555555555554, -0.02555555555555558, -0.07666666666666669, -0.1277777777777778,
           -0.1788888888888889, -0.23, -0.248, -0.25694117647058823, -0.26588235294117646, -0.27482352941176474,
           -0.2835294117647059, -0.29191176470588237, -0.3002941176470588, -0.3086764705882353, -0.3170588235294118,
           -0.32544117647058823, -0.33382352941176474, -0.3419852941176471, -0.34952941176470587, -0.35707352941176473,
           -0.36461764705882355, -0.37216176470588236, -0.37970588235294117, -0.38725000000000004, -0.39476470588235296,
           -0.4020294117647059, -0.4092941176470588, -0.41655882352941176, -0.4233823529411765, -0.4298088235294118,
           -0.43623529411764705, -0.44266176470588237, -0.4485588235294118, -0.4544264705882353, -0.46029411764705885,
           -0.4659558823529412, -0.471264705882353, -0.47657352941176473, -0.4818823529411765, -0.4865294117647059, -0.491,
           -0.4954705882352941, -0.49976470588235294, -0.5033970588235294, -0.5070294117647058, -0.5106617647058823,
           -0.5137647058823529, -0.5165588235294117, -0.5193529411764706, -0.5221029411764706, -0.5240588235294118,
           -0.526014705882353, -0.5279705882352942, -0.5295294117647059, -0.5306470588235295, -0.531764705882353,
           -0.5328823529411765, -0.5327500000000001, -0.5324705882352941, -0.5321911764705882, -0.5311176470588236,
           -0.5283235294117647, -0.5255294117647059, -0.5227352941176471, -0.5182941176470588, -0.513264705882353,
           -0.508235294117647, -0.5029411764705882, -0.496235294117647, -0.4895294117647058, -0.48282352941176465,
           -0.47352941176470587, -0.4623529411764706, -0.45117647058823535, -0.44, -0.44, -0.44, -0.44, -0.44, -0.44, -0.44,
           -0.44, -0.44, -0.44, -0.44, -0.44, -0.42444444444444446, -0.4088888888888889, -0.3933333333333333, -0.37777777777777777,
           -0.3622222222222222, -0.3466666666666667, -0.33111111111111113, -0.31555555555555553, -0.3, -0.3, -0.3, -0.3, -0.3,
           -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.23333333333333334, -0.16666666666666666, -0.09999999999999998,
           -0.033333333333333326, 0.033333333333333326, 0.10000000000000003, 0.16666666666666669, 0.23333333333333334, 0.3,
           0.3, 0.44, 0.44]
    
    return ew1,ew2


