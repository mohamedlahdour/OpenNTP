#!usr/bin/python3
# -*- coding: utf-8 -*-  
import json
import numpy 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as mpatch
from matplotlib.patches import Ellipse, Arc
import mpl_toolkits.mplot3d.art3d as art3d
import math 
import matplotlib.pyplot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
 This function is used to read 1D geometry from the input file
 and draw it in a figure using Matplotlib
"""
def visualization1D(self):
    Geometry_type = open('app/link/script01.py', "r" ).read()
    if Geometry_type in ['Slab Geometry']:
        filename = open('app/link/script.dir', "r" ).read()
        fmm_id = []
        assembly = []
        pin = []
        regmat = []
        nom    = []
        width_x = []
        width_y = []
        with open(filename) as json_data:
            data = json.load(json_data)
            nmat = data['data']['parameter']['Total number of materials']
            npc = data['data']['parameter']['Total number of pin cells']  
            na  = data['data']['parameter']['Total number of assemblies'] 
            core = data['data']['parameter']['Core']

            for i in range(na):
                assembly.append(data['data']['Assemblies'][i]['assembly'])

            for i in range(npc):
                pin.append(data['data']['PinCells'][i]['mat_fill'])

            Height = 1
            width = len(core)*len(pin[0])*len(assembly[0])
            fmm_id =numpy.zeros((width),dtype='i')
            i=0
            for j in range(len(core)):
                for k in range(len(assembly[0])):
                    for m in range(len(pin[0])):
                        fmm_id[i] = pin[assembly[core[j]-1][k]-1][m]
                        i+=1
            nx = len(core)
            ny = 1
           
            nxx = len(assembly[0])
            nyy = 1
            NX = nxx*nx
            NY = ny*nyy
            width_x.append(data['data']['PinCells'][0]['width'])
            width_y = [2]
            nx =  len(width_x[0])*NX
            ny =  1
            xcm  = width_x[0]*NX
            ycm  = width_y[0]*NY
            for j in range(nmat):
                nom.append(data['data']['materials'][j]['name'])
            for i in range(npc):
                regmat.append(data['data']['PinCells'][i]['mat_fill'])
            fig, ax = plt.subplots(figsize=(5, 4), dpi=100)
            widthx = [0]
            widthy = [0]
            somx = [0]
            somy = [0]
            for i in range(nx):
                widthx.append(xcm[i])
                somx.append(sum(widthx))

            widthy = [2.0]
            somy = [2.0]
            ycm = [2.0]
            rectangles = []
            red_patch = []
            mx = somx[0]
            my = somy[0]
            
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.77, box.height])
            for i in range(nx):
                for j in range(ny):
                    rectangles.append(mpatch.Rectangle((somx[i],somy[j]), xcm[i], ycm[j],linewidth=0.0,edgecolor='k'))
            colr = ['#ff6666','#ffcc99', '#99ff99', '#66b3ff','#c2c2f0','#ffb3e6', 
                             '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
            for i in range(nmat):
                red_patch.append(mpatches.Patch(color=colr[i], label=nom[i]))
            
            n=0
            for r in rectangles:
                ax.add_artist(r) 
                r.set_facecolor(color=colr[int(fmm_id[n])-1])

                rx, ry = r.get_xy()
                cx = rx + r.get_width()/2.0
                cy = ry + r.get_height()/2.0
                n+=1

            ax.set_ylim((1,5))
            ax.set_xlim((min(somx), max(somx)))
            #ax.set_xticklabels([])
            ax.set_yticklabels([]) 
            ax.set_xlabel('X [cm]')
            #ax.set_ylabel('Y [cm]')
            ax.set_title('Color by Materials') 
            clb = plt.legend(handles=red_patch,loc='center left',title="Materials",
                         fontsize='small',bbox_to_anchor=(1, 0.5))
            plt.show()
    else:
        filename = open('app/link/script.dir', "r" ).read()
        name    = [] 
        ray     = [] 
        with open(filename) as json_data:
            data = json.load(json_data)
            nmat = data['data']['parameter']['Total number of materials']
            npc = data['data']['parameter']['Total number of pin cells'] 
        for j in range(nmat):
            name.append(data['data']['materials'][j]['name'])

        for i in range(npc):
            ray.append(data['data']['PinCells'][i]['ray'][0])

        fig, ax = plt.subplots(figsize=(5, 4), dpi=100)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        red_patch = []
        levels = [0.0]+ray
        xlist = np.linspace(-max(ray)-2, max(ray)+2, 100)
        ylist = np.linspace(-max(ray)-2, max(ray)+2, 100)
        X, Y = np.meshgrid(xlist, ylist)
        Z = np.sqrt(X ** 2 + Y ** 2 )
        contour = plt.contour(X, Y, Z, levels, colors='k')
        plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
        plt.axes().set_aspect("equal")
        colr = ['#ff6666','#ffcc99', '#99ff99', '#66b3ff','#c2c2f0','#ffb3e6', 
                   '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
        c = ('#ff6666','#ffcc99', '#99ff99', '#66b3ff','#c2c2f0','#ffb3e6', 
                   '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6')
        contour_filled = plt.contourf(X, Y, Z, levels, colors=c)
        for i in range(nmat):
            red_patch.append(mpatch.Patch(color=colr[i], label=name[i]))
        plt.title('Color by Materials')
        plt.xlabel('R [cm]')
        plt.ylabel('R [cm]')
        plt.legend(handles=red_patch, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
 This function is used to read 2D geometry from the input file
 and draw it in a figure using Matplotlib
"""
def visualization2D(self):
    filename = open('app/link/script.dir', "r" ).read()
    fmm_id =[[]]
    assembly = []
    pin = []
    regmat = []
    nom    = []
    width_x = []
    width_y = []
    with open(filename) as json_data:
        data = json.load(json_data)
        nmat = data['data']['parameter']['Total number of materials']
        npc = data['data']['parameter']['Total number of pin cells']  # Number of Pin Cell
        na  = data['data']['parameter']['Total number of assemblies'] 
        core = data['data']['parameter']['Core']

        for i in range(na):
            assembly.append(data['data']['Assemblies'][i]['assembly'])

        for i in range(npc):
            pin.append(data['data']['PinCells'][i]['mat_fill'])
        Height = len(core[0])*len(pin[0][0])*len(assembly[0])
        width = len(core[0])*len(pin[0][0])*len(assembly[0])
        fmm_id =numpy.zeros(((Height),(width)),dtype='i')

        i1=0
        for j in range(len(core[0])):
            for k in range(len(assembly[0])):
                for m in range(len(pin[0])):
                    i2=0
                    for i in range(len(core[0])):
                        for l in range(len(assembly[0])):
                            for n in range(len(pin[0])):
                                fmm_id[i1][i2] = pin[ assembly[ core[j][i] -1][k][l] -1][m][n]
                                i2+=1
                    i1+=1
        nx = len(np.amax(data['data']['parameter']['Core'],axis=0))
        ny = len(np.amax(data['data']['parameter']['Core'],axis=1))
            
        nxx = len(assembly[0])
        nyy = len(assembly[0])
        NX = nxx*nx
        NY = ny*nyy
        width_x.append(data['data']['PinCells'][0]['width_x'])
        width_y.append(data['data']['PinCells'][0]['width_y'])
        nx =  len(width_x[0])*NX
        ny =  len(width_y[0])*NY
        xcm  = width_x[0]*NX
        ycm  = width_y[0]*NY
        for j in range(nmat):
            nom.append(data['data']['materials'][j]['nom'])
        for i in range(npc):
            regmat.append(data['data']['PinCells'][i]['mat_fill'])
        fig, ax = plt.subplots()
        widthx = [0]
        widthy = [0]
        somx = [0]
        somy = [0]
        for i in range(nx):
            widthx.append(xcm[i])
            somx.append(sum(widthx))
        for j in range(ny):
            widthy.append(ycm[j])
            somy.append(sum(widthy))
        rectangles = []
        red_patch = []
        mx = somx[0]
        my = somy[0]

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.77, box.height])
        for i in range(nx):
            for j in range(ny):
                rectangles.append(mpatch.Rectangle((somx[i],somy[j]), xcm[i], ycm[j],linewidth=0.0,edgecolor='k'))
        colr = ['#ff6666','#ffcc99', '#99ff99', '#66b3ff','#c2c2f0','#ffb3e6', 
                   '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
        for i in range(nmat):
            red_patch.append(mpatches.Patch(color=colr[i], label=nom[i]))
        n=nx-1
        m=0
        for r in rectangles:
            ax.add_artist(r) 
            r.set_facecolor(color=colr[int(fmm_id[n][m])-1])

            rx, ry = r.get_xy()
            cx = rx + r.get_width()/2.0
            cy = ry + r.get_height()/2.0
            n-=1
            if (n<0):
                n=nx-1
                m+=1
        ax.set_ylim((min(somy), max(somy)))
        ax.set_xlim((min(somx), max(somx)))
        #ax.set_xticklabels([])
        #ax.set_yticklabels([]) 
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Y [cm]')
        ax.set_title('Color by Materials') 
        clb = plt.legend(handles=red_patch, loc='center left', title="Materials", fontsize='small', bbox_to_anchor=(1, 0.5))
        plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
 This function is used to plot the Level Symmetric Gaussian Quadrature Sets
"""
def P_Ordinate(self):
    file_name = 'app/Output/Ordinates.h'
    positions = np.loadtxt(file_name)
    ss = len(positions)
    red_patch = []
    if (ss == 4):
        x = np.array( [positions[0]] )
        y = np.array( [positions[1]] )
        z = np.array( [positions[2]] )
        w = np.array( [positions[3]] )
    else:
        x = positions[:,0]
        y = positions[:,1]
        z = positions[:,2]
        w = positions[:,3]
    ww = []
    i=1
    for a in w:
        if a not in ww:
            ww.append(a)
            red_patch.append(mpatch.Patch(label="W%s= %s" %(i,a)))
            i+=1

    w = 50. * w / w.max()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d([0, 1])
    ax.set_ylim3d([0, 1])
    ax.set_zlim3d([0, 1])
    ax.set_autoscale_on(False)
    n = 100
    zr = []
    for i in range(0, len(x)):
        add = True
        for j in range(0, len(zr)):
            if (zr[j][0] == z[i]):
                add = False
        if (add==True):
            zr.append( [z[i], ((x[i]**2.)+(y[i]**2.))**0.5] )
    circle = []
    #box = Bbox(((-1,-1),(1,1)))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    for i in range(0, len(zr)):
        r = zr[i][1]
        d = 2*r
        zed = zr[i][0]
        circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="white", alpha=0.2))
        ax.add_patch(circle[-1])
        art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="x")
        circle[-1].set_clip_on(True)
    
        circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="white", alpha=0.2))
        ax.add_patch(circle[-1])
        art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="y")
        circle[-1].set_clip_on(True)
    
        circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="white", alpha=0.2))
        ax.add_patch(circle[-1])
        art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="z")
        circle[-1].set_clip_on(True)
    for i in range(0, len(x)):
        custom = ( x[i], y[i], z[i] )
        ax.scatter(x[i],y[i],z[i],color=custom,s=w[i]*5,zorder=2)
    ax.view_init(elev=20., azim=25)
    ax.set_xlabel(u'\u00B5')
    ax.set_ylabel(u'\u03B7') 
    ax.set_zlabel(u'\u03BE')
    if (len(positions) == 4):
        tt = '2'
    else:
        delta = 4 +4*8*len(positions)
        tt = str(int((-1+ math.sqrt(delta))/2))
    plt.title(r'$S_{'+''.join([i for i in tt[0:]])+'}$')
    #fig.savefig(file_name+".png", dpi=100, facecolor='white')
    plt.legend(handles=red_patch,handlelength=False, handletextpad=0, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
This function allowing to plot the pin power distribution for a 1d or 2D geometry
"""
def powerpf(self):
    data = np.loadtxt('app/Output/PinPower.h')
    max_columns = len(data[0]) - 1
    max_rows = len(data)
    x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
    y = [data[0][colnum + 1] for colnum in range(max_columns)]
    z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    x,y = np.meshgrid(x, y)

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('X [cm]')                         
    ax.set_ylabel('Y [cm]') 
    im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
         origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
    clb=fig.colorbar(im, format='%.1E')
    clb.set_label('Normalized Pin Power Distribution')
    plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
This function allowing to plot the power distribution for a 1d or 2D geometry
"""
def power(self):
    data = np.loadtxt('app/Output/PowerDistribution.h')
    max_columns = len(data[0]) - 1
    max_rows = len(data)
    x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
    y = [data[0][colnum + 1] for colnum in range(max_columns)]
    z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    x,y = np.meshgrid(x, y)
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('X [cm]')                         
    ax.set_ylabel('Y [cm]') 
    im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
             origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
    clb = fig.colorbar(im, format='%.1E')
    #ax.set_xticklabels([])
    #ax.set_yticklabels([])
    #clb.ax.set_title('Peaking Factor Distribution ', loc='left', fontsize=9)
    #plt.xticks(rotation=45)
    #plt.yticks(rotation=45)
    #ax.xaxis.tick_top()
    clb.set_label('Power Distribution')
    plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
This function allowing to plot the neutron flux for a 1d or 2D geometry
"""
def plot(self):
    filename = open('app/link/script.dir', "r" ).read()
    with open(filename) as json_data:
         data = json.load(json_data)
         ng = data['data']['parameter']['Total number of energy groups']         
    data = []
    M00 = open('app/link/script00.py', "r" ).read()
    if M00 == 'CP1D':
        data = np.loadtxt('app/Output/FLUX_CP.H')
    elif M00 == 'SN1D':
        data = np.loadtxt('app/Output/FLUX_SN.H')
        #N = data['data']['parameter']['Number of Angular Discretization']
    elif M00 == 'MOC1D':
        data = np.loadtxt('app/Output/FLUX_MOC.H')
    elif M00 == 'SN2D':
        data = np.loadtxt('app/Output/FG1')
    else:
        QMessageBox.warning(self, "Warning", "select the calculation method.")
    if int(len(data)) >= 0:           
        if M00 != 'SN2D':
            fig, ax = plt.subplots()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
            matrix = []  
            for line in data:
                matrix.append(line)
            max_columns = len(matrix[0]) - 1
            max_rows = len(matrix)
            x = [matrix[rownum][0] for rownum in range(max_rows)]
            y = [[matrix[rownum][colnum + 1] for rownum in range(max_rows)] for colnum in range(max_columns)]
            p = [0]*max_columns
            for i in range(max_columns):
                key = 'Group', int(i+1)
                p[i] = plt.plot(x,y[i] ,label="Group %s" %(max_columns-i),linewidth=1)
            if M00 == 'CP1D':
                plt.title("CP Method")
            elif M00 == 'SN1D':
                plt.title("SN Method") 
            elif M00 == 'MOC1D':
                plt.title("MOC")  
            plt.xlabel('Distance [cm]')
            plt.ylabel('Normalized Flux')
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
        elif M00 == 'SN2D':     
            for i in range(ng):   
                data = np.loadtxt('app/Output/FG'+str(i+1))
                max_columns = len(data[0]) - 1
                max_rows = len(data)
                x = [data[rownum+1][0] for rownum in range(max_rows-1)] 
                y = [data[0][colnum + 1] for colnum in range(max_columns)]
                z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows-1)] for colnum in range(max_columns)]
                x = np.array(x)
                y = np.array(y)
                z = np.array(z)
                x,y = np.meshgrid(x, y, sparse=True)
                fig = matplotlib.pyplot.figure()

                ax = fig.add_subplot(111)
                ax.set_xlabel('X [cm]')                         
                ax.set_ylabel('Y [cm]') 
                ax.set_title('Energy Group '+str(i+1))
                #ax.set_title(r'$S_{'+''.join(str(N))+'}$, Energy Gr '+str(ng))
                im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
                origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
                clb = fig.colorbar(mappable=im, format='%.1E')
                clb.set_label('Normalized scalar flux')
                plt.show()
        else:
            QMessageBox.warning(self, "Warning", "Select More than a Fine Number of Meshes")

