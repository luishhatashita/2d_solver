import math
import numpy as np
import matplotlib
font  = {'family'       : 'serif',
         'weight'       : 'normal',
         'size'         :  14,}
xtick = {'direction'    : 'in',
         'minor.visible': True,
         'top'          : True,}
ytick = {'direction'    : 'in',
         'minor.visible': True,
         'right'        : True,}
text  = {'usetex'       : True,}
matplotlib.rc('font' , **font )
matplotlib.rc('xtick', **xtick)
matplotlib.rc('ytick', **ytick)
matplotlib.rc('text' , **text )
import matplotlib.pyplot as plt

if __name__ == '__main__':
    #nx, ny, nhc = 33, 25, 1
    nx, ny, nhc = 65, 49, 1
    nrest = 1000 
    its = [nrest*i for i in range(3)]
    #its = [nrest*i for i in range(11)]

    x   = np.fromfile(f"./grid/g{nx}x{ny}u.bin", dtype=np.double)
    x   = x.reshape((nx,ny,2))
    x_h = np.fromfile(f"./grid/g{nx}x{ny}u_whc.bin", dtype=np.double)
    x_h = x_h.reshape((nx+2*nhc,ny+2*nhc,2))

    # Original grid nodes:
    fig, ax = plt.subplots(figsize=(8,6))
    for i in range(nx):
        ax.plot(x[i,:,0], x[i,:,1], color='k', marker='.', markersize=5)
    for j in range(ny):
        ax.plot(x[:,j,0], x[:,j,1], color='k', marker='.', markersize=5)
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u.png", dpi=300)
    plt.close(fig)

    # Grid nodes with halo cells:
    fig, ax = plt.subplots(figsize=(8,6))
    for i in range(nx+2*nhc):
        ax.plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
    for j in range(ny+2*nhc):
        ax.plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
    for i in range(nx):
        ax.plot(x[i,:,0], x[i,:,1], color='k', marker='.', markersize=5)
    for j in range(ny):
        ax.plot(x[:,j,0], x[:,j,1], color='k', marker='.', markersize=5)
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u_whc.png", dpi=300)
    plt.close(fig)

    # Grid nodes with cell centers and face centroids:
    xc  = np.fromfile(f"./grid/g{nx}x{ny}u_xc.bin", dtype=np.double)
    xc  = xc.reshape((nx+2*nhc-1,ny+2*nhc-1,2))
    xu  = np.fromfile(f"./grid/g{nx}x{ny}u_xu.bin", dtype=np.double)
    xu  = xu.reshape((nx+2*nhc,  ny+2*nhc-1,2))
    xv  = np.fromfile(f"./grid/g{nx}x{ny}u_xv.bin", dtype=np.double)
    xv  = xv.reshape((nx+2*nhc-1,ny+2*nhc,  2))
    fig, ax = plt.subplots(figsize=(12,8))
    for i in range(nx+2*nhc):
        ax.plot(x_h[i,:,0], x_h[i,:,1], 
                color='r',
                alpha=0.25, linewidth=0.75
        )
    for j in range(ny+2*nhc):
        ax.plot(x_h[:,j,0], x_h[:,j,1], 
                color='r',
                alpha=0.25, linewidth=0.75
        )
    ax.scatter(xc[:,:,0].ravel(), xc[:,:,1].ravel(), color='r' , s=2)
    ax.scatter(xu[:,:,0].ravel(), xu[:,:,1].ravel(), color='g' , s=2)
    ax.scatter(xv[:,:,0].ravel(), xv[:,:,1].ravel(), color='b' , s=2)
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u_c.png", dpi=300)
    plt.close(fig)

    # Projected cell face areas and volumes along lines:
    su  = np.fromfile(f"./grid/g{nx}x{ny}u_su.bin", dtype=np.double)
    su  = su.reshape((nx+2*nhc,ny+2*nhc-1,2))
    su  = 0.5*(su[1:,:,:]+su[:-1,:,:])
    sv  = np.fromfile(f"./grid/g{nx}x{ny}u_sv.bin", dtype=np.double)
    sv  = sv.reshape((nx+2*nhc-1,ny+2*nhc,2))
    sv  = 0.5*(sv[:,1:,:]+sv[:,:-1,:])
    v   = np.fromfile(f"./grid/g{nx}x{ny}u_v.bin", dtype=np.double)
    v   = v.reshape((nx+2*nhc-1,ny+2*nhc-1))
    fig, axs = plt.subplots(1, 2, figsize=(12,5))
    axs[0].plot(xc[1:-1,1,0], su[1:-1,1,0], linestyle='dashed', color='tab:blue')
    axs[0].plot(xc[1:-1,1,0], su[1:-1,1,1], linestyle='dotted', color='tab:blue')
    axs[0].plot(xc[1:-1,1,0], sv[1:-1,1,0], linestyle='dashed', color='tab:orange')
    axs[0].plot(xc[1:-1,1,0], sv[1:-1,1,1], linestyle='dotted', color='tab:orange')
    axs[0].plot(xc[1:-1,1,0],  v[1:-1,1]  , linestyle='solid',  color='tab:green')
    axs[0].set(
        xlabel = r'$x$ [m]',
        ylabel = r'$S_{\xi_x, \xi_y}, S_{\eta_x, \eta_y}, V \; [{\rm m^2, m^3}]$',
    )
    ix = int(nx/2)
    axs[1].plot(su[ix,1:-1,0], xc[ix,1:-1,1], linestyle='dashed', color='tab:blue')
    axs[1].plot(su[ix,1:-1,1], xc[ix,1:-1,1], linestyle='dotted', color='tab:blue')
    axs[1].plot(sv[ix,1:-1,0], xc[ix,1:-1,1], linestyle='dashed', color='tab:orange')
    axs[1].plot(sv[ix,1:-1,1], xc[ix,1:-1,1], linestyle='dotted', color='tab:orange')
    axs[1].plot( v[ix,1:-1]  , xc[ix,1:-1,1], linestyle='solid',  color='tab:green')
    axs[1].set(
        xlabel = r'$S_{\xi_x, \xi_y}, S_{\eta_x, \eta_y}, V \; [{\rm m^2, m^3}]$',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u_suv_v_lines.png", dpi=300)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(12,8), sharey='row', sharex='col')
    cs1 = axs[0,0].contourf(xc[1:-1,1:-1,0], xc[1:-1,1:-1,1], su[1:-1,1:-1,0])
    cb1 = fig.colorbar(cs1, ax=axs[0,0])
    cb1.ax.set_ylabel(r'$S_{\xi_x} \; [{\rm m^2}]$')
    axs[0,0].set(
        #xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
        #title  = r'$S_{\xi_x}$',
    )
    cs2 = axs[0,1].contourf(xc[1:-1,1:-1,0], xc[1:-1,1:-1,1], su[1:-1,1:-1,1])
    cb2 = fig.colorbar(cs2, ax=axs[0,1])
    cb2.ax.set_ylabel(r'$S_{\xi_y} \; [{\rm m^2}]$')
    axs[0,1].set(
        xlabel = r'$x$ [m]',
        #ylabel = r'$y$ [m]',
        #title  = r'$S_{\xi_y}$',
    )
    cs3 = axs[1,0].contourf(xc[1:-1,1:-1,0], xc[1:-1,1:-1,1], sv[1:-1,1:-1,0])
    cb3 = fig.colorbar(cs3, ax=axs[1,0])
    cb3.ax.set_ylabel(r'$S_{\eta_x} \; [{\rm m^2}]$')
    axs[1,0].set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
        #title  = r'$S_{\eta_x}$',
    )
    cs4 = axs[1,1].contourf(xc[1:-1,1:-1,0], xc[1:-1,1:-1,1], sv[1:-1,1:-1,1])
    cb4 = fig.colorbar(cs4, ax=axs[1,1])
    cb4.ax.set_ylabel(r'$S_{\eta_y} \; [{\rm m^2}]$')
    axs[1,1].set(
        xlabel = r'$x$ [m]',
        #ylabel = r'$y$ [m]',
        #title  = r'$S_{\eta_y}$',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u_suv_v.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8,6))
    cs = ax.contourf(xc[1:-1,1:-1,0], xc[1:-1,1:-1,1], v[1:-1,1:-1])
    cb = fig.colorbar(cs, ax=ax)
    cb.ax.set_ylabel(r'$V \; [{\rm m^3}]$')
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
        #title  = r'$S_{\eta_y}$',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u_v.png", dpi=300)
    plt.close(fig)


    # Vectors: primitive and conservative 
    for it in its:
        fig, axs = plt.subplots(2, 2, figsize=(12,8), sharey='row', sharex='col')
        Qv  = np.fromfile(f"./out/rest/Qv{it:05d}.bin", dtype=np.double)
        Qv  = Qv.reshape((nx+2*nhc-1, ny+2*nhc-1, 4))
        cs1 = axs[0,0].pcolormesh(x_h[:,:,0], x_h[:,:,1], Qv[:,:,0])
        #cs1 = axs[0,0].pcolormesh(xc[:,:,0], xc[:,:,1], Qv[:,:,0])
        #for i in range(nx+2*nhc):
        #    axs[0,0].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[0,0].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb1 = fig.colorbar(cs1, ax=axs[0,0])
        cb1.ax.set_ylabel(r'$p \; [{\rm Pa}]$')
        axs[0,0].set(
            #xlabel = r'$x$ [m]',
            ylabel = r'$y$ [m]',
            #title  = r'$S_{\xi_x}$',
        )
        cs2 = axs[0,1].pcolormesh(x_h[:,:,0], x_h[:,:,1], Qv[:,:,1])
        #cs2 = axs[0,1].pcolormesh(xc[:,:,0], xc[:,:,1], Qv[:,:,1])
        #for i in range(nx+2*nhc):
        #    axs[0,1].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[0,1].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb2 = fig.colorbar(cs2, ax=axs[0,1])
        cb2.ax.set_ylabel(r'$u \; [{\rm m/s}]$')
        axs[0,1].set(
            xlabel = r'$x$ [m]',
            #ylabel = r'$y$ [m]',
            #title  = r'$S_{\xi_y}$',
        )
        cs3 = axs[1,0].pcolormesh(x_h[:,:,0], x_h[:,:,1], Qv[:,:,2])
        #cs3 = axs[1,0].pcolormesh(xc[:,:,0], xc[:,:,1], Qv[:,:,2])
        #for i in range(nx+2*nhc):
        #    axs[1,0].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[1,0].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb3 = fig.colorbar(cs3, ax=axs[1,0])
        cb3.ax.set_ylabel(r'$v \; [{\rm m/s}]$')
        axs[1,0].set(
            xlabel = r'$x$ [m]',
            ylabel = r'$y$ [m]',
            #title  = r'$S_{\eta_x}$',
        )
        cs4 = axs[1,1].pcolormesh(x_h[:,:,0], x_h[:,:,1], Qv[:,:,3])
        #cs4 = axs[1,1].pcolormesh(xc[:,:,0], xc[:,:,1], Qv[:,:,3])
        #for i in range(nx+2*nhc):
        #    axs[1,1].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[1,1].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb4 = fig.colorbar(cs4, ax=axs[1,1])
        cb4.ax.set_ylabel(r'$T \; [{\rm K}]$')
        axs[1,1].set(
            xlabel = r'$x$ [m]',
            #ylabel = r'$y$ [m]',
            #title  = r'$S_{\eta_y}$',
        )
        fig.tight_layout()
        fig.savefig(f"./out/post/Qv{it:05d}.png", dpi=300)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8,6))
        u_mag = np.sqrt(Qv[:,:,1]*Qv[:,:,1] + Qv[:,:,2]*Qv[:,:,2])
        for i in range(nx+2*nhc):
            ax.plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        for j in range(ny+2*nhc):
            ax.plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        ax.quiver(xc[:,:,0], xc[:,:,1], Qv[:,:,1]/u_mag, Qv[:,:,2]/u_mag)
        ax.set(
            #xlabel = r'$x$ [m]',
            ylabel = r'$y$ [m]',
        )
        fig.tight_layout()
        fig.savefig(f"./out/post/v{it:05d}.png", dpi=300)
        plt.close(fig)

        fig, axs = plt.subplots(2, 2, figsize=(12,8), sharey='row', sharex='col')
        Q  = np.fromfile(f"./out/rest/Q{it:05d}.bin", dtype=np.double)
        Q  = Q.reshape((nx+2*nhc-1, ny+2*nhc-1, 4))
        cs1 = axs[0,0].pcolormesh(x_h[:,:,0], x_h[:,:,1], Q[:,:,0])
        #cs1 = axs[0,0].pcolormesh(xc[:,:,0], xc[:,:,1], Q[:,:,0])
        #for i in range(nx+2*nhc):
        #    axs[0,0].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[0,0].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb1 = fig.colorbar(cs1, ax=axs[0,0])
        cb1.ax.set_ylabel(r'$\rho \; [{\rm kg/m^3}]$')
        axs[0,0].set(
            #xlabel = r'$x$ [m]',
            ylabel = r'$y$ [m]',
            #title  = r'$S_{\xi_x}$',
        )
        cs2 = axs[0,1].pcolormesh(x_h[:,:,0], x_h[:,:,1], Q[:,:,1])
        #cs2 = axs[0,1].pcolormesh(xc[:,:,0], xc[:,:,1], Q[:,:,1])
        #for i in range(nx+2*nhc):
        #    axs[0,1].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[0,1].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb2 = fig.colorbar(cs2, ax=axs[0,1])
        cb2.ax.set_ylabel(r'$\rho u \; [{\rm kg/m^2 s}]$')
        axs[0,1].set(
            xlabel = r'$x$ [m]',
            #ylabel = r'$y$ [m]',
            #title  = r'$S_{\xi_y}$',
        )
        cs3 = axs[1,0].pcolormesh(x_h[:,:,0], x_h[:,:,1], Q[:,:,2])
        #cs3 = axs[1,0].pcolormesh(xc[:,:,0], xc[:,:,1], Q[:,:,2])
        #for i in range(nx+2*nhc):
        #    axs[1,0].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[1,0].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb3 = fig.colorbar(cs3, ax=axs[1,0])
        cb3.ax.set_ylabel(r'$\rho v \; [{\rm kg/m^2 s}]$')
        axs[1,0].set(
            xlabel = r'$x$ [m]',
            ylabel = r'$y$ [m]',
            #title  = r'$S_{\eta_x}$',
        )
        cs4 = axs[1,1].pcolormesh(x_h[:,:,0], x_h[:,:,1], Q[:,:,3])
        #cs4 = axs[1,1].pcolormesh(xc[:,:,0], xc[:,:,1], Q[:,:,3])
        #for i in range(nx+2*nhc):
        #    axs[1,1].plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
        #for j in range(ny+2*nhc):
        #    axs[1,1].plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
        cb4 = fig.colorbar(cs4, ax=axs[1,1])
        cb4.ax.set_ylabel(r'$\rho e_t \; [{\rm J/m^3}]$')
        axs[1,1].set(
            xlabel = r'$x$ [m]',
            #ylabel = r'$y$ [m]',
            #title  = r'$S_{\eta_y}$',
        )
        fig.tight_layout()
        fig.savefig(f"./out/post/Q{it:05d}.png", dpi=300)
        plt.close(fig)
