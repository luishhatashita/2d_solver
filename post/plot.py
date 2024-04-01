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
    x   = np.fromfile("./g65x49u.bin", dtype=np.double)
    x   = x.reshape((65,49,2))
    x_h = np.fromfile("./g65x49u_whc.bin", dtype=np.double)
    x_h = x_h.reshape((67,51,2))

    # Original grid nodes:
    fig, ax = plt.subplots(figsize=(8,6))
    for i in range(65):
        ax.plot(x[i,:,0], x[i,:,1], color='k', marker='.', markersize=5)
    for j in range(49):
        ax.plot(x[:,j,0], x[:,j,1], color='k', marker='.', markersize=5)
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig("./g65x49u.png", dpi=300)
    plt.close(fig)

    # Grid nodes with halo cells:
    fig, ax = plt.subplots(figsize=(8,6))
    for i in range(67):
        ax.plot(x_h[i,:,0], x_h[i,:,1], color='r', marker='.', markersize=5)
    for j in range(51):
        ax.plot(x_h[:,j,0], x_h[:,j,1], color='r', marker='.', markersize=5)
    for i in range(65):
        ax.plot(x[i,:,0], x[i,:,1], color='k', marker='.', markersize=5)
    for j in range(49):
        ax.plot(x[:,j,0], x[:,j,1], color='k', marker='.', markersize=5)
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig("./g65x49u_whc.png", dpi=300)
    plt.close(fig)

    # Grid nodes with cell centers and face centroids:
    xc  = np.fromfile("./g65x49u_xc.bin", dtype=np.double)
    xc  = xc.reshape((66,50,2))
    xu  = np.fromfile("./g65x49u_xu.bin", dtype=np.double)
    xu  = xu.reshape((67,50,2))
    xv  = np.fromfile("./g65x49u_xv.bin", dtype=np.double)
    xv  = xv.reshape((66,51,2))
    fig, ax = plt.subplots(figsize=(12,8))
    for i in range(67):
        ax.plot(x_h[i,:,0], x_h[i,:,1], 
                color='r',
                alpha=0.25, linewidth=0.75
        )
    for j in range(51):
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
    fig.savefig("./g65x49u_c.png", dpi=300)
    plt.close(fig)

    # Projected cell face areas and volumes along lines:
    su  = np.fromfile("./g65x49u_su.bin", dtype=np.double)
    su  = su.reshape((67,50,2))
    su  = 0.5*(su[1:,:,:]+su[:-1,:,:])
    sv  = np.fromfile("./g65x49u_sv.bin", dtype=np.double)
    sv  = sv.reshape((66,51,2))
    sv  = 0.5*(sv[:,1:,:]+sv[:,:-1,:])
    v   = np.fromfile("./g65x49u_v.bin", dtype=np.double)
    v   = v.reshape((66,50))
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
    axs[1].plot(su[30,1:-1,0], xc[30,1:-1,1], linestyle='dashed', color='tab:blue')
    axs[1].plot(su[30,1:-1,1], xc[30,1:-1,1], linestyle='dotted', color='tab:blue')
    axs[1].plot(sv[30,1:-1,0], xc[30,1:-1,1], linestyle='dashed', color='tab:orange')
    axs[1].plot(sv[30,1:-1,1], xc[30,1:-1,1], linestyle='dotted', color='tab:orange')
    axs[1].plot( v[30,1:-1]  , xc[30,1:-1,1], linestyle='solid',  color='tab:green')
    axs[1].set(
        xlabel = r'$S_{\xi_x, \xi_y}, S_{\eta_x, \eta_y}, V \; [{\rm m^2, m^3}]$',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig("./g65x49u_suv_v_lines.png", dpi=300)
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
    fig.savefig("./g65x49u_suv_v.png", dpi=300)
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
    fig.savefig("./g65x49u_v.png", dpi=300)
    plt.close(fig)
