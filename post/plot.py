import math
import json
import numpy as np
import matplotlib
font  = {'family'       : 'serif',
         'weight'       : 'normal',
         'size'         :  12,}
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

def plot_grid_figures(x, x_h, xc, xu, xv):
    # Original grid nodes:
    fig, ax = plt.subplots(figsize=(4,3))
    for i in range(nx):
        ax.plot(x[i,:,0], x[i,:,1], color='k', marker='.', markersize=5)
    for j in range(ny):
        ax.plot(x[:,j,0], x[:,j,1], color='k', marker='.', markersize=5)
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./grid/g{nx}x{ny}u.eps", format='eps')
    fig.savefig(f"./grid/g{nx}x{ny}u.png", dpi=300)
    plt.close(fig)

    # Grid nodes with halo cells:
    fig, ax = plt.subplots(figsize=(4,3))
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
    fig.savefig(f"./grid/g{nx}x{ny}u_whc.eps", format='eps')
    fig.savefig(f"./grid/g{nx}x{ny}u_whc.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4,3))
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
    fig.savefig(f"./grid/g{nx}x{ny}u_c.eps", format='eps')
    fig.savefig(f"./grid/g{nx}x{ny}u_c.png", dpi=300)
    plt.close(fig)

def plot_grid_metrics(xc, su, sv, v):
    fig, axs = plt.subplots(1, 2, figsize=(8,3))
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
    fig.savefig(f"./grid/g{nx}x{ny}u_suv_v_lines.eps", format='eps')
    #fig.savefig(f"./grid/g{nx}x{ny}u_suv_v_lines.png", dpi=300)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(8,6), sharey='row', sharex='col')
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

    fig, ax = plt.subplots(figsize=(4,3))
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

def plot_primitives(x_h, Qv, it):
    cmap_linear = matplotlib.cm.jet
    cmap_diverg = matplotlib.cm.seismic
    pref = 101325.0
    Tref = 300.0
    uref = 694.4
    fig, axs = plt.subplots(2, 2, figsize=(8,6), sharey='row', sharex='col')
    max_abs_pres = np.abs(Qv[1:-1,1:-1,0] - pref).max()
    norm_pres = matplotlib.colors.Normalize(vmin=(pref-max_abs_pres), vmax=(pref+max_abs_pres))
    cs1 = axs[0,0].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Qv[1:-1,1:-1,0], 
                              cmap=cmap_diverg, norm=norm_pres)
    axs[0,0].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb1 = fig.colorbar(cs1, ax=axs[0,0])
    cb1.ax.set_ylabel(r'$p \; [{\rm Pa}]$')
    axs[0,0].set(
        #xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
        #title  = r'$S_{\xi_x}$',
    )
    max_abs_u = np.abs(Qv[1:-1,1:-1,1] - uref).max()
    norm_u    = matplotlib.colors.Normalize(vmin=(uref-max_abs_u), vmax=(uref+max_abs_u))
    cs2 = axs[0,1].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Qv[1:-1,1:-1,1],
                              cmap=cmap_diverg, norm=norm_u)
    axs[0,1].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb2 = fig.colorbar(cs2, ax=axs[0,1])
    cb2.ax.set_ylabel(r'$u \; [{\rm m/s}]$')
    axs[0,1].set(
        #xlabel = r'$x$ [m]',
        #ylabel = r'$y$ [m]',
        #title  = r'$S_{\xi_y}$',
    )
    max_abs_v = np.abs(Qv[1:-1,1:-1,2]).max()
    norm_v    = matplotlib.colors.Normalize(vmin=(-max_abs_v), vmax=(max_abs_v))
    cs3 = axs[1,0].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Qv[1:-1,1:-1,2],
                              cmap=cmap_diverg, norm=norm_v)
    axs[1,0].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb3 = fig.colorbar(cs3, ax=axs[1,0])
    cb3.ax.set_ylabel(r'$v \; [{\rm m/s}]$')
    axs[1,0].set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
        #title  = r'$S_{\eta_x}$',
    )
    max_abs_T = np.abs(Qv[1:-1,1:-1,3] - Tref).max()
    norm_T    = matplotlib.colors.Normalize(vmin=(Tref-max_abs_T), vmax=(Tref+max_abs_T))
    cs4 = axs[1,1].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Qv[1:-1,1:-1,3],
                              cmap=cmap_diverg, norm=norm_T)
    axs[1,1].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
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

    fig, ax = plt.subplots(figsize=(4,3))
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

def plot_conservatives(x_h, Q, it):
    cmap_linear = matplotlib.cm.jet
    fig, axs = plt.subplots(2, 2, figsize=(8,6), sharey='row', sharex='col')
    cs1 = axs[0,0].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Q[1:-1,1:-1,0],
                              cmap=cmap_linear)
    axs[0,0].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb1 = fig.colorbar(cs1, ax=axs[0,0])
    cb1.ax.set_ylabel(r'$\rho \; [{\rm kg/m^3}]$')
    axs[0,0].set(
        #xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    cs2 = axs[0,1].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Q[1:-1,1:-1,1],
                              cmap=cmap_linear)
    axs[0,1].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb2 = fig.colorbar(cs2, ax=axs[0,1])
    cb2.ax.set_ylabel(r'$\rho u \; [{\rm kg/m^2 s}]$')
    axs[0,1].set(
        xlabel = r'$x$ [m]',
        #ylabel = r'$y$ [m]',
    )
    cs3 = axs[1,0].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Q[1:-1,1:-1,2],
                              cmap=cmap_linear)
    axs[1,0].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb3 = fig.colorbar(cs3, ax=axs[1,0])
    cb3.ax.set_ylabel(r'$\rho v \; [{\rm kg/m^2 s}]$')
    axs[1,0].set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    cs4 = axs[1,1].pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], Q[1:-1,1:-1,3],
                              cmap=cmap_linear)
    axs[1,1].plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    cb4 = fig.colorbar(cs4, ax=axs[1,1])
    cb4.ax.set_ylabel(r'$\rho e_t \; [{\rm J/m^3}]$')
    axs[1,1].set(
        xlabel = r'$x$ [m]',
        #ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./out/post/Q{it:05d}.png", dpi=300)
    plt.close(fig)

def plot_schlieren(x_h, su, sv, v, Q, it):
    cmap = matplotlib.cm.binary
    fig, ax = plt.subplots(figsize=(4,3))
    epx = su[:,:,0]/v[:,:]
    epy = su[:,:,1]/v[:,:]
    etx = sv[:,:,0]/v[:,:]
    ety = sv[:,:,1]/v[:,:]
    rho = Q[:,:,0]
    drhodx = 0.5*(rho[2:,1:-1]-rho[:-2,1:-1])*epx[1:-1,1:-1] \
           + 0.5*(rho[1:-1,2:]-rho[1:-1,:-2])*etx[1:-1,1:-1]
    drhody = 0.5*(rho[2:,1:-1]-rho[:-2,1:-1])*epy[1:-1,1:-1] \
           + 0.5*(rho[1:-1,2:]-rho[1:-1,:-2])*ety[1:-1,1:-1]
    drho_mag = np.sqrt(drhodx*drhodx + drhody*drhody) 
    #ax.pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], drhodx[:,:], cmap=cmap)
    ax.pcolormesh(x_h[1:-1,1:-1,0], x_h[1:-1,1:-1,1], drho_mag[:,:], cmap=cmap)
    ax.plot(x_h[1:-1,1,0], x_h[1:-1,1,1], color='black')
    ax.set(
        xlabel = r'$x$ [m]',
        ylabel = r'$y$ [m]',
    )
    fig.tight_layout()
    fig.savefig(f"./out/post/grad_rho_{it:05d}.png", dpi=300)
    plt.close(fig)

    return drho_mag

def plot_error_norms(nit, L2, Linfty):
    its = [i+1 for i in range(nit)]
    fig, axs = plt.subplots(2, 2, figsize=(8,6), sharey='row', sharex='col')
    axs[0,0].plot(its, L2[:,0]    , linestyle='solid',  color='black')
    axs[0,0].plot(its, Linfty[:,0], linestyle='dashed', color='black')
    axs[0,0].set(
        #xlabel = r'Iterations',
        ylabel = r'$L_2, L_\infty$ [-]',
        yscale = 'log',
        title  = r'$\rho/\rho_{ref}$',
    )
    axs[0,1].plot(its, L2[:,1]    , linestyle='solid',  color='black', label=r'$L_2$')
    axs[0,1].plot(its, Linfty[:,1], linestyle='dashed', color='black', label=r'$L_\infty$')
    axs[0,1].set(
        #xlabel = r'Iterations',
        #ylabel = r'$L_2, L_\infty$ [-]',
        yscale = 'log',
        title  = r'$\rho u/\rho_{ref} u_{ref}$',
    )
    axs[0,1].legend()
    axs[1,0].plot(its, L2[:,2]    , linestyle='solid',  color='black')
    axs[1,0].plot(its, Linfty[:,2], linestyle='dashed', color='black')
    axs[1,0].set(
        xlabel = r'Iterations',
        ylabel = r'$L_2, L_\infty$ [-]',
        yscale = 'log',
        title  = r'$\rho v/\rho_{ref} v_{ref}$',
    )
    axs[1,1].plot(its, L2[:,3]    , linestyle='solid',  color='black')
    axs[1,1].plot(its, Linfty[:,3], linestyle='dashed', color='black')
    axs[1,1].set(
        xlabel = r'Iterations',
        #ylabel = r'$L_2, L_\infty$ [-]',
        yscale = 'log',
        title  = r'$\rho e_t/\rho_{ref} e_{t,ref}$',
    )
    fig.tight_layout()
    fig.savefig(f"./out/post/error_norms.eps", format='eps')
    fig.savefig(f"./out/post/error_norms.png", dpi=300)
    plt.close(fig)

def get_validation_data(xc, drho_mag, Qv, Q):
    idx_x1 = np.argmax(drho_mag[:, 10]) 
    idx_x2 = np.argmax(drho_mag[:, 5]) 
    idx_x3 = np.nanargmax(Qv[:, 1, 0])-1
    #idx_x3 = np.argmax(drho_mag[:, 0]) 
    #idx_x3 = np.argmax(drho_mag[15:25, 0]) 
    data = {}
    angle = np.arctan2(xc[1+idx_x1, 10, 1] - xc[1+idx_x2, 5, 1],
                       xc[1+idx_x1, 10, 0] - xc[1+idx_x2, 5, 0])
    data['angle'] = angle*180/math.pi
    p_ref = 101325.0
    T_ref = 300.0
    u_ref = 694.4
    rho_ref = p_ref/(280.0*T_ref)
    p0_ref = p_ref + 0.5*rho_ref*u_ref*u_ref
    data['rho2_rho1'] = Q[1+idx_x3, 1, 0]/rho_ref
    data['p2_p1'] = Qv[1+idx_x3, 1, 0]/p_ref
    data['T2_T1'] = Qv[1+idx_x3, 1, 3]/T_ref
    c = math.sqrt(1.4*280*Qv[1+idx_x3, 0, 3])
    data['M2'] = Qv[1+idx_x3, 1, 1]/c
    p02 = Qv[1+idx_x3, 1, 0] \
        + 0.5*Q[1+idx_x3, 1, 0]*(Qv[1+idx_x3, 1, 1]*Qv[1+idx_x3, 1, 1]
                                +Qv[1+idx_x3, 1, 2]*Qv[1+idx_x3, 1, 2])
    data['p02_p01'] = p02/p0_ref
    return data

def plot_validation(its, validation):
    index = {0: 'angle', 1: 'rho2_rho1', 2: 'p2_p1', 3: 'T2_T1', 4: 'M2', 5: 'p02_p01'}
    valid = {0: 39.31, 1: 1.458, 2: 1.707, 3: 1.170, 4: 1.641, 5: 0.9846}
    label = {0: r'Wave Angle [$^\circ$]', 
             1: r'$\rho_2/\rho_1$ [-]', 
             2: r'$p_2/p_1$ [-]', 
             3: r'T$_2$/T$_1$ [-]', 
             4: r'M2 [-]',
             5: r'$p_{02}/p_{01}$ [-]',}
    for i in range(6):
        data = []
        for it in its:
            data.append(validation[f'{it}'][index[i]])
        fig, ax = plt.subplots(figsize=(4,3))
        ax.plot(its, data, 'ks', label='Present work')
        ax.plot(its, np.ones_like(its)*valid[i], 'k', label='Analytic solution')
        if i != 0:
            ax.set(
                xlabel = 'Iterations',
                ylabel = label[i],
                ylim   = (0.0, 2.0),
            )
        else:
            ax.set(
                xlabel = 'Iterations',
                ylabel = label[i],
                #ylim   = (0.0, 2.0),
            )
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"./out/post/{index[i]}.eps", format='eps')
        fig.savefig(f"./out/post/{index[i]}.png", dpi=300)
        plt.close(fig)


if __name__ == '__main__':
    #nx, ny, nhc = 33, 25, 1
    nx, ny, nhc = 65, 49, 1
    nrest = 125
    trest = 20 
    its = [nrest*i for i in range(trest+1)]
    #its = [nrest*i for i in range(11)]

    # Grid nodes:
    x   = np.fromfile(f"./grid/g{nx}x{ny}u.bin", dtype=np.double)
    x   = x.reshape((nx,ny,2))
    x_h = np.fromfile(f"./grid/g{nx}x{ny}u_whc.bin", dtype=np.double)
    x_h = x_h.reshape((nx+2*nhc,ny+2*nhc,2))

    # Cell centers and face centroids:
    xc  = np.fromfile(f"./grid/g{nx}x{ny}u_xc.bin", dtype=np.double)
    xc  = xc.reshape((nx+2*nhc-1,ny+2*nhc-1,2))
    xu  = np.fromfile(f"./grid/g{nx}x{ny}u_xu.bin", dtype=np.double)
    xu  = xu.reshape((nx+2*nhc,  ny+2*nhc-1,2))
    xv  = np.fromfile(f"./grid/g{nx}x{ny}u_xv.bin", dtype=np.double)
    xv  = xv.reshape((nx+2*nhc-1,ny+2*nhc,  2))

    #plot_grid_figures(x, x_h, xc, xu, xv)
    #print("--Grid figures ok!")

    # Projected cell face areas and volumes along lines:
    su  = np.fromfile(f"./grid/g{nx}x{ny}u_su.bin", dtype=np.double)
    su  = su.reshape((nx+2*nhc,ny+2*nhc-1,2))
    su  = 0.5*(su[1:,:,:]+su[:-1,:,:])
    sv  = np.fromfile(f"./grid/g{nx}x{ny}u_sv.bin", dtype=np.double)
    sv  = sv.reshape((nx+2*nhc-1,ny+2*nhc,2))
    sv  = 0.5*(sv[:,1:,:]+sv[:,:-1,:])
    v   = np.fromfile(f"./grid/g{nx}x{ny}u_v.bin", dtype=np.double)
    v   = v.reshape((nx+2*nhc-1,ny+2*nhc-1))

    #plot_grid_metrics(xc, su, sv, v)
    #print("--Grid metrics ok!")

    # Vectors: primitive and conservative 
    validation = {}
    for it in its:
        print(f"--Figures for iteration {it}...")
        Qv  = np.fromfile(f"./out/rest/Qv{it:05d}.bin", dtype=np.double)
        Qv  = Qv.reshape((nx+2*nhc-1, ny+2*nhc-1, 4))

        #plot_primitives(x_h, Qv, it)

        Q  = np.fromfile(f"./out/rest/Q{it:05d}.bin", dtype=np.double)
        Q  = Q.reshape((nx+2*nhc-1, ny+2*nhc-1, 4))

        #plot_conservatives(x_h, Q, it)

        drho_mag = plot_schlieren(x_h, su, sv, v, Q, it)

        data = get_validation_data(xc, drho_mag, Qv, Q)
        validation[f'{it}'] = data
        print(validation[f'{it}'])

    plot_validation(its, validation)
    j = json.dumps(validation, indent=4)
    with open('./out/post/validation.json', 'w') as f:
        print(j, file=f)

    L2     = np.fromfile(f"./out/rest/L2.bin", dtype=np.double)
    L2     = L2.reshape((nrest*trest, 4))
    Linfty = np.fromfile(f"./out/rest/Linfty.bin", dtype=np.double)
    Linfty = Linfty.reshape((nrest*trest, 4))

    plot_error_norms(nrest*trest, L2, Linfty)
    print("--Error norms ok!")

