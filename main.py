import numpy as np
import matplotlib.pyplot as plt
from acados_settings import setup


def main(x0):
    N = 25
    tf=2.0


    ocp_solver, integrator, track = setup(x0, N, tf)

    nx = ocp_solver.acados_ocp.dims.nx
    nu = ocp_solver.acados_ocp.dims.nu

    Nsim = 700

    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))

    simX[0,:] = x0

    t_prep        = np.zeros(Nsim)      
    t_feedback    = np.zeros(Nsim)     
    for i in range(Nsim):

        ocp_solver.options_set('rti_phase', 1)
        status = ocp_solver.solve()
        t_prep[i] = ocp_solver.get_stats('time_tot')
        ocp_solver.set(0,"lbx", simX[i,:])
        ocp_solver.set(0, "ubx", simX[i, :])

        ocp_solver.options_set('rti_phase',2)
        status=ocp_solver.solve()
        t_feedback[i] = ocp_solver.get_stats('time_tot')
        simU[i,:] = ocp_solver.get(0,"u")

        simX[i+1,:] = integrator.simulate(x=simX[i,:], u=simU[i,:])

    t_prep_ms     = 1e3 * t_prep
    t_feedback_ms = 1e3 * t_feedback

    print('\n==== acados timing (ms) =====')
    print(f'Preparation  –  min {t_prep_ms.min():7.3f}   '
        f'mean {t_prep_ms.mean():7.3f}   max {t_prep_ms.max():7.3f}')
    print(f'Feedback     –  min {t_feedback_ms.min():7.3f}   '
        f'mean {t_feedback_ms.mean():7.3f}   max {t_feedback_ms.max():7.3f}')
    print('==============================')
    fig, ax = plt.subplots()
    xses, yses = track.get_track_plot_params()
    ax.plot(simX[:,0], simX[:,1], '-o',color='r', label='closed-loop path')
    ax.plot(xses,yses)
    ax.axis('equal'); ax.grid(True)
    ax.legend()
    ax.set_title(f"MPC Closed‑Loop Simulation ({Nsim} steps)")
    ax.set_xlabel("x [m]"); plt.ylabel("y [m]")
    plt.show()


if __name__ == '__main__':
    main(x0 = np.array([0,0,0,0]))

