from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosModel
import numpy as np
import casadi as ca
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
from track import Track

track = Track()

def export_diffdrive_model():
    model_name  = "basic_diffdrive_model"
    x1 = ca.MX.sym('x1')
    y1 = ca.MX.sym('y1')
    psi = ca.MX.sym('psi')
    theta= ca.MX.sym('theta') #track progress variable

    x = ca.vertcat(x1,y1,psi,theta)

    v = ca.MX.sym('v')
    omega = ca.MX.sym('omega')
    var = ca.MX.sym('var') 


    u = ca.vertcat(v,omega,var)

    x1_dot = ca.MX.sym('x_dot')
    y1_dot = ca.MX.sym('y_dot')
    psi_dot = ca.MX.sym('psi_dot')
    theta_dot = ca.MX.sym('theta_dot')
    
    x_dot = ca.vertcat(x1_dot,
                    y1_dot,
                    psi_dot,
                    theta_dot)
    
    f_expl = ca.vertcat(v*ca.cos(psi),
                     v*ca.sin(psi),
                     omega,
                     var)

    f_impl = x_dot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = x_dot
    model.u = u
    model.name = model_name

    model.t_label = "$t$ [s]"
    model.x_labels = ["$x$", "$y$", "$\\theta$"]
    model.u_labels = ["$v$", "$\\omega$"]


    return model





def setup(x0, N, tf):
    ocp = AcadosOcp()
    model = export_diffdrive_model()
    ocp.model = model
    ocp.solver_options.tf = tf
    ocp.solver_options.N_horizon = N

    q_c = 1e4
    q_l = 3e4
    q_theta = 1e2
    
    R_v = 1
    R_omega = 1
    R_var = 10
    
    R = np.diag([R_v, R_omega,R_var])

    x = ocp.model.x[0]
    y = ocp.model.x[1]
    psi = ocp.model.x[2]
    theta=ocp.model.x[3]

    u = ocp.model.u


    x_full, y_full,theta_full = track.create_grid(density=1000)

    x_lin = ca.interpolant('x', 'linear', [theta_full], x_full)
    y_lin = ca.interpolant('y', 'linear', [theta_full], y_full)

    x_ref = x_lin(theta)                
    y_ref = y_lin(theta)

    dx_dtheta = ca.gradient(x_ref, theta)   
    dy_dtheta = ca.gradient(y_ref, theta)
    phi       = ca.atan2(dy_dtheta, dx_dtheta)

    


    e_c = ca.sin(phi)*(x-x_ref) - ca.cos(phi)*(y-y_ref)
    e_l = -ca.cos(phi)*(x-x_ref) - ca.sin(phi)*(y-y_ref)
    
    ocp.model.cost_expr_ext_cost = (q_c * e_c**2 + q_l * e_l**2 + u[0]**2*R_v + u[1]**2*R_omega + u[2]**2*R_var - q_theta * theta)

    ocp.cost.cost_type = 'EXTERNAL'
    v_max = 1
    omega_max = 13

    ocp.constraints.lbu = np.array([-v_max/2, -omega_max, 0])
    ocp.constraints.ubu = np.array([v_max, omega_max, 1])
    ocp.constraints.idxbu = np.array([0, 1,2])

    
    ocp.constraints.x0 =  x0

    ocp.constraints.constr_type = 'BGH'
    ocp.constraints.constr_type_0 = 'BGH'
    ocp.constraints.constr_type_e = 'BGH'

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    
    
    ocp_solver = AcadosOcpSolver(ocp)
    acados_integrator = AcadosSimSolver(ocp)
    

    return ocp_solver, acados_integrator

def main(x0):
    N = 25
    tf=2.0


    ocp_solver, integrator = setup(x0, N, tf)

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

