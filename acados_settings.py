from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver
import numpy as np
import casadi as ca
from track import Track
from diffdrive_model import export_diffdrive_model


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
    
    x = ocp.model.x[0]
    y = ocp.model.x[1]
    psi = ocp.model.x[2]
    theta=ocp.model.x[3]

    u = ocp.model.u

    track = Track()
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
    

    return ocp_solver, acados_integrator, track