# MPCC For Differential Drive Robots (Acados)

This project implements a real-time model preditive contouring controller (MPCC) for a differential drive robot to follow a 2D spline-based track. The general controller architecture is based on [this paper](https://arxiv.org/pdf/1711.07300).

Note that this module does **not** support dynamic path generations at runtime, this is because no full LUT is used. Symbolic track expressions are baked into the optimization problem when the code is generated. Hence it's mainly aimed at racing applications where a fixed track is to be implemented.

The main challenge was figuring out how to parametrize the track in a way that works symbolically with the CasADi/Acados framework.

## Track Parametrization:

1. **Track Generation**
   
   A small number of control points define the shape of the track. Using spline interpolation, a smooth path is generated.

2. **Arclength Parametrization**

    The total arclength of the track is computed numerically. This gives us the parameter s, which acts like a continuous progress variable along the track.

3. **Mapping $x(s)$, $y(s)$**
   
   NumPy’s CubicSpline is used to map $x$ and $y$ coordinates to the arclength array s. This results in two objects, $cs_x$ and $cs_y$, which can return position values given an arclength input.

4. **Discretization into $\theta_p$ Grid**

   A uniformly spaced LUT of arclength values called $\theta_p$ is defined, ranging from $0$ to the total track length $(L)$. Each value of $\theta_p$ is mapped to $(x, y)$ using the spline functions.

5. **Symbolic Track Lookup in CasADi**
   
    CasADi offers an ```interpolant()``` function that allows you to symbolically interpolate values from a lookup table. This means $\theta$ can be used as a symbolic decision variable and also obtain smooth $x(\theta)$ and $y(\theta)$ inside the optimization problem — as well as support differentiation.

    This feature does require the use of MX types in CasADi, so both the states and dynamics must be defined using using MX symbolics (**SX will not work**).

