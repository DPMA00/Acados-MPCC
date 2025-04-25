# MPCC For Differential Drive Robots (Acados)

This project implements a real-time model preditive contouring controller (MPCC) for a differential drive robot to follow a 2D spline-based track. The general controller architecture is based on [this paper](https://arxiv.org/pdf/1711.07300). 

The main challenge was figuring out how to parametrize the track in a way that works symbolically with the CasADi/Acados framework.

## Track Parametrization: What I Did

1. **Track Generation**
   
   A small number of control points define the shape of the track. Using spline interpolation, a smooth path is generated.

2. **Arclength Parametrization**

    The total arclength of the track is computed numerically. This gives us the parameter s, which acts like a continuous progress variable along the track.

3. **Mapping $x(s)$, $y(s)$**
   
   I used NumPy’s CubicSpline to map $x$ and $y$ coordinates to the arclength array s. This results in two objects, $cs_x$ and $cs_y$, which can return position values given an arclength input.

4. **Discretization into $\theta_p$ Grid**

   I defined a uniformly spaced lookup table (LUT) of arclength values called $\theta_p$, ranging from $0$ to the total track length $(L)$. Each value of $\theta_p$ is mapped to $(x, y)$ using the spline functions.

5. **Symbolic Track Lookup in CasADi**
   
    Here’s the interesting part: CasADi offers an ```interpolate()``` function that allows you to symbolically interpolate values from a lookup table. This means I can use $\theta$ as a symbolic decision variable and obtain smooth $x(\theta)$ and $y(\theta)$ inside the optimization problem — and even take derivatives.

    There is a catch, though: this feature requires MX graph types in CasADi, so both the states and dynamics must be defined using MX symbols (not SX). Once that's in place, it will work.