## Introduction
This project mainly studies convex optimization methods for solving diffusion network inference problems. This project is based on the paper named  "Feature-Enhanced Probabilistic Models for Diffusion Network Inference" and "Uncovering the Temporal Dynamics of Diffusion Networks". Instead of using mentioned convex optimation methods (matlab tool cvx and l-bfgs-b), we proposed a new method: projected newton's mothod with bisection and interpolant used in line searching. 

## Algorithm
As the hessian matrix for the problem is easy to solve and the problem formulation has box constraints, we dicided to try projected newton's method. Instead of using pure newton's method with step size = 1, we choose step size t ourselves.

# Bisection
We use bisection in the line searching steps. We are here to minimize f(x-t*h) (we call it phi(t)) with respect to the step size t. This problem is indeed piecewise convex since x-t*h must be maintained as non-negative. Thus, it naturally brings us with 2 cases: (1) the minimizer is at one of the points where x-t*h = 0 (we call them breakpoints) (2) the minimizer is between two consecutive breakpoints. 

Intuitively, in the first case, we just check on every breakpoint if it is a local minimizer and in the second case, we just go through every posssible intervel ([0,first breakpoint], [first breakpoint, second breakpoint], ... , [last breakpoint, a]) and do bisection on them (since phi is convex on these intervals) to find local minimizers. However, we can make the algorithm more efficiently by leveraging on the properties: phi_prime(0) <= 0 and we can find a point a > the largest breakpoint where phi_prime(a) >= 0 (cuz phi_prime(inf) > 0) and maintain an invariance "phi_prime_left <= 0 and phi_prime_right >= 0". This gives us an optimized algorithm:

1. We first do bisection on the indexes of the breakpoints. we add 0 to the front of the breakpoints and a to the end, and we let start_index = 0 and end_index = len(breakpoints). We get the mid_index and check phi_prime(left_midbreakpoint) and phi_prime(right_midbreakpoint).
2. If phi_prime(midbreakpoint) is inf or nan, set end_index = mid_index.
3. If phi_prime(left_midbreakpoint) <= 0 and phi_prime(right_midbreakpoint) >= 0, midbreakpoint is already a local minimizer.
4. Else, either (1) phi_prime(left_midbreakpoint) > 0 or (2) phi_prime(right_midbreakpoint) < 0. If (1), we set end_index = mid_index. If (2), we set start_index = mid_index.
5. The stopping condition is start_index + 1 = end_index. Then we do bisection on the interval between these 2 consecutive breakpoints. 

# Interpolant
Although the above algorithm is good enough to find a suitable step size, we optimize it by utilizing an intepolant function \tilde\phi which can be uniquely defined by the original \phi and has a simpler formula to find its local minimizer. The step size t will then be set to the minimizer of the interpolant.

The motivation for these forms of the interpolant is that the ϕ(t) is the sum of a linear function x(t − t1) + y for some coefficients x, y and logarithmic terms of the form log(t − tl) or log(tr − t) with negative coefficients. The interpolant models just one of these terms.

The latest code has six cases for how to select bp. Use the following notation in this note: t1 := bp_begin_j, t2 := bp_end_j, δ := t2 − t1.




