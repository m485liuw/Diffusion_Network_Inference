# Introduction
This project mainly studies convex optimization methods for solving diffusion network inference problems. This project is based on the paper named  "Feature-Enhanced Probabilistic Models for Diffusion Network Inference" and "Uncovering the Temporal Dynamics of Diffusion Networks". Instead of using mentioned convex optimation methods (matlab tool cvx and l-bfgs-b), we proposed a new method: projected newton's mothod with bisection and interpolant used. 

# Details
As the hessian matrix for the problem is easy to solve and the problem formulation has box constraints, we dicided to try projected newton's method. Instead of using pure newton's method with step size = 1, we choose step sizes by bisection. We are here to minimize f(x-t*h) (we call it phi(t)) with respect to the step size t. This problem is indeed piecewise convex since x-t*h must be maintained as non-negative. Thus, it naturally brings us with 2 cases: (1) the minimizer is at one of the points where x-t*h = 0 (we call them breakpoints) (2) the minimizer is between two consecutive breakpoints. We can prove phi_prime(0) <= 0 and we can find a point > the largest breakpoint where phi_prime of it >= 0, so we can maintain an invariance phi_prime_left <= 0 and phi_prime_right >= 0.





