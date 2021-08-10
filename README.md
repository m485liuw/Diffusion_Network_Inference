## Introduction
This project mainly studies convex optimization methods for solving diffusion network inference problems. This project is based on the paper named  "Feature-Enhanced Probabilistic Models for Diffusion Network Inference" and "Uncovering the Temporal Dynamics of Diffusion Networks". Instead of using mentioned convex optimation methods (matlab tool cvx and l-bfgs-b), we proposed a new method: projected newton's mothod with bisection and interpolant used in line searching. 

## Algorithm
As the hessian matrix for the problem is easy to solve and the problem formulation has box constraints, we dicided to try projected newton's method. Instead of using pure newton's method with step size = 1, we choose step size t ourselves.

# Bisection
We use bisection in the line searching steps. We are here to minimize f(x-t*h) (we call it phi(t)) with respect to the step size t. This problem is indeed piecewise convex since x-t*h must be maintained as non-negative. Thus, it naturally brings us with 2 cases: (1) the minimizer is at one of the points where x-t*h = 0 (we call them breakpoints) (2) the minimizer is between two consecutive breakpoints. 

Intuitively, in the first case, we just check on every breakpoint if it is a local minimizer and in the second case, we just go through every posssible intervel ([0,first breakpoint], [first breakpoint, second breakpoint], ... , [last breakpoint, a]) and do bisection on them (since phi is convex on these intervals) to find local minimizers. However, we can make the algorithm more efficiently by leveraging on the properties: phi_prime(0) <= 0 and we can find a point a > the largest breakpoint where phi_prime(a) >= 0 (cuz phi_prime(inf) > 0) and maintain an invariance "phi_prime_left <= 0 and phi_prime_right >= 0". This gives us an optimized algorithm:

1. We first do bisection on the indexes of the breakpoints. we add 0 to the front of the breakpoints and a to the end, and we let start_index = 0 and end_index = len(breakpoints). 
2. Repeat until start_index + 1 = end_index:
   Get the mid_index and check phi_prime(left_midbreakpoint) and phi_prime(right_midbreakpoint).
   If phi_prime(midbreakpoint) is inf or nan, set end_index = mid_index.
   If phi_prime(left_midbreakpoint) <= 0 and phi_prime(right_midbreakpoint) >= 0, midbreakpoint is already a local minimizer (record it as special_bp).
   Else, either (1) phi_prime(left_midbreakpoint) > 0 or (2) phi_prime(right_midbreakpoint) < 0. If (1), we set end_index = mid_index. If (2), we set start_index = mid_index.
3.  If special_bp is selected, set t as special_bp. Else, we do bisection on the interval between these 2 consecutive breakpoints. 

# Interpolant
Although the above algorithm is good enough to find a suitable step size, we optimize it by utilizing an interpolant function \tilde\phi which can be uniquely defined by the original \phi and has a simpler formula to find its local minimizer. The step size t will then be set to the minimizer of the interpolant.

Like the algorithm above, we do the same from step 1 to step 3 except in the else condition in step 3, we further divided into 5 cases:

1. if end_index = 2 and end_breakpoint > 1, then set t as 1 as as in pure newton's method.

From now on, we notate start_breakpoint as t1 and end_breakpoint as t2. t2-t1 as δ.

2. if end_breakpoint is inf, then an interpolant ˜ϕ of the form ˜ϕ(t) = x(t − t1) + y + z ln(t2 − t) is selected, where x, y, z are chosen so that ˜ϕ(t1) = ϕ(t1), ˜ϕ′(t1) = ϕ′(t1), and ˜ϕ′′(t1) = ϕ′′(t1).
3. if (ϕ(t2) − ϕ(t1)) / δ < (ϕ′(t1) + ϕ′(t2)) / 2 − tol, then an interpolant of the form ˜ϕ(t) = x(t − t1) + y + z ln(t2 − t + w) is used. The coefficients x, y, z, w are chosen so that ˜ϕ(t1) = ϕ(t1), ˜ϕ(t2) = ϕ(t2), ˜ϕ′(t1) = ϕ′(t1), ˜ϕ′(t2) = ϕ′(t2).
4. (ϕ(t2) − ϕ(t1)) / δ > (ϕ′(t1) + ϕ′(t2)) / 2 + tol, then an interpolant of the form ˜ϕ(t) = x(t − t1) + y + z ln(t − t1 + w) is used. The coefficients satisfy the same equations as the previous item.
5. Else, bisection is used.

One can prove that in Cases 3–5, the interpolant is uniquely determined by the data, it is convex (i.e., z ≤ 0), and there is a simple formula to find its minimizer. Determining the interpolant, i.e., finding the coefficients x, y, z, w, is accomplished by reducing the problem to finding the root of a univariate function of ¯ w := w/δ and then solving this function with bisection. Note that this bisection does not involve evaluating the original f but only involves a simple univariate function of ¯ w and so is very fast.

The motivation for these forms of the interpolant is that the ϕ(t) is the sum of a linear function x(t − t1) + y for some coefficients x, y and logarithmic terms of the form log(t − tl) or log(tr − t) with negative coefficients. The interpolant models just one of these terms.






