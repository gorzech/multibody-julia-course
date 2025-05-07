### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 0cd41f17-71d9-4ebc-868e-39741e37b13c
md"""
## Solving Nonlinear Equations

How to find a solution to an equation in the general form $f(x)=0$. For example, $e^{-x} \sin x = \cos x$, thus $f(x) = e^{-x} \sin x - \cos x$. We want to find an approximate solution rather than an exact one. This problem is also called root finding. Main areas:

* Implicit ODE solvers
* Optimization – finding the minimum of a function
* Kinematic analysis of mechanisms
"""

# ╔═╡ 40f4936b-5ec4-4c98-9f48-cc41054cfd80
md"""
## Brute Force Method

Based on the collection of points $(x, f(x))$ – used e.g., in plotting. Simply go through all points and check when function crosses axis. This requires a lot of computation, therefore the name: brute force. But this should be easy to implement.
"""

# ╔═╡ ce5a4e83-1a4b-4955-a6f3-67806b29b6cc
md"""
## Brute Force Method – Algorithm

We have a set of points $(x_i, y_i)$, $y_i = f(x_i)$ $i = 0, \ldots, n$, and $x_0 < \ldots < x_n$. Curve crosses x axis when $y_i$ and $y_{i+1}$ have different signs: $y_i y_{i+1} \leq 0$. That is $f(x) = 0$ is in $[x_i, x_{i+1}]$. Assuming linear variation in interval, the root is located at $x = x_i - \frac{(x_{i+1} - x_i)}{(y_{i+1} - y_i)} y_i$. Now, take a look at the implementation:

* `brute_force_root_finder.m`
* `demo_brute_force_root_finder.m`
"""

# ╔═╡ f02d3fcd-826b-430e-992a-767d4bcdbb11
md"""
## Hands-on – Preallocation

This command reallocates roots vector at each occurrence `roots = [roots; root];`. Matlab warns this slows down the program. Can we avoid this? Propose and implement better solution. Benchmark with `timeit` to test if preallocation helps. Try to test with a function with a larger number of roots. One solution for preallocation can be found in `brute_force_root_finder_preallocate.m`.
"""

# ╔═╡ e03a4810-10ee-46a9-af7f-5e64dfe250a0
md"""
## Newton-Raphson’s Method

Famous and widely used in practice. NR method is fast but does not guarantee that a solution will be found. Idea: construct a series of linear equations (as we know how to solve them) and expect they will bring us to the solution of a nonlinear one. We will test function $f(x) = x^2 - 9$ on interval $[0, 1000]$.

* `naive_Newton.m`
* `demo_naive_Newton.m`
"""

# ╔═╡ 24edbd93-df12-456a-9ca2-e02fded15963
md"""
## Newton-Raphson’s Method

 $x_0$ is our initial guess. Next, we approximate the function at this point with a tangent line $\tilde{f}(x) = f(x_0) + f'(x_0)(x - x_0)$. Tangent slope is $f'(x_0)$ and it touches $f(x)$ curve at $x_0$. Tangent crosses the x-axis at $x = x_0 - \frac{f(x_0)}{f'(x_0)}$. This is new candidate point $x_1$.
"""

# ╔═╡ 05574013-0aa9-469d-839a-5b5d9ed8ef81
md"""
## Newton-Raphson’s Method

Starting with $x_0$ iterate $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$. For $n = 0, 1, 2, \ldots$. Until $|f(x_n)| < \varepsilon$ and $\varepsilon$ is a small number. Look at implementation – we are not using an array for $x_n$.
"""

# ╔═╡ 4c10d2c9-8b50-43ca-ad59-cf12839a86cd
md"""
## Test NR With Hyperbolic Tangent

Use $f(x) = \tanh x$ function (with derivative $f'(x) = 1 - \tanh^2 x$). Use initial conditions $x_0 = 1.08$ and $x_0 = 1.09$.
"""

# ╔═╡ dbd44a94-de01-49e9-b08c-602753b1a6c7
md"""
## Hands-on – Why Newton’s Method Fails?

Solve $\tanh x = 0$ by Newton’s method and study the intermediate details of the algorithm. Start with $x_0 = 1.08$. Plot the tangent in each iteration of Newton’s method. Then repeat the calculations and the plotting when $x_0 = 1.09$. Explain what you observe.

`plot_Newton_tanh.m` (in debug mode).
"""

# ╔═╡ a8e9887c-3190-41d4-83c3-cd1ee295bc09
md"""
## Robust Newton-Raphson

Newton-Raphson’s method can diverge. In such a case, our naïve implementation may run forever. Thus, we should limit the number of evaluations. Check if numerical values are numbers and are finite. Call function $f(x)$ only once at each call. Once again, let’s check the implementation. Newton’s method works better if $x_0$ is close to a real solution.

* `Newton.m`
* `demo_Newtons_method.m`
"""

# ╔═╡ 80ee2155-30b7-433b-ab17-de7dc2f73be8
md"""
## Symbolic Computations for Derivative

Newton-Raphson requires function derivative $f'(x)$. In some cases, this may be troublesome. But often, we can use symbolic computations to compute it. File: `symbolic_differentiation.m`

```matlab
syms x;               % define x as a mathematical symbol
f_expr = x^2 - 9;     % symbolic expression for f(x)
dfdx_expr = diff(f_expr);    % compute f’(x) symbolically
% Turn f_expr and dfdx_expr into plain Matlab functions
f = matlabFunction(f_expr);
dfdx = matlabFunction(dfdx_expr);
dfdx(5); % will print 10
```
"""

# ╔═╡ 3690c1e7-493a-4d4f-bcef-89db29a369c1
md"""
## Solving Multiple Nonlinear Algebraic Equations

- Often, we must solve multiple nonlinear algebraic equations at once.
- Newton-Raphson’s method is suitable to be extended for this case.
- Suppose we have $n$ nonlinear equations of $n$ variables:

$$\begin{align*} f_1(x_1, x_2, \ldots, x_n) &= 0 \\ f_2(x_1, x_2, \ldots, x_n) &= 0 \\ \vdots &= \vdots \\ f_n(x_1, x_2, \ldots, x_n) &= 0 \end{align*}$$

By introducing vector notation $\mathbf{f} = [f_1, f_2, \ldots, f_n]^T$, $\mathbf{x} = [x_1, x_2, \ldots, x_n]^T$.
"""

# ╔═╡ ac009506-61fe-4dfd-94b5-d0e6d53e8b7e
md"""
## Nomenclature for Matrix Calculus

- Matrices are boldface uppercase letters, e.g., $\mathbf{A}$, $\mathbf{J}$. 
- Column vectors (matrices) and algebraic vectors are boldface lowercase letters, e.g., $\mathbf{x}$, $\mathbf{b}$. 
- Scalars are, as usual, in lightface letters, e.g., $x$, $y$, $f$. 
-  $\mathbf{a}^T$ is a row vector (array). 
-  $\mathbf{0}$ is the zero matrix, $\mathbf{I}$ is the identity matrix.
"""

# ╔═╡ bda51d4d-4e02-456f-988b-3d067702ea3f
md"""
## Solving Multiple Nonlinear Algebraic Equations

We can write $f(x) = 0$. Having specific example as: 

$$\begin{align*} x^2 &= y - x \cos(\pi x) \\ yx + e^{-y} &= x^{-1} \end{align*}$$ 

With vector notation we can write $\mathbf{x} = \begin{bmatrix} x \\ y \end{bmatrix}$, $f(x) = \begin{bmatrix} x^2 - y + x \cos(\pi x) \\ yx + e^{-y} - x^{-1} \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$.
"""

# ╔═╡ fcaf9506-93e0-49bf-9b01-04a09f930671
md"""
## Taylor Expansion

- We want to follow an idea from one variable function: linearize it and find a root of that linearized function. 
- For vector function $\mathbf{f}(x)$ we can use first two terms of a Taylor series expansion. 
- Having $f$ and its derivative at point $x_i$ we can approximate it at some point $x_{i+1}$ by the two terms of a Taylor series around $x_i$: 

$$f(x_{i+1}) \approx f(x_i) + \mathbf{f(x_i)}(x_{i+1} - x_i)$$ 

- Good, but what is $\mathbf{\nabla f(x_i)}$?
"""

# ╔═╡ f883bc04-8e42-4077-9518-e5de695a41b2
md"""
## Jacobian Matrix

$\mathbf{\nabla f} = \frac{\partial f}{\partial x}$ is the matrix of all partial derivatives of $f$. 

$\mathbf{\nabla f}$ is often denoted by $\mathbf{J}$ and called Jacobian matrix: 

$$\mathbf{J} = \mathbf{\nabla f} = \frac{\partial f}{\partial x} = \begin{bmatrix} \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \ldots & \frac{\partial f_1}{\partial x_n} \\ \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \ldots & \frac{\partial f_2}{\partial x_n} \\ \vdots & \vdots & \ddots & \vdots \\ \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \ldots & \frac{\partial f_n}{\partial x_n} \end{bmatrix}$$ 

In general Jacobian matrix does not have to be square.
"""

# ╔═╡ da8211e3-578d-485c-bbd0-bafc9c80602f
md"""
## Jacobian Matrix

- For our example 

 $\mathbf{x} = \begin{bmatrix} x \\ y \end{bmatrix}$, $f(x) = \begin{bmatrix} x^2 - y + x \cos(\pi x) \\ yx + e^{-y} - x^{-1} \end{bmatrix}$. 

- Jacobian matrix is as follows: 

$$\mathbf{J} = \mathbf{\nabla f} = \begin{bmatrix} \frac{\partial f_1}{\partial x} & \frac{\partial f_1}{\partial y} \\ \frac{\partial f_2}{\partial x} & \frac{\partial f_2}{\partial y} \end{bmatrix} = \begin{bmatrix} 2x + \cos(\pi x) - \pi x \sin(\pi x) & -1 \\ y + x^{-2} & x - e^{-y} \end{bmatrix}$$
"""

# ╔═╡ ba859435-b135-4b60-b78b-5295ba7b5bf7
md"""
## Newton-Raphson’s Method

The $i$-th step of Newton-Raphson iteration consists of two steps:

1. Solve the linear system $\mathbf{J(x_i)}\delta = -f(x_i)$ with respect to $\delta$.
2. Set $x_{i+1} = x_i + \delta$.

Linear systems are usually solved using a variant of Gaussian elimination. Matlab uses the well-known LAPACK package for this purpose. In Matlab, use the backslash operator (mldivide) to solve $Ax = b$ by $x = A \backslash b$. Never-ever use the inverse – it is slower and less accurate. The Jacobian is often a sparse matrix (i.e., contains mostly zero elements). In such cases, special methods may be more efficient (often iterative).
"""

# ╔═╡ 70524fd5-c20d-420b-9fac-13a96a45c321
md"""
## Hands-on – Simple Mechanism

Use the Newton-Raphson method to find $\phi_2$ and $d$ when $\phi_1 = 30^\circ$.
Constraints are as follows: $$\Phi(u) = \begin{bmatrix} b \cos \phi_1 + a \cos \phi_2 - d \ b \sin \phi_1 + a \sin \phi_2 - r \end{bmatrix} = 0$$
Jacobian matrix is: $$\Phi_u = \begin{bmatrix} -a \sin \phi_2 & -1 \ a \cos \phi_2 & 0 \end{bmatrix} = 0$$
Source: P. E. Nikravesh, Computer-aided analysis.
"""

# ╔═╡ 41ccdb5e-34b2-429f-95db-4c863e6060b5
md"""
## Assignment 2 of 2

Having $a = OA$, $b = AB$, and $d = OB$: $$\mathbf{x} = \begin{bmatrix} \theta \ d \end{bmatrix} $$ $$ f(x) = \begin{bmatrix} a \cos \phi + b \cos \theta - d \ a \sin \phi - b \sin \theta \end{bmatrix}$$

Solve $f(x) = 0$ for $x$ for given $\phi$. Write a program that solves $f(x,t) = 0$ for $x$ using Newton-Raphson’s method. Use $\dot{f}(x,t) = 0$ to solve for $\dot{x}$. Loop for $t = \text{linspace}(0,1,101)$.

Source: P. E. Nikravesh, Computer-aided analysis.
"""

# ╔═╡ 542529ea-75f1-4dcc-b408-8d0c71c73ed0
md"""
## Assignment 2 of 2

Solve the problem from the previous slide using $\phi = \frac{\pi}{6} + \omega t$. Create plots of $t$ versus angle $\theta$, displacement $d$, and their time derivatives. Write a brief report (up to 2 A4 pages with 12 pt font) about your solution. Submit the report in PDF. Attach sources as a zip archive. Remember about proper code indent, comments, and meaningful variable and file names.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═0cd41f17-71d9-4ebc-868e-39741e37b13c
# ╠═40f4936b-5ec4-4c98-9f48-cc41054cfd80
# ╠═ce5a4e83-1a4b-4955-a6f3-67806b29b6cc
# ╠═f02d3fcd-826b-430e-992a-767d4bcdbb11
# ╠═e03a4810-10ee-46a9-af7f-5e64dfe250a0
# ╠═24edbd93-df12-456a-9ca2-e02fded15963
# ╠═05574013-0aa9-469d-839a-5b5d9ed8ef81
# ╠═4c10d2c9-8b50-43ca-ad59-cf12839a86cd
# ╠═dbd44a94-de01-49e9-b08c-602753b1a6c7
# ╠═a8e9887c-3190-41d4-83c3-cd1ee295bc09
# ╠═80ee2155-30b7-433b-ab17-de7dc2f73be8
# ╠═3690c1e7-493a-4d4f-bcef-89db29a369c1
# ╠═ac009506-61fe-4dfd-94b5-d0e6d53e8b7e
# ╠═bda51d4d-4e02-456f-988b-3d067702ea3f
# ╠═fcaf9506-93e0-49bf-9b01-04a09f930671
# ╠═f883bc04-8e42-4077-9518-e5de695a41b2
# ╠═da8211e3-578d-485c-bbd0-bafc9c80602f
# ╠═ba859435-b135-4b60-b78b-5295ba7b5bf7
# ╠═70524fd5-c20d-420b-9fac-13a96a45c321
# ╠═41ccdb5e-34b2-429f-95db-4c863e6060b5
# ╠═542529ea-75f1-4dcc-b408-8d0c71c73ed0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
