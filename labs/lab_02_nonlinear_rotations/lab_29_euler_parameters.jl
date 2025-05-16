### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 9db429e6-3198-11f0-12cb-87e8ec62e217
md"""
## Task: Local Frame Vectors and Orientation Matrix

You are given two vectors $\boldsymbol{a}$ and $\boldsymbol{b}$ defined along the positive $x$ and $y$ axes of a local frame. Their global components are:

$$\boldsymbol{a} =
\begin{bmatrix}
0.1107 \\
0.3924 \\
1.1286
\end{bmatrix}, \quad
\boldsymbol{b} =
\begin{bmatrix}
-1.9450 \\
1.5330 \\
-0.3422
\end{bmatrix}$$

Perform the following:

---

### (a) Orthogonality Check

Test whether the vectors $\boldsymbol{a}$ and $\boldsymbol{b}$ are orthogonal.

- What operation allows you to test this?
- What result would confirm orthogonality?

---

### (b) Construct Local Frame Unit Vectors

Determine the global components of the unit vectors: $\boldsymbol{u}_x$, $\boldsymbol{u}_y$, and $\boldsymbol{u}_z$.

---

### (c) Construct the Orientation Matrix $\boldsymbol{A}$

- Use the three unit vectors to build the $3 \times 3$ orientation matrix $\boldsymbol{A}$.
- This matrix maps local coordinates in the $xyz$ frame to global coordinates.

"""

# ╔═╡ 9f25e059-19d5-4a5a-be34-4926a9479ee4
md"""
## Task: Euler Parameters and Orientation Matrix

Determine the four Euler parameters for the transformation matrix $\boldsymbol{A}$ from above.

Test -matrix A for orthogonality, using the identity $\boldsymbol{A}^\top\boldsymbol{A} = \boldsymbol{I}$.
"""

# ╔═╡ 78bad3e1-056b-46f6-879f-e79d2d155601
md"""
## Task: Euler Parameters with reversed signs

Show that if the signs of all four Euler parameters are reversed, i.e., if $\boldsymbol{p} \rightarrow -\boldsymbol{p}$, then the transformation matrix $\boldsymbol{A}$ is not affected.
"""

# ╔═╡ 895ff4dc-8968-4694-b4c4-a55beac3900f
md"""
## Task: Vector Transformations Using Local Coordinates and Euler Parameters

Points $B$ and $C$ are defined in the local coordinate systems of two different rigid bodies:

- Point $B$ on body $i$ has local coordinates:  
  $$\boldsymbol{s}_i^B = [1.0,\ 1.0,\ -0.5]^T$$

- Point $C$ on body $j$ has local coordinates:  
  $$\boldsymbol{s}_j^C = [-2.0,\ 1.5,\ -1.0]^T$$

The origins of bodies $i$ and $j$ in global coordinates are:  
$$\boldsymbol{r}_i = [-1.2,\ 0.4,\ 3.1]^T, \quad
\boldsymbol{r}_j = [0.4,\ 4.5,\ 0.5]^T$$

Euler parameters (unit quaternions) describing the orientation of the bodies are:  
$$\boldsymbol{p}_i = [0.343,\ -0.564,\ 0.604,\ 0.447]^T, \quad
\boldsymbol{p}_j = [0.270,\ 0.732,\ -0.331,\ 0.531]^T$$

---

### Your Tasks

**(a)** Compute the **global coordinates** of points $B_i$ and $C_j$.

**(b)** Determine the **global components** of vector $\boldsymbol{d} = \overrightarrow{B_i C_j}$.

**(c)** Express the vector $\boldsymbol{d}$ in the **local frame** of body $i$  
(i.e., compute local components in the $x_i y_i z_i$ frame).

**(d)** Express the same vector in the **local frame** of body $j$  
(i.e., components in the $x_j y_j z_j$ frame).

---

Use quaternion-to-rotation matrix conversion and appropriate frame transformations to solve each step.
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
# ╟─9db429e6-3198-11f0-12cb-87e8ec62e217
# ╟─9f25e059-19d5-4a5a-be34-4926a9479ee4
# ╟─78bad3e1-056b-46f6-879f-e79d2d155601
# ╟─895ff4dc-8968-4694-b4c4-a55beac3900f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
