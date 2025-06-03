### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ bb4550ab-de3f-4890-9fdd-4e006a2d698d
begin
	# using DifferentialEquations
	# import PlotlyJS
	# using Plots
	# plotlyjs()
	# using LaTeXStrings
	using StaticArrays
	using PlutoUI
	using LinearAlgebra
	using Printf
	using Base64
	function embed_image(path_to_image::AbstractString; height::Integer=240, type::String="png")
	    img_bytes = read(path_to_image)
	    # Build a data URL: "data:image/{type};base64,{...}"
	    data_url = "data:image/$(type);base64," * base64encode(img_bytes)
    	attrs = (:height => height,)
	    PlutoUI.Resource(data_url, MIME("image/$(type)"), attrs)
	end
	using ForwardDiff
end;

# ╔═╡ b8f0fade-e453-42e3-9ecb-803d41328d10
md"""
# Programming Multibody Systems in Julia

### Lecture 9: Spatial Rigid Body Constraints Derivatives

Grzegorz Orzechowski
"""

# ╔═╡ 3d62a7d0-930d-461c-a4ef-ee38019d49b0
md"""
## Spatial Holonomic Constraints in Multibody Systems

- **Holonomic constraint**: $g(q,t)=0$ restricts the configuration space of a multibody system.

- **Levels of enforcement**  
  1. **Velocity**   $A(q,t)\,\dot q + b(q,t)=0$, $A=\partial g/\partial q$  
  2. **Acceleration**   $A\,\ddot q + \Gamma = 0$, $\Gamma = \dot A\,\dot q + \dot b$

- **Today’s focus**
  - Derive $A(q)$ and $\Gamma$ for:
    - orthogonality, point-on-line, coincident-point,
      simple, and driving constraints.
  - Express each Jacobian two ways  
    1. w.r.t. unit quaternions $q$ (configuration level)  
    2. w.r.t. body-frame angular velocity $\omega$ (velocity level)  
  - Convert between them via  
    $$\omega = 2\,G(q)\,\dot q \quad\Longleftrightarrow\quad
      \dot q = \tfrac12\,E(q)\,\omega$$
    where $G(q)$ (or $E(q)$) maps quaternion rates to angular velocity.

- **Goal**: clear recipes for $A$, its coordinate-space form, and  
  the right-hand side $\Gamma$—ready for implementation in Julia.

"""

# ╔═╡ aa843098-61a4-4702-9a8e-bfc55a81084d
md"""

## Quaternion ↔ Angular Velocity — both directions

**World-frame convention**

$$G(q)=
\begin{pmatrix}
-q_1 & -q_2 & -q_3\\
 q_0 & -q_3 & \;\;q_2\\
 q_3 & \;\;q_0 & -q_1\\
-q_2 & \;\;q_1 & \;\;q_0
\end{pmatrix}, \qquad
\dot q=\tfrac12\,G(q)\,\omega_W, \qquad
\omega_W = 2\,G(q)^{\!\top}\dot q.$$

**Body-frame convention**

$$E(q)=
\begin{pmatrix}
-q_1 & -q_2 & -q_3\\
 q_0 & \;\;q_3 & -q_2\\
-q_3 & \;\;q_0 & \;\;q_1\\
 q_2 & -q_1 & \;\;q_0
\end{pmatrix}, \qquad
\dot q=\tfrac12\,E(q)\,\omega_B, \qquad
\omega_B = 2\,E(q)^{\!\top}\dot q.$$

**Jacobian conversion**

If a velocity-level constraint is written as  

$$A_\omega\,\omega_B + b = 0,$$

then the coordinate-space Jacobian acting on $\dot q$ is  

$$A_q = A_\omega\,\bigl(2\,E(q)^{\!\top}\bigr), \qquad
A_q\,\dot q + b = 0.$$

(The same replacement with $2\,G(q)$ applies when $\omega_W$ is used.)

"""

# ╔═╡ 3c41f096-a7c3-4b93-ac8d-1f34808bdc25
md"""
## Why the $G(q)$ / $E(q)$ maps matter

- **Different layers need different coordinates**  
  *Integrators* evolve quaternions ($q,\dot q$) on $S^3$,  
  *constraint rows* are usually written in terms of angular velocity $\omega$ (a 3-vector).  
  We need a **bridge** to keep both views consistent.

- **Clean separation of roles**  
  $q$ — the configuration variable for Lie-group RK4  
  $\omega$ — the physically intuitive “spin” that appears in joint constraints and Kane/DAE forms.

- **Jacobian reuse**  
  1. Derive $A_\omega$ once from rigid-body geometry (easy, 3×3 blocks).  
  2. Convert to $A_q$ when assembling the full DAE:  
     $A_q = A_\omega\,(2\,G(q))$ (or $2\,E(q)$ for body-frame $\omega$).  
  3. Same $G(q)$ also converts measured $\dot q$ back to $\omega$ for logging.

- **Numerical stability**  
  Working in $\omega$ avoids quaternion-norm constraints at the velocity level,  
  while the Lie-group integrator keeps $q$ on $S^3$ at the configuration level.

*Bottom line*: $G(q)$ / $E(q)$ let us write concise constraint Jacobians in angular-velocity form and seamlessly plug them into a quaternion-based simulation without ever violating either representation.
"""

# ╔═╡ 560d04a2-4fc5-4ca8-86a5-d23000fefcc0
md"""
## How to map the $\Gamma$ term to coordinate space (body-frame ω)

Start with the acceleration-level constraint written in **twist space**

$$A_\omega\,\alpha_B + A_v\,a + \Gamma_\omega \;=\; r(t),
\qquad
\alpha_B,\omega_B \text{ in body frame.}$$

### 1. Express angular acceleration via quaternion coordinates  

Body-frame mapping  
$$\dot q \;=\; \tfrac12\,E(q)\,\omega_B, 
\qquad
\omega_B \;=\; 2\,E(q)^{\!\top}\dot q.$$

Differentiate once:

$$\alpha_B \;=\; 2\,E(q)^{\!\top}\ddot q 
              + 2\,\dot E(q)^{\!\top}\dot q.$$

### 2. Substitute α into the twist-space equation  

$$A_\omega\bigl(2E^{\!\top}\ddot q + 2\dot E^{\!\top}\dot q\bigr)
+ A_v\,a + \Gamma_\omega = r(t).$$

### 3. Collect coordinate-space terms  

- **Coordinate Jacobian**

  $$A_q \;=\; A_\omega\,(2\,E(q)^{\!\top}).$$

- **Coordinate-space Γ**

  $$\Gamma_q \;=\; \Gamma_\omega \;+\; A_\omega\,
      \bigl(2\,\dot E(q)^{\!\top}\dot q\bigr).$$

Resulting equation in coordinate space:

$$A_q\,\ddot q + A_v\,a + \Gamma_q \;=\; r(t).$$

### Practical notes

*  $E(q)$ is the $4\times3$ body-frame conversion matrix shown earlier.  
*  $\dot E(q)$ can be formed analytically or with automatic differentiation; it depends bilinearly on $q$ and $\dot q$.  
* Exactly the same recipe works for world-frame ω by replacing $E(q)$ with $G(q)$.  
* In code, compute `Γ_q = Γ_ω + A_ω * (2 * dEdt' * q̇)` once per step; reuse `A_q` already formed for the velocity level.

This completes the bridge: **$(A_\omega, \Gamma_\omega) \;\longrightarrow\; (A_q, \Gamma_q)$** so the same constraint can be enforced in Lie-group, quaternion-based integrators without loss of physical meaning.
"""

# ╔═╡ 8f9b2378-0a3b-4f9e-a835-6a9f1ea7df3b
md"""
## Elementary Implicit Constraints

The constraints of numerous joints can be formulated using the four elementary holonomic constraints specified below (Nikravesh [73, 74]).  They are expressed using the connection vector $\boldsymbol{l}$ between two body points $P_i$ and $P_j$, whose body-fixed positions are $\boldsymbol{c}_i$, $\boldsymbol{c}_j$ and whose global positions are $\boldsymbol{r}_i$, $\boldsymbol{r}_j$.

"""

# ╔═╡ 8633ca2c-ab27-4d71-b2a4-9de6b24068a0
md"""

## Type I: Coinciding Points (Ball Joint)

**Description:** A coincident point constraint requires that two points (one on each body) remain coincident in space. This is essentially a spherical joint (ball-and-socket joint) with 3 degrees of constraint. For instance, if points $P_i$ and $P_j$ coincide, there are 3 scalar constraints, $b_{ij}=3$.

$(embed_image("../../assets/figures/lecture_32_constraint_point.png", type="png", height=200))

### Position Level

The connection vector is zero:

$$\boldsymbol{g}_{ij} \equiv \boldsymbol{l} \equiv \boldsymbol{r}_j + \boldsymbol{c}_j - \boldsymbol{r}_i - \boldsymbol{c}_i = \boldsymbol{0}.$$

"""

# ╔═╡ 4c3754c8-1369-410e-a0d8-e22db2299bdb
md"""

### Velocity Level

Define

$$\dot{\boldsymbol{l}}
= \dot{\boldsymbol{r}}_j + \dot{\boldsymbol{c}}_j - \dot{\boldsymbol{r}}_i - \dot{\boldsymbol{c}}_i,
\quad
\dot{\boldsymbol{r}}_i = \boldsymbol{v}_i,
\quad
\dot{\boldsymbol{r}}_j = \boldsymbol{v}_j,$$

and note

$$\dot{\boldsymbol{c}}_i = \tilde{\boldsymbol{\omega}}_i\,\boldsymbol{c}_i = -\tilde{\boldsymbol{c}}_i\,\boldsymbol{\omega}_i,
\quad
\dot{\boldsymbol{c}}_j = \tilde{\boldsymbol{\omega}}_j\,\boldsymbol{c}_j = -\tilde{\boldsymbol{c}}_j\,\boldsymbol{\omega}_j.$$

Then the time-derivative of the constraint becomes

$$\dot{\boldsymbol{g}}_{ij}
= \dot{\boldsymbol{l}}
= \boldsymbol{v}_j - \tilde{\boldsymbol{c}}_j\,\boldsymbol{\omega}_j - \boldsymbol{v}_i + \tilde{\boldsymbol{c}}_i\,\boldsymbol{\omega}_i
= \boldsymbol{0}.$$

In block form:

$$\dot{\boldsymbol{g}}_{ij}
=
\bigl[\tilde{\boldsymbol{c}}_i \;\; -\boldsymbol{E}\bigr]
\begin{bmatrix}\boldsymbol{\omega}_i\\\boldsymbol{v}_i\end{bmatrix}
+
\bigl[-\,\tilde{\boldsymbol{c}}_j \;\; \boldsymbol{E}\bigr]
\begin{bmatrix}\boldsymbol{\omega}_j\\\boldsymbol{v}_j\end{bmatrix}
= \boldsymbol{0},$$

or more compactly

$$\boldsymbol{G}_i
\begin{bmatrix}\boldsymbol{\omega}_i\\\boldsymbol{v}_i\end{bmatrix}
+
\boldsymbol{G}_j
\begin{bmatrix}\boldsymbol{\omega}_j\\\boldsymbol{v}_j\end{bmatrix}
= \boldsymbol{0}.$$

"""

# ╔═╡ 2510042b-3ce0-4287-b5f2-ed0158f3decb
md"""

### Acceleration Level

Taking another time derivative, with $\dot{\boldsymbol{\omega}}_i=\boldsymbol{\alpha}_i$, $\dot{\boldsymbol{v}}_i=\boldsymbol{a}_i$, etc., gives

$$\ddot{\boldsymbol{g}}_{ij}
= \ddot{\boldsymbol{l}}
= \dot{\boldsymbol{v}}_j 
- \tilde{\boldsymbol{c}}_j\,\dot{\boldsymbol{\omega}}_j 
- \dot{\tilde{\boldsymbol{c}}}_j\,\boldsymbol{\omega}_j 
- \dot{\boldsymbol{v}}_i 
+ \tilde{\boldsymbol{c}}_i\,\dot{\boldsymbol{\omega}}_i
+ \dot{\tilde{\boldsymbol{c}}}_i\,\boldsymbol{\omega}_i 
= \boldsymbol{0}.$$

Defining the $(3\times1)$ “Coriolis” vector

$$\bar{\boldsymbol{\gamma}}_{ij}
\equiv
\dot{\tilde{\boldsymbol{c}}}_i\,\boldsymbol{\omega}_i
- \dot{\tilde{\boldsymbol{c}}}_j\,\boldsymbol{\omega}_j
=
\tilde{\boldsymbol{\omega}}_j\,\tilde{\boldsymbol{\omega}}_j\,\boldsymbol{c}_j
- \tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_i\,\boldsymbol{c}_i,$$

one obtains the final acceleration-level constraint.

---

  

**Note:** The coincident point constraint yields **3 equations** and removes 3 relative translational degrees of freedom between the bodies (allowing free relative rotation if no other constraints added). It is often combined with an orthogonality constraint to form a universal joint, or with a full orientation lock to form a welded joint.
"""

# ╔═╡ 28b5d088-3f17-40cd-acd5-cf85417c199c
md"""
## Type II: Constant Projection -- Type 2 Perpendicularity


The projection $l_{0}$ of the connection vector $\boldsymbol l$ of the points $P_i$ and $P_j$ onto the unit vector $\boldsymbol e_i$ on body $i$ is constant.

$(embed_image("../../assets/figures/lecture_32_constraint_projection.png", type="png", height=200))

"""

# ╔═╡ e1402bdf-815a-4544-a843-c23ef52cfd32
md"""

### Position Level

The scalar constraint reads  

$$g_{ij} \;\equiv\;\boldsymbol e_i^T\,\boldsymbol l \;-\; l_{0} \;=\; 0,
\qquad
\boldsymbol l \;=\;\boldsymbol r_j + \boldsymbol c_j - \boldsymbol r_i - \boldsymbol c_i.$$

"""

# ╔═╡ 8a4318e2-c69e-4491-acfa-18957513c8ae
md"""

### Velocity Level

Taking the time derivative of the position constraint gives  

$$\dot g_{ij}
\;\equiv\;
\dot{\boldsymbol e}_i^T\,\boldsymbol l
\;+\;
\boldsymbol e_i^T\,\dot{\boldsymbol l}
\;=\;0.$$

Here  

$$\dot{\boldsymbol e}_i \;=\;\tilde{\boldsymbol\omega}_i\,\boldsymbol e_i
\;=\;-\,\tilde{\boldsymbol e}_i\,\boldsymbol\omega_i,$$

and $\dot{\boldsymbol l}$ is as in the coinciding‐points case.  One can show this becomes  

$$\dot g_{ij}
\;=\;
\bigl[\,
\boldsymbol e_i^T\,(\tilde{\boldsymbol c}_i + \tilde{\boldsymbol l})
\;-\;\boldsymbol e_i^T
\,\bigr]
\begin{bmatrix}\boldsymbol\omega_i\\\boldsymbol v_i\end{bmatrix}
\;+\;
\bigl[\,-\,\boldsymbol e_i^T\,\tilde{\boldsymbol c}_j
\;\;\;\boldsymbol e_i^T\,\bigr]
\begin{bmatrix}\boldsymbol\omega_j\\\boldsymbol v_j\end{bmatrix}
\;=\;0,$$

or compactly  

$$\boldsymbol G_i\begin{bmatrix}\boldsymbol\omega_i\\\boldsymbol v_i\end{bmatrix}
\;+\;
\boldsymbol G_j\begin{bmatrix}\boldsymbol\omega_j\\\boldsymbol v_j\end{bmatrix}
\;=\;0.$$

"""

# ╔═╡ 1af6ee11-fd22-45cb-a0b2-a86cf220f734
md"""

### Acceleration Level

Differentiating again yields  

$$\ddot g_{ij}
\;\equiv\;
\ddot{\boldsymbol e}_i^T\,\boldsymbol l
\;+\;2\,\dot{\boldsymbol e}_i^T\,\dot{\boldsymbol l}
\;+\;\boldsymbol e_i^T\,\ddot{\boldsymbol l}
\;=\;0.$$

With  

$$\ddot{\boldsymbol e}_i
\;=\;
\tilde{\boldsymbol\omega}_i\,\dot{\boldsymbol e}_i
\;+\;\widetilde{\dot{\boldsymbol e}}_i\,\boldsymbol\omega_i,$$  

one obtains the acceleration‐level constraint in the form 

$$\boldsymbol G_i\begin{bmatrix}\dot{\boldsymbol\omega}_i\\\boldsymbol a_i\end{bmatrix}
\;+\;
\boldsymbol G_j\begin{bmatrix}\dot{\boldsymbol\omega}_j\\\boldsymbol a_j\end{bmatrix}
\;=\;
\bar{\bar{\boldsymbol\gamma}}_{ij},$$  

where the $(3\times1)$ vector of non‐homogeneous terms is  

$$\bar{\bar{\boldsymbol\gamma}}_{ij}
\;=\;
\boldsymbol l^T\,\widetilde{\dot{\boldsymbol e}}_i
\;+\;
\boldsymbol e_i^T\bigl(\tilde{\boldsymbol c}_i\,\boldsymbol\omega_i
\;-\;\tilde{\boldsymbol c}_j\,\boldsymbol\omega_j\bigr)
\;+\;2\,\dot{\boldsymbol e}_i^T\,\dot{\boldsymbol l}.$$

"""

# ╔═╡ f32c2302-6c1d-4b62-a935-9f07e6a1b7ca
md"""
## Type III: Constant Angle -- Orthogonality Constraint Generalization

The angle $\varphi$ between the unit vectors $\boldsymbol e_i$ on body $i$ and $\boldsymbol e_j$ on body $j$ is constant.

$(embed_image("../../assets/figures/lecture_32_constraint_constant_angle.png", type="png", height=200))

"""

# ╔═╡ 3d2ccdd8-06f7-4da8-bae3-8b368c7e8f5d
md"""

### Position Level

The scalar position constraint is  

$$g_{ij} \;\equiv\;\boldsymbol e_i^T\,\boldsymbol e_j \;-\;\cos\varphi \;=\;0.$$

"""

# ╔═╡ c4ee6262-53d8-4229-94c4-9205387742e8
md"""

### Velocity Level

Differentiating in time gives  

$$\dot g_{ij}
\;\equiv\;
\dot{\boldsymbol e}_i^T\,\boldsymbol e_j
\;+\;
\boldsymbol e_i^T\,\dot{\boldsymbol e}_j
\;=\;0.$$

With  

$$\dot{\boldsymbol e}_i = \tilde{\boldsymbol\omega}_i\,\boldsymbol e_i = -\,\tilde{\boldsymbol e}_i\,\boldsymbol\omega_i,
\quad
\dot{\boldsymbol e}_j = \tilde{\boldsymbol\omega}_j\,\boldsymbol e_j = -\,\tilde{\boldsymbol e}_j\,\boldsymbol\omega_j,$$

one arrives at the form  

$$\dot g_{ij}
\;=\;
\bigl[\;\boldsymbol e_i^T\,\tilde{\boldsymbol e}_j\;\;\;0^T\bigr]
\begin{bmatrix}\boldsymbol\omega_i\\\boldsymbol v_i\end{bmatrix}
\;+\;
\bigl[\;-\,\boldsymbol e_i^T\,\tilde{\boldsymbol e}_j\;\;\;0^T\bigr]
\begin{bmatrix}\boldsymbol\omega_j\\\boldsymbol v_j\end{bmatrix}
\;=\;0,$$

or compactly  

$$\boldsymbol G_i
\begin{bmatrix}\boldsymbol\omega_i\\\boldsymbol v_i\end{bmatrix}
\;+\;
\boldsymbol G_j
\begin{bmatrix}\boldsymbol\omega_j\\\boldsymbol v_j\end{bmatrix}
\;=\;0.$$

"""

# ╔═╡ 58a6879b-aabc-4f52-8619-fa9f10c19cea
md"""

### Acceleration Level

Another time‐derivative yields  

$$\ddot g_{ij}
\;\equiv\;
\ddot{\boldsymbol e}_i^T\,\boldsymbol e_j
\;+\;
\boldsymbol e_i^T\,\ddot{\boldsymbol e}_j
\;+\;
2\,\dot{\boldsymbol e}_i^T\,\dot{\boldsymbol e}_j
\;=\;0.$$

With  

$$\ddot{\boldsymbol e}_i = \tilde{\boldsymbol\omega}_i\,\dot{\boldsymbol e}_i + \widetilde{\dot{\boldsymbol e}}_i\,\boldsymbol\omega_i,
\quad
\ddot{\boldsymbol e}_j = \tilde{\boldsymbol\omega}_j\,\dot{\boldsymbol e}_j + \widetilde{\dot{\boldsymbol e}}_j\,\boldsymbol\omega_j,$$

this can be written in the form  

$$\boldsymbol G_i
\begin{bmatrix}\dot{\boldsymbol\omega}_i\\\boldsymbol a_i\end{bmatrix}
\;+\;
\boldsymbol G_j
\begin{bmatrix}\dot{\boldsymbol\omega}_j\\\boldsymbol a_j\end{bmatrix}
\;=\;\bar{\bar{\boldsymbol\gamma}}_{ij},$$

where the non‐homogeneous term is  

$$\bar{\bar{\boldsymbol\gamma}}_{ij}
\;=\;
\boldsymbol e_i^T\,\tilde{\boldsymbol\omega}_i\,\dot{\boldsymbol e}_i
\;+\;
\boldsymbol e_i^T\,\tilde{\boldsymbol\omega}_j\,\dot{\boldsymbol e}_j
\;+\;2\,\dot{\boldsymbol e}_i^T\,\dot{\boldsymbol e}_j
\;-\;2\,\boldsymbol e_i^T\,\tilde{\boldsymbol\omega}_j\,\boldsymbol\omega_j\;\boldsymbol e_j.$$

"""

# ╔═╡ 864382a1-e451-4e60-9e2f-8baef0c1f615
md"""

## Type IV: Simple Coordinate Constraint

**Description**  
Fix point $P$ on body $A$ to a constant global height $c$ along $\hat z$.

$$\Phi_s(q_A,X_A)=P_z-c=0,\qquad
P=X_A+R_A(q_A)\,r_{AP}.$$

---

### Velocity level  

$$\dot\Phi_s=\dot P_z=0.$$

Jacobian row acting on $\bigl[\;\omega_A^{\top}\;v_A^{\top}\bigr]^{\!\top}$:

$$A_\omega=
\begin{bmatrix}
(r_{AP}^G\times\hat z)^{\!\top}\ \ & 0\;0\;1
\end{bmatrix},
\qquad
A_\omega
\begin{bmatrix}\omega_A\\ v_A\end{bmatrix}=0.$$

---

### Acceleration level  

$$\ddot\Phi_s=\ddot P_z=0\;\;\Longrightarrow\;\;
A_\omega
\begin{bmatrix}\alpha_A\\ a_A\end{bmatrix}
+\Gamma = 0,$$

with the **same** $A_\omega$ and

$$\Gamma=(\omega_A\times(\omega_A\times r_{AP}^G))\!\cdot\hat z.$$

---

### Coordinate-space Jacobian  

Use $\dot q_A=\tfrac12\,E(q_A)\,\omega_A$:

$$A_q=
A_\omega
\begin{bmatrix}
2\,E(q_A)\\[2pt] I_3
\end{bmatrix},
\qquad
A_q
\begin{bmatrix}\dot q_A\\ v_A\end{bmatrix}=0.$$

---

**Interpretation**  
- Removes **one** translational DOF.  
- Rotational term $(r_{AP}^G\times\hat z)^{\!\top}\omega_A$ captures the vertical speed induced by body spin.  
- $\Gamma$ contains the familiar centripetal term, ensuring no vertical acceleration of $P$.
"""

# ╔═╡ 67d4d9ec-d011-41ce-a879-a21755daf0ca
md"""

## Type V: Driving Coordinate Constraint

**Description**  
Lock a single coordinate to a **time-law** $f(t)$ (rheonomic).  
Example: prescribe the height of point $P$ on body $A$:

$$\Phi_d(q_A,X_A,t)=P_z-f(t)=0.$$

---

### Velocity level  

Differentiate once:

$$\dot\Phi_d=\dot P_z-\dot f(t)=0.$$

Jacobian row on $\bigl[\;\omega_A^{\top}\;v_A^{\top}\bigr]^{\!\top}$:

$$A_\omega=
\begin{bmatrix}
\;(r_{AP}^G \times \hat z)^{\!\top} & 0\;0\;1
\end{bmatrix},
\qquad
A_\omega
\begin{bmatrix}\omega_A\\ v_A\end{bmatrix}
= \dot f(t).$$

---

### Acceleration level  

Second derivative:

$$\ddot\Phi_d=\ddot P_z-\ddot f(t)=0
\;\;\Longrightarrow\;\;
A_\omega
\begin{bmatrix}\alpha_A\\ a_A\end{bmatrix}
+ \Gamma = \ddot f(t),$$

with the same $A_\omega$ as above and  
$\Gamma=(r_{AP}^G\times\hat z)^{\!\top}(\omega_A\times r_{AP}^G)$.

---

### Coordinate-space form  

Convert to quaternion rates via  
$\dot q_A=\tfrac12\,E(q_A)\,\omega_A$:

$$A_q=
A_\omega
\begin{bmatrix}
2\,E(q_A) \\[2pt] I_3
\end{bmatrix},\qquad
A_q
\begin{bmatrix}\dot q_A\\ v_A\end{bmatrix}
= \dot f(t).$$

---

**Key points**

- Same geometric Jacobian as the simple constraint.  
- Time-law appears as **known RHS**: $\dot f(t)$ (velocity), $\ddot f(t)$ (acceleration).  
- Works for any prescribed distance, angle, or slider stroke by updating $f(t)$ and the direction vectors.

"""

# ╔═╡ fa578742-2dc4-439a-993b-a99ced9e84fa
md"""

## Implicit Constraints of the Revolute Joint

With the help of the four elementary constraints, the implicit constraints of numerous joint types can be formulated.  As an example, the $b_{ij} = 5$ implicit constraints of the revolute joint are set up.

$(embed_image("../../assets/figures/lecture_32_constraint_revolute.png", type="png", height=200))
"""

# ╔═╡ 98cf5350-ff12-47f6-a1ca-d72b7d68708b
md"""
### Position Level

- **Type I (Coinciding Points)**  
  The points $P_i$ on body $i$ and $P_j$ on body $j$ lie on the same axis of rotation, so  
  $$\boldsymbol{r}_j + \boldsymbol{c}_j \;-\; \boldsymbol{r}_i \;-\; \boldsymbol{c}_i \;=\; \boldsymbol{0},$$
  which contributes **3 constraints**.

- **Type III (Constant Angle)**  
  Let $\boldsymbol{u}_j$ be the unit vector along the rotation axis on body $j$. On body $i$, pick two orthonormal unit normals $\boldsymbol{m}_i,\boldsymbol{n}_i$ perpendicular to the joint axis. We require each of these normals to remain perpendicular to $\boldsymbol{u}_j$. Thus:
  $$\boldsymbol{m}_i^T\,\boldsymbol{u}_j \;=\; 0,$$
  $$\boldsymbol{n}_i^T\,\boldsymbol{u}_j \;=\; 0.$$
  Each of these is **1 constraint**, so altogether Type III adds **2 constraints**.

Putting them together, the $b_{ij}=5$ position‐level constraints for the revolute joint are  
$$\begin{bmatrix}
\underbrace{\boldsymbol{r}_j + \boldsymbol{c}_j - \boldsymbol{r}_i - \boldsymbol{c}_i}_{\text{Type I (3 equations)}}\\
\underbrace{\boldsymbol{m}_i^T\,\boldsymbol{u}_j}_{\text{Type III (1 equation)}}\\
\underbrace{\boldsymbol{n}_i^T\,\boldsymbol{u}_j}_{\text{Type III (1 equation)}}
\end{bmatrix}
\;=\;
\begin{bmatrix}
\boldsymbol{0}\\
0\\
0
\end{bmatrix}.$$

"""

# ╔═╡ e812ca78-78bd-4c82-9446-a9828eb3577f
md"""

### Velocity Level

At the velocity level, each of the five constraints yields one row in the combined Jacobian blocks $\boldsymbol G_i$ and $\boldsymbol G_j$. In block‐matrix form (cf. Eq. (9.15)), we write  
$$\dot{\boldsymbol{g}}_{ij}
\;\equiv\;
\underbrace{
  \begin{bmatrix}
    \tilde{\boldsymbol{c}}_i & -\boldsymbol{E} \\
    \boldsymbol{m}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T \\
    \boldsymbol{n}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T
  \end{bmatrix}
}_{\;=\;\boldsymbol G_i\;}
\begin{bmatrix}\boldsymbol{\omega}_i\\\boldsymbol{v}_i\end{bmatrix}
\;+\;
\underbrace{
  \begin{bmatrix}
    -\,\tilde{\boldsymbol{c}}_j & \boldsymbol{E} \\
    -\,\boldsymbol{m}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T \\
    -\,\boldsymbol{n}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T
  \end{bmatrix}
}_{\;=\;\boldsymbol G_j\;}
\begin{bmatrix}\boldsymbol{\omega}_j\\\boldsymbol{v}_j\end{bmatrix}
\;=\;\boldsymbol{0}.$$

Here:  
-  $\tilde{\boldsymbol{c}}_i$ and $\tilde{\boldsymbol{c}}_j$ are the skew‐symmetric matrices for vectors $\boldsymbol c_i,\boldsymbol c_j$.  
-  $\tilde{\boldsymbol{u}}_j$ is the skew‐matrix of the joint‐axis vector $\boldsymbol u_j$.  
-  $\boldsymbol{E}$ is the $3\times3$ identity.  
- The first block row (Type I) has $\tilde{\boldsymbol{c}}_i,\,-\tilde{\boldsymbol{c}}_j$ and $\boldsymbol E$ accordingly.  
- The next two rows (Type III) involve $\boldsymbol m_i^T\,\tilde{\boldsymbol u}_j$ and $\boldsymbol n_i^T\,\tilde{\boldsymbol u}_j$.

Thus,  
$$\boldsymbol G_i = 
\begin{bmatrix}
\tilde{\boldsymbol{c}}_i & -\boldsymbol E\\
\boldsymbol{m}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T\\
\boldsymbol{n}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T
\end{bmatrix},
\qquad
\boldsymbol G_j = 
\begin{bmatrix}
-\,\tilde{\boldsymbol{c}}_j & \boldsymbol E\\
-\,\boldsymbol{m}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T\\
-\,\boldsymbol{n}_i^T\,\tilde{\boldsymbol{u}}_j & \boldsymbol{0}^T
\end{bmatrix}.$$

"""

# ╔═╡ 9c0d0f7f-7466-421f-b35e-6f0a6e2694c4
md"""

### Acceleration Level

At the acceleration level (cf. Eq. (9.16)), we have 

$$\ddot{\boldsymbol{g}}_{ij}
=
\boldsymbol G_i
\begin{bmatrix}\dot{\boldsymbol{\omega}}_i\\\boldsymbol{a}_i\end{bmatrix}
+
\boldsymbol G_j
\begin{bmatrix}\dot{\boldsymbol{\omega}}_j\\\boldsymbol{a}_j\end{bmatrix}
\;=\;
\bar{\bar{\boldsymbol\gamma}}_{ij},$$

where the \((5\times1)\) non‐homogeneous term is constructed from the Type I Coriolis term (Eq. (9.22)) and the two Type III contributions (Eq. (9.36)). In block form:  

$$\bar{\bar{\boldsymbol\gamma}}_{ij}
=
\begin{bmatrix}
\tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_j\,\boldsymbol{c}_j \;-\; \tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_i\,\boldsymbol{c}_i \\[0.5em]
\boldsymbol{m}_i^T\,\bigl(\tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_i + \tilde{\boldsymbol{\omega}}_j\,\tilde{\boldsymbol{\omega}}_j - 2\,\tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_j\bigr)\,\boldsymbol{u}_j \\[0.5em]
\boldsymbol{n}_i^T\,\bigl(\tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_i + \tilde{\boldsymbol{\omega}}_j\,\tilde{\boldsymbol{\omega}}_j - 2\,\tilde{\boldsymbol{\omega}}_i\,\tilde{\boldsymbol{\omega}}_j\bigr)\,\boldsymbol{u}_j
\end{bmatrix}.$$

"""

# ╔═╡ 6c248888-faa8-437a-b375-2079c0795ed6
md"""

## Conclusion

- **Two equivalent Jacobians**
  -  $A_q$ — acts on quaternion‐rates $\dot q$ (configuration view).  
  -  $A_\omega$ — acts on angular velocity $\omega$ (physical twist view).  

- **Bridge** 

  $$A_q = A_\omega\,[\,2\,E(q)\,], \qquad
    \dot q = \tfrac12\,E(q)\,\omega$$

  (or $G(q)$ for world–frame $\omega$).

- **Workflow**
  1. Write holonomic constraint $g(q,X)=0$.  
  2. Differentiate once ⇒ $A_\omega\,\omega + A_v\,v + b = 0$.  
  3. Insert $2\,E(q)$ when you need the same row for quaternion integrators.  
  4. Acceleration level: $A_\omega\,\alpha + A_v\,a + \Gamma = 0$ with  
     $\Gamma = \dot A_\omega\,\omega + \dot A_v\,v + \dot b$.

- **Why keep both forms?**
  -  $A_q$ plugs into Lie-group RK4 / shooting methods on $S^3$.  
  -  $A_\omega$ aligns with textbook $\bigl[J\,\dot q = 0,\;J\,\ddot q+\dot J\,\dot q=0\bigr]$ DAE formulations and is more intuitive for joint design.

Understanding and switching between these two representations lets you derive constraints once and reuse them consistently across *kinematics*, *dynamics*, and *control* in quaternion-based multibody simulations.
"""

# ╔═╡ 23942a6b-f7b0-4646-989e-c576823b3437
md"""
## Julia Example for Ball Joint
"""

# ╔═╡ c0ced6cc-048d-4ff8-b414-e3755c3b8f3a
# ========== small helpers =====================================================
tilde(v::SVector{3}) = @SMatrix [ 0.0 -v[3]  v[2];
                                  v[3] 0.0  -v[1];
                                 -v[2] v[1]  0.0 ]

# ╔═╡ b428ec1e-e56e-47ac-b723-b033219b93b0
rotmat(q::SVector{4}) = begin
    q0,q1,q2,q3 = q
    @SMatrix [
        1-2(q2^2+q3^2)   2(q1*q2-q0*q3)    2(q1*q3+q0*q2);
        2(q1*q2+q0*q3)   1-2(q1^2+q3^2)    2(q2*q3-q0*q1);
        2(q1*q3-q0*q2)   2(q2*q3+q0*q1)    1-2(q1^2+q2^2)
    ]
end

# ╔═╡ 9c44685e-27e4-42e6-b362-9cf9eb9a45a5
# body-frame conversion E(q)  (4×3)  and its t-derivative
function Emat(q::SVector{4})
    q0,q1,q2,q3 = q
    @SMatrix [
        -q1  -q2  -q3;
         q0   q3  -q2;
        -q3   q0   q1;
         q2  -q1   q0
    ]
end

# ╔═╡ b448ee2f-465f-4d16-a62a-3bc2b8209153
Edot(q,q̇) = Emat(q̇)   # linear dependence on q

# ╔═╡ ad867600-289d-4452-acf2-e76c5885a05e
# ========== light container ===================================================
struct BodyState
    q  ::SVector{4}   # unit quaternion (scalar first)
    X  ::SVector{3}   # world position
    ωB ::SVector{3}   # body-frame angular velocity
    v  ::SVector{3}   # world linear velocity
end

# ╔═╡ f0a55f46-b7e6-43d1-a63f-eeec1f6d774e
# ========== internal pre-compute (shared by all functions) ====================
function _prep(bA::BodyState,rA::SVector{3},
               bB::BodyState,rB::SVector{3})
    RA = rotmat(bA.q);  RB = rotmat(bB.q)
    cA = RA*rA;         cB = RB*rB
    P_A = bA.X + cA;    P_B = bB.X + cB
    ωA_W = bA.ωB;    ωB_W = bB.ωB
    (;RA,RB,cA,cB,P_A,P_B,ωA_W,ωB_W)
end

# ╔═╡ 88c5f68d-85e2-41fc-a068-862da3708a5e
# ========== 1) constraint vector g ============================================
g_ball_joint(bA,rA,bB,rB) = _prep(bA,rA,bB,rB).P_B -
                             _prep(bA,rA,bB,rB).P_A

# ╔═╡ de3e3be0-359e-4ceb-94e2-50ab1bf3b051
# ========== 2) twist-space Jacobian Aω ========================================
function Aω_ball_joint(bA,rA,bB,rB)
    C = _prep(bA,rA,bB,rB)
    AωA =  tilde(C.cA)
    AωB =  tilde(C.cB)
    hcat(AωA, -I(3),  -AωB,  I(3))   # size 3×12
end

# ╔═╡ 50c5f06a-be83-45bd-af13-7ec7670432b1
# ========== 3) Γ in twist space ===============================================
function Γω_ball_joint(bA,rA,bB,rB)
    C = _prep(bA,rA,bB,rB)
    (C.ωA_W × (C.ωA_W × C.cA)) -
    (C.ωB_W × (C.ωB_W × C.cB))
end

# ╔═╡ a856b5d2-2d58-4324-a925-0ed884ba8545
# ========== 4) coordinate-space Jacobian Aq ===================================
function Aq_ball_joint(bA,rA,bB,rB)
    C = _prep(bA,rA,bB,rB)
    EA, EB = Emat(bA.q), Emat(bB.q)
    AωA =  tilde(C.cA)
    AωB = -tilde(C.cB)
    Aq = hcat(AωA*(2*EA'), -I(3),
          AωB*(2*EB'),  I(3))          # size 3×14
end

# ╔═╡ ad14cde5-752a-4b3b-a7b1-c58b0ac1d731
# ========== 5) Γ in coordinate space =========================================
function Γq_ball_joint(bA,rA,bB,rB)
    C = _prep(bA,rA,bB,rB)
    EA, EB = Emat(bA.q), Emat(bB.q)
    # quaternion rates
    q̇A = 0.5*EA*bA.ωB
    q̇B = 0.5*EB*bB.ωB
    # extra conversion term
    extraA = 2 .* (Edot(bA.q,q̇A)' * q̇A)
    extraB = 2 .* (Edot(bB.q,q̇B)' * q̇B)
    Γω_ball_joint(bA,rA,bB,rB) +
    (tilde(C.cA)*extraA + (-tilde(C.cB))*extraB)
end


# ╔═╡ 1a62acef-e214-44fb-9e1a-749412be82f6
# ========== demo ==============================================================
bodyA = BodyState(
    SVector(1.0,0,0,0),  SVector(0,0,0),
    SVector(0,0,0),      SVector(0,0,0)
)

# ╔═╡ 78883d92-3f12-42f0-8129-fb21615dc76c
bodyB = BodyState(
    SVector(cosd(15),0,sind(15),0),   # 30° about x
    SVector(1,0,0),
    SVector(0,1,0),               # body-frame spin around +y
    SVector(0,0,0)
)

# ╔═╡ 0852e9e1-a866-45c2-9dec-1a9896ed3013
EA, EB = Emat(bodyA.q), Emat(bodyB.q)

# ╔═╡ d98bd951-2924-436c-b335-479bf250f7ae
EB' * EB

# ╔═╡ 4e73e349-fb21-4acf-9b6f-3438994b7f88
EB * EB' + bodyB.q * bodyB.q'

# ╔═╡ 09a233c5-1dd6-4806-9372-889d9dadbe22
rA = @SVector [0.0,0.0,0.0]

# ╔═╡ 1f6ed89f-4d01-486a-91ee-117d0528bc07
rB = @SVector [1.0,0.0,0.0]

# ╔═╡ 61d37bbc-21b5-4c7f-aa9e-6e9884e5bbad
println("g   = ", g_ball_joint(bodyA,rA, bodyB,rB))

# ╔═╡ 0270587f-ad25-4dd0-9f87-2414b613364d
println("Aω  =\n", Aω_ball_joint(bodyA,rA, bodyB,rB))

# ╔═╡ eec92226-1cba-4bab-bb6f-9c961832d074
println("Γω  = ", Γω_ball_joint(bodyA,rA, bodyB,rB))

# ╔═╡ e7fa14a6-efcd-4ad2-a526-3b6b9868e26c
println("Aq  =\n", Aq_ball_joint(bodyA,rA, bodyB,rB))

# ╔═╡ 37804ea6-aaec-4700-985d-f852ebe0f110
println("Γq  = ", Γq_ball_joint(bodyA,rA, bodyB,rB))

# ╔═╡ 44f7be93-2ab5-412b-86a4-db7c2855504a
md"""
## Now, let's use ForwardDiff.jl
"""

# ╔═╡ 2bfaab7a-7e4f-4df2-84dc-b95f2155ea46
#---------------------------------------------------------------- pack/unpack ---
pack(b::BodyState) = vcat(b.q, b.X)              # 7 numbers / body

# ╔═╡ c820b7e8-4913-45a9-bb92-38200d7810ea
pack(bA, bB)   = vcat(pack(bA), pack(bB))    # 14

# ╔═╡ ba3965ee-93de-44ff-bb41-8665fc4e6521
function unpack(p, templA, templB)
    qA = SVector{4}(p[1:4]);   XA = SVector{3}(p[5:7])
    qB = SVector{4}(p[8:11]);  XB = SVector{3}(p[12:14])
    BodyState(qA,XA,templA.ωB,templA.v),
    BodyState(qB,XB,templB.ωB,templB.v)
end

# ╔═╡ 1427f884-857a-4c5a-87b8-78908d8f53b8
# --- quick StaticArrays block-diagonal ----------------------------------------
function blkdiag(mats::SMatrix...)
    rows = sum(size(M,1) for M in mats)
    cols = sum(size(M,2) for M in mats)
    out  = zeros(MMatrix{rows,cols,Float64})  # ← mutable
    r = c = 1
    for M in mats
        out[r:r+size(M,1)-1, c:c+size(M,2)-1] .= M
        r += size(M,1);  c += size(M,2)
    end
    SMatrix(out)                              # convert back to immutable
end

# ╔═╡ e6b81f40-f242-4915-b1ef-49afdd932288
#------------------------------------------------- Aq via ForwardDiff.jacobian --
function Aq_AD(bA, rA, bB, rB)
    p0   = pack(bA, bB)
    gfun = pvec -> begin
        bA′,bB′ = unpack(pvec,bA,bB)
        g_ball_joint(bA′,rA,bB′,rB)
    end
    ForwardDiff.jacobian(gfun, p0) |> SMatrix{3,14}
end

# ╔═╡ d227a20f-3639-431a-832a-50d84189c883
function Aω_AD(bA, rA, bB, rB)
    EA, EB = Emat(bA.q), Emat(bB.q)
	Aq = Aq_AD(bA, rA, bB, rB)
    T = blkdiag(0.5*EA, SMatrix{3,3}(I(3)), 0.5*EB, SMatrix{3,3}(I(3)))   # maps [ω_B;v]→[q̇;v]
    Aq * T                                      # 3×12
end

# ╔═╡ d31117c7-d025-42f3-8d6e-73e7c28b1904
#-------------------------------- Γq via Hessian contraction  (immutable-safe) --
function Γq_AD(bA,rA,bB,rB)
    p0  = pack(bA,bB)

    dqA = 0.5 * Emat(bA.q) * bA.ωB
    dqB = 0.5 * Emat(bB.q) * bB.ωB
    ṗ   = vcat(dqA, bA.v, dqB, bB.v)             # 14-vector (SVector)

    gfun = p -> begin
        bA′,bB′ = unpack(p,bA,bB)
        g_ball_joint(bA′,rA,bB′,rB)
    end

    Γ = SVector{3,Float64}( ntuple(k -> begin
            Hk = ForwardDiff.hessian(pp -> gfun(pp)[k], p0)
            dot(ṗ,  Hk * ṗ)                      # ṗᵀ Hk ṗ
        end, 3) )
    return Γ
end

# ╔═╡ 0fd04dfb-7008-4dc4-a395-21a50e6e6622
println("Aq_AD  =\n", Aq_AD(bodyA, rA, bodyB, rB))

# ╔═╡ c20e7849-ae12-4d47-947c-4dfe3309d86d
Aq_ball_joint(bodyA,rA, bodyB,rB) - Aq_AD(bodyA, rA, bodyB, rB)

# ╔═╡ 44d57b3d-4404-4d6a-a0fe-78548c392fd0
println("Aω_AD  =\n", Aω_AD(bodyA,rA, bodyB,rB))

# ╔═╡ bb0bf79d-3c4a-421f-aa83-b080efeb199d
Aω_ball_joint(bodyA,rA, bodyB,rB) - Aω_AD(bodyA,rA, bodyB,rB)

# ╔═╡ c5373ff8-42f2-4230-a575-26ee93eb7753
println("Γq  = ", Γq_AD(bodyA, rA, bodyB, rB))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Base64 = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
ForwardDiff = "~0.10.38"
PlutoUI = "~0.7.62"
StaticArrays = "~1.9.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "8f6bbb71f0e3ceda94f87040cf321bd1686786f7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "d3de2694b52a01ce61a036f18ea9c0f61c4a9230"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.62"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─b8f0fade-e453-42e3-9ecb-803d41328d10
# ╠═bb4550ab-de3f-4890-9fdd-4e006a2d698d
# ╟─3d62a7d0-930d-461c-a4ef-ee38019d49b0
# ╟─aa843098-61a4-4702-9a8e-bfc55a81084d
# ╟─3c41f096-a7c3-4b93-ac8d-1f34808bdc25
# ╟─560d04a2-4fc5-4ca8-86a5-d23000fefcc0
# ╟─8f9b2378-0a3b-4f9e-a835-6a9f1ea7df3b
# ╟─8633ca2c-ab27-4d71-b2a4-9de6b24068a0
# ╟─4c3754c8-1369-410e-a0d8-e22db2299bdb
# ╟─2510042b-3ce0-4287-b5f2-ed0158f3decb
# ╟─28b5d088-3f17-40cd-acd5-cf85417c199c
# ╟─e1402bdf-815a-4544-a843-c23ef52cfd32
# ╟─8a4318e2-c69e-4491-acfa-18957513c8ae
# ╟─1af6ee11-fd22-45cb-a0b2-a86cf220f734
# ╟─f32c2302-6c1d-4b62-a935-9f07e6a1b7ca
# ╟─3d2ccdd8-06f7-4da8-bae3-8b368c7e8f5d
# ╟─c4ee6262-53d8-4229-94c4-9205387742e8
# ╟─58a6879b-aabc-4f52-8619-fa9f10c19cea
# ╟─864382a1-e451-4e60-9e2f-8baef0c1f615
# ╟─67d4d9ec-d011-41ce-a879-a21755daf0ca
# ╟─fa578742-2dc4-439a-993b-a99ced9e84fa
# ╟─98cf5350-ff12-47f6-a1ca-d72b7d68708b
# ╟─e812ca78-78bd-4c82-9446-a9828eb3577f
# ╟─9c0d0f7f-7466-421f-b35e-6f0a6e2694c4
# ╟─6c248888-faa8-437a-b375-2079c0795ed6
# ╟─23942a6b-f7b0-4646-989e-c576823b3437
# ╠═c0ced6cc-048d-4ff8-b414-e3755c3b8f3a
# ╠═b428ec1e-e56e-47ac-b723-b033219b93b0
# ╠═9c44685e-27e4-42e6-b362-9cf9eb9a45a5
# ╠═b448ee2f-465f-4d16-a62a-3bc2b8209153
# ╠═ad867600-289d-4452-acf2-e76c5885a05e
# ╠═f0a55f46-b7e6-43d1-a63f-eeec1f6d774e
# ╠═88c5f68d-85e2-41fc-a068-862da3708a5e
# ╠═de3e3be0-359e-4ceb-94e2-50ab1bf3b051
# ╠═50c5f06a-be83-45bd-af13-7ec7670432b1
# ╠═a856b5d2-2d58-4324-a925-0ed884ba8545
# ╠═ad14cde5-752a-4b3b-a7b1-c58b0ac1d731
# ╠═1a62acef-e214-44fb-9e1a-749412be82f6
# ╠═78883d92-3f12-42f0-8129-fb21615dc76c
# ╠═0852e9e1-a866-45c2-9dec-1a9896ed3013
# ╠═d98bd951-2924-436c-b335-479bf250f7ae
# ╠═4e73e349-fb21-4acf-9b6f-3438994b7f88
# ╠═09a233c5-1dd6-4806-9372-889d9dadbe22
# ╠═1f6ed89f-4d01-486a-91ee-117d0528bc07
# ╠═61d37bbc-21b5-4c7f-aa9e-6e9884e5bbad
# ╠═0270587f-ad25-4dd0-9f87-2414b613364d
# ╠═eec92226-1cba-4bab-bb6f-9c961832d074
# ╠═e7fa14a6-efcd-4ad2-a526-3b6b9868e26c
# ╠═37804ea6-aaec-4700-985d-f852ebe0f110
# ╟─44f7be93-2ab5-412b-86a4-db7c2855504a
# ╠═2bfaab7a-7e4f-4df2-84dc-b95f2155ea46
# ╠═c820b7e8-4913-45a9-bb92-38200d7810ea
# ╠═ba3965ee-93de-44ff-bb41-8665fc4e6521
# ╠═1427f884-857a-4c5a-87b8-78908d8f53b8
# ╠═e6b81f40-f242-4915-b1ef-49afdd932288
# ╠═d227a20f-3639-431a-832a-50d84189c883
# ╠═d31117c7-d025-42f3-8d6e-73e7c28b1904
# ╠═0fd04dfb-7008-4dc4-a395-21a50e6e6622
# ╠═c20e7849-ae12-4d47-947c-4dfe3309d86d
# ╠═44d57b3d-4404-4d6a-a0fe-78548c392fd0
# ╠═bb0bf79d-3c4a-421f-aa83-b080efeb199d
# ╠═c5373ff8-42f2-4230-a575-26ee93eb7753
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
