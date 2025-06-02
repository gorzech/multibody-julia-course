### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ de59c8ec-b5e1-4adf-bacd-c10258eef778
begin
	# using DifferentialEquations
	# import PlotlyJS
	# using Plots
	# plotlyjs()
	# using LaTeXStrings
	# using StaticArrays
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
	# using ForwardDiff
end;

# ╔═╡ 0df19b60-3f9c-11f0-2860-5b7ac0586d16
md"""
# Programming Multibody Systems in Julia

### Lecture 10: Spatial Rigid Multibody Dynamics

Grzegorz Orzechowski
"""

# ╔═╡ 3ec42369-8374-422d-8366-506a2e843787
md"""
## Dynamics – Equations of Motion
Particle: $$\boldsymbol{f} = m \ddot{\boldsymbol{r}}$$

System of particles: $$\boldsymbol{m} \ddot{\boldsymbol{r}} = \boldsymbol{f}$$

Rigid body: $$\boldsymbol{M} \ddot{\boldsymbol{q}} = \boldsymbol{f}$$

Where $\boldsymbol{M}$ is the mass matrix, $\boldsymbol{q}$ are coordinates, and $\boldsymbol{f}$ is the vector of generalized forces.
"""

# ╔═╡ 942e3726-a441-4774-b54b-d4d1ea7ecafe
md"""

## Orientation Representation: Unit Quaternions (Euler Parameters)
- Use **unit quaternions** (4 parameters) to represent 3D orientation without singularities. A unit quaternion $q = [q_0,\,q_1,\,q_2,\,q_3]$ (with $\lVert q\rVert = 1$) corresponds to a rotation in $\mathbb{R}^3$.
- Quaternions avoid gimbal lock by using redundancy (one extra parameter) and an implicit unit-length constraint. The three independent parameters (Euler *principal* parameters) cover the full rotation group $\mathrm{SO}(3)$.
- The unit quaternion evolves with the body’s angular velocity $\boldsymbol{\omega}$ via the kinematic equation:
  $$\dot{q} \;=\; \frac{1}{2}\;q \otimes (0,\,\boldsymbol{\omega}),$$
  (treating $\boldsymbol{\omega}=[\omega_x,\omega_y,\omega_z]$ as a pure quaternion). This couples the quaternion’s rate to the angular velocity.
- **Normalization:** During numerical integration, $q$ should be re-normalized periodically (or enforced via an algebraic constraint) to maintain $\lVert q\rVert = 1$. This prevents drift from the unit hypersphere.

"""

# ╔═╡ f9a0bb53-035b-4dfc-9a42-db3068eccf17
md"""

## Translational Coordinates and Spatial Velocities
- Use **Cartesian coordinates** $\mathbf{r} = [x,\,y,\,z]^T$ for a body’s position. Linear motion is simply the time derivative: $\dot{\mathbf{r}} = \mathbf{v}$, where $\mathbf{v}$ is the **linear velocity**.
- A free rigid body in space has 6 degrees of freedom (3 translational, 3 rotational). We describe rotational motion with the **angular velocity** vector $\boldsymbol{\omega}$ in $\mathbb{R}^3$. Together, $(\mathbf{v},\,\boldsymbol{\omega})$ form a **spatial velocity** (twist) describing the body’s motion.
- **Linear vs.\ Angular Velocity:** The linear velocity refers to the center‐of‐mass motion (translational), while angular velocity describes the rate of change of orientation (rotational) about the body’s center of mass. These are independent components of the state.
- Kinematically, $\ddot{\mathbf{r}} = \mathbf{a}$ (acceleration) and $\dot{\boldsymbol{\omega}} = \boldsymbol{\alpha}$ (angular acceleration) describe changes in $\mathbf{v}$ and $\boldsymbol{\omega}$ respectively. Both will appear in the equations of motion.

"""

# ╔═╡ 0ff028b1-8431-492e-8a77-403f67c9b212
md"""
- **Single spatial rigid body** (origin at center of mass):
  
  In velocity‐level coordinates, let  
  $$\mathbf{r} = \begin{bmatrix}x \\ y \\ z\end{bmatrix}, 
    \qquad
    \boldsymbol{\omega} = \begin{bmatrix}\omega_x \\ \omega_y \\ \omega_z\end{bmatrix},$$
  be the linear and angular velocities, respectively.  Then the Newton–Euler equations can be written in block‐matrix form as

  $$\underbrace{
      \begin{bmatrix}
        m\,\mathbf{I}_{3\times3} & \mathbf{0}_{3\times3} \\[6pt]
        \mathbf{0}_{3\times3}      & \mathbf{I}_{3\times3}
      \end{bmatrix}
    }_{\mathbf{M}}
    \,
    \underbrace{
      \begin{bmatrix}
        \ddot{\mathbf{r}} \\[4pt]
        \dot{\boldsymbol{\omega}}
      \end{bmatrix}
    }_{\displaystyle 
      \begin{matrix}\text{linear\,acceleration}\\[-3pt]\text{angular\,acceleration}\end{matrix}
    }
    \;=\;
    \underbrace{
      \begin{bmatrix}
        \mathbf{F}_{\mathrm{ext}} \\[4pt]
        \mathbf{M}_{\mathrm{ext}}
      \end{bmatrix}
    }_{\displaystyle 
      \begin{matrix}\text{net external}\\[-3pt]\text{force on COM}\end{matrix}
    }\,.$$

  - Here $m$ is the mass, $\mathbf{I}$ is the $3\times3$ inertia tensor about the COM (in body coordinates), and  
    $\mathbf{I}_{3\times3}$ is the $3\times3$ identity.  
  - $\mathbf{F}_{\mathrm{ext}} \in \mathbb{R}^3$ is the total external force (e.g.\ gravity, contact),  
    $\mathbf{M}_{\mathrm{ext}} \in \mathbb{R}^3$ is the total external moment (torque) **about** the COM.

"""

# ╔═╡ fe72f0b0-4a2c-4dc5-b5e6-8f6e8af815b0
md"""

## Equations of Motion: Unconstrained Rigid Body
- **Newton–Euler Equations:** For a rigid body of mass $m$ and inertia tensor $\mathbf{I}$ (about its center of mass), the equations of motion split into translation and rotation:
  - **Translational:**
    $$m\,\ddot{\mathbf{r}} \;=\; \mathbf{F}_{\mathrm{ext}},$$
    where $\mathbf{F}_{\mathrm{ext}}$ is the sum of external forces (e.g.\ gravity) on the body’s center of mass.
  - **Rotational:**
    $$\mathbf{I}\,\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times (\mathbf{I}\,\boldsymbol{\omega}) \;=\; \mathbf{M}_{\mathrm{ext}}.$$
    This is Euler’s rotation equation, relating the angular acceleration to the applied external moments $\mathbf{M}_{\mathrm{ext}}$ (torques) and the gyroscopic term $\boldsymbol{\omega} \times (\mathbf{I}\,\boldsymbol{\omega})$.
- The orientation quaternion $q$ enters indirectly: $\mathbf{I}$ can be taken constant in body coordinates (principal axes), while the rotation matrix $R(q)$ (from body to world) updates to track orientation. The quaternion evolves as described earlier, and $R = R(q)$ is used to map vectors between frames.
- **Energy and Momentum:** In absence of external forces/torques, linear momentum $m\,\mathbf{v}$ and angular momentum $\mathbf{L}=\mathbf{I}\,\boldsymbol{\omega}$ are conserved. Rotation about a principal axis with no external torque results in steady rotation (free rigid body motion or Euler’s top).

"""

# ╔═╡ 43ebdfe9-b0ce-4656-9592-d66d52570b8b
md"""
## Dynamics – Equations of Motion

- For spatial rigid‐body dynamics,  

  $$\mathbf{M}\,\ddot{\mathbf{q}}
    \;+\;
    \begin{bmatrix}
      \mathbf{0}_{3\times 1} \\[6pt]
      \,\boldsymbol{\omega} \times \bigl(\mathbf{J}\,\boldsymbol{\omega}\bigr)\,
    \end{bmatrix}
    \;=\;
    \begin{bmatrix}
      f_x \\[4pt]
      f_y \\[4pt]
      f_z \\[4pt]
      M_x \\[4pt]
      M_y \\[4pt]
      M_z
    \end{bmatrix},$$

  where

  $$\mathbf{M}
    \;=\;
    \begin{bmatrix}
      m\,\mathbf{I}_{3\times 3} & \mathbf{0}_{3\times 3} \\[6pt]
      \mathbf{0}_{3\times 3}     & \mathbf{J}
    \end{bmatrix},
    \qquad
    \dot{\mathbf{q}} \;=\;
    \begin{bmatrix}
      \dot x \\[4pt]
      \dot y \\[4pt]
      \dot z \\[4pt]
      \omega_x \\[4pt]
      \omega_y \\[4pt]
      \omega_z
    \end{bmatrix}.$$

  - Here:
    -  $m$ is the mass of the body.
    -  $\mathbf{I}_{3\times 3}$ is the $3\times 3$ identity matrix.
    -  $\mathbf{J}\in\mathbb{R}^{3\times 3}$ is the full inertia tensor about the center of mass (expressed in body coordinates).
    -  $(f_x,\,f_y,\,f_z)$ is the net external force on the center of mass.
    -  $(M_x,\,M_y,\,M_z)$ is the net external moment (torque) about the center of mass.
    - The term $\boldsymbol{\omega}\times(\mathbf{J}\,\boldsymbol{\omega})$ is the gyroscopic (Euler) term for rotational dynamics.
- To make it work, place the origin at the body’s center of mass!
"""

# ╔═╡ 5f8a633d-f1b4-406c-b732-30ff783cc040
md"""
## Basic Forces

- **Gravitational Force**  
  A uniform gravitational field $\mathbf{g}_{\mathrm{acc}}\in\mathbb{R}^3$ acting on a rigid body of mass $m$ produces a spatial generalized force (force + moment) about the center of mass:  

  $$\mathbf{f}^{\mathrm{gravity}}
    \;=\;
    \begin{bmatrix}
      m\,\mathbf{g}_{\mathrm{acc}} \\[6pt]
      \mathbf{r}_{\mathrm{com}}\times \bigl(m\,\mathbf{g}_{\mathrm{acc}}\bigr)
    \end{bmatrix}.$$

  - If the origin of body coordinates is placed at the center of mass ($\mathbf{r}_{\mathrm{com}}=\mathbf{0}$), then  

    $$\mathbf{f}^{\mathrm{gravity}}
      \;=\;
      \begin{bmatrix}
        m\,\mathbf{g}_{\mathrm{acc}} \\[4pt]
        \mathbf{0}
      \end{bmatrix}.$$

"""

# ╔═╡ 03c8a932-8f2a-4725-827d-4eee7e0f59b7
md"""

- **Single Force**  
  A force $\mathbf{f}_i = [\,f_x,\;f_y,\;f_z\,]^T$ applied at a point $P$ with position vector $\mathbf{r}_P$ (in world coordinates) yields the spatial generalized force:  

  $$\mathbf{f}^{\mathrm{single}}
    \;=\;
    \begin{bmatrix}
      \,f_x \\[4pt]
      \,f_y \\[4pt]
      \,f_z \\[6pt]
      \mathbf{r}_P \times \mathbf{f}_i
    \end{bmatrix}
    \;=\;
    \begin{bmatrix}
      f_x \\[4pt]
      f_y \\[4pt]
      f_z \\[6pt]
      r_{P,y}\,f_z \;-\; r_{P,z}\,f_y \\[4pt]
      r_{P,z}\,f_x \;-\; r_{P,x}\,f_z \\[4pt]
      r_{P,x}\,f_y \;-\; r_{P,y}\,f_x
    \end{bmatrix}.$$

"""

# ╔═╡ 96cf32e2-23f1-44b4-b7bd-fee1a0734522
md"""

- **Single Torque (Moment)**  
  A pure moment (torque) $\boldsymbol{n} = [\,n_x,\;n_y,\;n_z\,]^T$ about the center of mass introduces no net force but only a moment:  

  $$\mathbf{f}^{\mathrm{torque}}
    \;=\;
    \begin{bmatrix}
      0 \\[4pt]
      0 \\[4pt]
      0 \\[6pt]
      n_x \\[4pt]
      n_y \\[4pt]
      n_z
    \end{bmatrix}.$$

- **Summary of Spatial “Basic Forces”**  
  In block‐vector form, any external action can be written as  

  $$\mathbf{f}
    \;=\;
    \begin{bmatrix}
      \mathbf{F}_{\mathrm{net}} \\[6pt]
      \mathbf{M}_{\mathrm{net}}
    \end{bmatrix}
    \;=\;
    \begin{bmatrix}
      \sum_i \mathbf{f}_i \;+\; m\,\mathbf{g}_{\mathrm{acc}} \\[6pt]
      \sum_i \bigl(\mathbf{r}_i \times \mathbf{f}_i \bigr)
      \;+\;\sum_j \boldsymbol{n}_j 
      \;+\;\mathbf{r}_{\mathrm{com}}\times \bigl(m\,\mathbf{g}_{\mathrm{acc}}\bigr)
    \end{bmatrix}.$$

  - Here $\mathbf{f}_i$ are the applied forces at points $\mathbf{r}_i$, and $\boldsymbol{n}_j$ are any pure torques about the center of mass.  
  - If the body origin is at the COM, then the gravity‐induced moment $\mathbf{r}_{\mathrm{com}}\times (m\,\mathbf{g}_{\mathrm{acc}})$ vanishes.
"""

# ╔═╡ 436cec7e-3e97-48da-b0e5-f296ba6eb14e
md"""
## Equations of Motion

- **System of Unconstrained Bodies**  
  For a collection of $n$ rigid bodies in 3D (each with 6 DOF), the block‐diagonal mass‐inertia matrix and acceleration vector satisfy  
  
$$\begin{bmatrix}
      \mathbf{M}_{1} & \mathbf{0}        & \cdots & \mathbf{0}        \\[6pt]
      \mathbf{0}     & \mathbf{M}_{2}    & \cdots & \mathbf{0}        \\[6pt]
      \vdots         & \vdots            & \ddots & \vdots            \\[6pt]
      \mathbf{0}     & \mathbf{0}        & \cdots & \mathbf{M}_{n}
    \end{bmatrix}
    \,
    \begin{bmatrix}
      \ddot{\mathbf{q}}_{1} \\[4pt]
      \ddot{\mathbf{q}}_{2} \\[4pt]
      \vdots                \\[4pt]
      \ddot{\mathbf{q}}_{n}
    \end{bmatrix}
    \;+\;
    \begin{bmatrix}
      \displaystyle
      \begin{bmatrix}
        \mathbf{0}_{3} \\[4pt]
        \boldsymbol{\omega}_{1}\times\bigl(\mathbf{J}_{1}\,\boldsymbol{\omega}_{1}\bigr)
      \end{bmatrix} \\[12pt]
      \displaystyle
      \begin{bmatrix}
        \mathbf{0}_{3} \\[4pt]
        \boldsymbol{\omega}_{2}\times\bigl(\mathbf{J}_{2}\,\boldsymbol{\omega}_{2}\bigr)
      \end{bmatrix} \\[12pt]
      \vdots \\[12pt]
      \displaystyle
      \begin{bmatrix}
        \mathbf{0}_{3} \\[4pt]
        \boldsymbol{\omega}_{n}\times\bigl(\mathbf{J}_{n}\,\boldsymbol{\omega}_{n}\bigr)
      \end{bmatrix}
    \end{bmatrix}
    \;=\;
    \begin{bmatrix}
      \mathbf{f}_{1} \\[4pt]
      \mathbf{f}_{2} \\[4pt]
      \vdots         \\[4pt]
      \mathbf{f}_{n}
    \end{bmatrix}
    .$$

"""

# ╔═╡ c3a1b76e-cee8-4036-85c3-ced87079d14a
md"""

- For each body $i$:  
    -  $\displaystyle \mathbf{M}_{i} \;=\; 
      \begin{bmatrix}
        m_{i}\,\mathbf{I}_{3\times 3} & \mathbf{0}_{3\times 3} \\[6pt]
        \mathbf{0}_{3\times 3}         & \mathbf{J}_{i}
      \end{bmatrix},$  
      where $m_{i}$ is the mass and $\mathbf{J}_{i}\in\mathbb{R}^{3\times 3}$ is the full inertia tensor about the center of mass (in body‐fixed coordinates).  
    -  $\displaystyle \dot{\mathbf{q}}_{i} = \begin{bmatrix}\dot{\mathbf{r}}_{i} \\ \boldsymbol{\omega}_{i}\end{bmatrix}$, with $\dot{\mathbf{r}}_{i}=[\dot{x}_{i},\,\dot{y}_{i},\,\dot{z}_{i}]^T$ and $\boldsymbol{\omega}_{i}=[\omega_{x,i},\,\omega_{y,i},\,\omega_{z,i}]^T$.  
    -  $\displaystyle \mathbf{f}_{i} = \begin{bmatrix}\mathbf{F}_{i} \\ \mathbf{M}_{i}\end{bmatrix}$, where $\mathbf{F}_{i}\in\mathbb{R}^3$ is the net external force and $\mathbf{M}_{i}\in\mathbb{R}^3$ is the net external torque about the COM.  
    - The term $\boldsymbol{\omega}_{i}\times\bigl(\mathbf{J}_{i}\,\boldsymbol{\omega}_{i}\bigr)$ is the **gyroscopic (Euler) term** that arises in the rotational dynamics.  


"""

# ╔═╡ 8f2cf452-ba97-49a6-a722-2a7c7a4ae56d
md"""

- **Compact Notation**  
  Stacking all bodies’ coordinates and forces gives  
  $$\mathbf{M}\,\ddot{\mathbf{q}} \;+\; \mathbf{G}(\boldsymbol{\omega}) \;=\; \mathbf{f}, 
    \quad
    \mathbf{M} = \mathrm{diag}\bigl(\mathbf{M}_{1},\,\mathbf{M}_{2},\,\dots,\,\mathbf{M}_{n}\bigr),
    \quad
    \mathbf{q} = \begin{bmatrix}\mathbf{q}_{1} \\[4pt] \mathbf{q}_{2} \\[2pt] \vdots \\[2pt] \mathbf{q}_{n}\end{bmatrix},
    \quad
    \mathbf{f} = \begin{bmatrix}\mathbf{f}_{1} \\[4pt] \mathbf{f}_{2} \\[2pt] \vdots \\[2pt] \mathbf{f}_{n}\end{bmatrix},$$  
  where  
  $$\mathbf{G}(\boldsymbol{\omega})
    \;=\;
    \begin{bmatrix}
      \displaystyle
      \begin{bmatrix}
        \mathbf{0}_{3} \\[4pt]
        \boldsymbol{\omega}_{1}\times(\mathbf{J}_{1}\,\boldsymbol{\omega}_{1})
      \end{bmatrix} \\[12pt]
      \displaystyle
      \begin{bmatrix}
        \mathbf{0}_{3} \\[4pt]
        \boldsymbol{\omega}_{2}\times(\mathbf{J}_{2}\,\boldsymbol{\omega}_{2})
      \end{bmatrix} \\[12pt]
      \vdots \\[12pt]
      \displaystyle
      \begin{bmatrix}
        \mathbf{0}_{3} \\[4pt]
        \boldsymbol{\omega}_{n}\times(\mathbf{J}_{n}\,\boldsymbol{\omega}_{n})
      \end{bmatrix}
    \end{bmatrix}.$$

"""

# ╔═╡ ddbdff62-1f24-454f-981b-c5892f991623
md"""
- **Remarks**  
  1. In an **unconstrained** system, bodies do not share joint constraints; they interact only through external forces (e.g., springs or contacts), which appear inside each $\mathbf{f}_{i}$.  
  2. If you have “springs + dampers + actuators” between bodies $i$ and $j$, then each contributes a force $\mathbf{F}_{ij}$ at the contact point; the moment in $\mathbf{f}_{i}$ is $\mathbf{r}_{ij}\times \mathbf{F}_{ij}$, and in $\mathbf{f}_{j}$ the equal-and-opposite $-\mathbf{F}_{ij}$, etc.  
  3. To write the **full ODE** for body $i$, you include the gyroscopic term explicitly:  
     $$\mathbf{M}_{i}\,\ddot{\mathbf{q}}_{i}
       \;+\;
       \begin{bmatrix}
         \mathbf{0}_{3} \\[4pt]
         \boldsymbol{\omega}_{i}\times\bigl(\mathbf{J}_{i}\,\boldsymbol{\omega}_{i}\bigr)
       \end{bmatrix}
       \;=\;
       \mathbf{f}_{i}
       .$$  
"""

# ╔═╡ 834471f2-792d-4f28-b129-ba94bd53d8b2
md"""

## Constrained Multibody Dynamics (Holonomic Constraints)
- In multibody systems, **constraints** (joint connections, fixed distances, etc.) impose algebraic relations $\boldsymbol{\phi}(q) = \mathbf{0}$ on the coordinates. For example, a pin joint fixes relative position, or a fixed orientation constraint between bodies.
- We enforce holonomic constraints using **Lagrange multipliers**. The equations of motion become a Differential-Algebraic System:

  $$\mathbf{M}(q)\,\ddot{q} + \mathbf{h}(q,\,\dot{q}) = \mathbf{Q}_{\mathrm{ext}} + \mathbf{J}^T(q)\,\boldsymbol{\lambda},$$

  along with the constraint equations $\boldsymbol{\phi}(q)=\mathbf{0}$. Here $\mathbf{J} = \frac{\partial \boldsymbol{\phi}}{\partial q}$ is the constraint Jacobian, and $\boldsymbol{\lambda}$ are Lagrange multipliers yielding constraint forces.
-  $\mathbf{h}(q,\,\dot{q})$ encompasses Coriolis, centrifugal, and gravity terms (if using a Lagrangian formulation). In a Newton–Euler formulation, the unconstrained equations for each body are augmented by constraint forces $\mathbf{J}^T \boldsymbol{\lambda}$ that enforce $\boldsymbol{\phi}(q)=\mathbf{0}$.
- The full system is a mix of differential and algebraic equations (DAE). One can solve this by solving for $\ddot{q}$ and $\boldsymbol{\lambda}$ simultaneously at each time step (index-3 DAE), or by reducing the coordinates (constraint reduction). The Lagrange multiplier approach keeps coordinates in Cartesian/quaternion form for simplicity.

"""

# ╔═╡ 25046354-8633-474d-bd3c-643a31ca3e7a
md"""
## System of Constrained Bodies

- **Equations of Motion with Lagrange Multipliers**  

  $$\mathbf{M}\,\ddot{\mathbf{q}} \;+\; \mathbf{C}_{q}^T\,\boldsymbol{\lambda} \;=\; \mathbf{f}, 
    \quad
    \mathbf{C}(\mathbf{q}) \;=\; \mathbf{0}.$$  

  -  $\boldsymbol{\lambda}$ is the vector of unknown Lagrange multipliers.  
  - Note that $\mathbf{C}_{q}^T\,\boldsymbol{\lambda}$ represents the generalized constraint forces.  
  - Directly solving this DAE is not straightforward; specialized techniques (e.g., penalty methods or stabilization) are required for numerical integration.

- **Augmented System via Constraint Differentiation**  
  Differentiate the holonomic constraint $\mathbf{C}(\mathbf{q}) = \mathbf{0}$ twice with respect to time:  

  $$\ddot{\mathbf{C}}
    \;=\;
    \mathbf{C}_{q}\,\ddot{\mathbf{q}} \;+\; \dot{\mathbf{C}}_{q}\,\dot{\mathbf{q}}
    \;=\;
    \mathbf{0}
    \;\;\Longrightarrow\;\;
    \mathbf{C}_{q}\,\ddot{\mathbf{q}} \;-\; \mathbf{g} \;=\; \mathbf{0},$$  

  where  

  $$\mathbf{g} \;=\; -\,\dot{\mathbf{C}}_{q}\,\dot{\mathbf{q}}
    \quad\text{(collects velocity‐dependent terms).}$$  

  Combining this with the equations of motion yields the block‐linear system:  

  $$\begin{bmatrix}
      \mathbf{M}       & \mathbf{C}_{q}^T \\[6pt]
      \mathbf{C}_{q}   & \mathbf{0}
    \end{bmatrix}
    \begin{bmatrix}
      \ddot{\mathbf{q}} \\[4pt]
      \boldsymbol{\lambda}
    \end{bmatrix}
    \;=\;
    \begin{bmatrix}
      \mathbf{f} \\[4pt]
      \mathbf{g}
    \end{bmatrix}.$$  

  - This system simultaneously solves for accelerations $\ddot{\mathbf{q}}$ and multipliers $\boldsymbol{\lambda}$.  
  - It is an index‐3 DAE, which can be numerically unstable if integrated directly. Use constraint stabilization (e.g.\ Baumgarte’s method) or coordinate reduction to obtain a stable integration scheme.  
"""

# ╔═╡ 788552d6-4235-4d85-b0db-fe3aa5ce46c5
md"""

## Constraint Stabilization (Baumgarte’s Method)
- **Drift Problem:** Numerical integration of DAEs can lead to slight violation of the constraints $\boldsymbol{\phi}(q)=\mathbf{0}$ (constraint drift). Small position errors can accumulate and cause constraint forces to become inaccurate.
- **Baumgarte Stabilization:** To combat drift, Baumgarte’s method adds damping and stiffness feedback terms to the constraint equations. Instead of enforcing $\ddot{\boldsymbol{\phi}}(q,\,\dot{q})=\mathbf{0}$ exactly, one imposes a **damped** form:

  $$\ddot{\boldsymbol{\phi}} \;+\; 2\alpha\,\dot{\boldsymbol{\phi}} \;+\; \beta^2\,\boldsymbol{\phi} \;=\; \mathbf{0},$$

  where $\alpha$ and $\beta$ are positive stabilization parameters. This is essentially a PD controller on the constraint: $\boldsymbol{\phi}$ is driven to zero like a critically damped oscillator (often one chooses $\beta=\alpha$ for critical damping).
- The modified acceleration-level constraints lead to a slightly altered set of equations for $\ddot{q}$ and $\boldsymbol{\lambda}$, which ensure that any small violations in $\boldsymbol{\phi}$ or $\dot{\boldsymbol{\phi}}$ are corrected over time. Baumgarte stabilization thus improves numerical stability at the cost of making the constraints “soft” (they behave like stiff springs if $\alpha,\beta$ are large).
- **Tuning:** The parameters $(\alpha,\,\beta)$ must be chosen with care—too large can make the system stiff or unstable, too small and the constraints drift slowly. When chosen well, the system’s constraint errors are quickly damped out without disturbing the physical motion significantly.

"""

# ╔═╡ 3c09fd43-8aa7-4d01-93f0-e14a67918172
md"""
# System of Constrained Bodies

- Finally, we can solve the following system for accelerations:
$$\begin{bmatrix}
    \mathbf{M} & \mathbf{C}_{q}^T \\[6pt]
    \mathbf{C}_{q} & \mathbf{0}
  \end{bmatrix}
  \begin{bmatrix}
    \ddot{\mathbf{q}} \\[4pt]
    \boldsymbol{\lambda}
  \end{bmatrix}
  \;=\;
  \begin{bmatrix}
    \mathbf{f} \\[4pt]
    g \;-\; 2\,\alpha\,\dot{C} \;-\; \beta^{2}\,C
  \end{bmatrix}$$

-  $\mathbf{M}$ – diagonal mass‐inertia matrix  
-  $C,\,\mathbf{C}_{q},\,\dot{C},\,g$ are as defined in kinematic analysis  
-  $\alpha$ and $\beta$ are constant (Baumgarte stabilization parameters)  
-  $\mathbf{f}$ is the generalized force vector  
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Base64 = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
PlutoUI = "~0.7.62"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "987fcbd7bed4a2524a7d452cee1f6c2ad944a8a7"

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

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

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

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

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
# ╟─0df19b60-3f9c-11f0-2860-5b7ac0586d16
# ╠═de59c8ec-b5e1-4adf-bacd-c10258eef778
# ╟─3ec42369-8374-422d-8366-506a2e843787
# ╟─942e3726-a441-4774-b54b-d4d1ea7ecafe
# ╟─f9a0bb53-035b-4dfc-9a42-db3068eccf17
# ╠═0ff028b1-8431-492e-8a77-403f67c9b212
# ╟─fe72f0b0-4a2c-4dc5-b5e6-8f6e8af815b0
# ╟─43ebdfe9-b0ce-4656-9592-d66d52570b8b
# ╟─5f8a633d-f1b4-406c-b732-30ff783cc040
# ╟─03c8a932-8f2a-4725-827d-4eee7e0f59b7
# ╟─96cf32e2-23f1-44b4-b7bd-fee1a0734522
# ╟─436cec7e-3e97-48da-b0e5-f296ba6eb14e
# ╟─c3a1b76e-cee8-4036-85c3-ced87079d14a
# ╟─8f2cf452-ba97-49a6-a722-2a7c7a4ae56d
# ╟─ddbdff62-1f24-454f-981b-c5892f991623
# ╟─834471f2-792d-4f28-b129-ba94bd53d8b2
# ╟─25046354-8633-474d-bd3c-643a31ca3e7a
# ╟─788552d6-4235-4d85-b0db-fe3aa5ce46c5
# ╟─3c09fd43-8aa7-4d01-93f0-e14a67918172
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
