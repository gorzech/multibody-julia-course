### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a871cd04-2cbb-11f0-0511-456e02603a5c
begin
	using Plots
	using LinearAlgebra
	using LaTeXStrings
	using PlutoUI
	using Test
	using Random
	using Base64
	function embed_image(path_to_image::AbstractString; height::Integer=240, type::String="png")
	    img_bytes = read(path_to_image)
	    # Build a data URL: "data:image/{type};base64,{...}"
	    data_url = "data:image/$(type);base64," * base64encode(img_bytes)
    	attrs = (:height => height,)
	    PlutoUI.Resource(data_url, MIME("image/$(type)"), attrs)
	end
end;

# ╔═╡ 40679803-a7b6-49a7-9a75-5d4d141420b2
using Printf

# ╔═╡ 7ba54843-ed1b-4233-8586-8dabe56b8c88
md"""
# Programming Multibody Systems in Julia

### Lab 2: Numerical Integration

Grzegorz Orzechowski
"""

# ╔═╡ 42485ba0-09cb-4a84-bf4f-35c24df8fbce
md"""
## Pluto Lab Workflow – Teams‑Integrated

**Flow**

1. **Kick‑off mini‑demo (~5 min)**  
  - Download today's Pluto notebook from GitHub repository.
    - `lab_12_numerical_integration.jl`
    - `lab_12_bonus_ode.jl`
  - Run first cell, outline goals & any new widgets  

2. **Self‑paced Blocks**  
  - Students tackle tasks from current block in Pluto  
  - Questions via **Teams Chat** / raise‑hand  
  - Each ticks status columns in the shared **Excel “Progress Board”**  

3. **Pit‑stop recap (~5 min)**  
  - When ~70 % of sheet cells are green, screen‑share model answer  
  - One student explains their approach; address FAQs  

4. **Repeat 2-3 for all blocks**

5. **Wrap‑up (5 min)**  
  - Students add short reflection at the end of this notebook
  - **Export HTML/PDF → upload to Teams › Assignments › Lab 12 - Integration**

---

### Teams Assets

- **GitHub** – today’s notebook
- **Excel Progress Board** – columns for task, colour‑coded  
- **Chat / Meeting** – live Q&A, screen‑shares at pit‑stops  
- **Assignments** – collects exported notebooks for evidence  

"""

# ╔═╡ 92f106b4-3183-4948-b306-6edb67dade4f
md"""
## Note for next week classes

Install [Visual Studi Code](https://code.visualstudio.com/) and related Julia extensions.

"""

# ╔═╡ 73c5828c-1287-45fb-a065-bdfcfafbef40
md"""
## Demo: Writing Tests

- Write simple test script to test the function add

"""

# ╔═╡ ad8c62ba-8ada-4e43-9c4d-9b1c7457a696
add(a, b) = a + b

# ╔═╡ b7366a38-2ba6-4c59-9fbb-b736baf5edb7
md"""
- Write two tests for 1+1=2 and 0.1+0.2=0.3 in single `@testset`
"""

# ╔═╡ 4bd134bf-013f-48da-85aa-eda3b963e705
@testset "Test add function on simple examples" begin
	@test add(1, 1) == 2
	@test add(0.1, 0.2) ≈ 0.3 atol = 1e-14
	@test isapprox(add(0.1, 0.2), 0.3, atol = 1e-14)
end

# ╔═╡ aaf17a1e-caad-41d1-b992-20600c0f87b4
@printf("%.17f\n", 0.1)

# ╔═╡ fde8b52a-ad59-4458-ab30-61012ae17050
@printf("%.17f\n", 0.2)

# ╔═╡ 9fcc30e7-fb7e-4975-9bc5-bdf596545759
@printf("%.17f\n", 0.1 + 0.2)

# ╔═╡ b4cc1bc0-34f7-4e1d-bf75-c89e4daee3b5
@printf("%.17f\n", 0.3)

# ╔═╡ fcab19f8-0703-42ea-a771-6aded71efef2
md"""
## Computing Integrals

- Good introduction to numerical methods
- You cannot solve most integrals analytically
- Well-known concept
- Exact mathematics vs numerical mathematics
  - For integration we do not have to bother too much with rounding errors
- Focus on how to write error-free code
"""

# ╔═╡ b77aeac2-8a9a-46e2-9251-62342018354a
md"""
## Integration

$$\int_{a}^{b} f(x) dx = F(b) - F(a), \text{ where } f(x) = \frac{dF}{dx}$$

- How to find anti-derivative $F(x)$?
- Above provide exact solution, but can be challenging or even impossible
- Numerical methods provide approximation only, but solution is straightforward
- Should we care, that solution might not be exact?
"""

# ╔═╡ d8a987dd-a332-4c55-b4b9-eca4f2b1491c
md"""
## Basic Ideas

$$\int_{a}^{b} f(x) dx$$

- Most methods split interval $[a,b]$ into smaller ones
- Usually, points are evenly distributed
- Thus, the integration points are 

$$x_i = a + ih, , i = 0, 1, \ldots, n, , h = \frac{b - a}{n}$$

- and 

$$\int_{a}^{b} f(x) , dx = \sum_{i=0}^{n-1} \int_{x_{i-1}}^{x_i} f(x) , dx$$
"""

# ╔═╡ 5e5885ab-93eb-48cf-94d1-4c744356ef60
md"""
- We end up with number of small intervals $[x_i, x_{i+1}]$
- Over small interval f can be approximated with simple function
    - Constant
    - Linear
    - Parabola
- Varying interval length (and approximation function) we can adjust accuracy
"""

# ╔═╡ d9ba5f01-dedd-4543-9594-1e135732de90
md"""
## Computational Example

- It is advantageous to use as test example integral with known solution
- Car velocity is given as $v(t) = 3t^2 e^{t^3}$
- What is the distance after one second?

$$\int_{0}^{1} v(t) , dt = V(1) - V(0) \approx 1.718$$

- where anti-derivative is given as $V(t) = e^{t^3} - 1$
"""

# ╔═╡ cfa76d39-7114-422f-8098-81e8c6379cf4
"""
plot_3t2exp_t3(; t_min=0.0, t_max=1.0, n=200, nmarkers=6) → Plot

Plot y=3 t^2 exp(t^3) over [t_min,t_max] with
• solid white line
• yellow fill down to y=0
• red‐edged, white‐filled markers at nmarkers evenly‐spaced t’s
"""
function plot_3t2exp_t3(; t_min=0.0, t_max=1.0, n::Int=200, nmarkers::Int=6)
	# high‐res curve
	t = range(t_min, t_max, length=n)
	y = 3 .* t.^2 .* exp.(t.^3)

	# base plot: white curve on black
	plt = plot(
	  t, y;
	  linewidth            = 2,
	  xlabel               = L"t",
	  ylabel               = L"y(t) = 3t^2e^{t^3}",
	  legend               = false,
	)

	# yellow fill down to y=0
	plot!(plt, t, y; fillrange=0, fillcolor=:yellow, alpha=1.0)

	return plt
end;

# ╔═╡ f47504cf-2a32-4484-ac67-29d5d07797e9
plot_3t2exp_t3()

# ╔═╡ 8c101a8a-15f7-4ce1-b2ec-0118426e4494
md"""
## Composite Trapezoidal Rule

 - In this example we have intervals of length $h=0.25$ s. We can approximate integral as 

$$\int_0^1 v(t) dt \approx h \frac{(v(0)+v(0.25))}{2} + h \frac{(v(0.25)+v(0.5))}{2}$$
$$+ h \frac{(v(0.5)+v(0.75))}{2} + h \frac{(v(0.75)+v(1))}{2} = 1.9227$$ 

Compared with 1.718, it differs by about 12%. 

Let use more trapezoids!

"""

# ╔═╡ 64bbb77a-8f89-4ca8-b284-0523c7f201e9
"""
    plot_trapz_3t2exp_t3(n_traps; t_min=0.0, t_max=1.0, n_fine=200)

Plot y=3 t² e^{t³} over [t_min,t_max], shade the composite trapezoids
in yellow, and draw black edges around each trapezoid.
"""
function plot_trapz_3t2exp_t3(n_traps::Int; t_min::Float64=0.0,
                                       t_max::Float64=1.0,
                                       n_fine::Int=200)

    # 1) smooth curve
    t_fine = range(t_min, t_max; length=n_fine)
    y_fine = 3 .* t_fine.^2 .* exp.(t_fine.^3)

    # 2) nodes for trapezoids
    xs = range(t_min, t_max; length=n_traps+1)
    ys = 3 .* xs.^2 .* exp.(xs.^3)

    # 3) Base plot: white curve on black
    plt = plot(
        t_fine, y_fine;
        linewidth  = 2,
        xlabel     = L"t",
        ylabel     = L"y(t)=3t^2e^{t^3}",
        legend     = false,
        framestyle = :box,
    )

    # 4) Draw and outline each trapezoid
    for i in 1:n_traps
        x0, x1 = xs[i],   xs[i+1]
        y0, y1 = ys[i],   ys[i+1]
        # polygon coords for the trapezoid
        xpoly = [x0, x0, x1, x1]
        ypoly = [0.0, y0,  y1,  0.0]
        # filled yellow with black border
        plot!(
            plt,
            xpoly, ypoly;
            seriestype = :shape,
            fillcolor  = :yellow,
            alpha      = 0.8,
            linecolor  = :black,
            linewidth  = 1.5,
            label      = ""
        )
    end

    # 5) Overlay the true curve and markers
    plot!(plt, t_fine, y_fine; color=:blue, linewidth=2)
    scatter!(
        plt,
        xs, ys;
        marker            = :circle,
        markersize        = 6,
        markerstrokecolor = :red,
        markercolor       = :white,
        label             = ""
    )

    return plt
end; 

# ╔═╡ 7cff897b-b18d-406c-98e2-427d84638e44
@bind n_trapezoids_bind Slider(2:15, default=4, show_value=true)

# ╔═╡ 92760dd1-60d0-41a5-836f-967c2c188b26
plot_trapz_3t2exp_t3(n_trapezoids_bind)

# ╔═╡ 377084d5-ea28-4d6a-bb84-d948d7481335
md"""
## General Algorithm

$$\int_a^b f(x) dx = \int_{x_0}^{x_1} f(x) dx + \int_{x_1}^{x_2} f(x) dx + \ldots + \int_{x_{n-1}}^{x_n} f(x) dx$$ 
$$\approx h \frac{f(x_0) + f(x_1)}{2} + h \frac{f(x_1) + f(x_2)}{2} + \ldots + h \frac{f(x_{n-1}) + f(x_n)}{2}$$ 
$$= h/2 [f(x_0) + 2f(x_1) + 2f(x_2) + \ldots + 2f(x_{n-1}) + f(x_n)]$$ 
$$= h \left[ \frac{f(x_0)}{2} + \sum_{i=1}^{n-1} f(x_i) + \frac{f(x_n)}{2} \right]$$

"""

# ╔═╡ 02f9ece7-d540-4c4e-9510-b3063533ec3e
md"""
## D.1 Trapezoidal Integration

$$\int_a^b f(x) , dx \approx h \left[ \frac{f(x_0)}{2} + \sum_{i=1}^{n-1} f(x_i) + \frac{f(x_n)}{2} \right]$$ 

Write function `trapezoidal(f, a, b, n)`. Return the approximation of the integral. Test your implementation on the function $v(t) = 3t^2 e^{t^3}$ with anti-derivative $V(t) = e^{t^3} - 1$
"""

# ╔═╡ 30ed1925-12f8-4c0e-aa94-38b18142a875
function trapezoidal(f, a, b, n)
	h = (b - a) / n
	integral = 0.5(f(a) + f(b))
	for i = 1:n-1
		integral += f(a + i * h)
	end
	# integral += sum(f.(a .+ (1:n-1) .* h))
	return integral * h
end

# ╔═╡ bbecb121-f67f-4b79-a644-2549c4ebd427
v(t) = 3t^2 * exp(t^3)

# ╔═╡ 8ba1a37c-ecd7-4866-bd7a-c4eea406de1c
V(t) = exp(t^3) - 1

# ╔═╡ bcb543e7-46b6-4bec-b8c2-d98af9ea03ed
trapezoidal(v, 0, 2, 10_000)

# ╔═╡ 7d1afb22-7df0-4188-a8ee-02225d54fcf2
V(2) - V(0)

# ╔═╡ 01dc37bd-2111-4d5a-b5f9-93aafa073ebf
@testset "Compare trapezoidal with exact calculations" begin
	v(t) = 3t^2 * exp(t^3)
	V(t) = exp(t^3) - 1
	a = -1
	b = 3
	result = trapezoidal(v, a, b, 10_000)
	expected = V(b) - V(a)
	@test result ≈ expected rtol = 1e-4
end

# ╔═╡ fb137bca-4905-427a-9330-7e983a89ef9c
md"""
## Composite Midpoint Method

Approximate area with rectangles. This method is even more accurate than the trapezoid. Rectangle height is measured in the middle of the interval (therefore its name). 

$$\int_0^1 v(t) dt \approx h v\left(\frac{0+0.25}{2}\right) + h v\left(\frac{0.25+0.5}{2}\right) $$
$$+ h v\left(\frac{0.5+0.75}{2}\right) + h v\left(\frac{0.75+1}{2}\right) = 1.6190$$ 
Compared with 1.718, it differs by about 6%
"""

# ╔═╡ 74276d46-bcd3-4cbf-befb-06675fe47d0e
"""
composite_midpoint_plot(n; t_min=0.0, t_max=1.0, n_fine=200)

Plot y(t)=3 t² e^{t³} on [t_min,t_max] (white curve) and overlay the
composite midpoint‐rule rectangles (yellow fill, black border),
with markers at each midpoint.
"""
function composite_midpoint_plot(n::Int;
								  t_min::Float64=0.0,
								  t_max::Float64=1.0,
								  n_fine::Int=200)

	# Base function
	f(t) = 3*t^2 * exp(t^3)

	# 1) smooth curve for plotting
	t_fine = range(t_min, t_max; length=n_fine)
	y_fine = f.(t_fine)

	# 2) partition nodes and midpoints
	xs = range(t_min, t_max; length=n+1)
	ms = (xs[1:end-1] .+ xs[2:end]) ./ 2    # midpoints
	hs = f.(ms)                             # heights at midpoints

	# 3) plot the smooth curve
	plt = plot(
		t_fine, y_fine;
		linewidth  = 2,
		xlabel     = L"t",
		ylabel     = L"y(t)=3t^2e^{t^3}",
		legend     = false
	)

	# 4) draw each midpoint‐rule rectangle
	for i in 1:n
		x0, x1 = xs[i], xs[i+1]
		h       = hs[i]
		xpoly   = [x0, x0, x1, x1]
		ypoly   = [0.0, h,   h,   0.0]
		plot!(
			plt,
			xpoly, ypoly;
			seriestype = :shape,
			fillcolor  = :yellow,
			linecolor  = :black,
			alpha      = 0.8,
			label      = ""
		)
	end

	# 5) re‐overlay the white curve on top of rectangles
	plot!(plt, t_fine, y_fine; color=:blue, linewidth=2)

	# 6) mark the midpoints
	scatter!(
		plt,
		ms, hs;
		marker            = :circle,
		markersize        = 6,
		markercolor       = :white,
		markerstrokecolor = :red,
		label             = ""
	)

	return plt
end;

# ╔═╡ 73ebd107-39d5-4963-a866-20e341d5e481
@bind n_rectangles_bind Slider(2:15, default=4, show_value=true)

# ╔═╡ 535c6122-3abc-4092-9886-f5d7c1a517a1
composite_midpoint_plot(n_rectangles_bind)

# ╔═╡ 42ce2dda-f4a9-4d07-be60-4674d732ffa3
md"""
## General Algorithm for the Midpoint Rule

$$\int_a^b f(x) dx = \int_{x_0}^{x_1} f(x) dx + \int_{x_1}^{x_2} f(x) dx + \ldots + \int_{x_{n-1}}^{x_n} f(x) dx$$ 

$$\approx h f\left(\frac{x_0 + x_1}{2}\right) + h f\left(\frac{x_1 + x_2}{2}\right) + \ldots + h f\left(\frac{x_{n-1} + x_n}{2}\right)$$ 

$$= h \left[ f\left(\frac{x_0 + x_1}{2}\right) + f\left(\frac{x_1 + x_2}{2}\right) + \ldots + f\left(\frac{x_{n-1} + x_n}{2}\right) \right]$$ 

$$= h \sum_{i=0}^{n-1} f(x_i), \quad x_i = \left(a + \frac{h}{2}\right) + ih$$
"""

# ╔═╡ e9175143-62ad-4071-8811-2518cbd91a0b
md"""
## D.2 Midpoint Integration

$$\int_a^b f(x) , dx \approx h \sum_{i=0}^{n-1} f(x_i), \quad x_i = \left(a + \frac{h}{2}\right) + ih$$ Write function `midpoint(f, a, b, n)`. Return the approximation of the integral. Check your implementation on the function $v(t) = 3t^2 e^{t^3}$ with anti-derivative $V(t) = e^{t^3} - 1$
"""

# ╔═╡ 307d5cc7-56cf-4e67-97ac-9af4e674dfac
function midpoint(f, a, b, n)
	h = (b - a) / n
	x_i = [a + 0.5h + i*h for i in 0:n-1]
	h * sum(f.(x_i))
end

# ╔═╡ 8e2ffd33-aad6-4693-b967-f681fcbdfea9
@testset "Compare midpoint with exact calculations" begin
	v(t) = 3t^2 * exp(t^3)
	V(t) = exp(t^3) - 1
	a = -1
	b = 3
	result = midpoint(v, a, b, 10_000)
	expected = V(b) - V(a)
	@test result ≈ expected rtol = 1e-4
end

# ╔═╡ c5a3a9b4-565e-4633-a096-6b72744e255f
md"""
## D.3 Compare Methods

Write a program to compare trapezoidal and midpoint methods. Use as example $e^{-y^2}$ on the interval from 0 to 2. Print results for $n=2, 2^2,\ldots,2^{20}$. Remember about nice formatting.
"""

# ╔═╡ fad9614a-30b3-447b-adb5-e0190ba12add
md"""
## Remark About Integration

There are many numerical integration methods:

- Simpson’s rule
- Gauss quadrature

$$\int_a^b f(x) , dx \approx \sum_{i=0}^{n-1} w_i f(x_i)$$

They differ in the way weights $w_i$ and points $x_i$ are chosen. Higher accuracy can be achieved by optimizing the location of $x_i$ 2.
"""

# ╔═╡ 27ea89c8-a090-4f57-aa92-da59d175eb3c
md"""
## Testing

We have done two types of tests:

- Comparison with exact solution (error decreases as $n$ increases)
- Check if the integral value stabilizes as $n$ grows

Unfortunately, those are very weak tests for software checking. There are still types of errors that those tests might not catch. Usually, the longer we use the procedure, the more we can trust it 3.
"""

# ╔═╡ 7dba23c3-a80c-40b4-94ff-a90bf9ebc7ca
md"""
## Proper Test Procedure

There are three serious ways to test the implementation of numerical methods via unit tests:

- Comparing with hand-computed results usually using problem with few operations – small $n$
- Solving a problem without numerical errors. E.g., trapezoidal rule must be exact for linear functions. Produced error should be zero (to machine precision)
- Demonstrate correct convergence rates. A strong test when we can compute exact errors, is to see how fast the error goes to zero as $n$ grows. In the trapezoidal and midpoint rules it is known that the error depends on $n$ as $n^{-2}$ as $n \to \infty$ 4.
"""

# ╔═╡ d12adf01-c32b-4f91-8b27-8072a6dc2eb1
md"""
## E.1 Hand-computed Results Test

We have already done this (and test the result looks the same). We can even simplify this for two trapezoids for $\int_0^1 v(t) , dt$, $v(t)=3t^2 e^{t^3}$

$$\frac{1}{2} h[v(0)+v(0.5)]+\frac{1}{2} h[v(0.5)+v(1)]=2.463642041244344$$

when $h=0.5$. Use above equation to write first test in that file. Verify if your test works.
"""

# ╔═╡ ab083dd7-d805-4174-8460-1fd296d8e41f
@testset "Hand-computed results for 2 trapezoids" begin
	expected = 2.463642041244344
	v(t) = 3t^2*ℯ^(t^3)
	result = trapezoidal(v, 0, 1, 2)
	@test expected ≈ result atol = 1e-14
end

# ╔═╡ d7fa4c50-5781-4ad9-9aee-dd5d0461ed77
md"""
## E.2 Problem Without Numerical Errors

Usually, numerical results contain unknown approximation errors. Approximation error vanishes – may be present in simple mathematical problems. E.g., the trapezoidal rule is exact for integration of linear functions. Specific test can be $\int_{1.2}^{4.4} (6x-4) , dx$. We try to avoid special numbers like 0 and 1.
"""

# ╔═╡ bdf30f39-2665-49fd-b1a0-67c5679ccc76
md"""
Use $\int_{1.2}^{4.4} (6x-4) , dx$ to write second test. Perform computation using three different values of $n$. Function anti-derivative is $F(x)=3x^2-4x$. Verify if your test works.
"""

# ╔═╡ 81a01d9e-e71b-4482-bf24-7d7fb6b40904
@testset "Trapezoidal - problem with no numerical errors" begin
	a = 1.2
	b = 4.4
	f(x) = 6x - 4
	F(x) = 3x^2 - 4x
	expected = F(b) - F(a)

	n = 2
	result = trapezoidal(f, a, b, n)
	@test expected ≈ result atol = 1e-14

	n = 7
	result = trapezoidal(f, a, b, n)
	@test expected ≈ result atol = 1e-14

	n = 385
	result = trapezoidal(f, a, b, n)
	@test expected ≈ result atol = 1e-14
end

# ╔═╡ d8c2e883-4ed1-4fe9-a148-0b969ed24a51
md"""
## Correct Convergence Rates

- Approximation errors are usually unknown, but we often may assume a certain asymptotic behavior of the error. 
- E.g., for trapezoidal rule when we double $n$ error is reduced by factor of about 4 – error convergence to zero as $n^{-2}$ – the convergence rate is 2. 
- Error of numerical integration methods usually converge to zero as $n^{-p}$. 
- Assume the error depends on $n$ as $E=Cn^r$ where $C$ is unknown constant and $r$ is the convergence rate.
"""

# ╔═╡ 1675513d-e387-4a5b-b08f-b4df538047a6
md"""

Consider $q$ experiments (numerical) with $n:n_1,n_2,\ldots,n_q$ and corresponding errors $E_1,E_2,\ldots,E_q$. 

For two consecutive experiments, $i$ and $i-1$, we have the error model:

$$E_i=Cn_i^r$$ 

$$E_{i-1}=Cn_{i-1}^r$$

$$\Rightarrow r_i=\frac{\ln(E_i/E_{i-1})}{\ln(n_i/n_{i-1})}$$

 $i$ index for $r$ was added as value of $r$ varies with $i$. Hopefully, $r_i$ approaches correct convergence rate as $n$ increases and $i \to q$.
"""

# ╔═╡ 20979253-7451-495b-8629-0b2443d79375
md"""
## E.4 Convergence Rate Test

The approximation errors in the trapezoidal rule are proportional to $n^{-2}$ for $i=1,2,\ldots,q$

$$n_i=2^i$$

Compute integral with $n_i$ intervals. Compute the error $E_i$. Estimate $r_i=\frac{\ln(E_i/E_{i-1})}{\ln(n_i/n_{i-1})}$ if $i>1$.
"""

# ╔═╡ f67f2677-c22a-476c-bd6c-adcdd198df10
md"""
##  F.1 Adaptive integration

Suppose we want to use the trapezoidal or midpoint method to compute an integral $\int_a^b f(x)dx$ with an error less than a prescribed tolerance $\epsilon$. What is the appropriate size of $n$? To answer this question, we may enter an iterative procedure where we compare the results produced by $n$ and $2n$ intervals, and if the difference is smaller than $\epsilon$, the value corresponding to $2n$ is returned. Otherwise, we double $n$ and repeat the procedure.

_Hint_: It may be a good idea to organize your code so that the function `adaptive_integration` can be used easily in future programs you write.

a) Write a function `adaptive_integration(f, a, b, tol, method="midpoint")` that implements the idea above (`tol` corresponds to the tolerance $\epsilon$, and `method` can be midpoint or trapezoidal).

b) Test the method on $\int_0^2 x^2 dx$ and $\int_0^2 \sqrt{x} dx$ for $\epsilon = 10^{-1}, 10^{-10}$ and write out the exact error.

c) Make a plot of $n$ versus $\epsilon \in [10^{-1}, 10^{-10}]$ for $\int_0^2 \sqrt{x} dx$. Use logarithmic scale for $\epsilon$.

_Remarks_: The type of method explored in this exercise is called adaptive because it tries to adapt the value of $n$ to meet a given error criterion. The true error can very seldom be computed (since we do not know the exact answer to the computational problem), so one has to find other indicators of error such as when here where changes in integral value as number of intervals doubled taken reflect error.

"""

# ╔═╡ d1334a48-646d-46f1-9055-fbad2821c699
function adaptive_integration(f, a, b, tol, method="midpoint")
	g = method == "midpoint" ? midpoint : trapezoidal

	n = 1
	value_n = g(f, a, b, n)
	n *= 2
	value_2n = g(f, a, b, n)
	while abs(value_2n - value_n) > tol
		value_n = value_2n
		n *= 2
		value_2n = g(f, a, b, n)
	end
	return value_2n, n
end

# ╔═╡ 3c873718-5942-4a38-94b5-e810e01ec965
v1(x) = x^2

# ╔═╡ 54e90d28-a8b5-4751-9db3-b02fa9c1b622
V1(x) = (1/3) * x^3

# ╔═╡ 3cbef440-daa4-4b30-a3b2-47ea7b438ba4
value1, n1 = adaptive_integration(v1, 0, 2, 1e-1, "midpoint")

# ╔═╡ 4f423781-62f0-48cf-96eb-c3fa2b9e7b8b
actual_value1 = V1(2) - V1(0)

# ╔═╡ e4907b58-d295-4434-91e7-3499784657e2
value1a, n1a = adaptive_integration(v1, 0, 2, 1e-10, "midpoint")

# ╔═╡ 435075b9-d18d-4879-ade7-4b4d0cfbc9f7
error1a = abs(actual_value1 - value1a)

# ╔═╡ 9578a1d6-4d4b-410b-aaeb-12ad621eff9e
tol_values = 10.0.^(-1:-0.5:-10)

# ╔═╡ 953fe954-1164-4d28-99dc-b2dde07be796
n_values = zeros(Int, length(tol_values))

# ╔═╡ dfc23663-ab5a-4dce-b2ce-5d775b3f73a7
v2(x) = √x

# ╔═╡ ea7e8b8a-1cbb-4d38-859d-47d404ad74c6
for (i, tv) in enumerate(tol_values)
	_, ni = adaptive_integration(v2, 0, 2, tv, "midpoint")
	n_values[i] = ni
end

# ╔═╡ 6f02ffe6-fc9f-4bff-9e70-8e2834e91cff
plot(tol_values, n_values, xscale=:log10, yscale=:log10)

# ╔═╡ 839e5694-7cf3-45b9-ac1a-860c1b2ab555
md"""
##  F.2 Revisit the fit of sines to a function

This is a continuation of Exercise 2.18. The task is to approximate a given function $f(t)$ on $[-\pi, \pi]$ by a sum of sines,

$$S_N(t) = \sum_{n=1}^{N} b_n \sin(nt).$$

We are now interested in computing the unknown coefficients $b_n$ such that $S_N(t)$ is in some sense the best approximation to $f(t)$. One common way of doing this is to first set up a general expression for the approximation error, measured by "summing up" the squared deviation of $S_N$ from $f$:

$$E = \int_{-\pi}^{\pi} (S_N(t)-f(t))^2 dt.$$

We may view $E$ as a function of $b_1, \ldots, b_N$. Minimizing $E$ with respect to $b_1, \ldots, b_N$ will give us a best approximation, in the sense that we adjust $b_1, \ldots, b_N$ such that $S_N$ deviates as little as possible from $f$.

Minimization of a function of $N$ variables, $E(b_1, \ldots, b_N)$ is mathematically performed by requiring all the partial derivatives to be zero:

$$\frac{\partial E}{\partial b_1} = 0, \frac{\partial E}{\partial b_2} = 0, \ldots, \frac{\partial E}{\partial b_N} = 0.$$

a) Compute the partial derivative $\frac{\partial E}{\partial b_1}$ and generalize to the arbitrary case $\frac{\partial E}{\partial b_N}$, $1 \leq n \leq N$.

b) Show that 

$$b_n = \frac{1}{\pi} \int_{-\pi}^{\pi} f(t) \sin(nt) dt.$$

c) Write a function `integrate_coeffs(f, N, M)` that computes $b_1, \ldots, b_N$ by numerical integration using $M$ intervals in the trapezoidal rule.

d) A remarkable property of the trapezoidal rule is that it is exact for integrals $\int_{-\pi}^{\pi} \sin(nt) dt$ (when subintervals are of equal size). Use this property to create a cell `test_integrate_coeff` to verify the implementation of `integrate_coeffs`.

e) Implement the choice $f(t)= \frac{1}{\pi} t$ as a Matlab function $f(t)$ and call `integrate_coeffs(f, 3, 100)` to see what optimal choice of $b_1, b_2, b_3$.

f) Make a function `plot_approx(f, N, M, filename)` where you plot $f(t)$ together with best approximation $S_N$ computed above using $M$ intervals for numerical integration.

g) Run `plot_approx(f, N, M, filename)` for $f(t)=\frac{1}{\pi} t$ for $N=3, 6, 12, 24$. Observe how approximation improves.

h) Run `plot_approx` for $f(t)= e^{-(t-\pi)}$ and $N=100$. Observe a fundamental problem: regardless of $N$, $S_N(-\pi)=0$, not $e^{2\pi}\approx 535$.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Base64 = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
LaTeXStrings = "~1.4.0"
Plots = "~1.40.13"
PlutoUI = "~0.7.23"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "c99d3ffc4fd019629e3b3aef6096c47bc1899e28"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "64e15186f0aa277e174aa81798f7eb8598e0157e"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "7ffa4049937aeba2e5e1242274dc052b0362157a"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.14"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "98fc192b4e4b938775ecd276ce88f539bcec358e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.14+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "f93655dc73d7a0b4a368e3c0bce296ae035ad76e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.16"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

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

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

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

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd10d2cc78d34c0e2a3a36420ab607b611debfbb"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.7"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

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

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

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

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9216a80ff3682833ac4b733caa8c00390620ba5d"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.0+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3b31172c032a1def20c98dae3f2cdc9d10e3b561"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "809ba625a00c605f8d00cd2a9ae19ce34fc24d68"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.13"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5152abbdab6488d5eec6a01029ca6697dff4ec8f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.23"

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

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "b81c5035922cc89c2d9523afc6c54be512411466"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.5"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "85c7811eddec9e7f22615371c3cc81a504c508ee"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+2"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5db3e9d307d32baba7067b13fc7b5aa6edd4a19a"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.36.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0ba42241cb6809f1a278d0bcb976e0483c3f1f2d"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+1"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "068dfe202b0a05b8332f1e8e6b4080684b9c7700"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.47+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "c950ae0a3577aec97bfccf3381f66666bc416729"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.8.1+0"
"""

# ╔═╡ Cell order:
# ╟─7ba54843-ed1b-4233-8586-8dabe56b8c88
# ╠═a871cd04-2cbb-11f0-0511-456e02603a5c
# ╟─42485ba0-09cb-4a84-bf4f-35c24df8fbce
# ╟─92f106b4-3183-4948-b306-6edb67dade4f
# ╟─73c5828c-1287-45fb-a065-bdfcfafbef40
# ╠═ad8c62ba-8ada-4e43-9c4d-9b1c7457a696
# ╟─b7366a38-2ba6-4c59-9fbb-b736baf5edb7
# ╠═4bd134bf-013f-48da-85aa-eda3b963e705
# ╠═40679803-a7b6-49a7-9a75-5d4d141420b2
# ╠═aaf17a1e-caad-41d1-b992-20600c0f87b4
# ╠═fde8b52a-ad59-4458-ab30-61012ae17050
# ╠═9fcc30e7-fb7e-4975-9bc5-bdf596545759
# ╠═b4cc1bc0-34f7-4e1d-bf75-c89e4daee3b5
# ╟─fcab19f8-0703-42ea-a771-6aded71efef2
# ╟─b77aeac2-8a9a-46e2-9251-62342018354a
# ╟─d8a987dd-a332-4c55-b4b9-eca4f2b1491c
# ╟─5e5885ab-93eb-48cf-94d1-4c744356ef60
# ╟─d9ba5f01-dedd-4543-9594-1e135732de90
# ╟─cfa76d39-7114-422f-8098-81e8c6379cf4
# ╟─f47504cf-2a32-4484-ac67-29d5d07797e9
# ╟─8c101a8a-15f7-4ce1-b2ec-0118426e4494
# ╟─64bbb77a-8f89-4ca8-b284-0523c7f201e9
# ╠═7cff897b-b18d-406c-98e2-427d84638e44
# ╟─92760dd1-60d0-41a5-836f-967c2c188b26
# ╟─377084d5-ea28-4d6a-bb84-d948d7481335
# ╟─02f9ece7-d540-4c4e-9510-b3063533ec3e
# ╠═30ed1925-12f8-4c0e-aa94-38b18142a875
# ╠═bbecb121-f67f-4b79-a644-2549c4ebd427
# ╠═8ba1a37c-ecd7-4866-bd7a-c4eea406de1c
# ╠═bcb543e7-46b6-4bec-b8c2-d98af9ea03ed
# ╠═7d1afb22-7df0-4188-a8ee-02225d54fcf2
# ╠═01dc37bd-2111-4d5a-b5f9-93aafa073ebf
# ╟─fb137bca-4905-427a-9330-7e983a89ef9c
# ╟─74276d46-bcd3-4cbf-befb-06675fe47d0e
# ╠═73ebd107-39d5-4963-a866-20e341d5e481
# ╠═535c6122-3abc-4092-9886-f5d7c1a517a1
# ╟─42ce2dda-f4a9-4d07-be60-4674d732ffa3
# ╟─e9175143-62ad-4071-8811-2518cbd91a0b
# ╠═307d5cc7-56cf-4e67-97ac-9af4e674dfac
# ╠═8e2ffd33-aad6-4693-b967-f681fcbdfea9
# ╟─c5a3a9b4-565e-4633-a096-6b72744e255f
# ╟─fad9614a-30b3-447b-adb5-e0190ba12add
# ╟─27ea89c8-a090-4f57-aa92-da59d175eb3c
# ╟─7dba23c3-a80c-40b4-94ff-a90bf9ebc7ca
# ╟─d12adf01-c32b-4f91-8b27-8072a6dc2eb1
# ╠═ab083dd7-d805-4174-8460-1fd296d8e41f
# ╟─d7fa4c50-5781-4ad9-9aee-dd5d0461ed77
# ╟─bdf30f39-2665-49fd-b1a0-67c5679ccc76
# ╠═81a01d9e-e71b-4482-bf24-7d7fb6b40904
# ╟─d8c2e883-4ed1-4fe9-a148-0b969ed24a51
# ╟─1675513d-e387-4a5b-b08f-b4df538047a6
# ╟─20979253-7451-495b-8629-0b2443d79375
# ╟─f67f2677-c22a-476c-bd6c-adcdd198df10
# ╠═d1334a48-646d-46f1-9055-fbad2821c699
# ╠═3c873718-5942-4a38-94b5-e810e01ec965
# ╠═54e90d28-a8b5-4751-9db3-b02fa9c1b622
# ╠═3cbef440-daa4-4b30-a3b2-47ea7b438ba4
# ╠═4f423781-62f0-48cf-96eb-c3fa2b9e7b8b
# ╠═e4907b58-d295-4434-91e7-3499784657e2
# ╠═435075b9-d18d-4879-ade7-4b4d0cfbc9f7
# ╠═9578a1d6-4d4b-410b-aaeb-12ad621eff9e
# ╠═953fe954-1164-4d28-99dc-b2dde07be796
# ╠═dfc23663-ab5a-4dce-b2ce-5d775b3f73a7
# ╠═ea7e8b8a-1cbb-4d38-859d-47d404ad74c6
# ╠═6f02ffe6-fc9f-4bff-9e70-8e2834e91cff
# ╟─839e5694-7cf3-45b9-ac1a-860c1b2ab555
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
