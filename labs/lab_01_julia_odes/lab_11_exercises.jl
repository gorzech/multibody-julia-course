### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ a9775808-2742-11f0-0831-8ff037918d61
md"""
# Programming Multibody Systems in Julia

### Lab 1: Julia Basics

Grzegorz Orzechowski
"""

# ╔═╡ ca4be4a5-8d18-43b6-b6f7-4dec7f35090f
md"""
## Volume of a Cube

- Write a program that computes the volume $V$ of a cube with sides of length $L=4$ cm and prints the result to the screen.
- Both $𝑉$ and $𝐿$ should be defined as separate variables in the program.
- Save the notebook and confirm that the correct result is printed.
"""

# ╔═╡ 8f1e7ed7-1e21-42f1-b111-b33128279e8d
# Put your solution here. The example solution is hidden below.



# ╔═╡ 675863d7-9deb-452e-845c-63edf399885a
begin
	L = 4
	V = L * L
	println("Volume of the cube with side L = $L cm is V = $V cm².")
end

# ╔═╡ 1d0b6de8-da11-4784-8172-a45957753e9f
md"""
## Volumes of Three Cubes

We are interested in the volume $𝑉$ of a cube with length $𝐿$: $𝑉=𝐿^3$, computed for three different values of $𝐿$.

1. Use the `range` function to compute three values of $𝐿$, equally spaced on the interval $[1, 3]$.
2. Caompte $V$ of `V = L.^3`. 
3. Make a plot of V versus L.
"""

# ╔═╡ ccb234af-b5f1-4ce7-94e1-f31172c312e5
md"""
## Compute π example

- Consider two computational schemes for the number 𝜋:
  - by Leibniz (1646–1716): $𝜋=8∑_{𝑘=0}^∞\frac{1}{(4𝑘+1)(4𝑘+3)}$
  - and by Euler (1707–1783): $𝜋=\sqrt{6∑_{𝑘=1}^∞\frac{1}{𝑘^2} }$
- Write a program that takes $𝑁$ as input from the user and plots the error development with both schemes as the number of iterations approaches $𝑁$. 
- Your program should also print out the final error achieved with both schemes, i.e., when the number of terms is $𝑁$. 
- Run the program with $𝑁=100$ and explain briefly what the graphs show.
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
# ╟─a9775808-2742-11f0-0831-8ff037918d61
# ╟─ca4be4a5-8d18-43b6-b6f7-4dec7f35090f
# ╠═8f1e7ed7-1e21-42f1-b111-b33128279e8d
# ╟─675863d7-9deb-452e-845c-63edf399885a
# ╟─1d0b6de8-da11-4784-8172-a45957753e9f
# ╠═ccb234af-b5f1-4ce7-94e1-f31172c312e5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
