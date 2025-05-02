### A Pluto.jl notebook ###
# v0.20.6

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
## Hands-on – Volume of a Cube

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
	println("Volume of the cube with side L = ",L," cm is V = ",V," cm².")
end

# ╔═╡ 1d0b6de8-da11-4784-8172-a45957753e9f


# ╔═╡ Cell order:
# ╟─a9775808-2742-11f0-0831-8ff037918d61
# ╟─ca4be4a5-8d18-43b6-b6f7-4dec7f35090f
# ╠═8f1e7ed7-1e21-42f1-b111-b33128279e8d
# ╟─675863d7-9deb-452e-845c-63edf399885a
# ╠═1d0b6de8-da11-4784-8172-a45957753e9f
