### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ 441cc450-405b-11f0-32d4-b533235d9307
begin
	using LinearAlgebra
	using ForwardDiff
	using SparseArrays
	using Test
	using PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ 61976551-986d-4297-8f15-dad98120346a
md"""
# Programming Multibody Systems in Julia

### Lab 3.1: Constraints and their derivatives

Following the nomenclature and equations from: E. J. Haug, _Computer Aided Kinematics and Dynamics of Mechanical Systems: Basic Methods._ Boston: Allyn and Bacon, 1989.
 

Grzegorz Orzechowski
"""

# ╔═╡ 9425a84f-29af-4b74-815e-6f1c0876bd35
md"""
## Helper functions
"""

# ╔═╡ f180ffcc-79e3-4554-a4d8-d7172dd6c133
# ========== small helpers =====================================================
skew(a) = [ 0.0 -a[3]  a[2];
          a[3] 0.0  -a[1];
          -a[2] a[1]  0.0 ]

# ╔═╡ 714a2372-2df1-422a-912b-a77223c022da
rotmat(p) = begin
    e0,e1,e2,e3 = p
    [
        2(e0^2 + e1^2) - 1    2(e1 * e2 - e0 * e3)    2(e1 * e3 + e0 * e2);
        2(e1 * e2 + e0 * e3)  2(e0^2 + e2^2) - 1      2(e2 * e3 - e0 * e1);
        2(e1 * e3 - e0 * e2)  2(e2 * e3 + e0 * e1)    2(e0^2 + e3^2) - 1
    ]
end

# ╔═╡ ab0159aa-f6e5-4ca5-b2d4-f9e666e3902a
# body-frame conversion E(p)  (4×3)  and its t-derivative
function Emat(p)
    e0,e1,e2,e3 = p
    [
        -e1  -e2  -e3;
         e0   e3  -e2;
        -e3   e0   e1;
         e2  -e1   e0
    ]'
end

# ╔═╡ ca7454e1-843d-488e-bda0-e1ca4e95236c
# body-frame conversion E(p)  (4×3)  and its t-derivative
function Gmat(p)
    e0,e1,e2,e3 = p
    [
        -e1  e0  e3 -e2;
		-e2 -e3  e0  e1;
		-e3  e2 -e1  e0
    ]
end

# ╔═╡ b3e08dc3-a0aa-440e-a530-5ef5568694e9
Gdot(p,ṗ) = Gmat(ṗ)   # linear dependence on p 

# ╔═╡ 7b2a97cb-223a-4603-bd89-d813ff8a831b
md"""
## Body state
"""

# ╔═╡ 831493a6-a56f-4f8c-9862-ce331a615e32
# ========== light container ===================================================
struct BodyState{T<:Number}
    p  :: Vector{T}   # unit quaternion (scalar first)
    r  :: Vector{T}   # world position
    ω  :: Vector{T}   # body‐frame angular velocity
    v  :: Vector{T}   # world linear velocity
end

# ╔═╡ daae7f68-1249-409d-915e-53ba6aadbc2c
#---------------------------------------------------------------- pack/unpack ---
pack(b::BodyState) = vcat(b.p, b.r)              # 7 numbers / body

# ╔═╡ 173d35d7-ba79-4faa-bdf8-86a964f822d1
pack(bA, bB)   = vcat(pack(bA), pack(bB))    # 14

# ╔═╡ 7a1ba4ee-801b-4902-bf48-ce20521fead1
function unpack(p::Vector{T}, templA, templB) where T
    qA = p[1:4];   rA = p[5:7]
    qB = p[8:11];  rB = p[12:14]
    BodyState{T}(qA, rA, templA.ω, templA.v),
    BodyState{T}(qB, rB, templB.ω, templB.v)
end

# ╔═╡ a9fb78b9-061e-417a-a89b-e41c818761d1
md"""
## Transformations and AD functions

"""

# ╔═╡ a26861f4-72c4-422f-8665-eda6af5ff32c
function Aω_to_Aq(Aω, pA::Vector, pB::Vector)
	GA, GB = Gmat(pA), Gmat(pB)
	T = blockdiag(sparse(2GA), sparse(I, 3, 3), sparse(2GB), sparse(I, 3, 3))
	Aq = Aω * T
end

# ╔═╡ 5fbe31e5-961f-466b-809f-0f0b8b14ac41
Aω_to_Aq(Aω, bA::BodyState, bB::BodyState) = Aω_to_Aq(Aω, bA.p, bB.p)

# ╔═╡ df32ac4d-da15-44a2-80f3-d8c072fbd04d
#------------------------------------------------- Aq via ForwardDiff.jacobian --
function Aq_AD(cfun::Function, bA, bB)
    q0   = pack(bA, bB)
    _jfun = pvec -> begin
        _bA, _bB = unpack(pvec, bA, bB)
        cfun(_bA, _bB)
    end
    ForwardDiff.jacobian(_jfun, q0)
end

# ╔═╡ 2efb9c7c-eea1-4546-a23c-9cf74a3f07d0
function Aq_to_Aω(Aq, pA::Vector, pB::Vector)
    GA, GB = Gmat(pA), Gmat(pB)
	T = blockdiag(sparse(0.5GA'), sparse(I, 3, 3), sparse(0.5GB'), sparse(I, 3, 3))
   # maps [ω_B;v]→[ṗ;v]
    Aq * T                                      # 3×12
end

# ╔═╡ ea90acb0-59e0-466a-9ff2-8af6e40474b0
Aq_to_Aω(Aq, bA::BodyState, bB::BodyState) = Aq_to_Aω(Aq, bA.p, bB.p)

# ╔═╡ f5e1ff86-a36e-489e-9c4c-6d1cdb5e83de
ω_to_ṗ(ω, p) = 0.5 * Gmat(p)' * ω

# ╔═╡ 1c35f174-7419-49f4-bb8f-6e5dabdec1d3
md"""
## Body precompute
"""

# ╔═╡ 0176290a-1ff8-4c45-90a4-dde1e160dbd3
# ========== internal pre-compute (shared by all functions) ====================
function _prep(bA::BodyState, sA′,
               bB::BodyState, sB′)
    RA = rotmat(bA.p);  RB = rotmat(bB.p)
    sA = RA*sA′;         sB = RB*sB′
    rA = bA.r;          rB = bB.r
    ωA = bA.ω;          ωB = bB.ω
    (;RA,RB,sA,sB,rA,rB,ωA,ωB)
end

# ╔═╡ 402a6ddb-003d-4979-9135-794c5cf07196
# ========== demo ==============================================================
bodyA = BodyState{Float64}(
    [0.6425754631219991,
 -0.08032193289024989,
  0.24096579867074963,
  0.7228973960122489],  [0.1,0.2,0.3],
    [0,0,0],     [0,0,0]
)

# ╔═╡ 83917cb0-f9e6-4a35-a487-c2097011c8a9
bodyB = BodyState{Float64}(
    [cosd(15),0,sind(15),0],   # 30° about x
    [1,-0.7,0.4],
    [0,1,0],               # body-frame spin around +y
    [0,0,0]
)

# ╔═╡ d9bda04d-d218-42ed-8517-ea4af5f1a102
sA′ = [0.3,-0.1,0.2]

# ╔═╡ 7fb1909f-2f0b-40cc-b010-50193d47c531
sB′ = [1.0,0.7,-0.2]

# ╔═╡ 0eaea5fd-a169-40ea-8078-40e0d06f8f81
aA′ = [-0.3,0.2,-0.1]

# ╔═╡ b6a78e8b-16ae-40b4-9e43-8e9fdd6771c3
md"""
## Ball joint
"""

# ╔═╡ 45998a24-667f-41fa-b2b7-b131eefacc19
md"""
### Constraint equation
"""

# ╔═╡ b8c9ae41-100b-454a-bc54-d4fd4fbd6dd0
# ========== 1) constraint vector g ============================================
function g_ball_joint(bA, sA′, bB, sB′)
	C = _prep(bA, sA′, bB, sB′)
	C.rB + C.sB - C.rA - C.sA
end

# ╔═╡ e92263af-7da1-45b0-b92b-fdda7a52acc4
@testset "Ball joint call works" begin
	g = @test_nowarn g_ball_joint(bodyA, sA′, bodyB, sB′)
	@test length(g) == 3
	@test all(isfinite.(g))
end

# ╔═╡ ce9f7614-7543-4711-b07d-000977ab4379
md"""
### Jacobians
"""

# ╔═╡ 20b30699-a2a5-4818-939d-33a48e242122
# ========== 2) twist-space Jacobian Aω ========================================
function Aω_ball_joint(bA, sA′, bB, sB′)
    C = _prep(bA, sA′, bB, sB′)
    hcat(C.RA * skew(sA′), -I(3),  -C.RB * skew(sB′),  I(3))   # size 3×12
end

# ╔═╡ 24818e33-398b-4328-9806-a7c0f84965dd
# ========== 4) coordinate-space Jacobian Aq ===================================
function Aq_ball_joint(bA, sA′, bB, sB′)
	Aω = Aω_ball_joint(bA, sA′, bB, sB′)
	Aω_to_Aq(Aω, bA, bB)   # size 3×14
end

# ╔═╡ 925aac96-4940-4bf0-ae0e-1eb7493f1224
#------------------------------------------------- Aq via ForwardDiff.jacobian --
function Aq_ball_joint_AD(bA, sA′, bB, sB′)
    gfun(_bA, _bB) = g_ball_joint(_bA, sA′, _bB, sB′)
    Aq_AD(gfun, bA, bB)
end

# ╔═╡ 44bf0eed-9f27-470a-8c6d-2eba65856f1d
function Aω_ball_joint_AD(bA, sA′, bB, sB′)
	Aq = Aq_ball_joint_AD(bA, sA′, bB, sB′)
    Aq_to_Aω(Aq, bA, bB)                                     # 3×12
end

# ╔═╡ 93a7d201-6179-408b-86e9-0a26b87a753b
@testset "Compare ball joint Aω with its AD version" begin
	Aω = Aω_ball_joint(bodyA, sA′, bodyB, sB′)
	Aω_AD = Aω_ball_joint_AD(bodyA, sA′, bodyB, sB′)
	@test Aω ≈ Aω_AD
end

# ╔═╡ d846fa4f-b28c-46c1-b848-3fdf06364d2f
@testset "Compare ball joint Aq*q̇ with its AD version" begin
	Aq = Aq_ball_joint(bodyA, sA′, bodyB, sB′)
	Aq_AD = Aq_ball_joint_AD(bodyA, sA′, bodyB, sB′)
	q̇ = vcat(
		ω_to_ṗ(bodyA.ω, bodyA.p),
		bodyA.v,
		ω_to_ṗ(bodyB.ω, bodyB.p),
		bodyB.v
	)
	@test Aq * q̇ ≈ Aq_AD * q̇
end

# ╔═╡ 7d5d4ccf-ccf4-43fd-bc6e-02cd0a1154d7
md"""
## Orthogonality Constraint -- Type 1 Perpendicularity
"""

# ╔═╡ f60f0875-4600-4823-823b-05edde788c70
md"""
### Constraint equation
"""

# ╔═╡ 36f81f5c-2b28-4224-ab57-bb3bb11a32a3
# ========== 1) constraint vector g ============================================
function g_orthogonality(bA, sA′, bB, sB′)
	C = _prep(bA, sA′, bB, sB′)
	[C.sA ⋅ C.sB]
end

# ╔═╡ 428eae10-32c8-4d3a-ad13-dd46086b5db7
@testset "Ball joint call works" begin
	g = @test_nowarn g_orthogonality(bodyA, sA′, bodyB, sB′)
	@test length(g) == 1
	@test all(isfinite.(g))
end

# ╔═╡ 213fbd77-6f61-4270-8fbc-bc7b3d06d06d
md"""
### Jacobians
"""

# ╔═╡ bf2f22d1-2750-491b-911e-eb9f8bd12dc8
# ========== 2) twist-space Jacobian Aω ========================================
function Aω_orthogonality(bA, sA′, bB, sB′)
    C = _prep(bA, sA′, bB, sB′)
    hcat(-C.sB' * C.RA * skew(sA′), zeros(1, 3),  -C.sA' * C.RB * skew(sB′),  zeros(1, 3))   # size 1×12
end

# ╔═╡ 6978bc40-b9ef-46c7-af29-4d0bf9ae1be4
# ========== 4) coordinate-space Jacobian Aq ===================================
function Aq_orthogonality(bA, sA′, bB, sB′)
	Aω = Aω_orthogonality(bA, sA′, bB, sB′)
	Aω_to_Aq(Aω, bA, bB)   # size 1×14
end

# ╔═╡ 612cb996-7223-40e3-afc6-bad022f9a9fa
#------------------------------------------------- Aq via ForwardDiff.jacobian --
function Aq_orthogonality_AD(bA, sA′, bB, sB′)
    gfun(_bA, _bB) = g_orthogonality(_bA, sA′, _bB, sB′)
    Aq_AD(gfun, bA, bB)
end

# ╔═╡ 08e0e485-3d0d-4b26-853e-44c4d8d565f8
function Aω_orthogonality_AD(bA, sA′, bB, sB′)
	Aq = Aq_orthogonality_AD(bA, sA′, bB, sB′)
    Aq_to_Aω(Aq, bA, bB)                                     # 1×12
end

# ╔═╡ 13564452-ed30-490a-ae44-e78b15ce218c
@testset "Compare orthogonality Aω with its AD version" begin
	Aω = Aω_orthogonality(bodyA, sA′, bodyB, sB′)
	Aω_AD = Aω_orthogonality_AD(bodyA, sA′, bodyB, sB′)
	@test Aω ≈ Aω_AD
end

# ╔═╡ cd8f504c-53ba-4ee9-9e89-d1bab4afe944
@testset "Compare orthogonality Aq*q̇ with its AD version" begin
	Aq = Aq_orthogonality(bodyA, sA′, bodyB, sB′)
	Aq_AD = Aq_orthogonality_AD(bodyA, sA′, bodyB, sB′)
	q̇ = vcat(
		ω_to_ṗ(bodyA.ω, bodyA.p),
		bodyA.v,
		ω_to_ṗ(bodyB.ω, bodyB.p),
		bodyB.v
	)
	@test Aq * q̇ ≈ Aq_AD * q̇
end

# ╔═╡ aa93c255-c0c9-45d1-9c3c-8e46415f972e
md"""
## Constant Projection -- Type 2 Perpendicularity
"""

# ╔═╡ 8e4a4613-b299-40a9-b3e6-e3097d84bc56
md"""
### Constraint equation
"""

# ╔═╡ 8c465897-efa7-4698-be78-f9a804156224
# ========== 1) constraint vector g ============================================
function g_orthogonality_2(bA, aA′, sA′, bB, sB′)
	C = _prep(bA, sA′, bB, sB′)
	aA = C.RA * aA′
	d = C.rB + C.sB - C.rA - C.sA
	[aA ⋅ d]
end

# ╔═╡ 8e8e75ff-f749-4881-a0a1-2bb4d07ad11d
@testset "Ball joint call works" begin
	g = @test_nowarn g_orthogonality_2(bodyA, aA′, sA′, bodyB, sB′)
	@test length(g) == 1
	@test all(isfinite.(g))
end

# ╔═╡ 418b2863-e8df-479e-8515-bcab806773bd
md"""
### Jacobians
"""

# ╔═╡ f099af94-0bd9-4fee-b8cf-0760ca0888cf
# ========== 2) twist-space Jacobian Aω ========================================
function Aω_orthogonality_2(bA, aA′, sA′, bB, sB′)
    C = _prep(bA, sA′, bB, sB′)
	aA = C.RA * aA′
	d = C.rB + C.sB - C.rA - C.sA
    collect(hcat(
		aA′' * skew(sA′) - d' * C.RA * skew(aA′), 
		-aA',  
		-aA' * C.RB * skew(sB′), 
		aA'
	))   # size 1×12
end

# ╔═╡ 682a8df8-52dc-4ed1-8210-231693c19268
# ========== 4) coordinate-space Jacobian Aq ===================================
function Aq_orthogonality_2(bA, aA′, sA′, bB, sB′)
	Aω = Aω_orthogonality_2(bA, aA′, sA′, bB, sB′)
	Aω_to_Aq(Aω, bA, bB)   # size 1×14
end

# ╔═╡ ad59816c-4f1a-47cb-b0ac-50b22f007164
#------------------------------------------------- Aq via ForwardDiff.jacobian --
function Aq_orthogonality_2_AD(bA, aA′, sA′, bB, sB′)
    gfun(_bA, _bB) = g_orthogonality_2(_bA, aA′, sA′, _bB, sB′)
    Aq_AD(gfun, bA, bB)
end

# ╔═╡ 9d1581f1-bac6-43b1-8eb8-4690408a5f9c
function Aω_orthogonality_2_AD(bA, aA′, sA′, bB, sB′)
	Aq = Aq_orthogonality_2_AD(bA, aA′, sA′, bB, sB′)
    Aq_to_Aω(Aq, bA, bB)                                     # 1×12
end

# ╔═╡ 35453b85-1a91-4253-91ed-773c2e096a03
@testset "Compare orthogonality_2 Aω with its AD version" begin
	Aω = Aω_orthogonality_2(bodyA, aA′, sA′, bodyB, sB′)
	Aω_AD = Aω_orthogonality_2_AD(bodyA, aA′, sA′, bodyB, sB′)
	@test Aω ≈ Aω_AD
end

# ╔═╡ 8c0a6edf-e556-4aa8-9103-a239f9a20225
@testset "Compare orthogonality_2 Aq*q̇ with its AD version" begin
	Aq = Aq_orthogonality_2(bodyA, aA′, sA′, bodyB, sB′)
	Aq_AD = Aq_orthogonality_2_AD(bodyA, aA′, sA′, bodyB, sB′)
	q̇ = vcat(
		ω_to_ṗ(bodyA.ω, bodyA.p),
		bodyA.v,
		ω_to_ṗ(bodyB.ω, bodyB.p),
		bodyB.v
	)
	@test length(Aq_AD * q̇) == 1
	@test Aq * q̇ ≈ Aq_AD * q̇
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
ForwardDiff = "~0.10.38"
PlutoUI = "~0.7.62"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "005369054785afef08ba47bee3a26c358e59603b"

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

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
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

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

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
# ╟─61976551-986d-4297-8f15-dad98120346a
# ╠═441cc450-405b-11f0-32d4-b533235d9307
# ╟─9425a84f-29af-4b74-815e-6f1c0876bd35
# ╟─f180ffcc-79e3-4554-a4d8-d7172dd6c133
# ╟─714a2372-2df1-422a-912b-a77223c022da
# ╠═ab0159aa-f6e5-4ca5-b2d4-f9e666e3902a
# ╠═ca7454e1-843d-488e-bda0-e1ca4e95236c
# ╠═b3e08dc3-a0aa-440e-a530-5ef5568694e9
# ╟─7b2a97cb-223a-4603-bd89-d813ff8a831b
# ╠═831493a6-a56f-4f8c-9862-ce331a615e32
# ╠═daae7f68-1249-409d-915e-53ba6aadbc2c
# ╠═173d35d7-ba79-4faa-bdf8-86a964f822d1
# ╠═7a1ba4ee-801b-4902-bf48-ce20521fead1
# ╟─a9fb78b9-061e-417a-a89b-e41c818761d1
# ╠═a26861f4-72c4-422f-8665-eda6af5ff32c
# ╠═5fbe31e5-961f-466b-809f-0f0b8b14ac41
# ╠═df32ac4d-da15-44a2-80f3-d8c072fbd04d
# ╠═2efb9c7c-eea1-4546-a23c-9cf74a3f07d0
# ╠═ea90acb0-59e0-466a-9ff2-8af6e40474b0
# ╠═f5e1ff86-a36e-489e-9c4c-6d1cdb5e83de
# ╟─1c35f174-7419-49f4-bb8f-6e5dabdec1d3
# ╠═0176290a-1ff8-4c45-90a4-dde1e160dbd3
# ╠═402a6ddb-003d-4979-9135-794c5cf07196
# ╠═83917cb0-f9e6-4a35-a487-c2097011c8a9
# ╠═d9bda04d-d218-42ed-8517-ea4af5f1a102
# ╠═7fb1909f-2f0b-40cc-b010-50193d47c531
# ╠═0eaea5fd-a169-40ea-8078-40e0d06f8f81
# ╟─b6a78e8b-16ae-40b4-9e43-8e9fdd6771c3
# ╟─45998a24-667f-41fa-b2b7-b131eefacc19
# ╠═b8c9ae41-100b-454a-bc54-d4fd4fbd6dd0
# ╠═e92263af-7da1-45b0-b92b-fdda7a52acc4
# ╟─ce9f7614-7543-4711-b07d-000977ab4379
# ╠═20b30699-a2a5-4818-939d-33a48e242122
# ╠═24818e33-398b-4328-9806-a7c0f84965dd
# ╠═925aac96-4940-4bf0-ae0e-1eb7493f1224
# ╠═44bf0eed-9f27-470a-8c6d-2eba65856f1d
# ╠═93a7d201-6179-408b-86e9-0a26b87a753b
# ╠═d846fa4f-b28c-46c1-b848-3fdf06364d2f
# ╟─7d5d4ccf-ccf4-43fd-bc6e-02cd0a1154d7
# ╟─f60f0875-4600-4823-823b-05edde788c70
# ╠═36f81f5c-2b28-4224-ab57-bb3bb11a32a3
# ╠═428eae10-32c8-4d3a-ad13-dd46086b5db7
# ╟─213fbd77-6f61-4270-8fbc-bc7b3d06d06d
# ╠═bf2f22d1-2750-491b-911e-eb9f8bd12dc8
# ╠═6978bc40-b9ef-46c7-af29-4d0bf9ae1be4
# ╠═612cb996-7223-40e3-afc6-bad022f9a9fa
# ╠═08e0e485-3d0d-4b26-853e-44c4d8d565f8
# ╠═13564452-ed30-490a-ae44-e78b15ce218c
# ╠═cd8f504c-53ba-4ee9-9e89-d1bab4afe944
# ╟─aa93c255-c0c9-45d1-9c3c-8e46415f972e
# ╟─8e4a4613-b299-40a9-b3e6-e3097d84bc56
# ╠═8c465897-efa7-4698-be78-f9a804156224
# ╠═8e8e75ff-f749-4881-a0a1-2bb4d07ad11d
# ╠═418b2863-e8df-479e-8515-bcab806773bd
# ╠═f099af94-0bd9-4fee-b8cf-0760ca0888cf
# ╠═682a8df8-52dc-4ed1-8210-231693c19268
# ╠═ad59816c-4f1a-47cb-b0ac-50b22f007164
# ╠═9d1581f1-bac6-43b1-8eb8-4690408a5f9c
# ╠═35453b85-1a91-4253-91ed-773c2e096a03
# ╠═8c0a6edf-e556-4aa8-9103-a239f9a20225
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
