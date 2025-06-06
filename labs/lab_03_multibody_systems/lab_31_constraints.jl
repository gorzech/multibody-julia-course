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

# ╔═╡ 42fd79b4-ca40-491a-aa7d-1bcd02bbe258
using FiniteDiff

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

# ╔═╡ 0176290a-1ff8-4c45-90a4-dde1e160dbd3
# ========== internal pre-compute (shared by all functions) ====================
function _prep(bA::BodyState, rA,
               bB::BodyState, rB)
    RA = rotmat(bA.p);  RB = rotmat(bB.p)
    cA = RA*rA;         cB = RB*rB
    rA = bA.r;          rB = bB.r
    ωA = bA.ω;          ωB = bB.ω
    (;RA,RB,cA,cB,rA,rB,ωA,ωB)
end

# ╔═╡ b8c9ae41-100b-454a-bc54-d4fd4fbd6dd0
# ========== 1) constraint vector g ============================================
function g_ball_joint(bA, rA, bB, rB)
	C = _prep(bA, rA, bB, rB)
	C.rB + C.cB - C.rA - C.cA
end

# ╔═╡ 20b30699-a2a5-4818-939d-33a48e242122
# ========== 2) twist-space Jacobian Aω ========================================
function Aω_ball_joint(bA, rA, bB, rB)
    C = _prep(bA, rA, bB, rB)
    hcat(C.RA * skew(rA), -I(3),  -C.RB * skew(rB),  I(3))   # size 3×12
end

# ╔═╡ 24818e33-398b-4328-9806-a7c0f84965dd
# ========== 4) coordinate-space Jacobian Aq ===================================
function Aq_ball_joint(bA, rA, bB, rB)
    GA, GB = Gmat(bA.p), Gmat(bB.p)
	Aω = Aω_ball_joint(bA, rA, bB, rB)
	T = blockdiag(sparse(2GA .+ bA.p'), sparse(I, 3, 3), sparse(2GB), sparse(I, 3, 3))
	Aq = Aω * T   # size 3×14
end

# ╔═╡ 925aac96-4940-4bf0-ae0e-1eb7493f1224
#------------------------------------------------- Aq via ForwardDiff.jacobian --
function Aq_AD(bA, rA, bB, rB)
    p0   = pack(bA, bB)
    gfun = pvec -> begin
        bA′, bB′ = unpack(pvec, bA, bB)
        g_ball_joint(bA′, rA, bB′, rB)
    end
    ForwardDiff.jacobian(gfun, p0)
	# FiniteDiff.finite_difference_jacobian(gfun, p0)
end

# ╔═╡ 44bf0eed-9f27-470a-8c6d-2eba65856f1d
function Aω_AD(bA, rA, bB, rB)
    GA, GB = Gmat(bA.p), Gmat(bB.p)
	Aq = Aq_AD(bA, rA, bB, rB)
	T = blockdiag(sparse(0.5GA'), sparse(I, 3, 3), sparse(0.5GB'), sparse(I, 3, 3))
   # maps [ω_B;v]→[ṗ;v]
    Aq * T                                      # 3×12
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

# ╔═╡ d2f7efd5-c1e2-49bb-abb0-7f5c660542b5
GA, GB = Gmat(bodyA.p), Gmat(bodyB.p)

# ╔═╡ 0f583fd5-833d-4a30-8831-13e3cadae215
T1 = blockdiag(sparse(0.5GA'), sparse(I, 3, 3), sparse(0.5GB'), sparse(I, 3, 3))

# ╔═╡ 8b81855f-500a-433d-8b19-ff4fe66806d0
T2 = blockdiag(sparse(2GA), sparse(I, 3, 3), sparse(2GB), sparse(I, 3, 3))

# ╔═╡ d9bda04d-d218-42ed-8517-ea4af5f1a102
rA = [0.3,-0.1,0.2]

# ╔═╡ 7fb1909f-2f0b-40cc-b010-50193d47c531
rB = [1.0,0.7,-0.2]

# ╔═╡ 7c030370-0356-438f-974b-52314018eeb8
Aω = Aω_ball_joint(bodyA, rA, bodyB, rB)

# ╔═╡ 25f353f6-6e27-4ea4-a69f-9bef135e744f
Ap = Aω * T2

# ╔═╡ dc436914-a57d-4c3f-a75b-1388f942f8a9
Aω - Ap * T1

# ╔═╡ 9d029d01-ad30-486b-ae8c-1f11d6f6b615
Ap2 = Aq_AD(bodyA, rA, bodyB, rB)

# ╔═╡ 6443ac96-2df7-428c-9f27-9610d8b155f4
Aω2 = Ap2 * T1

# ╔═╡ f0dc7e08-a62f-4f5f-be9a-a3d7ded5295c
Ap2 - Aω2 * T2

# ╔═╡ e2ceb5b3-c460-4cb0-85c1-fa2d1ed5d2e7
md"""
g = 
"""

# ╔═╡ 55fa8d56-0933-42c1-abff-f47659dc6fb8
g_ball_joint(bodyA, rA, bodyB, rB)

# ╔═╡ bcc69f84-5c3e-4b91-9cf3-3c25aefe8931
md"Aω  ="

# ╔═╡ b51ac72e-60ab-4c93-b77c-6c0b789112dc
md"Aq  ="

# ╔═╡ 2a392a53-d288-4314-bf50-f0fa0c2e5d68
Aq_ball_joint(bodyA, rA, bodyB, rB)

# ╔═╡ c49a1969-231c-4de3-9a46-6f4bfb16cae0
Aq_AD(bodyA, rA, bodyB, rB)

# ╔═╡ 5a821300-4a0d-45e5-9a4c-8752684e182c
Aq_ball_joint(bodyA, rA, bodyB, rB) - Aq_AD(bodyA, rA, bodyB, rB)

# ╔═╡ 5741e629-c94c-4722-b316-8c56df8de428
md"""
Aω_AD  =
"""

# ╔═╡ 24306d5f-a731-4778-b43b-78c173c1ec91
Aω_AD(bodyA, rA, bodyB, rB)

# ╔═╡ 2fde2821-3734-4e14-adf0-7e5aca958868
md"""
Aω\_ball\_joint = 
"""

# ╔═╡ 6a72636b-a262-4fdb-bb84-5846209f2f46
Aω_ball_joint(bodyA, rA, bodyB, rB)

# ╔═╡ 0367e53f-5980-4dfd-b130-f1540382d1f8
Aω_ball_joint(bodyA, rA, bodyB, rB) - Aω_AD(bodyA, rA, bodyB, rB)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
FiniteDiff = "~2.27.0"
ForwardDiff = "~0.10.38"
PlutoUI = "~0.7.62"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.5"
manifest_format = "2.0"
project_hash = "0fa2669321de02e4e7e736374728addf49006448"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

    [deps.Adapt.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "bebb10cd3f0796dd1429ba61e43990ba391186e9"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

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

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "f089ab1f834470c525562030c8cfde4025d5e915"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.27.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

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
# ╠═2941e478-38fb-4fb4-8cc9-58a867657825
# ╠═71d305d6-929a-458e-9a07-0e8cb9a2d865
# ╟─7b2a97cb-223a-4603-bd89-d813ff8a831b
# ╠═831493a6-a56f-4f8c-9862-ce331a615e32
# ╠═daae7f68-1249-409d-915e-53ba6aadbc2c
# ╠═173d35d7-ba79-4faa-bdf8-86a964f822d1
# ╠═7a1ba4ee-801b-4902-bf48-ce20521fead1
# ╠═0176290a-1ff8-4c45-90a4-dde1e160dbd3
# ╠═b8c9ae41-100b-454a-bc54-d4fd4fbd6dd0
# ╠═20b30699-a2a5-4818-939d-33a48e242122
# ╠═24818e33-398b-4328-9806-a7c0f84965dd
# ╠═42fd79b4-ca40-491a-aa7d-1bcd02bbe258
# ╠═925aac96-4940-4bf0-ae0e-1eb7493f1224
# ╠═44bf0eed-9f27-470a-8c6d-2eba65856f1d
# ╠═d2f7efd5-c1e2-49bb-abb0-7f5c660542b5
# ╠═0f583fd5-833d-4a30-8831-13e3cadae215
# ╠═8b81855f-500a-433d-8b19-ff4fe66806d0
# ╠═7c030370-0356-438f-974b-52314018eeb8
# ╠═25f353f6-6e27-4ea4-a69f-9bef135e744f
# ╠═dc436914-a57d-4c3f-a75b-1388f942f8a9
# ╠═9d029d01-ad30-486b-ae8c-1f11d6f6b615
# ╠═6443ac96-2df7-428c-9f27-9610d8b155f4
# ╠═f0dc7e08-a62f-4f5f-be9a-a3d7ded5295c
# ╠═402a6ddb-003d-4979-9135-794c5cf07196
# ╠═83917cb0-f9e6-4a35-a487-c2097011c8a9
# ╠═d9bda04d-d218-42ed-8517-ea4af5f1a102
# ╠═7fb1909f-2f0b-40cc-b010-50193d47c531
# ╟─e2ceb5b3-c460-4cb0-85c1-fa2d1ed5d2e7
# ╠═55fa8d56-0933-42c1-abff-f47659dc6fb8
# ╟─bcc69f84-5c3e-4b91-9cf3-3c25aefe8931
# ╟─b51ac72e-60ab-4c93-b77c-6c0b789112dc
# ╠═2a392a53-d288-4314-bf50-f0fa0c2e5d68
# ╠═c49a1969-231c-4de3-9a46-6f4bfb16cae0
# ╠═5a821300-4a0d-45e5-9a4c-8752684e182c
# ╟─5741e629-c94c-4722-b316-8c56df8de428
# ╠═24306d5f-a731-4778-b43b-78c173c1ec91
# ╟─2fde2821-3734-4e14-adf0-7e5aca958868
# ╠═6a72636b-a262-4fdb-bb84-5846209f2f46
# ╠═0367e53f-5980-4dfd-b130-f1540382d1f8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
