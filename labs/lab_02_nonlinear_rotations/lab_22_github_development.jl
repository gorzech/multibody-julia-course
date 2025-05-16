### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ d76e7ad0-324f-11f0-2d4d-43b6632ef29e
md"""
# **Spatial Multibody Package Sprint ­– Session Plan (Updated)**  
*Pluto + Teams + GitHub; no coding required between meetings*

---

## **Timeline & Goals**

| Date | Time | Focus | “Done” for the block |
|------|------|-------|-----------------------|
| **Today** | **1 h** | **Design Huddle** | Core API & repo skeleton agreed, sprint board with 6 starter tickets |
| **Tomorrow** | 09:00-12:00 | **Sprint-01 A** – *Spatial Kinematics core* | `Body`, `Joint`, `assemble(model)`, Euler-parameter step green |
| | 13:00-16:00 | **Sprint-01 B** – *Constraints & Viz* | `constraints()` + `Cq()` via AD; NR loop; 2-body demo plot; tag **`v0.1-kinematics`** |
| **(2-week gap)** | — | **NO required coding** | Students may ignore repo |
| **+2 weeks** | 3 h | **Integration Lab** | All PRs merged; tests green; `M(q)` + gravity forces; tag **`v0.1.5-prepDynamics`** |
| **Final Lab** | 09:00-12:00 | **Sprint-02 A** – *Dynamics kernel* | RHS `q̈ = M⁻¹F`; Euler-Cromer on falling block example |
| | 13:00-16:00 | **Sprint-02 B** – *Polish & Release* | Pluto demo (3-body crank); docs build; tag **`v1.0-courseRelease`**; lightning talks |

---

## **Design Huddle (1 h today)**

1. **Vision brief** – scope ➜ *Euler params, 6-DOF bodies, AD Jacobian*.  
2. **Break-out decisions (15 min)** – type hierarchy, units, `StaticArrays` vs `Vector`.  
3. **Roles & repo tour** – Kinematics lead, Dynamics lead, Docs lead, CI sheriff.  
4. **Board populate** – 6 tickets into **Sprint-01** column.  
5. **Definition-of-Done for Sprint-01** locked (next slide).

---

## **Sprint-01 Definition of Done**

- ✅ `make_system` returns an empty system to pass around  
- ✅ `constraints()` & `jacobian()` via `ForwardDiff` pass tests  
- ✅ Pluto notebook reproduces 2-body demo plot  
- ✅ PR merged with green CI & formatter  
- ✅ Each student has ≥ 1 merged PR

---

## **Tomorrow – Flow**

### **Sprint-01 A (09:00-12:00)**  
*Stand-up → code → rolling PR reviews → micro-retro*

### **Sprint-01 B (13:00-16:00)**  
*Quick sync → code → final merges → tag `v0.1-kinematics` → 1-line reflection cell*

---

## **Integration Lab (+2 weeks, 3 h)**  

1. Green-build gate (fix tests).  
2. Add mass matrix `M(q)` & gravity; minimal docs.  
3. Tag `v0.1.5-prepDynamics`.

---

## **Final Lab – Sprint-02**

| Block | Tasks | Outcome |
|-------|-------|---------|
| **A (09:00-12:00)** | Implement `accelerations()` + Euler-Cromer step; validate falling block | Dynamics kernel green |
| **B (13:00-16:00)** | Docs deploy, Pluto demo, tag `v1.0`; 5-min team talks | Package release & showcase |

*Pass rule:* ≥ 1 merged PR **and** tests green **and** reflection note.

---

## **Logistics, Tools & “What’s a PR?”**

| Tool | Purpose |
|------|---------|
| **GitHub Pull Request (PR)** | *Branch → commit → **PR** → peer review → merge*.  PR ≤ 200 LOC, 1 reviewer ≠ author, CI green = merge. |
| **Teams** | Channel **#stand-up** (daily 3-liner), Files (Pluto notebooks), Excel progress board. |
| **Pluto** | Live coding & slide export. |
| **CI Actions** | Run `Pkg.test()` + `JuliaFormatter`. |

Stay inside the live blocks → keep PRs small → merge when green → **ship a real spatial multibody package in four sessions**.
"""

# ╔═╡ 893d7002-5a27-437f-8221-a5b1d2142995
md"""
### How we use Git & GitHub in this course

`TL;DR  =  branch ➜ commit ➜ push ➜ PR ➜ review ➜ merge`


|**Step**|**Command / Click**|**What it Does**|**Tips**|
|---|---|---|---|
|**1. Clone the repo**|git clone https://github.com/gorzech/MultibodyDynamicsLite.jl|Copies the remote repo to your laptop/lab PC|First time only|
|**2. Create a feature branch**|git checkout -b kinematics/step_euler|Keeps your work isolated|Name = subsystem / short task|
|**3. Code & commit**|git add src/kinematics.jl git commit -m "feat: Euler-param propagation"|Saves a coherent change-set|_Small_ commits (≤ 200 LOC)|
|**4. Push your branch**|git push --set-upstream origin kinematics/step_euler|Uploads branch to GitHub||
|**5. Open a Pull Request (PR)**|GitHub → “Compare & Pull request”|Asks for code review|Fill the template checklist|
|**6. Fix CI / review comments**|edit → commit → push (GitHub auto-updates the PR)|CI = tests + formatter|Green ✔ means passing|
|**7. Merge**|Click **“Merge Pull Request”** after approval & green CI|Your code lands on main|Delete branch when prompted|
|**8. Update local main**|git checkout main git pull|Keeps your copy fresh|Do this before starting the next task|

**Ground rules**

1. **One PR ≤ 200 LOC** (easy to review).
    
2. **1 reviewer ≠ author** (peer learning).
    
3. **CI must be green** (tests & formatter) before merge.
    
4. **No direct pushes to main.**
    
"""

# ╔═╡ c504abc4-0863-4b95-bfca-eca01bd14406
md"""
Cheat-sheet:

git status – what’s modified? 

git log --oneline --graph – branch history 

git switch - – last branch
"""

# ╔═╡ 90769042-4aa9-4401-8569-22b6e6cef99c
md"""
## Lab-Machine / Tool Checklist

|**Tool**|**Required version / note**|
|---|---|
|**Julia**|1.10 (pre-installed in lab image)|
|**Git**|v2.3+ with credential helper enabled|
|**VS Code**|Julia extension (for students who prefer editors)|
|**Pluto.jl**|Bundled via Pkg.add("Pluto") or preloaded|
|**ForwardDiff.jl**|In project [deps] – already added via the students’ PR|
|**StaticArrays.jl**|Ditto|
|**Browser**|Edge/Chrome/Firefox – for Pluto & GitHub|
|**Teams desktop app**|Chat + progress board tab|

_No additional licences or secrets are needed.  The CI runs entirely on free GitHub Actions minutes._
"""

# ╔═╡ bf242fe3-940a-4724-bb19-9b4636aeecc8
md"""
# **“Design Huddle” — 60-minute micro-workshop**  
*Kick-starts the MultibodyDynamicsLite.jl package project*

---

## 0. Prep (instructor, 5 min before start)

| Item | Action |
|------|--------|
| Teams meeting | Start/record; share screen on **Outline slide** |
| Repo | Ensure empty skeleton pushed to **main** |
| Shared doc | Open `DESIGN.md` in Teams (Word or OneNote live) |
| Project board | Open **Projects › Sprint-01** tab (columns: *Backlog · Sprint-01 · Review · Done*) |

---

## 1. Minute-by-minute plan

| Clock | Activity | Facilitator notes | Output artefact |
|------:|----------|-------------------|-----------------|
| **00 – 05** | **Welcome & scope** | *Why spatial MBS?* MVP = rigid bodies, Euler params, revolute & spherical ... joints, AD Jacobian | Everyone aligned |
| **05 – 10** | **Repo tour** | Show `src/`, `test/`, `.github/`, CI badge | Students know where to poke |
| **10 – 25** | **Break-out pairs** — *three design questions* | Split into 3 rooms (*)<br>1. **Data types**: `Body6D`, `Joint`, `Model`?<br>2. **Units**: plain Float vs Unitful?<br>3. **Linear algebra backend**: `StaticArrays` vs `Vector` | Each pair writes bullet answer in **DESIGN.md** |
| **25 – 35** | **Reconvene, decision log** | Read bullets aloud, quick thumbs-up/-down → pick default; unresolved gets **`design` Issue** | Updated **DESIGN.md** (+ Issues if needed) |
| **35 – 40** | **Role & branching policy** | • Kinematics lead, Dynamics lead, Docs lead, CI sheriff<br>• Branch naming: `subsystem/topic`<br>• 1 PR ≤ 200 LOC, needs 1 review | Names written in DESIGN.md |
| **40 – 55** | **Populate Sprint-01 board** | On projector create 6 starter Issues:<br>1. `feat: Body6D`<br>2. `feat: step_euler!`<br>3. `feat: assemble(model)`<br>4. `feat: constraints(revolute)`<br>5. `feat: jacobian via AD`<br>6. `docs: Getting-Started`<br>Drag them into **Sprint-01** column and self-assign | Visible Kanban board |
| **55 – 58** | **Define “Done” for Sprint-01** | Type into board description:<br>• Euler step passes test• `constraints` + `Cq` green<br>• Pluto 2-body demo runs<br>• ≥ 1 merged PR per student | DoD text block |
| **58 – 60** | **Logistics & adjourn** | Remind: **stand-up tomorrow**; clone repo tonight if possible | Everyone clear on next steps |

* (*) If only 5 students, keep discussion in plenary instead of breakout.

---

## 2. Facilitation checklist

- ☐ **Timer visible** (phone or Teams timer)  
- ☐ Keep **DESIGN.md** cursor at top so edits are obvious  
- ☐ For each decision write **R** (resolved) or **Q** (open question) tag in doc  
- ☐ After board is populated ask “Does every card have an owner?”  
- ☐ Screenshot final board & post in Teams channel

---

## 3. Artifacts leaving the room

| File / link | Where | Purpose |
|-------------|-------|---------|
| `DESIGN.md` | Repo root | Authoritative design record |
| Project board snapshot | Teams thread | Task list for tomorrow |
| Meeting recording | Teams → Files | Catch-up for absentees |

---

*Now the rails are laid; tomorrow’s two three-hour sprints can go straight into coding.*  
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
# ╟─d76e7ad0-324f-11f0-2d4d-43b6632ef29e
# ╟─893d7002-5a27-437f-8221-a5b1d2142995
# ╟─c504abc4-0863-4b95-bfca-eca01bd14406
# ╟─90769042-4aa9-4401-8569-22b6e6cef99c
# ╠═bf242fe3-940a-4724-bb19-9b4636aeecc8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
