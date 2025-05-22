### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ f09aec17-c235-4c1c-bf32-c9d9cd345feb
md"""
# Programming Multibody Systems in Julia

### Lecture 7: Julia Package Development Guide

Grzegorz Orzechowski
"""

# ╔═╡ 58a2b3d4-2990-4599-b148-0432826b801b
md"""
## Introduction

Welcome! In this Pluto.jl slide deck, we’ll explore how to **get started with Julia packages** – from creating your own package to testing and debugging it. These slides are geared towards Julia beginners (with some numerical methods background) who have used Pluto notebooks and are ready to organize code into packages. We will cover the following topics:

- **Creating a New Package** with `Pkg.generate` and managing it with `Pkg.dev`
    
- **Revise.jl** for live coding and automatic reloading of code changes
    
- **Best Practices for Structs** – defining types and dealing with changes
    
- **Modules & Exports** – organizing code and public API conventions
    
- **Testing** your package with Julia’s `Test.jl` framework
    
- **Debugging** Julia code (especially using VS Code’s tools)
    
- **Using VS Code** with the Julia plugin for development productivity
    
- **StaticArrays.jl** for efficient small fixed-size matrices/vectors
    
- **Toy Project Example** – a mini project (Newton’s method solver) that ties these ideas together, plus an alternative project idea (2-mass spring system)
    

"""

# ╔═╡ f5887bab-576a-4991-a077-58413e7f3a5f
md"""
## Creating a New Package (Pkg.generate & Pkg.dev)

  

When starting a Julia project that you intend to reuse or share, it’s best to create it as a **package**. Julia’s package manager Pkg can scaffold a new package and help manage its environment:

1. **Generate a package template:** In the Julia Pkg REPL (press `]` in the REPL), use the generate command. For example:
    

```julia
(@v1.x) pkg> generate MyPackage
  Generating  project MyPackage:
    MyPackage/Project.toml
    MyPackage/src/MyPackage.jl
```

This creates a new folder **MyPackage/** with a `Project.toml` (package metadata) and a starter source file `src/MyPackage.jl`.
    
2. **Activate and develop your package:** Navigate into the new package directory and activate it as an environment. You can then use `Pkg.dev` to develop the package in your main environment. For example, from the Pkg REPL:
    

```julia
(@v1.x) pkg> activate MyPackage
(MyPackage) pkg> dev .
```

Using `dev` registers your package in the active environment by path, so you can using `MyPackage` in that environment without needing to publish it. The `dev` command tells Julia to track the local package directory (so changes are reflected). After `dev`, your global environment’s `Project.toml` will list MyPackage with a path, meaning it’s using your local code .
    
3. **Add dependencies:** Within the package’s environment (e.g. (`MyPackage`) prompt), use `Pkg.add("DependencyName")` to add any packages your code needs. This updates MyPackage’s `Project.toml` `[deps]` and Manifest. For example, if your package uses `StaticArrays`, do `pkg> add StaticArrays`. These deps will be isolated to MyPackage’s environment.
    
4. **Write your package code** in `src/MyPackage.jl` (more on module setup on the next slides). After generating, that file contains a placeholder `module MyPackage ... end`. Fill it with your types, functions, and export statements.
    
5. **Using the package:** With the package developed (or after `dev`), you can use it in Julia just like any installed package. For example, in the Julia REPL (with the global `env` active), `using MyPackage` will load your code (and precompile if needed). If Julia can’t find your package, ensure you ran `Pkg.dev` in the right environment or that the package is activated.
    

  

**Why use generate and dev?** This workflow ensures your package has its own environment and UUID, and `dev` lets you work on it locally while still being able to load it in other projects. It keeps dependencies tidy and makes it easy to eventually publish your package.

"""

# ╔═╡ 833bfc7b-b5dd-44ea-8d18-c98c3f53525d
md"""
## Revise.jl for Live Coding and Reloading

  

One of the biggest boosts to productivity in Julia development is **Revise.jl**, which automatically reloads code changes into a running session. This means you can edit your package code and **see changes without restarting Julia**.

- **What Revise does:** Revise watches the files of your modules and updates function definitions on the fly whenever you save changes. For example, if you edit a function in `MyPackage/src/MyPackage.jl`, Revise will update the method in your running session so that calling it reflects the new code. This enables a workflow where you run your package, test something, tweak the code, and re-run without restarting Julia each time.
    
- **Using Revise:** Simply do `using Revise` at the start of your Julia session (before loading your package). If you’re working in VS Code, note that the Julia extension **loads Revise by default** in the REPL, so you often get this functionality automatically. In a Pluto notebook, you might not need Revise (since Pluto re-reacts to code changes), but for developing packages in a REPL or VS Code, it’s invaluable.
    
- **Limitations:** Revise can handle redefinitions of functions and most code changes, but **it cannot redefine types (structs) or constants** in a live session. This is actually a limitation of Julia itself: once a `struct` type is defined, you cannot change its field types or add fields without restarting. Revise will notify you if it detects an impossible change (e.g., modifying a struct). In such cases, you’ll have to stop and restart your Julia session to incorporate the change.
    
- **Best practices with Revise:** Load Revise at the very beginning (`using Revise`) so it tracks any subsequent `using MyPackage` or `include` calls. Edit your package source files with your favorite editor; when you save, Revise takes care of updating the code in the running session. This gives a live-coding feel even though Julia is compiled. Many Julia developers keep a session open for days using Revise to avoid recompiling everything repeatedly.
    
- **Note:** Since you can’t change types on the fly, plan your data structures ahead (or see next slide for some strategies). But for function logic, Revise covers you – you can iteratively improve and test your code quickly.
    

"""

# ╔═╡ 5481a651-d77b-4966-911f-3ff2bbdf050a
md"""
## Best Practices: Defining and Changing struct Types
 

Defining your own structs (custom types) is central to Julia, but it comes with a caveat: you **cannot redefine a struct in the same session**. Here are some tips for working with structs during development:

- **Immutability vs Mutability:** By default, a `struct` in Julia is immutable (fields can’t be changed after construction, which is good for performance). You can use `mutable struct` if you need to change field values of an instance. However, **neither can be redefined completely without restarting** the session. Design your struct with the fields it needs from the start.
    
- **Plan your type design:** Think about the data and operations you’ll need. If you realize you forgot a field or used a wrong type, you’ll have to either make a new type or restart and redefine. For example, if you defined `struct Particle; x::Float64; v::Float64; end` but you need an acceleration field, you can’t just add it on the fly – you must make a new definition (or restart Julia). Revise will throw an error if it detects a struct definition change.
    
- **Workarounds for development:** During early prototyping, some developers use tricks to avoid restarts. One strategy is to use **placeholder types** like NamedTuples or dictionaries for quick prototyping of data fields, then convert to a formal struct when things stabilize. Another common trick is to define your struct with a temporary name and update it (e.g., define `StructName1`, then later redefine as `StructName2` and update any references), essentially versioning the type name. This is a bit tedious, and often it’s simpler to restart your Julia session after a major type change. The key point is that changing a type’s definition is a big deal in Julia’s compiled world.
    
- **Define minimal fields:** Only store what you absolutely need in the struct’s fields (for instance, avoid including redundant data that can be computed from other fields). This will reduce how often you feel the need to change the type later. Additional computed properties can often be provided via functions instead of stored fields.
    
- **Use constructors for flexibility:** You can provide custom inner or outer constructors for your struct to accept different input forms. This way, if the way you construct the object changes, you might adjust a constructor without changing the struct fields themselves.
    
- **Mutable vs Immutable performance:** Immutable structs are generally faster and better for most purposes (they can be stack-allocated and inlined by the compiler). Use `mutable struct` only if you truly need to modify the fields (e.g., for algorithmic reasons or performance in a tight loop where changing an object in-place is needed). Even with a `mutable struct`, you **still cannot add or remove fields** at runtime; “mutable” only means you can change field values, not the type’s definition.

Remember, this limitation is only about the development phase. Once your types are settled, you likely won’t be redefining them anyway. It’s a small upfront design consideration for a big gain in performance and clarity in Julia code.

"""

# ╔═╡ c5aa9dd2-aacd-4e39-aceb-693679b03ed8
md"""
## Module and Export Conventions

Julia packages are just modules (namespaces) with some extra file structure. Here’s how to organize your package module and use exports effectively:

- **Module setup:** In `src/MyPackage.jl`, you’ll typically have:
    

```julia
module MyPackage
# ... import or using statements for dependencies ...

export foo, Bar  # names to expose publicly

# internal includes (if you split code into multiple files)
include("utilities.jl")
include("core.jl")

# define your functions, types, etc.
struct Bar
    x::Float64
    y::Float64
end

foo(b::Bar) = 2b.x + b.y  # example function

end # module
```

- When you do `using MyPackage`, Julia will execute this file and bring the module into scope. The `module ... end` block defines a global namespace for your package.
    
- **Exports:** The `export` keyword inside a module declares which names (functions, types, constants) should be exported, i.e., made available to users of your package without needing the `MyPackage`. prefix. For example, if you `export foo, Bar`, a user can just do `using MyPackage; foo(x)` instead of `MyPackage.foo(x)`. **Only export your public API** – the functions and types you expect users to call directly. Keep internal helper functions or intermediate names un-exported (users can still access them via `MyPackage.func` if needed, but they won’t clutter the namespace on `using MyPackage`).
    
- **Including multiple source files:** As your package grows, you might split code into multiple files for organization (as hinted with `include("utilities.jl")`, etc.). These included files should not have their own module declarations (they exist within the main module). It’s common to have one top-level module (named after your package) and include other `.jl` files that contribute to that module. For example, you might have `src/MyPackage.jl` that includes `src/somefeature.jl`, etc., all under `module MyPackage ... end`.
    
- **Module naming:** By convention, the module name matches the package name and is in CamelCase (e.g., package `MyPackage.jl` contains module `MyPackage`). Also, source file names typically match the module (with the `.jl` extension).
    
- **Project structure recap:** A Julia package has a layout like:
    

```
MyPackage/       (directory)
  Project.toml   (package metadata: name, UUID, deps, compat)
  src/
      MyPackage.jl   (main source file defining module MyPackage)
  test/
      runtests.jl    (test script, if you have tests)
```

- The `[deps]` in `Project.toml` lists package dependencies; when you `Pkg.add` things in the package’s env, this file updates. The UUID in `Project.toml` uniquely identifies your package (generated by `Pkg.generate`). Julia uses it to avoid name conflicts.
    
- **Import vs using:** Inside your module, if you rely on other packages, you can use `using SomeDependency` to bring its exports in, or `import SomeDependency` to bring it in without exports. Use whichever is appropriate (if you only need a few names, you might do `import X: name1, name2`). Be mindful of extending functions from other modules (use `import` for the function before extending via new methods).

In summary, treat the module as the encapsulation of your package’s code. Use exports to provide a clean interface, and keep everything else internal to avoid namespace pollution for users.

"""

# ╔═╡ 9fb6a51c-f183-49f4-b584-083d2391c631
md"""

## Testing Your Package with Test.jl

Testing is crucial, even for small projects. Julia’s built-in **Test** standard library makes it easy to write and run tests:

- **Test environment setup:** Ensure that `Test` is listed as a dependency in your package. When you generated the package, a `test` directory with a `runtests.jl` might be created (depending on Julia version). If not, create `test/runtests.jl`. Also `add Test` to your package’s environment (e.g., `Pkg.add("Test")` in the package environment) – in modern Julia, generating a package often includes `Test` by default in the `Project.toml` `[extras]` and creates a test file.
    
- **Writing tests:** In `test/runtests.jl`, start by `using MyPackage` and `using Test`. Then you can define test sets. For example:
    

```julia
# test/runtests.jl
using MyPackage, Test

@testset "MyPackage functionality" begin
    # Simple example tests:
    @test myfunction(2) == 4    # expected output
    @test myfunction(0) == 0

    # You can group further using nested @testset if needed
    @testset "Edge cases" begin
        @test_throws DomainError myfunction(-1)  # example of testing an error
    end
end
```

- Here we use `@test` for checking conditions (each will report pass/fail) and `@testset` to group related tests with a label. There are also macros like `@test_throws` (to verify that calling a function throws a certain exception), `@testsets` can be nested, etc.
    
- **Running tests:** To run your package’s tests, you have a few options:
    
    - **Pkg REPL:** Activate the package (or just be in the main environment if you did `dev`) and do `pkg> test MyPackage`. This will load `test/runtests.jl` and execute it. You will see a summary of test results (e.g., how many passed) .
        
    - **Programmatically:** In a Julia session, you can do using `Pkg; Pkg.test("MyPackage")`.
        
    - **VS Code test explorer:** If you’re using VS Code, the Julia extension might show a “test” icon or allow you to run tests from the UI (it discovers tests by looking at the `test` folder).
        
    - **From Pluto:** Pluto is not typically used to run package tests, but you could in theory include `using Pkg; Pkg.test("MyPackage")` in a Pluto cell (though it’s more common to run tests via REPL or CI).
        
    
- **Continuous Integration (CI):** If your students eventually put code on GitHub, they can set up GitHub Actions to run `Pkg.test` on each push, which is great practice. For now, running tests locally is fine.
    
- **Why test?** Writing tests ensures your package works as expected and helps catch regressions when you modify code. It also provides usage examples for anyone reading the tests. Encourage students to test both typical cases and edge cases (e.g., what if input is negative, zero, NaN, etc., depending on context).
    

  

Remember to include `using Test` in your test script. The Test stdlib provides many tools for writing robust tests (like nearly-equal comparisons with ≈, or @testset setup/teardown functionality), which they can explore as needed .

"""

# ╔═╡ a5560aa7-e985-4504-9534-8c71fae525c4
md"""

## Debugging in VS Code

  

Debugging Julia code can be done with print statements (@show or println debugging), but an interactive debugger is often more efficient for complex issues. **VS Code’s Julia extension** includes a powerful debugger interface. Here’s how to use it:

- **Setting up:** Open your project folder in Visual Studio Code with the Julia extension installed. Make sure the extension is using the correct Julia executable (check the Julia extension settings if necessary). Typically, you’ll start a REPL and include your package using Revise for iterative work. But for step-by-step debugging, you use the debugger.
    
- **Launching the debugger:** You can debug a Julia file by opening it and clicking the **Run and Debug** button (or pressing F5). If you have a Main.jl or some script to run, that could be your entry. For package code, you might write a small script or use the test script to drive the debugger. When you start debugging, VS Code opens a **Julia Debug** console and runs your code in debug mode.
    
- **Breakpoints:** Set breakpoints by clicking in the left margin next to a line number in your source code. A red dot will appear indicating the breakpoint. When running under the debugger, the execution will pause at that line before executing it . You can set multiple breakpoints, including conditional breakpoints (right-click to add a condition).
    
- **Stepping through code:** Once paused at a breakpoint, you’ll see a toolbar with commands: _Continue_, _Step Over_, _Step Into_, _Step Out_, etc. _Step Into_ goes into function calls, _Step Over_ executes the current line without diving into it, and _Continue_ resumes until the next breakpoint (or program end). As you step, you can observe the flow of control.
    
- **Inspecting variables:** When paused, the **Variables** panel in VS Code shows local variables and their values at the current scope. You can expand complex structures to see their fields. This is extremely useful to check if values are what you expect at a given point in execution.
    
- **REPL integration:** While paused, you can also use the **Debug Console** to execute Julia commands in the current debug context. For instance, you can call a function or check the value of an expression to probe deeper into state.
    
- **Performance note:** The debugger works by interpreting code (using the lower-level Julia interpreter). This is slower than normal execution. You might notice that code runs much slower in debug mode. That’s normal – use it for logical debugging, not for timing/performance checks. There is a option called “Compiled Mode” where you can exclude certain functions from interpretation to speed things up if needed.
    
- **Alternative: Debugger.jl in REPL:** Another way is to use the Debugger.jl interface via REPL (for example, using the @enter macro on a function call, or Juno’s debugger if anyone uses Juno/Atom). However, for beginners, the VS Code graphical debugger is more intuitive.
    

  

Encourage students to use the debugger to step through their code, especially if they are getting unexpected results. Seeing exactly what each part of the code is doing (and the state of variables) is invaluable for learning and bug-fixing.

"""

# ╔═╡ 16482ac0-83ec-4253-ac08-2c3a182ba86a
md"""

## Using VS Code with Julia

  

VS Code isn’t just for debugging – it’s a **full-featured development environment for Julia**. Some tips on using VS Code effectively with Julia:

- **Install the Julia extension:** This is provided by Julia Computing/JuliaLang. Once installed, it adds Julia syntax highlighting, error highlighting, integrated REPL, debugger, and more. Make sure to configure the extension’s Julia path to your Julia installation if it doesn’t auto-detect.
    
- **Integrated REPL:** You can launch a Julia REPL inside VS Code (for example, by pressing Ctrl+Shift+P and selecting “Julia: Start REPL” or simply use the Julia Editor run button). The REPL will appear in a terminal panel. One great feature is inline execution: you can send code from the editor to the REPL (e.g., pressing Alt+Enter will execute the current line or selection in the REPL). This is very handy for testing functions in your package as you write them. With Revise loaded, if you edit and save your package file, the REPL (with using MyPackage) will automatically update.
    
- **Code intelligence:** VS Code provides autocompletion, signature help, and documentation on hover for Julia code. It uses the Language Server. This means as you type, you might get suggestions for function names or see docstrings. Encourage students to use this to discover functions or recall argument orders.
    
- **Linting:** The Julia extension will underline syntax errors or potential issues in your code as you edit, which can catch mistakes early (like undefined variables, type mismatches, etc.). Pay attention to those squiggly lines!
    
- **Plot pane and tables:** If your code produces plots (using packages like Plots.jl) or displays DataFrames, VS Code can show these in dedicated panes. It’s useful for visualizing results without leaving the environment.
    
- **Workspace view:** There is a **Julia Workspace** view (in the Side Panel) that shows global variables and their values in the active Julia session. This is similar to what you may have seen in MATLAB or other IDEs. You can inspect arrays, etc., from there. (Note: With recent updates, the workspace is being revamped, but the concept remains.)
    
- **Version control integration:** If the students use git (which is a good idea even for small projects), VS Code will highlight changes, allow commits, etc., all within the editor. This is not Julia-specific, but a good development practice.
    
- **Using Pluto and VS Code together:** While Pluto notebooks are great for exploration and teaching, VS Code is often better for structuring code into a package. It’s perfectly fine to use Pluto for trying things out (or for interactive visualization) and then move stable code into a package that you edit in VS Code. VS Code can even preview Pluto notebooks or run them, but that’s optional. The key is to know when to transition from notebook-style to package development: usually when code needs to be reused, tested, or grown in complexity.
    

  

In summary, VS Code with the Julia plugin provides a robust environment for coding, much like an IDE for other languages, and can greatly help with organizing projects and improving productivity. Given that the Julia extension even loads Revise.jl automatically and has nice debug capabilities, it complements the workflow of package development perfectly.

"""

# ╔═╡ 99b50fd7-bf78-4e8b-b5f8-ce8d7291fc87
md"""

## StaticArrays.jl for Small Fixed-Size Matrices/Vectors

  

When working with linear algebra or geometry in Julia, one common scenario is dealing with **small vectors or matrices** – for example, 3D coordinates, rotation matrices, small Jacobians, etc. **StaticArrays.jl** is a package that provides a specialized array type for this use case, offering major performance benefits for small sizes.

- **What are StaticArrays?** They are fixed-size array types whose size is known at compile time (e.g., a 3-element vector or a 3x3 matrix). The StaticArrays package provides types like SVector{N,T} (immutable static vector of length N) and MVector{N,T} (mutable static vector), as well as SMatrix{R,C,T} for matrices, among others. Under the hood, a static array is essentially represented as a tuple in memory, which the compiler can fully unroll and optimize .
    
- **Why use them for small sizes?** For small arrays (a typical rule of thumb is up to 10×10 or so ), StaticArrays can be much faster than regular Array. There are a few reasons:
    
    1. _No heap allocations:_ StaticArrays are allocated on the stack (or even in CPU registers for very small ones) rather than the heap, so creating and destroying them is cheap (no garbage collection overhead) .
        
    2. _Compile-time unrolling:_ Operations on StaticArrays are unrolled at compile-time. For example, adding two SVector{3,Float64} will generate three separate add instructions – effectively no loop, just straight-line code for each element . This eliminates loop overhead and can enable the compiler to optimize aggressively (even vectorize operations via SIMD).
        
    3. _Inlined in structures:_ If you have a struct that holds a static vector, the whole thing can often be stored without pointers, making memory access very cache-friendly .
        
    
- **Using StaticArrays:** First, add the package (Pkg.add("StaticArrays")) and using StaticArrays. You can create a static array with the @SVector or @SMatrix macros for convenience, or by calling the type constructors. For example:
    

```julia
using StaticArrays
v = @SVector [1.0, 2.0, 3.0]         # SVector{3,Float64}
w = @SVector [2.0, 1.0, 0.0]         # another 3-element static vector
sum_v = v + w                        # yields SVector{3,Float64}([3.0, 3.0, 3.0])
M = @SMatrix [1 2; 3 4]              # 2x2 static matrix
detM = det(M)                        # StaticArrays defines linear algebra for SMatrix
```

- StaticArrays behave a lot like regular arrays, but you’ll notice they print with their dimension in the type. You can index into them (v[1] etc.) just like normal. Many functions from Julia’s LinearAlgebra work with them (det, inv, solve, etc.), or have specialized methods for static arrays.
    
- **When not to use StaticArrays:** If the array size might be large or not known at compile time, stick to regular Array. StaticArrays generate specialized code for each size, so using a static array of size 1000x1000 would lead to enormous compiled code (and heavy compilation time). They shine for small, fixed problems like 2D/3D physics, small matrices in algorithms, etc. A good rule: if you can count the size on your fingers, StaticArrays might help; if it’s data-sized (like 1000s of elements), use regular arrays.
    
- **Example scenario:** In a physics simulation with a few bodies, using SVector{3,Float64} for positions and velocities can drastically cut down on allocations and improve speed. Or for a Newton’s method Jacobian in 2 variables, using an SMatrix{2,2} could make each iteration faster by avoiding creating temporary heap arrays for Jacobians.
    


In our upcoming toy project example, we’ll highlight where StaticArrays could be applied for efficiency. But even if students don’t use StaticArrays immediately, it’s good to know about this tool as they start writing performance-sensitive Julia code involving small vector/matrix math.

"""

# ╔═╡ cbf450ac-911f-4ec2-82fc-a47237dec16b
md"""

## Toy Project Idea: Newton’s Method Solver

Now that we’ve covered the toolbox, let’s apply it in a small project. One suggested toy project is to implement a **basic nonlinear equation solver using Newton’s method**. This can solidify concepts like module structure, using external packages (for possibly linear algebra or StaticArrays), and writing tests. The idea:

- **What to build:** A Julia package (let’s call it NewtonSolver.jl) that provides a function to find roots of a nonlinear equation or system using Newton’s method. Newton’s method is an iterative algorithm to solve f(x) = 0 by linearizing f around a guess and iteratively improving that guess: $x_{new} = x - J(f)^{-1} * f(x)$ (where J is the Jacobian matrix for multivariate f, or derivative for single-variable case).
    
- **Scope:** For simplicity, the solver can assume the user provides the function $f$ and its Jacobian $J$. This avoids implementing automatic differentiation or finite differencing (which is possible but beyond a basic intro). We can handle both the scalar case (f: ℝ -> ℝ) and small vector case (f: ℝ^n -> ℝ^n). This is a good use case for StaticArrays if n is small, since each iteration might involve inverting a small Jacobian matrix.
    
- **Basic algorithm:**
    
    1. Start with an initial guess x0 (scalar or vector).
        
    2. Compute f(x0) and J(x0).
        
    3. Solve the linear system J(x0) * Δ = f(x0) for the step Δ (in scalar case, Δ = f(x0)/f’(x0)). In code, we can do Δ = J(x) \ f(x) – Julia’s \ will do the right thing for numbers or matrices.
        
    4. Update the guess: x1 = x0 - Δ.
        
    5. Check convergence: if ‖f(x1)‖ is below a tolerance, or the change ‖Δ‖ is small, we consider it converged.
        
    6. Otherwise, repeat from step 2 with x1. Also put a safeguard on maximum iterations.
        
    
- **Example problem:** A simple test could be finding √2 by solving f(x) = x^2 - 2 = 0. For a 2D example, one could solve a system like:
    
    f_1(x,y) = x^2 + y^2 - 1 = 0 \\ f_2(x,y) = x - y = 0,
    
    which has solutions (x=y=±1/√2). Newton’s method can find one of them given a decent initial guess.
    
- **Project structure:** We will create NewtonSolver.jl with a module NewtonSolver that exports a function (e.g., newton_solve). The package will have tests to verify the solver on known problems.
    
- **Why this project?:** It’s simple enough to do in a few lines of code, but non-trivial in that it uses linear algebra and demonstrates the benefit of Julia’s multiple dispatch (the same code can work for scalars or StaticVector inputs, for instance). It also gives a taste of numerical methods in Julia.
    
- **Alternate project idea:** Instead of Newton’s method, another option is a **minimal multibody system simulation** – e.g., two masses connected by a spring (and maybe to walls) with ODE integration. That would involve formulating differential equations m*x'' = -kx etc., and solving them over time (possibly using the DifferentialEquations.jl library or a simple integrator loop). This is a bit more involved, but a great exercise in structuring a physics simulation (defining a struct for the system state, maybe using StaticArrays for state vectors, and iterating an ODE solver). Depending on interest, students could attempt this after the Newton solver. For now, we’ll illustrate the Newton solver example.
    

"""

# ╔═╡ c0d16ea3-8965-42dd-8c5c-452656a80846
md"""

## Example Package: NewtonSolver.jl (Outline & Code)

  

Let’s outline our Newton solver package and show simplified code for it. This will tie together package generation, module definition, testing, and usage of some concepts above:

```
NewtonSolver/              ← Package directory
├── Project.toml           ← Contains name "NewtonSolver", UUID, deps like StaticArrays (if used), Test (in extras)
├── src/
│   └── NewtonSolver.jl    ← Module definition and implementation
└── test/
    └── runtests.jl        ← Tests for the package
```

Now, key parts of the source and test files:

```julia
# src/NewtonSolver.jl
module NewtonSolver
export newton_solve

\"\"\" Solve f(x) = 0 using Newton's method.
    `f` should map a number or vector to same-shaped output (Residual).
    `J` should give the Jacobian (derivative matrix) at a point.
    `x0` is the initial guess (Number or AbstractVector).
    Returns an approximate root x such that f(x) ≈ 0.
\"\"\"
function newton_solve(f, J, x0; tol::Real=1e-8, maxiter::Int=20)
    x = x0  # can be a scalar or vector
    for iter in 1:maxiter
        # Compute function and Jacobian
        local Fx = f(x)
        local Jx = J(x)
        # Solve Jx * Δ = Fx  (if scalar, this is just Δ = Fx/Jx)
        # In Julia, the "\" operator solves linear systems.
        Δ = Jx \ Fx  
        # Update estimate
        x = x - Δ
        # Check convergence (using norm for vector or abs for scalar)
        if (isa(Fx, Number) ? abs(Fx) : norm(Fx)) < tol
            @info "Newton converged in $iter iterations."
            return x
        end
    end
    @warn "NewtonSolver.newton_solve: reached max iterations, solution may not be converged."
    return x
end

end # module
```

A few notes on this implementation:

- We wrote it to accept both scalar and vector x. The use of Δ = Jx \ Fx works for both cases: if x is a number, Jx and Fx are numbers and “\” does division; if x is a vector of length _n_, Jx should be an n×n matrix and Fx a vector of length _n_, so “\” does a linear solve.
    
- We included a simple convergence check using the norm of the residual f(x). In a real solver, one might also check the size of the update Δ.
    
- We provide an informational log when it converges and a warning if it doesn’t within maxiter.
    
- This function is generic: it doesn’t assume specific types for x. It could work with x as a Vector{Float64}, or an SVector from StaticArrays, or even a scalar (in which case Jx and Fx are just Floats). Julia’s multiple dispatch and generic math make it so one method can handle all those cases.
    

  

Now, the testing file:

```julia
# test/runtests.jl
using NewtonSolver
using Test

# 1) Test scalar case (simple root finding for sqrt(2)):
@testset "Scalar Newton solver" begin
    # f(x) = x^2 - 2, root at sqrt(2)
    let 
        f(x) = x^2 - 2
        J(x) = 2x           # derivative of f
        root = newton_solve(f, J, 1.0)   # start from initial guess 1.0
        @test abs(root - sqrt(2)) < 1e-6
    end

    # Edge case: if initial guess is very far or zero derivative, etc.
    # (For brevity, not fully testing those here.)
end

# 2) Test 2D vector case (circle and line intersection):
@testset "Vector Newton solver" begin
    import StaticArrays  # use static vectors for convenience
    let 
        # Use SVector for x, and SMatrix for Jacobian for performance
        using StaticArrays
        function f_vec(u)
            x, y = u  # u is a tuple or static vector of (x,y)
            return SVector{2}(x^2 + y^2 - 1,
                               x - y)
        end
        function J_vec(u)
            x, y = u
            return SMatrix{2,2}(2x, 2y,
                                 1,   -1)
        end
        x0 = @SVector [0.5, 0.5]              # initial guess
        sol = newton_solve(f_vec, J_vec, x0)
        # The solution should satisfy x ≈ y and x^2 + y^2 ≈ 1 (so x ≈ ±0.707)
        @test abs(sol[1] - sol[2]) < 1e-8      # x ≈ y
        @test abs(sol[1]^2 + sol[2]^2 - 1) < 1e-8  # on the unit circle
    end
end
```

In the tests, we tried a scalar example and a 2D vector example:

- The scalar test checks that newton_solve finds approximately √2 for the equation x²–2=0.
    
- The vector test uses StaticArrays (SVector and SMatrix) to solve the two-equation system described. It verifies the solution properties instead of comparing to a known exact solution directly.
    

  

Students can run these tests with ] test NewtonSolver in the package environment. If all is well, they’ll see “Test Summary: | Pass …”.

  

This project demonstrates how to set up a Julia package, write a generic algorithm, and test it. It’s a foundation they can extend – for example:

- Add a fallback if no Jacobian is provided (perhaps use finite differences or ForwardDiff.jl to compute J).
    
- Improve convergence criteria or add a maximum step size.
    
- Package could be expanded to include other methods (secant method, etc.) or made into a small library.
    

  

Encourage students to experiment: break the code (e.g., give a bad initial guess) and use the debugger or printouts to see how it behaves, or try the alternative project (mass-spring ODE) for a different challenge. Happy coding!
"""

# ╔═╡ 2070d6c6-a89a-40eb-b39a-efdbea46402f
md"""
# Working with Julia Packages: Part 2 – Intermediate and Advanced Topics

This follow-up to the beginner-focused Pluto slides explores **intermediate and advanced topics** in Julia package development. These are especially useful as your packages grow or become shared with others. Topics include compatibility handling, documentation, CI, extensibility patterns, and performance considerations.

---

## 📦 1. Package Compatibility and `[compat]` Section

The `[compat]` section in `Project.toml` defines the version bounds for Julia and your dependencies. It ensures that your package won’t break if a dependency introduces breaking changes.

```toml
[compat]
julia = "1.10"
StaticArrays = "1"
```

**Tips:**

* Be specific but not overly restrictive.
* Use `Pkg.update()` and `Pkg.resolve()` to manage versions.
* Tools like `CompatHelper.jl` (GitHub bot) can help automate updates.

---

## 📚 2. Documentation with Documenter.jl

**Documenter.jl** builds beautiful documentation hosted on GitHub Pages.

### Basic setup:

1. Add Documenter to your `[extras]` and `[targets]`.
2. Create `docs/make.jl`:

   ```julia
   using Documenter, MyPackage

   makedocs(
       modules = [MyPackage],
       sitename = "MyPackage.jl",
       format = Documenter.HTML()
   )
   ```
3. Build docs locally with `julia --project=docs docs/make.jl`.
4. Setup GitHub Action to deploy on push.

**Write docstrings:**

```julia
\"\"\"
    newton_solve(f, J, x0; tol, maxiter)

Finds a root of `f(x)=0` using Newton's method.
\"\"\"
```

Use triple-quoted markdown for docstrings. Documenter can parse them into HTML.

---

## 🧩 3. Extensible APIs and Multiple Dispatch

Design your package so others can extend it. Some tips:

* Use abstract types and traits.
* Separate interface from implementation.

```julia
abstract type AbstractSolver end
struct NewtonSolver <: AbstractSolver end

solve(::NewtonSolver, f, J, x0) = ...
```

This lets users define new solvers by subtyping `AbstractSolver` and extending `solve`.

---

## 🧰 4. Advanced Debugging Tools

Besides VS Code:

* `Infiltrator.jl` – use `@infiltrate` to pause and inspect state.
* `Debugger.jl` – REPL step-through debugging (`@enter`, `@run`, etc.).

These can help debug deeper issues where print-debugging is insufficient.

---

## ✨ 5. Formatting, Linting, and Style

Use tools to maintain code quality:
- `JuliaFormatter.jl` to auto-format code:
  ```julia
  using JuliaFormatter
  format("src/"; indent=4)
````

* `JET.jl`, `Lint.jl` for catching bugs or anti-patterns.
* Use `aqua` for testing project hygiene.

---

## ✅ Summary

| Topic             | Purpose                                |
| ----------------- | -------------------------------------- |
| `[compat]`        | Dependency/version safety              |
| Documenter.jl     | Automatic docs with testing            |
| Traits/API design | Extensible package architecture        |
| StaticArrays      | Performance for small vectors          |
| Doctests          | Testable, live documentation           |
| Formatting        | Style consistency & error prevention   |

This concludes **Part 2** of our Julia package development guide. These tools and practices help you write high-quality, maintainable, and shareable packages.

Happy hacking! 🚀

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
# ╟─f09aec17-c235-4c1c-bf32-c9d9cd345feb
# ╟─58a2b3d4-2990-4599-b148-0432826b801b
# ╟─f5887bab-576a-4991-a077-58413e7f3a5f
# ╟─833bfc7b-b5dd-44ea-8d18-c98c3f53525d
# ╟─5481a651-d77b-4966-911f-3ff2bbdf050a
# ╟─c5aa9dd2-aacd-4e39-aceb-693679b03ed8
# ╟─9fb6a51c-f183-49f4-b584-083d2391c631
# ╟─a5560aa7-e985-4504-9534-8c71fae525c4
# ╟─16482ac0-83ec-4253-ac08-2c3a182ba86a
# ╟─99b50fd7-bf78-4e8b-b5f8-ce8d7291fc87
# ╟─cbf450ac-911f-4ec2-82fc-a47237dec16b
# ╟─c0d16ea3-8965-42dd-8c5c-452656a80846
# ╠═2070d6c6-a89a-40eb-b39a-efdbea46402f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
