# Environment Setup Instructions

## Step 1: Install Julia and juliaup

- Download and install the [juliaup](https://github.com/JuliaLang/juliaup) tool as described here: [https://julialang.org/install/](https://julialang.org/install/).

## Step 2: Run Julia

- Open terminal.
- Navigate to main course directory
- Run Julia by typing: `julia`
- Execute command `1+1` and verify results :)

## Step 3: Project Environment

- In julia type `]` (right bracket sign).
    You should change terminal mode to package manager. The command prompt should change from `julia>` to `(@1.11) pkg>` or similar.
- Type `activate .` to activate new project.
- Type `instantiate` to install all project's packages.

## Step 4: Run Pluto

- Press Backspace to return to default julia mode. Command prompt should change to `julia>`.
- Type:
    ```julia
    using Pluto
    Pluto.run()
    ```