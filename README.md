# PerfumeEngineering
This repository contains the necessary codes to reproduce the results of the manuscript https://chemrxiv.org/engage/chemrxiv/article-details/656703395bc9fcb5c9a3a476

```markdown
# Scent Diffusion Model

This repository contains code for modeling the diffusion of various scent molecules in air using Julia.

## Requirements

The following Julia packages are required:

- `SparseArrays`
- `BandedMatrices`
- `GCIdentifier`
- `Clapeyron`
- `Distributions`
- `Plots`
- `LaTeXStrings`
- `PyCall`
- `Conda`
- `ForwardDiff`
- `OrdinaryDiffEq`
- `UnPack`
- `LinearAlgebra`
- `QuasiMonteCarlo`
- `StatsPlots`

Additionally, the RDKit library is required and can be installed via pip:

```bash
pip install rdkit-pypi
```

## Installation

1. Clone the repository and navigate into it:

    ```bash
    git clone https://github.com/yourusername/scent-diffusion-model.git
    cd scent-diffusion-model
    ```

2. Activate the Julia environment and install dependencies:

    ```julia
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    ```

3. Install RDKit via PyCall:

    ```julia
    using PyCall
    run(`$(PyCall.python) -m pip install rdkit-pypi`)
    ```

## Usage

The main script sets up the model, selects ingredients, defines utility functions, and performs the diffusion simulation. Run the script using Julia:

```julia
include("main_script.jl")
```

For more details on the code, refer to the inline comments and function definitions within the script.

## License

This project is licensed under the MIT License.
```
