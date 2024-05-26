using Pkg
Pkg.activate(".")
using SparseArrays
using BandedMatrices


#Mixture Properties

using GCIdentifier
using Clapeyron
using Distributions, Plots
gr()
using LaTeXStrings
include("set_size.jl")
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=1.5, 
    framestyle=:box, 
    label=nothing, 
    grid=false,
    titlefontsize = 10
)


using PyCall
using Conda
run(`$(PyCall.python) -m pip install rdkit-pypi`)
Chem = pyimport("rdkit.Chem")
mol_descritprs = pyimport("rdkit.Chem.rdMolDescriptors")
mol_loops = pyimport("rdkit.Chem.rdmolops")
mol_weight = pyimport("rdkit.Chem.Descriptors")

#-----------

# Ingredient universe selection (it can be more, but I will not go into it)
limonene = get_groups_from_smiles("CC1=CCC(CC1)C(=C)C", UNIFACGroups)
pinene = get_groups_from_smiles("CC1=CCC2CC1C2(C)C", UNIFACGroups)
geraniol = get_groups_from_smiles("CC(=CCCC(=CCO)C)C", UNIFACGroups)
linalool = get_groups_from_smiles("CC(=CCCC(C)(C=C)O)C", UNIFACGroups)
l_acetate = get_groups_from_smiles("CC(=CCCC(C)(C=C)OC(=O)C)C", UNIFACGroups)
vanilin = get_groups_from_smiles("COC1=C(C=CC(=C1)C=O)O", UNIFACGroups)
galaxolide = get_groups_from_smiles("CC1COCC2=CC3=C(C=C12)C(C(C3(C)C)C)(C)C", UNIFACGroups)
tonalide = get_groups_from_smiles("CC1CC(C2=C(C1(C)C)C=C(C(=C2)C(=O)C)C)(C)C", UNIFACGroups)
etanol =  get_groups_from_smiles("CCO", UNIFACGroups)
water = get_groups_from_smiles("O", UNIFACGroups)


all_smiles = ["CC1=CCC(CC1)C(=C)C", "CC1=CCC2CC1C2(C)C", "CC(=CCCC(=CCO)C)C","CC(=CCCC(C)(C=C)O)C",
"CC(=CCCC(C)(C=C)OC(=O)C)C", "COC1=C(C=CC(=C1)C=O)O", "CC1COCC2=CC3=C(C=C12)C(C(C3(C)C)C)(C)C",
"CC1CC(C2=C(C1(C)C)C=C(C(=C2)C(=O)C)C)(C)C", "CCO", "O"]
#


map_idxs = Dict(1 => 1, 2 => 4, 3 => 2, 4 => 3, 5 => 9, 6 => 5, 7 => 6)

#------ utility
function composition(molecule)
    # Check that there is a valid molecule
    if !isnothing(molecule)

        # Add hydrogen atoms--RDKit excludes them by default
        molecule_with_Hs = mol_loops.AddHs(molecule)
        comp = Dict()

        # Get atom counts
        for atom in molecule_with_Hs.GetAtoms()
            atomic_num = atom.GetAtomicNum()
            if haskey(comp, atomic_num)
                comp[atomic_num] += 1
            else
                comp[atomic_num] = 1
            end
        end

#=         # If charged, add charge as "atomic number" 0
        charge = getFormalCharge(molecule_with_Hs)
        if charge != 0
            comp[0] = charge
        end =#

        return comp
    end
end

vv = composition(Chem.MolFromSmiles("CC1=CCC(CC1)C(=C)C")) #Testing function 2
vv = composition(Chem.MolFromSmiles("CC1=CCC2CC1C2(C)C")) #Testing function 2

#---- end of utility




function Di_air(component_smile; P = 1.0 , T = 296.0, DiffusionVolumes = Dict(6 => 16.5, 1 => 1.98, 8 => 5.48, 7 => 5.69, 17 => 19.5, 16 => 17.0))
    #Dict(6 => 16.5, 1 => 1.98, 8 => 5.48, 7 => 5.69, 17 => 19.5, 16 => 17.0)
    #Dict(6 => 15.9, 1 => 2.31, 8 => 6.11, 7 => 5.69, 17 => 19.5, 16 => 17.0)

    #Denominator
    AromaticRingsNum = mol_descritprs.CalcNumAromaticRings(Chem.MolFromSmiles(component_smile))
    HeteroRingsNum = mol_descritprs.CalcNumHeterocycles(Chem.MolFromSmiles(component_smile))
    Composition = composition(Chem.MolFromSmiles(component_smile))

    if component_smile == "O"
        AVolume = 12.7
    else
        AVolume = 0.0 + (AromaticRingsNum + HeteroRingsNum)*-18.30
    end

    if component_smile == "O"
        nothing
    else
    for key in keys(Composition)
        if haskey(DiffusionVolumes, key)
           AVolume += Composition[key]*DiffusionVolumes[key]
        else
            AVolume += Composition[key]*0.0
        end 

    end
end

    BVolume = 20.1 #Air

    AMolWeight = mol_weight.ExactMolWt(Chem.MolFromSmiles(component_smile))
    BMolWeight = 28.960

    Dₐᵦ = 1.00*10^-3*(T^1.75)*sqrt(1/AMolWeight + 1/BMolWeight)/(P*(AVolume^(1/3) + BVolume^(1/3))^2)*10^-4 #m2/s
    #Dₐᵦ = (1.43*10^-7)*(T^(1.5))*(1/AMolWeight + 1/BMolWeight)/(P*2.0*(AVolume^(1/3) + BVolume^(1/3))^2)*3600 #m2/h
    return Dₐᵦ, AMolWeight
end

Dair = [Di_air(i, T = 298.15, P = 1.0)[1] for i in all_smiles]
MolMass =  [Di_air(i, T = 298.15, P = 1.0)[2] for i in all_smiles]


#Vapor pressure Antoine coefficients
As = [6.81591, 8.64144, 10.9356, 7.23914757036692]
Bs = [2075.62, 3287.427, 4535.023, 1582.19525234045]
Cs = [-16.65, 0.0, 0.0, -48.854769331095]


#Activity coefficient
mixture = UNIFAC([limonene, geraniol, vanilin, etanol]; puremodel = BasicIdeal)
positions = [1, 3 , 6, 9]

#Properties (limonene, geraniol, vanilin, etanol)
D_quaternary = Dair[positions] #Diffusivity
PVap_quaternary = @. 10^(As - Bs/(298.15 + Cs))*10^3 #Vapor pressure
C_T = 101325.0/(Clapeyron.R̄ * 298.15) #mol/m3
A_lg = 0.071 #m^2
L = 2.0 #m


#--------- PDE setting
using ForwardDiff, OrdinaryDiffEq, UnPack, LinearAlgebra


#stencil
n = 40
xmax = L
h = xmax/n
domain = 0.0:h:n*h |> collect

∂² = BandedMatrix{Float64}(undef, (n + 1, n + 1), (1,1))
∂²[band(0)] .= -2
∂²[band(1)] .= 1
∂²[band(-1)] .= 1
∂²
∂² = ∂²[2:end - 1, 2:end - 1]
∂²

#Initial condition
x₀ = [0.12, 0.12, 0.06, 0.7] #Limonene, geraniol, vanillin, ethanol
N₀ = 1e-3 #mmol
u0 = zeros(size(∂², 2) + 1, 4) .+ 0.0
n₀x₀ = N₀*x₀ 
u0[end, :] = n₀x₀

cachey0 = [1.0; zeros(size(∂², 2) - 1)]

# Sampling Diffusivity and Antoine´s Equation
σ_Dair = 8.429798847332153e-7
Tsats_σ_eth = 0.1468180958847352
Tsats_σ = [1.0, 5.0, 5.0, Tsats_σ_eth]

sigmas = 4.0

# Saturation temperature pdfs
Tsat_pdfs = [Truncated(Normal(0.0, 1.0), -sigmas, sigmas), Truncated(Normal(0.0, 1.0), -sigmas, sigmas),
Truncated(Normal(0.0, 1.0), -sigmas, sigmas), Truncated(Normal(0.0, 1.0), -sigmas, sigmas)]

#Diffusivity pdfs
Dair_pdfs = [Truncated(Normal(0.0, 1.0), -sigmas, sigmas), Truncated(Normal(0.0, 1.0), -sigmas, sigmas),
Truncated(Normal(0.0, 1.0), -sigmas, sigmas), Truncated(Normal(0.0, 1.0), -sigmas, sigmas)]



using QuasiMonteCarlo, Distributions, Random
particles = 2_000
N = particles
Random.seed!(1234)
lb = zeros(8)
ub = ones(8)
s1 = QuasiMonteCarlo.sample(N, lb, ub, SobolSample())
tosample_d = [Dair_pdfs; Tsat_pdfs]
SamplesQuasi = map(x -> quantile.(tosample_d, s1[:, x]), 1:particles)
scatter(SamplesQuasi[1, :], SamplesQuasi[2, :])


#Cache for solutions
mixture_size = 4
variables = 3 + 1 # 3 positions + 1 liquid phase 
all_variables = mixture_size*variables
saveats = 100.0
#saveats = 1.0
tmax = 3600.0*5
#tmax = 60.0*10.0
Chains = ones(all_variables, Int(tmax/saveats + 1), particles);
idx_50 = argmin(abs.(domain[2:end - 1] .- 0.20))
idx_100 = argmin(abs.(domain[2:end - 1] .- 1.0))

function save_at_chains(ODEsolutionArray, i)
    # z = 0.05 cm
    c_1_entry =  (@view ODEsolutionArray[1, 1, :])*C_T
    c_2_entry =  (@view ODEsolutionArray[1, 2, :])*C_T
    c_3_entry =  (@view ODEsolutionArray[1, 3, :])*C_T
    c_4_entry =  (@view ODEsolutionArray[1, 4, :])*C_T

    Chains[1, :, i] = copy(c_1_entry)
    Chains[2, :, i] = copy(c_2_entry)
    Chains[3, :, i] = copy(c_3_entry)
    Chains[4, :, i] = copy(c_4_entry)

    # z ~ 50 cm
    c_1_mid =  (@view ODEsolutionArray[idx_50, 1, :])*C_T
    c_2_mid =  (@view ODEsolutionArray[idx_50, 2, :])*C_T
    c_3_mid =  (@view ODEsolutionArray[idx_50, 3, :])*C_T
    c_4_mid =  (@view ODEsolutionArray[idx_50, 4, :])*C_T

    Chains[5, :, i] = copy(c_1_mid)
    Chains[6, :, i] = copy(c_2_mid)
    Chains[7, :, i] = copy(c_3_mid)
    Chains[8, :, i] = copy(c_4_mid)

    #z ~ 100 cm
    c_1_end =  (@view ODEsolutionArray[idx_100, 1, :])*C_T
    c_2_end =  (@view ODEsolutionArray[idx_100, 2, :])*C_T
    c_3_end =  (@view ODEsolutionArray[idx_100, 3, :])*C_T
    c_4_end =  (@view ODEsolutionArray[idx_100, 4, :])*C_T

    Chains[9, :, i] = copy(c_1_end)
    Chains[10, :, i] = copy(c_2_end)
    Chains[11, :, i] = copy(c_3_end)
    Chains[12, :, i] = copy(c_4_end)


    #xᵢ 
    x₁ = (@view ODEsolutionArray[end, 1, :])./
    ( (@view ODEsolutionArray[end, 1, :]) + 
     (@view ODEsolutionArray[end, 2, :]) +
     (@view ODEsolutionArray[end, 3, :]) +
     (@view ODEsolutionArray[end, 4, :]))

    x₂ = (@view ODEsolutionArray[end, 2, :])./
     ( (@view ODEsolutionArray[end, 1, :]) + 
      (@view ODEsolutionArray[end, 2, :]) +
      (@view ODEsolutionArray[end, 3, :]) +
      (@view ODEsolutionArray[end, 4, :]) )

    x₃ = (@view ODEsolutionArray[end, 3, :])./
      ( (@view ODEsolutionArray[end, 1, :]) + 
       (@view ODEsolutionArray[end, 2, :]) +
       (@view ODEsolutionArray[end, 3, :]) +
       (@view ODEsolutionArray[end, 4, :]) )

    x₄ = (@view ODEsolutionArray[end, 4, :])./
       ( (@view ODEsolutionArray[end, 1, :]) + 
        (@view ODEsolutionArray[end, 2, :]) +
        (@view ODEsolutionArray[end, 3, :]) +
        (@view ODEsolutionArray[end, 4, :]) ) 
    
    Chains[13, :, i] = copy(x₁)
    Chains[14, :, i] = copy(x₂)
    Chains[15, :, i] = copy(x₃)
    Chains[16, :, i] = copy(x₄)
    

end

struct ScentDiffusionModel{T1, T2, T3, T4}
    ∂²::T1
    mixture::T2
    Diff::T3
    PVap::T4
end
    
function(f::ScentDiffusionModel)(du, u, p, t)

@unpack ∂², mixture, Diff, PVap = f

        y0 = activity_coefficient(mixture, 1.0, 298.15, [u[end, 1], u[end, 2], u[end, 3], u[end, 4]]).*10.0.^(As .- Bs./(Cs .+ (p[5:8].*Tsats_σ .+ 298.15)))*10^3/101325.0.*
        u[end, :]./(u[end, 1] + u[end, 2] + u[end, 3] + u[end, 4])

        for i in 1:4 
        # Gas phase balance
        du[1:end - 1, i] = (p[i]*σ_Dair + Diff[i]) * (∂²*u[1:end - 1, i] + y0[i]*cachey0)/h^2 #Last parts account for boundary value
        
        # Liquid phase balance
        du[end , i] = (p[i]*σ_Dair + Diff[i]) * A_lg * C_T * (-u[2, i] + 4*u[1, i] - 3*y0[i])/(2*h)
        #du[end , i] = Diff[i] * A_lg * C_T * (u[1, i] - y0[i])/h
        
        end
end

#Testing solution
rhs_k = ScentDiffusionModel(∂², mixture, D_quaternary, PVap_quaternary)
tspan_k = (0.0, 3600.0*5)
#tspan_k = (0.0, 60.0*10)
p = zeros(8)
problem_k = ODEProblem(rhs_k, u0, tspan_k, p)
@time solution_k = solve(problem_k, FBDF(autodiff = true),
 saveat = 100., abstol = 1e-8, reltol = 1e-8);

solarray_k = Array(solution_k)
plot(solution_k.t, solarray_k[idx_50, 4, :]*C_T)

samples = SamplesQuasi

#Every time you run this loop, you should initialize Chains again
for i in 1:size(samples, 1)

    #QuasiMonteCarlo samples
    p = [samples[i][1],samples[i][2],samples[i][3],samples[i][4],  
        samples[i][5], samples[i][6], samples[i][7], samples[i][8]]
 
    new_prob = remake(problem_k; p = p)
    @time solution = Array(solve(new_prob, FBDF(autodiff = true),
     abstol = 1e-8, reltol = 1e-8, saveat = saveats))
    #Saving at chains cache
    println(i)
    save_at_chains(solution, i)
end

using StatsPlots


#Change the index of Chains array if you want to visualize from another distance of the releasing source
#save_at_chains function has legends for the meaning of each index of the Chains array. 
eth_100 = errorline(0.0:100.0:3600.0*5, Chains[4, :, :], color = "red", errorstyle=:ribbon, errortype = :percentile, 
label=L"\mu", groupcolor = nothing, percentiles = [0.05, 99.5], secondarylinealpha = 0.7, xlabel = "t/s", ylabel = L"\textrm{C}/\textrm{mol}\cdot\textrm{m}^{-3}",
 title = "Ethanol (z = 0.05 m)", size = set_size("thesis"; fraction = 0.8).*72.0)

lim_100 = plot([],[])
 errorline!(lim_100, 0.0:100.0:3600.0*5, Chains[1, :, :], errorstyle=:ribbon, errortype = :percentile, color = "green",
label=L"\mu", groupcolor = nothing, percentiles = [0.05, 99.5], secondarylinealpha = 0.7, xlabel = "t/s", ylabel = L"\textrm{C}/\textrm{mol}\cdot\textrm{m}^{-3}",
 title = "Limonene (z = 0.05 m)", size = set_size("thesis"; fraction = 0.9).*72.0)

ger_100 = plot([],[])
plot!(ger_100, [],[])
errorline!(ger_100, 0.0:100.0:3600.0*5, Chains[2, :, :], errorstyle=:ribbon, errortype = :percentile, color = "green",
label=L"\mu", groupcolor = nothing, percentiles = [0.05, 99.5], secondarylinealpha = 0.7, xlabel = "t/s", ylabel = L"\textrm{C}/\textrm{mol}\cdot\textrm{m}^{-3}",
 title = "Geraniol (z = 0.05 m)", size = set_size("thesis"; fraction = 0.9).*72.0)


van_100 = plot([],[])
plot!(van_100, [],[])
plot!(van_100, [],[])
errorline!(van_100, 0.0:100.0:3600.0*5, Chains[3, :, :], errorstyle=:ribbon, errortype = :percentile, color = "green",
label=L"\mu", groupcolor = nothing, percentiles = [0.05, 99.5], secondarylinealpha = 0.7, xlabel = "t/s", ylabel = L"\textrm{C}/\textrm{mol}\cdot\textrm{m}^{-3}",
 title = "Vanillin (z = 0.05 m)", size = set_size("thesis"; fraction = 0.9).*72.0)

all_100 = plot(lim_100, ger_100, van_100, eth_100, layout=(2, 2),
 formatter = :scientific,
 size = (382*1.8, 236*1.8))