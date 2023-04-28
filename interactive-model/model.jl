### A Pluto.jl notebook ###
# v0.19.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 21d55edb-3665-4e65-9ad8-b7501f895bff
begin
	import Pkg
	Pkg.activate(Base.current_project())
	Pkg.instantiate()

	using Plots
	using Unitful
	using Measurements
	using PlutoUI
	using LaTeXStrings
	using DataFrames
	using CSV
	
	TableOfContents(title="SEMF Model")
end

# ╔═╡ 810eda6e-7ce6-11ed-2541-230735c5fffe
md"""
# SEMF Model Testing

Implements the SEMF formula and create an interactive model where one can adjust the strength of each coefficient in the equation.
"""

# ╔═╡ de4d1aaa-1751-4dc7-a90e-328020469db9
md"""
## SEMF Formulae
"""

# ╔═╡ 211c5d22-b5a3-4936-bfb6-a159f578acf8
struct Model
	a_V::Number 	# Volume
	a_S::Number 	# Surface
	a_C::Number 	# Coulomb
	a_A::Number 	# Asymmetry
	a_P::Number 	# Parity 
	k_P::Number 	# Parity
end

# ╔═╡ ca7ea201-30b8-4b49-9782-e0b3bd1a2242
begin
	# Model coefficients from select authors
	wapstra = Model(
		14.1u"MeV", 
		13u"MeV", 
		0.595u"MeV", 
		19u"MeV",
		33.5u"MeV",
		-0.75
	)
	rohlf = Model(
		15.75u"MeV",
		17.8u"MeV",
		0.711u"MeV",
		23.7u"MeV",
		11.18u"MeV",
		-0.5
	)
end

# ╔═╡ 3a8ca0f4-ce5c-4d45-94ad-1839b16612cd
begin
	function SEMF(N::Int, Z::Int, model::Model)
		A = N + Z

		# Parity sign
		sgn = (A % 2 == 1) ? 0 : ((Z % 2 == 1) ? -1 : 1)

		# Calculate binding energy
		E_B = (
			model.a_V*A 
			- model.a_S*A^(2/3) 
			- model.a_C*Z*(Z - 1)/A^(1/3) 
			- model.a_A*(N - Z)^2/A 
			+ model.a_P*sgn*A^model.k_P
		)
		return E_B
	end

	function SEMF(; N=nothing, Z=nothing, A=nothing, model::Model)
		# Parse arguments
		if count(isnothing, (N, Z, A)) != 1
			throw(ArgumentError("Exactly two of N, Z, and A should be passed."))
		end

		# Compute missing argument
		if isnothing(N) 
			N = A - Z
		elseif isnothing(Z) 
			Z = A - N
		elseif isnothing(A)
			A = N + Z
		end

		# Check if number of nucleons is valid
		if count(x -> x < 0, (N, Z, A)) > 0
			throw(DomainError("too many nucleons for the given A, or negative numbers passed."))
		end
		
		return SEMF(N, Z, model)
	end

	function calculate_SEMF(N_max::Int, Z_max::Int, model::Model)
		N = 0:N_max
		Z = 0:Z_max
		energy = Matrix(undef, length(Z), length(N))
		
		for N in N, Z in Z
			energy[Z+1, N+1] = SEMF(; N=N, Z=Z, model=model) / (N + Z)
		end
		energy[energy .< 0u"MeV"] .= 0u"MeV"
		
		return N, Z, energy
	end
end

# ╔═╡ c47c6a22-8f9c-468e-8e2f-5ba638204018
md"""
## Wapstra, A. H. (1958) Coefficients
"""

# ╔═╡ ddbc1b67-9bfb-45ea-b484-892ec475fb8a
begin
	N_length = 140
	Z_length = 90
	N, Z, energy = calculate_SEMF(N_length, Z_length, wapstra)
end

# ╔═╡ 22c020c4-83cc-4e35-afda-0da246fedfc3
let
	# For some reason the heatmap doesn't get displayed on the first execution of this cell. Run it again to make it show up.
	p = plot(
		heatmap(N, Z, getproperty.(energy, :val), c=:viridis),
		xlabel = "Number of neutrons " * L"N",
		ylabel = "Number of protons " * L"Z",
		colorbar_title = "Binding energy per nucleon (MeV)",
	)
	contour!(p, N, Z, getproperty.(energy, :val), c=cgrad(:binary, rev=true), levels=6, lw=2)
end

# ╔═╡ 5d0fc8ef-b39e-4c98-a2c0-541384bce152
md"""
## Custom Model w/ Sliders
"""

# ╔═╡ 9e20d133-ee03-42a6-bf00-67b972f2ec77
md"""
Custom model sliders (MeV). From left to right: $a_V$, $a_S$, $a_C$, $a_A$, $a_P$, $k_P$.

$(@bind a_V_slider Slider(0:0.1:20, default=rohlf.a_V.val)) 
$(@bind a_S_slider Slider(10:0.1:20, default=rohlf.a_S.val))  
$(@bind a_C_slider Slider(0:0.001:1, default=rohlf.a_C.val))  
$(@bind a_A_slider Slider(15:0.1:25, default=rohlf.a_A.val))  
$(@bind a_P_slider Slider(0:0.1:40, default=rohlf.a_P.val))  
$(@bind k_P_slider Slider(-1:0.01:0, default=rohlf.k_P))  

"""

# ╔═╡ 05119c1c-578a-4a82-adf7-f2a1cb296160
silder_model = Model(
	((a_V_slider, a_S_slider, a_C_slider, a_A_slider, a_P_slider) .* u"MeV")...,
	k_P_slider
)

# ╔═╡ ae7c6e8a-93c9-43ba-8df1-d7da86825c35
let
	N, Z, energy = calculate_SEMF(140, 90, silder_model)
	
	p = plot(
		heatmap(N, Z, getproperty.(energy, :val), c=:viridis),
		xlabel = "Number of neutrons " * L"N",
		ylabel = "Number of protons " * L"Z",
		colorbar_title = "Binding energy per nucleon (MeV)",
	)
	contour!(p, N, Z, getproperty.(energy, :val), c=cgrad(:binary, rev=true), levels=6, lw=2)
end

# ╔═╡ 3ccffeb5-6482-43f3-b3d8-b0cdbe7b618d
md"""
## Experimental Data

Test out plotting the experimental data from AME2020.
"""

# ╔═╡ f6cedfae-fb0e-4e31-bde2-2f2aece8ebcc
begin
	df = DataFrame(CSV.File(joinpath(@__DIR__, "../data/binding_energy_per_A.csv")))
	df[!, "E"] .= df[!, "E"] .* u"keV"
	df[!, "U(E)"] .= df[!, "U(E)"] .* u"keV"
	df = df[df[!, "known"] .== 1, :]
	df = df[df[!, "Z"] .< Z_length, :]
	df = df[df[!, "N"] .< N_length, :]
end

# ╔═╡ 500f98de-48be-4675-b907-37f64ded2aa1
begin
	energy_exp = Matrix{Any}(missing, Z_length+1, N_length+1)
	for row in eachrow(df)
		energy_exp[row.Z+1, row.N+1] = row.E
	end
	energy_exp = uconvert.(u"MeV", energy_exp)
end

# ╔═╡ 20451e36-afea-43ff-af9f-1a2d095c3440
let
	# Plot the experimental binding energies
	p = plot(
		heatmap(N, Z, energy_exp, c=:viridis),
		xlabel = "Number of neutrons " * L"N",
		ylabel = "Number of protons " * L"Z",
		colorbar_title = "Binding energy per nucleon (MeV)",
	)
	# contour!(p, N, Z, getproperty.(energy, :val), c=cgrad(:binary, rev=true), levels=6, lw=2)
end

# ╔═╡ cbfd4d8e-c499-473e-bec5-48f28d8f7cb1
let
	# Plot the difference in binding energy from the SEMF prediction.
	p = plot(
		heatmap(N, Z, uconvert.(u"keV", energy_exp .- energy), c=:viridis),
		xlabel = "Number of neutrons " * L"N",
		ylabel = "Number of protons " * L"Z",
		colorbar_title = "Binding energy per nucleon (keV)",
	)
	# contour!(p, N, Z, getproperty.(energy, :val), c=cgrad(:binary, rev=true), levels=6, lw=2)
end

# ╔═╡ Cell order:
# ╟─810eda6e-7ce6-11ed-2541-230735c5fffe
# ╠═21d55edb-3665-4e65-9ad8-b7501f895bff
# ╟─de4d1aaa-1751-4dc7-a90e-328020469db9
# ╠═211c5d22-b5a3-4936-bfb6-a159f578acf8
# ╠═ca7ea201-30b8-4b49-9782-e0b3bd1a2242
# ╠═3a8ca0f4-ce5c-4d45-94ad-1839b16612cd
# ╟─c47c6a22-8f9c-468e-8e2f-5ba638204018
# ╠═ddbc1b67-9bfb-45ea-b484-892ec475fb8a
# ╠═22c020c4-83cc-4e35-afda-0da246fedfc3
# ╟─5d0fc8ef-b39e-4c98-a2c0-541384bce152
# ╟─9e20d133-ee03-42a6-bf00-67b972f2ec77
# ╟─05119c1c-578a-4a82-adf7-f2a1cb296160
# ╠═ae7c6e8a-93c9-43ba-8df1-d7da86825c35
# ╟─3ccffeb5-6482-43f3-b3d8-b0cdbe7b618d
# ╠═f6cedfae-fb0e-4e31-bde2-2f2aece8ebcc
# ╠═500f98de-48be-4675-b907-37f64ded2aa1
# ╠═20451e36-afea-43ff-af9f-1a2d095c3440
# ╠═cbfd4d8e-c499-473e-bec5-48f28d8f7cb1
