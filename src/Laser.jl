# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
# %%
module Laser
# modules
using Parameters
using Plots
using PhysicalConstants
# export
export Pulse, Time, getEt, getAt, get_env
# unit conversion
include("unit.jl")
# Et => give E(t) and A(t) = -\int dt' E(t')
# At => give A(t) and E(t) = -dA/dt
abstract type FormType end
struct Et <: FormType end
struct At <: FormType end
const _formtype_dict = Dict{String, T where T <: FormType}("Et" => Et(), "At" => At())
struct Pulse{TE<:FormType, TF <: AbstractFloat,TS <: AbstractString,TI <: Integer}
	form::TE
	I::TF
	E::TF
	A::TF
	λ::TF
	ω::TF
	T::TF
	ϕ::TF
	env_type::TS
	ncyc::TI
	tL::TF
	t0::TF
	tR::TF
end
function Pulse(fint, wlen, cep, env_type, ncyc_env, form_type = "Et")
	@assert haskey(_formtype_dict, form_type)
	I = float(fint)
	E = i2e(fint)
	λ = float(wlen)
	ω = wlen2au(wlen)
	A = E / ω
	T = 2π / ω
	ϕ = deg2rad(float(cep))
	tL = 0.
	t0 = T * ncyc_env / 2
	tR = T * ncyc_env
	return Pulse(_formtype_dict[form_type], I, E, A, λ, ω, T, ϕ, env_type, ncyc_env, tL, t0, tR)
end
include("getEtAt_E.jl")
include("getEtAt_A.jl")
include("plotfunc.jl")
include("time.jl")
end
# %%
