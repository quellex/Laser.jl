# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
# vector potential
function getAt(t::TF, pulse::Pulse{At, TF}) where TF
	@unpack tL, tR, T, A, ω, ϕ = pulse
	if t < tL || tR < t
		return 0
	end
	return A*get_env(t, pulse)*sin(ω*t + ϕ)
end
function getAt(t::TVF, pulse::Pulse{At}) where {TVF <: AbstractVector{<:AbstractFloat}}
	At = zeros(eltype(t), size(t))
	for it in 1:length(t)
		At[it] = getAt(t[it], pulse)
	end
	return At
end
function getAt(t::TF, laser::Vector{Pulse{At, TF}}) where {TF <: AbstractFloat}
	At = zero(t)
	for ip in 1:length(laser)
		At += getAt(t, laser[ip])
	end
	return At
end
function getAt(t::TVF, laser::Vector{Pulse{At}}) where {TVF <: AbstractVector{<:AbstractFloat}}
	At = zeros(eltype(t), size(t))
	for it in 1:length(t)
		At[it] = getAt(t[it], laser)
	end
	return At
end
# electric field
# Et = -dA/dt
function getEt(t::TF, pulse::Pulse{At})where TF <: AbstractFloat
	@unpack tL, tR, T = pulse
	if t < tL || tR < t
		return 0
	end
	dt = T*1e-5
	Et =  (1/60dt)*(getAt(t + 3dt, pulse) - getAt(t - 3dt, pulse)) -
			(3/20dt)*(getAt(t + 2dt, pulse) - getAt(t - 2dt, pulse)) +
			(3/4dt)*(getAt(t + dt, pulse) - getAt(t - dt, pulse))
	return -Et
end
function getEt(t::TVF, pulse::Pulse{At}) where {TVF <: AbstractVector{<:AbstractFloat}}
	Et = zeros(eltype(t), size(t))
	for it in 1:length(t)
		Et[it] = getEt(t[it], pulse)
	end
	return Et
end
function getEt(t::TF, laser::Vector{Pulse{At}}) where {TF <: AbstractFloat}
	Et = zero(t)
	for ip in 1:length(laser)
		Et += getEt(t, laser[ip])
	end
	return Et
end
function getEt(t::TVF, laser::Vector{Pulse{At}}) where {TVF <: AbstractVector{<:AbstractFloat}}
	Et = zeros(eltype(t), size(t))
	for it in 1:length(t)
		Et[it] = getEt(t[it], laser)
	end
	return Et
end
