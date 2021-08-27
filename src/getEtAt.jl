# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
## envelope
function get_env(t::TF, pulse::Pulse{TL,TF})::TF where {TL,TF}
    @unpack tL, tR, T, env_type, ncyc = pulse
    if t <= tL || tR <= t
        return 0
    end
    if env_type == :sin2
        return sin((π * t) / (ncyc * T))^2
    elseif env_type == :sin4
        return sin((π * t) / (ncyc * T))^4
    else
        env_type == :flat || env_type == :delta
        return 1
    end
end
function get_env(t::TVF, pulse::Pulse)::TVF where {TVF <: AbstractVector{<:AbstractFloat}}
    env = zeros(eltype(t), size(t))
    for it in eachindex(t)
        env[it] = get_env(t[it], pulse)
    end
    return env
end
## electric field
function getEtAt(t::TF, pulse::Pulse{TL,TF})::Tuple{TF,TF} where {TL,TF}
    return (getEt(t, pulse), getAt(t, pulse))
end
function getEt(t::TF, pulse::Pulse{Et,TF})::TF where {TF}
    @unpack E, ω, ϕ = pulse
    return E * get_env(t, pulse) * sin(ω * t + ϕ)
end
function getEt(t::TF, pulse::Pulse{At,TF})::TF where {TF}
    @unpack tL, tR, T = pulse
    if t <= tL || tR <= t
        return zero(TF)
    end
    dt = T * 1e-5
    Et =
        (1 / 60dt) * (getAt(t + 3dt, pulse) - getAt(t - 3dt, pulse)) -
        (3 / 20dt) * (getAt(t + 2dt, pulse) - getAt(t - 2dt, pulse)) +
        (3 / 4dt) * (getAt(t + dt, pulse) - getAt(t - dt, pulse))
    return -Et
end
function getEt(t::TVF, pulse::Pulse)::TVF where {TVF <: AbstractVector{<:AbstractFloat}}
    Et = zeros(eltype(t), length(t))
    for it in eachindex(t)
        Et[it] = getEt(t[it], pulse)
    end
    return Et
end
## vector potential
function getAt(t::TF, pulse::Pulse{Et,TF})::TF where {TF}
    @unpack tL, tR, T = pulse
    if t <= tL || tR <= t
        return 0
    end
    nt = 1
    tprev = pulse.tprev[1]
    At = pulse.Atprev[1]
    if tprev > t
        tprev = 0
        At = 0
    end
    while (t - tprev) / nt > T / 100000
        nt *= 10
    end
    dt = (t - tprev) / nt
    for it = 1:nt
        t0 = tprev + (it - 1) * dt
        t1 = t0 + dt / 2
        t2 = t0 + dt
        At += (-dt / 6) * (getEt(t0, pulse) + 4getEt(t1, pulse) + getEt(t2, pulse))
    end
    pulse.tprev[1] = t
    pulse.Atprev[1] = At
    return At
end
function getAt(t::TF, pulse::Pulse{At,TF})::TF where {TF}
    @unpack tL, tR, T, A, ω, ϕ = pulse
    if t <= tL || tR <= t
        return 0
    end
    return A * get_env(t, pulse) * sin(ω * t + ϕ)
end
function getAt(t::TVF, pulse::Pulse)::TVF where {TVF <: AbstractVector{<:AbstractFloat}}
    At = zeros(eltype(t), size(t))
    for it = eachindex(t)
        At[it] = getAt(t[it], pulse)
    end
    return At
end
