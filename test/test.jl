# %%
using Plots
# %%
include("../src/Laser.jl")
# %%
function dfdt(f, dt)
    n = length(f)
    g = copy(f)
    g[1] = (f[2] - f[1]) / dt
    g[n] = (f[n] - f[n - 1]) / dt
    for i in 2:n - 1
        g[i] = (f[i + 1] - f[i - 1]) / 2dt
    end
    return g
end
 # %%
pulseEt = Laser.Pulse(1, 800, 0, :sin4, 10, :Et)
pulseAt = Laser.Pulse(1, 800, 0, :sin4, 10, :At)
pulse = pulseAt
time = Laser.Time(pulse, 1000, 10)
# %%
Et = Laser.getEt(time.t, pulse)
# %%
At = Laser.getAt(time.t, pulse)
# %%
dAdt = dfdt(At, time.dt)
# %%
p1 = plot(time.t, Et, label="Et", color=:black, linewidth=3.0)
p1 = plot!(time.t, -dAdt, label="-dAdt")
p2 = plot(time.t, At, label="At")
plot(p1, p2, layout=(2, 1))
# %%
