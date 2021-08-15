# %%
using Plots
# %%
include("../src/Laser.jl")
# %%
pulseEt = Laser.Pulse(1, 800, 0, :sin4, 10, :Et)
pulseAt = Laser.Pulse(1, 800, 0, :sin4, 10, :At)
pulse = pulseEt
time = Laser.Time(pulse, 1000, 10)
# %%
Et = Laser.getEt(time.t, pulseAt)
# %%
At = Laser.getAt(time.t, pulseAt)
# %%
p1 = plot(time.t, Et, label="Et")
# %%
p2 = plot(time.t, At, label="At")
# %%
plot(p1, p2, layout=(2, 1))
# %%
pulseEt.tprev[1]
