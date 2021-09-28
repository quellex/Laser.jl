# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
function get_fluence(time::Time, pulse::Pulse)
	env = get_env(time.t, pulse)
	dt = au2as(time.dt) * 1e-18 # dt in second
	flux = pulse.I * env.^2
	fluence = copy(flux)
	fluence[1] = flux[1]
	fluence[2] = (flux[1] + flux[2])
	for it in 3:length(fluence)
		fluence[it] = fluence[it - 1] + (flux[it - 1] + flux[it]) / 2
	end
	return fluence .* dt
end
