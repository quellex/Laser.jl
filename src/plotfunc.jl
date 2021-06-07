# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
function plot_Et!(p, t, pulse::Pulse; show_env=false)
	Et = getEt(t, pulse)
	plot!(p, t, Et, label="")
	if show_env
		envelope = pulse.E * get_env(t, pulse)
		plot!(p, t, envelope, label="")
	end
	return p
end
function plot_Et(t, pulse::Pulse; show_env=false)
	p = plot(xlabel="time[a.u.]", ylabel="electric field [a.u.]")
	plot_Et!(p, t, pulse, show_env=show_env)
	return p
end
function plot_Et(t, laser::Vector{Pulse}; show_all=true, show_env=false)
	p = plot(xlabel="time[a.u.]", ylabel="electric field [a.u.]")
	if show_all
		for pulse in laser
			plot_Et!(p, t, pulse, show_env=show_env)
		end
	end
	Et = getEt(t, laser)
	plot!(p, t, Et, label="total")
	return p
end
function plot_At!(p, t, pulse::Pulse)
	At = getAt(t, pulse)
	plot!(p, t, At, label="")
	return p
end
function plot_At(t, pulse::Pulse)
	p = plot(xlabel="time[a.u.]", ylabel="vector potential [a.u.]")
	plot_At!(p, t, pulse)
	return p
end
function plot_At(t, laser::Vector{Pulse}; show_all=true)
	p = plot(xlabel="time[a.u.]", ylabel="vector potential [a.u.]")
	if show_all
		for pulse in laser
			plot_At!(p, t, pulse)
		end
	end
	At = getAt(t, laser)
	plot!(p, t, At, label="total")
	return p
end
