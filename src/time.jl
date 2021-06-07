# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
struct Time{TF <: AbstractFloat,TVF <: AbstractVector{TF},TI <: Integer}
	nstep_cyc::TI
	dt::TF   # real time propagation
	dh::Complex{TF}   # imaginary time-propagation
	nstep_tot::TI # num. of time step
	tcyc::TVF # time step in optical cycle of the first laser
	t::TVF # time step in atomic unit
end
function Time(λ, nstep_cyc, ncyc_tot)
	T = 2π / wlen2au(λ)
	nstep_cyc = nstep_cyc
	dt = T / nstep_cyc
	dh = dt / im
	ncyc_tot = ncyc_tot
	nstep_tot = nstep_cyc * ncyc_tot
	tcyc = collect(range(0., ncyc_tot, length=nstep_tot))
	t = tcyc * T
	return Time(nstep_cyc, dt, dh, nstep_tot, tcyc, t)
end
