abstract type TimeEvoAlg end
abstract type RungKuttaAlg <: TimeEvoAlg end

struct RungKutta1 <: RungKuttaAlg end
struct RungKutta2 <: RungKuttaAlg end
struct RungKutta4 <: RungKuttaAlg end
struct KrylovExpm <: TimeEvoAlg 
	maxiter::Int
	tol::Float64
end
KrylovExpm(; maxiter::Int=10000000, tol::Real=1.0e-12) = KrylovExpm(maxiter, convert(Float64, tol))

timeevo(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, t::Real, dt::Real, alg::RungKuttaAlg) = timeevo(copy(state), h, t, dt, alg)

timeevo!(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, t::Real, dt::Real, alg::RungKutta1) = _first_order!(state, h, t, dt)
timeevo!(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, t::Real, dt::Real, alg::RungKutta2) = _rungekutta2order!(state, h, t, dt)
timeevo!(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, t::Real, dt::Real, alg::RungKutta4) = _rungekutta4order!(state, h, t, dt)

timeevo(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, tspan::Tuple{Real, Real}, dt::Real, alg::RungKuttaAlg) = timeevo!(copy(state), h, tspan, dt, alg)
function timeevo!(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, tspan::Tuple{Real, Real}, dt::Real, alg::RungKuttaAlg)
	ti, tf = tspan
	tdiff = tf - ti
	nsteps = round(Int, abs(tdiff / dt)) 	
	if nsteps == 0
		nsteps = 1
	end
	dt = tdiff / nsteps
	for i in 1:nsteps
		tj = ti + (i-1) * dt
		timeevo!(state, h, tj, dt, alg)
	end
	return state
end

# function timeevo(state::VecOrMat, h::Union{Hamiltonian, Lindbladian}, tspan::Tuple{Real, Real}, alg::KrylovExpm)
# 	ti, tf = tspan
# 	tdiff = tf - ti
# 	return timeevo(state, h, tdiff, alg)
# end

function timeevo(state::AbstractVector, h::Hamiltonian, dt::Real, alg::KrylovExpm)
	isconstant(h) || error("only constant Hamiltonian supported")
	tmp, info = exponentiate(h, -im*dt, state; tol=alg.tol, maxiter=alg.maxiter, ishermitian=true)
	(info.converged >= 1) || error("krylov subspace method fails to converge")
	return tmp
end
function timeevo(state::AbstractMatrix, h::Lindbladian, dt::Real, alg::KrylovExpm)
	isconstant(h) || error("only constant Lindbladian supported")
	tmp, info = exponentiate(h, dt, state; tol=alg.tol, maxiter=alg.maxiter, ishermitian=false)
	(info.converged >= 1) || error("krylov subspace method fails to converge")
	return tmp
end

function _first_order!(x, f, t, dt)
	# x .+= f(t + dt/2, x) .* dt
	axpy!(dt, f(t + dt/2, x), x)
end

function _rungekutta2order!(x, f, t, dt)
	k1 = f(t, x)
	k2 = f(t+dt, x + dt*k1)
	h = 0.5*dt
	axpy!(h, k1, x)
	axpy!(h, k2, x)
	# @. x += h * (k1 + k2)
end

# function _rungekutta4order!(x, f, t, dt)
# 	k1 = dt .* f(t, x)
# 	k2 = dt .* f(t+dt/2, x + k1 ./ 2)
# 	k3 = dt .* f(t+dt/2, x + k2 ./ 2)
# 	k4 = dt .* f(t+dt, x + k3)
# 	@. x += (k1/6 + k2/3 + k3/3 + k4/6)
# end

function _rungekutta4order!(x, f, t, dt)
	k1 = rmul!(f(t, x), dt)
	k2 = rmul!(f(t+dt/2, x + k1 / 2), dt)
	k3 = rmul!(f(t+dt/2, x + k2 / 2), dt)
	k4 = rmul!(f(t+dt, x + k3), dt)
	axpy!(1/6, k1, x)
	axpy!(1/3, k2, x)
	axpy!(1/3, k3, x)
	axpy!(1/6, k4, x)
end

