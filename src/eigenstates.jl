

function groundstate(h::Hamiltonian; maxiter::Int=100000)
	isconstant(h) || throw(ArgumentError("groundstate only defined for constant Hamiltonian"))
	d = size(h, 1)
	eigval, eigvec, info = eigsolve(h, randn(eltype(h), d), 1, :SR; ishermitian=true, maxiter=maxiter)	
	(info.converged>=1) || error("eigsolve fails to converge")
	return eigval[1], eigvec[1]
end

function steadystate(h::Lindbladian; maxiter::Int=100000)
	isconstant(h) || throw(ArgumentError("steadystate only defined for constant Lindbladian"))
	d = div(size(h, 1), 2)
	eigval, eigvec, info = eigsolve(h, randn(eltype(h), d, d), 1, EigSorter(abs; rev = false); ishermitian=false, maxiter=maxiter)
	(info.converged>=1) || error("eigsolve fails to converge")
	rho = eigvec[1]
	rho ./= tr(rho)
	return eigval[1], rho
end