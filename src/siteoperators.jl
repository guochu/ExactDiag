

# spin operators
spin_x() = [0. 1; 1 0]
spin_y() = [0. -im; im 0]
spin_z() = [-1. 0; 0 1]
spin_up() = [0. 0; 0 1]
spin_down() = [1. 0; 0 0]

# spin states 
spin_state_up() = [0., 1]
spin_state_down() = [1., 0]

pauli_matrices() = spin_x(), spin_y(), spin_z()

function spin_half_matrices()
	x, y, z = pauli_matrices()
	sp = real((x + im * y) / 2)
	return Dict("x"=>x, "y"=>y, "z"=>z, "+"=>sp, "-"=>sp', "↑"=>spin_up(), "↓"=>spin_down())
end



### bosonic operators
function bosonaoperator(; d::Int)
	(d <= 1) && error("d must be larger than 1")
	a = zeros(Float64, d, d)
	for i = 1:(d - 1)
		a[i, i+1] = sqrt(i)
	end
	return a
end
bosonadagoperator(; d::Int) = adjoint(bosonaoperator(d=d))
function bosondensityoperator(; d::Int) 
	a = bosonaoperator(d=d)
	return a' * a
end
function bosonoccupationoperator(n::Int; d::Int)
	(0 <= n <= d-1) || throw(BoundsError(0:d-1, n))
	r = zeros(d, d)
	r[n+1, n+1] = 1
	return r
end

function boson_matrices(; d::Int)
	a = bosonaoperator(d=d)
	adag = a'
	n = adag * a
	return Dict("a"=>a, "adag"=>adag, "n"=>n)
end

function fock_state(n::Int; d::Int)
	(0 <= n <= d-1) || throw(BoundsError(0:d-1, n))
	r = zeros(d)
	r[n+1] = 1
	return r
end