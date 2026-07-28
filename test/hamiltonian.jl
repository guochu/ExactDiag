println("------------------------------------")
println("|            Hamiltonian           |")
println("------------------------------------")

function randstate(L)
	m = randn(ComplexF64, L)
	return normalize!(m)
end

@testset "Hamiltonian definition and mv" begin
	
	d = 5
	x, y, z = pauli_matrices()
	p = boson_matrices(d=d)
	n, adag, a = p["n"], p["adag"], p["a"]

	# rabi Hamiltonian
	m1 = 0.5*kron(y, one(n)) + kron(one(z), n) 
	m2 = kron(x, a+adag)

	H1 = Hamiltonian(m1 + m2)
	H2 = Hamiltonian(m1, [TimeDependentOperator(m2, sin)])


	@test isconstant(H1)
	@test !isconstant(H2)

	@test eltype(H1) == ComplexF64
	@test eltype(H2) == ComplexF64

	@test size(H1, 1) == size(H1, 2) == 2 * d
	@test size(H2, 1) == size(H2, 2) == 2 * d


	@test ishermitian(matrix(H1))
	@test matrix(H1) ≈ matrix(H2(pi/2))

	v = rand(eltype(H1), size(H1, 1))
	@test H1(v) ≈ H2(pi/2)(v)
	@test H1(v) ≈ H2(pi/2, v)

	# time evolution
	tol = 1.0e-12

	psi = randstate(2*d)
	t = 0.2
	# exact evolution
	psi2 = timeevo(psi, H1, t, KrylovExpm())
	@test abs(norm(psi2)-1) < tol

	# 2 order rung kutta
	tol2 = 1.0e-3
	psi3 = timeevo(psi, -im*H2, (0, t), 0.01, RungKutta2())
	@test abs(norm(psi3)-1) < tol2

	psi4 = timeevo(psi, -im*H1, (0, t), 0.01, RungKutta2())
	@test norm(psi2 - psi4) < tol2


	# 4 order rung kutta
	tol2 = 1.0e-6
	psi3 = timeevo(psi, -im*H2, (0, t), 0.01, RungKutta4())
	@test abs(norm(psi3)-1) < tol2

	psi4 = timeevo(psi, -im*H1, (0, t), 0.01, RungKutta4())
	@test norm(psi2 - psi4) < tol2


end

