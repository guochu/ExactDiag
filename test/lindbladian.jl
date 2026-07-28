println("------------------------------------")
println("|            Lindbladian           |")
println("------------------------------------")


function randdm(L)
	m = randn(ComplexF64, L, L)
	m = m' * m
	m ./= tr(m)
	return m
end

@testset "Lindbladian definition and mv" begin

	d = 3
	s = spin_half_matrices()
	x, y, z, sp, sm = s["x"], s["y"], s["z"], s["+"], s["-"]
	p = boson_matrices(d=d)
	n, adag, a = p["n"], p["adag"], p["a"]

	Gamma = 0.7

	# rabi Hamiltonian
	m1 = 0.5*kron(z, one(n)) + kron(one(z), n) 
	m2 = kron(y, a+adag)

	H1 = Hamiltonian(m1+m2)
	H2 = Hamiltonian(m1, [TimeDependentOperator(m2, sin)])

	jumps = [kron(sp, one(n)), Gamma * kron(sm, one(n))]
	L1 = Lindbladian(H1, jumps)
	L2 = Lindbladian(H2, jumps)

	@test size(L1, 1) == size(L1, 2) == 4 * d^2
	@test size(L2, 1) == size(L2, 2) == 4 * d^2

	@test eltype(L1) == ComplexF64
	@test eltype(L2) == ComplexF64

	@test isconstant(L1)
	@test !isconstant(L2)

	rho = randdm(2 * d)

	tol = 1.0e-15
	rho1 = L1(rho)

	@test abs(tr(rho1)) < tol
	@test abs(tr(L2(0.4, rho))) < tol

	L1m = matrix(L1)
	@test size(L1m, 1) == size(L1m, 2) == 4 * d^2

	rho2 = reshape(L1m * reshape(rho, length(rho)), size(rho))
	@test norm(rho1 - rho2) < tol

	@test L1m ≈ matrix(L2(pi/2))
	@test rho1 ≈ L2(pi/2)(rho)

	# time evolution
	tol = 1.0e-12
	tol2 = 1.0e-6

	t = 0.2
	rho2 = timeevo(rho, L1, t, KrylovExpm())
	@test abs(tr(rho2)-1) < tol

	# 2 order rung kutta
	tol2 = 1.0e-3
	rho3 = timeevo(rho, L2, (0, t), 0.01, RungKutta2())
	@test abs(tr(rho3)-1) < tol2

	rho4 = timeevo(rho, L1, (0, t), 0.01, RungKutta2())
	@test norm(rho2 - rho4) < tol2


	# 4 order rung kutta
	tol2 = 1.0e-6
	rho3 = timeevo(rho, L2, (0, t), 0.01, RungKutta4())
	@test abs(tr(rho3)-1) < tol2

	rho4 = timeevo(rho, L1, (0, t), 0.01, RungKutta4())
	@test norm(rho2 - rho4) < tol2
end