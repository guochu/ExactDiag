
"""
    groundstate(h::Hamiltonian; maxiter=Defaults.maxiter_eig,
                tol=Defaults.maxiter_tol, krylovdim=Defaults.krylovdim)
        -> (eigval, eigvec)

常量 Hamiltonian 的基态（KrylovKit `eigsolve` 最小实部本征对）。
`maxiter` / `tol` / `krylovdim` 透传 eigsolve，缺省取 `ExactDiag.Defaults`。
"""
function groundstate(h::Hamiltonian; maxiter::Int=Defaults.maxiter_eig,
                     tol::Real=Defaults.maxiter_tol, krylovdim::Integer=Defaults.krylovdim)
	isconstant(h) || throw(ArgumentError("groundstate only defined for constant Hamiltonian"))
	d = size(h, 1)
	eigval, eigvec, info = eigsolve(h, randn(eltype(h), d), 1, :SR; ishermitian=true,
	                                maxiter=maxiter, tol=tol, krylovdim=krylovdim)
	(info.converged>=1) || error("eigsolve fails to converge")
	return eigval[1], eigvec[1]
end

"""
    steadystate(h::Lindbladian; maxiter=Defaults.maxiter_eig,
                tol=Defaults.maxiter_tol, krylovdim=Defaults.krylovdim)
        -> (eigval, rho)

常量 Lindbladian 的稳态（|λ| 最小的 Liouvillian 本征对，迹归一密度矩阵）。
`maxiter` / `tol` / `krylovdim` 透传 eigsolve，缺省取 `ExactDiag.Defaults`。
"""
function steadystate(h::Lindbladian; maxiter::Int=Defaults.maxiter_eig,
                     tol::Real=Defaults.maxiter_tol, krylovdim::Integer=Defaults.krylovdim)
    isconstant(h) || throw(ArgumentError("steadystate only defined for constant Lindbladian"))
    d = size(h.hc, 1)   # 物理维度（Lindbladian 作用在 d×d 密度矩阵上；size(h,1) = d²）
    eigval, eigvec, info = eigsolve(h, randn(eltype(h), d, d), 1, EigSorter(abs; rev = false);
                                    ishermitian=false, maxiter=maxiter, tol=tol,
                                    krylovdim=krylovdim)
    (info.converged>=1) || error("eigsolve fails to converge")
    rho = eigvec[1]
    rho ./= tr(rho)
    return eigval[1], rho
end
