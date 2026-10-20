"""
    ExactDiag.Defaults

包级数值缺省值（KrylovKit `eigsolve` 系列调用共用：`groundstate` / `steadystate`
的 maxiter、tol、krylovdim）。调用方可用 keyword 覆盖。
"""
module Defaults

export maxiter_eig, maxiter_tol, krylovdim

"eigsolve 最大迭代数（谱类求解的 maxiter）。"
const maxiter_eig::Int = 100000

"eigsolve 收敛容差（tol，与 KrylovKit.KrylovDefaults 一致）。"
const maxiter_tol::Float64 = 1e-12

"eigsolve Krylov 子空间维数（与 KrylovKit.KrylovDefaults 一致）。"
const krylovdim::Int = 30

end # module
