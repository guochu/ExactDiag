using Test, LinearAlgebra, Random

@testset "MCWF trajectory" begin
    Random.seed!(20261020)
    # 二能级自发辐射：系综 ⟨P_e(t)⟩ ≈ e^{−Γt}（ExactDiag 约定 C = √(Γ/2)σ⁻ ⇒ 衰率 Γ）
    Γ = 1.0
    sm = [0 1; 0 0] .* 1.0im                    # σ⁻：|e⟩→|g⟩（基序 |g⟩,|e⟩）
    L = Lindbladian(zeros(ComplexF64, 2, 2), [sqrt(Γ / 2) .* sm])
    times = collect(0:0.2:3.0)
    ntraj = 4000
    acc = zeros(length(times))
    for _ in 1:ntraj
        ψ0 = [0.0im, 1.0]                       # |e⟩
        traj = mcwf_trajectory(L, ψ0, times; dt = 0.005)
        acc .+= [abs2(tr[2]) for tr in traj]
    end
    mean_pe = acc ./ ntraj
    for (k, t) in enumerate(times)
        σ = sqrt(exp(-Γ * t) * (1 - exp(-Γ * t)) / ntraj)
        @test abs(mean_pe[k] - exp(-Γ * t)) < 5σ + 5e-3
    end

    # 无耗散极限（零耗散算符）：范数与态保持
    L0 = Lindbladian(zeros(ComplexF64, 2, 2), [zeros(ComplexF64, 2, 2)])
    traj = mcwf_trajectory(L0, [1.0im, 0.0], [0.0, 1.0]; dt = 0.01)
    @test norm(traj[2]) ≈ 1.0 atol = 1e-12
    @test abs2(traj[2][1]) ≈ 1.0 atol = 1e-12
end
