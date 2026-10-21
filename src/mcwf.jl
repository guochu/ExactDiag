# Monte Carlo wavefunction（量子轨迹）——单条轨迹积分器（IonChainSimulator 计划
# §2.2/§9.4：MCWF 作为通用求解能力回馈 ExactDiag；系综编排/统计在上层包）。
#
# 算法（Mølmer–Castin–Dalibard 一阶变体）：
#   1. 抽 p ~ U(0,1)；以无跳变有效方程 dψ/dt = hc·ψ 演化（hc = −iH − ΣC†C = L.hc，
#      即 Lindbladian 的反厄米部分——范数衰减 d‖ψ‖²/dt = −2Σ⟨C†C⟩ 与本包
#      dissipator 约定 dρ += 2Σ(CρC†−½{C†C,ρ}) 自洽：通道 j 跳变率 = 2⟨C†C⟩）；
#   2. 跳变事件在步末检验：若 ‖ψ‖² < p，按 w_j = ‖C_jψ‖² 选通道，
#      ψ ← C_jψ/‖C_jψ‖，重抽 p；否则不重归一（范数携带下一次跳变信息）；
#   3. 输出时刻返回归一化副本。
# 跳变时刻取步末（偏差 O(dt)）：dt 需满足 max(2⟨C†C⟩)·dt << 1。

"""
    mcwf_trajectory(L::Lindbladian, ψ0, times; dt, rng=default_rng()) -> Vector{Vector{ComplexF64}}

恒定 Lindbladian 的单条量子轨迹（Monte Carlo wavefunction），返回各输出时刻的
（归一化）态矢量。`times` 须严格递增且 `times[1]` 为初态时刻。
系综平均与方差估计由上层包编排（如 IonChainSimulator.mcwf_ensemble）。
"""
function mcwf_trajectory(L::Lindbladian, ψ0::AbstractVector,
                         times::AbstractVector{<:Real}; dt::Real,
                         rng::Random.AbstractRNG=Random.default_rng())
    isconstant(L) || throw(ArgumentError("mcwf_trajectory only supports constant Lindbladian"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    all(diff(times) .> 0) || throw(ArgumentError("times must be strictly increasing"))
    ψ = complex.(collect(ψ0))
    ψ ./= norm(ψ)
    out = Vector{Vector{ComplexF64}}(undef, length(times))
    out[1] = copy(ψ)
    p = rand(rng)
    for k in 2:length(times)
        nsteps = max(1, round(Int, (times[k] - times[k-1]) / dt))
        h = (times[k] - times[k-1]) / nsteps
        for s in 1:nsteps
            p = _mcwf_step!(ψ, L, h, p, rng)
        end
        nrm = norm(ψ)
        out[k] = ψ ./ nrm
    end
    return out
end

"一步无跳变 RK4 + 跳变判定；返回更新后的随机阈值 p。"
function _mcwf_step!(ψ::AbstractVector, L::Lindbladian, h::Real, p::Real,
                     rng::Random.AbstractRNG)
    hc = L.hc
    k1 = hc * ψ
    k2 = hc * (ψ .+ h / 2 .* k1)
    k3 = hc * (ψ .+ h / 2 .* k2)
    k4 = hc * (ψ .+ h .* k3)
    ψ .+= (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    n2 = real(dot(ψ, ψ))
    if n2 < p
        C = L.dissipators
        isempty(C) && error("norm decayed without dissipators (inconsistent Lindbladian)")
        w = [real(dot(Cj * ψ, Cj * ψ)) for Cj in C]
        j = _weighted_sample(w, rand(rng))
        ψ .= C[j] * ψ
        ψ ./= norm(ψ)
        return rand(rng)
    end
    return p
end

function _weighted_sample(w::AbstractVector{<:Real}, r::Real)
    tot = sum(w)
    tot > 0 || throw(ArgumentError("zero total jump weight"))
    acc = 0.0
    for (j, wj) in enumerate(w)
        acc += wj / tot
        r < acc && return j
    end
    return length(w)
end
