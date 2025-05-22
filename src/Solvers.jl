module Solvers

export Solver, LaxFriedrichs, HLL, HLLC, Roe, compute_flux

using ..ProblemDefinition
using ..Reconstruction

# 抽象求解器类型
abstract type Solver end

# 具体求解器类型
struct LaxFriedrichs <: Solver end
struct HLL <: Solver end
struct HLLC <: Solver end
struct Roe <: Solver end 


# 求解器主接口函数
"""
    compute_flux(solver::Solver, state, problem::ShockTubeProblem, recon::ReconstructionMethod) -> flux

计算数值通量。

- `solver`: 求解器类型实例（如 `LaxFriedrichs()`）
- `state`: 当前状态变量，包含 `ρ`、`u`、`p`、`E`，包括幽灵单元
- `problem`: 问题参数，包括 `γ` 等
- `recon`: 重构方法实例（如 `LinearReconstruction()` 或 `WENOReconstruction()`）

返回值：

- `flux`: 在每个界面上的通量数组，大小为 `(3, nx + 2)`，包括所有界面
"""
function compute_flux(solver::Solver, state, problem::ShockTubeProblem, recon::ReconstructionMethod)
    if solver isa LaxFriedrichs
        return compute_flux_lax_friedrichs(state, problem, recon)
    elseif solver isa HLL
        return compute_flux_hll(state, problem, recon)
    elseif solver isa HLLC 
        return compute_flux_hllc(state, problem, recon)
    elseif solver isa Roe 
        return compute_flux_roe(state, problem, recon)
    else
        error("Unknown solver type")
    end
end

# 实现具体的求解器

# Lax-Friedrichs 求解器
function compute_flux_lax_friedrichs(state, problem::ShockTubeProblem, recon::ReconstructionMethod)
    γ = problem.γ
    nx = length(state.ρ) - 2  # 去掉幽灵单元

    # 提取变量，包括幽灵单元
    ρ_ext = state.ρ
    u_ext = state.u
    p_ext = state.p
    # E_ext = state.E

    # 重构变量
    ρL, ρR = reconstruct(recon, ρ_ext)
    uL, uR = reconstruct(recon, u_ext)
    pL, pR = reconstruct(recon, p_ext)
    # EL, ER = reconstruct(recon, E_ext)

    # 初始化通量数组，大小为 (3, nx + 2)，包括所有界面
    flux = zeros(3, nx + 2)

    for i in 1:nx + 1  # 所有界面，包括幽灵单元之间的界面
        # 计算平均状态
        ρ_avg = 0.5 * (ρL[i] + ρR[i])
        u_avg = 0.5 * (uL[i] + uR[i])
        p_avg = 0.5 * (pL[i] + pR[i])
        E_avg = p_avg / (γ - 1) + 0.5 * ρ_avg * u_avg^2

        # 计算通量
        flux[1, i] = ρ_avg * u_avg
        flux[2, i] = ρ_avg * u_avg^2 + p_avg
        flux[3, i] = (E_avg + p_avg) * u_avg
    end

    return flux
end

# HLL 求解器
function compute_flux_hll(state, problem::ShockTubeProblem, recon::ReconstructionMethod)
    γ = problem.γ
    nx = length(state.ρ) - 2  # 去掉幽灵单元

    # 提取变量，包括幽灵单元
    ρ_ext = state.ρ
    u_ext = state.u
    p_ext = state.p
    # E_ext = state.E

    # 重构变量
    ρL, ρR = reconstruct(recon, ρ_ext)
    uL, uR = reconstruct(recon, u_ext)
    pL, pR = reconstruct(recon, p_ext)
    # EL, ER = reconstruct(recon, E_ext)
    EL = pL ./ (γ - 1) + 0.5 * ρL .* uL.^2
    ER = pR ./ (γ - 1) + 0.5 * ρR .* uR.^2

    # 初始化通量数组，大小为 (3, nx + 2)
    flux = zeros(3, nx + 2)

    for i in 1:nx + 1  # 所有界面
        # 左右守恒变量
        UL = [ρL[i], ρL[i] * uL[i], EL[i]]
        UR = [ρR[i], ρR[i] * uR[i], ER[i]]

        # 计算左、右侧声速
        aL = sqrt(γ * pL[i] / ρL[i])
        aR = sqrt(γ * pR[i] / ρR[i])

        # 计算波速
        sL = min(uL[i] - aL, uR[i] - aR)
        sR = max(uL[i] + aL, uR[i] + aR)

        # 计算左、右侧通量
        FL = [ρL[i] * uL[i],
              ρL[i] * uL[i]^2 + pL[i],
              (EL[i] + pL[i]) * uL[i]]
        FR = [ρR[i] * uR[i],
              ρR[i] * uR[i]^2 + pR[i],
              (ER[i] + pR[i]) * uR[i]]

        # 计算 HLL 通量
        if sL >= 0
            flux[:, i] = FL
        elseif sR <= 0
            flux[:, i] = FR
        else
            flux[:, i] = (sR * FL - sL * FR + sL * sR * (UR - UL)) / (sR - sL)
        end
    end

    return flux
end

# HLLC 求解器
function compute_flux_hllc(state, problem::ShockTubeProblem, recon::ReconstructionMethod)
    γ = problem.γ
    nx = length(state.ρ) - 2  # 去掉幽灵单元

    # 提取变量，包括幽灵单元
    ρ_ext = state.ρ
    u_ext = state.u
    p_ext = state.p

    # 重构变量
    ρL, ρR = reconstruct(recon, ρ_ext)
    uL, uR = reconstruct(recon, u_ext)
    pL, pR = reconstruct(recon, p_ext)

    # 计算重构后的总能量 EL, ER
    EL = pL ./ (γ - 1) .+ 0.5 .* ρL .* uL.^2
    ER = pR ./ (γ - 1) .+ 0.5 .* ρR .* uR.^2

    # 初始化通量数组
    flux = zeros(3, nx + 1)

    for i in 1:nx + 1  # 所有界面
        # 左侧状态
        ρL_i = ρL[i]
        uL_i = uL[i]
        pL_i = pL[i]
        EL_i = EL[i]

        # 右侧状态
        ρR_i = ρR[i]
        uR_i = uR[i]
        pR_i = pR[i]
        ER_i = ER[i]

        # 计算左、右侧的比焓
        HL_i = (EL_i + pL_i) / ρL_i
        HR_i = (ER_i + pR_i) / ρR_i

        # 计算左、右侧声速
        aL = sqrt(γ * pL_i / ρL_i)
        aR = sqrt(γ * pR_i / ρR_i)

        # 估计波速 sL 和 sR
        sL = min(uL_i - aL, uR_i - aR)
        sR = max(uL_i + aL, uR_i + aR)

        # 避免除零错误，添加一个小量 ε
        ε = 1e-6
        if sR - sL < ε
            sR += ε
            sL -= ε
        end

        # 计算中间波速 s_star
        numerator = pR_i - pL_i + ρL_i * uL_i * (sL - uL_i) - ρR_i * uR_i * (sR - uR_i)
        denominator = ρL_i * (sL - uL_i) - ρR_i * (sR - uR_i)
        s_star = numerator / (denominator + ε)  # 避免除零

        # 构建左、右侧守恒变量和通量
        UL_i = [ρL_i, ρL_i * uL_i, EL_i]
        UR_i = [ρR_i, ρR_i * uR_i, ER_i]

        FL_i = [ρL_i * uL_i,
                ρL_i * uL_i^2 + pL_i,
                (EL_i + pL_i) * uL_i]
        FR_i = [ρR_i * uR_i,
                ρR_i * uR_i^2 + pR_i,
                (ER_i + pR_i) * uR_i]

        if sL >= 0
            flux[:, i] = FL_i
        elseif sL <= 0 <= s_star
            # 计算左侧星区域的守恒变量
            ρ_star_L = ρL_i * (sL - uL_i) / (sL - s_star)
            m_star_L = ρ_star_L * s_star
            E_star_L = ρ_star_L * ( (EL_i / ρL_i) + (s_star - uL_i) * (s_star + pL_i / (ρL_i * (sL - uL_i))) )
            U_star_L = [ρ_star_L, m_star_L, E_star_L]
            flux[:, i] = FL_i + sL * (U_star_L - UL_i)
        elseif s_star <= 0 <= sR
            # 计算右侧星区域的守恒变量
            ρ_star_R = ρR_i * (sR - uR_i) / (sR - s_star)
            m_star_R = ρ_star_R * s_star
            E_star_R = ρ_star_R * ( (ER_i / ρR_i) + (s_star - uR_i) * (s_star + pR_i / (ρR_i * (sR - uR_i))) )
            U_star_R = [ρ_star_R, m_star_R, E_star_R]
            flux[:, i] = FR_i + sR * (U_star_R - UR_i)
        else  # sR <= 0
            flux[:, i] = FR_i
        end
    end

    return flux
end

# Roe 求解器
function compute_flux_roe(state, problem::ShockTubeProblem, recon::ReconstructionMethod)
    γ = problem.γ
    nx = length(state.ρ) - 2  # 去掉幽灵单元

    # 提取变量，包括幽灵单元
    ρ_ext = state.ρ
    u_ext = state.u
    p_ext = state.p

    # 重构变量
    ρL, ρR = reconstruct(recon, ρ_ext)
    uL, uR = reconstruct(recon, u_ext)
    pL, pR = reconstruct(recon, p_ext)

    # 计算总能量
    EL = pL ./ (γ - 1) .+ 0.5 .* ρL .* uL.^2
    ER = pR ./ (γ - 1) .+ 0.5 .* ρR .* uR.^2

    # 构建守恒变量 UL, UR
    UL = zeros(3, nx + 1)
    UR = zeros(3, nx + 1)
    UL[1, :] = ρL
    UL[2, :] = ρL .* uL
    UL[3, :] = EL
    UR[1, :] = ρR
    UR[2, :] = ρR .* uR
    UR[3, :] = ER

    # 初始化通量数组
    flux = zeros(3, nx + 1)

    for i in 1:nx + 1
        # 左右状态
        ρL_i = ρL[i]
        uL_i = uL[i]
        pL_i = pL[i]
        EL_i = EL[i]

        ρR_i = ρR[i]
        uR_i = uR[i]
        pR_i = pR[i]
        ER_i = ER[i]

        # 左右守恒变量
        UL_i = [ρL_i, ρL_i * uL_i, EL_i]
        UR_i = [ρR_i, ρR_i * uR_i, ER_i]

        # 左右物理通量
        FL_i = [ρL_i * uL_i,
                ρL_i * uL_i^2 + pL_i,
                (EL_i + pL_i) * uL_i]
        FR_i = [ρR_i * uR_i,
                ρR_i * uR_i^2 + pR_i,
                (ER_i + pR_i) * uR_i]

        # 计算 Roe 平均量
        sqrt_ρL = sqrt(ρL_i)
        sqrt_ρR = sqrt(ρR_i)
        ρ_avg = sqrt_ρL * sqrt_ρR
        u_avg = (sqrt_ρL * uL_i + sqrt_ρR * uR_i) / (sqrt_ρL + sqrt_ρR)
        HL = (EL_i + pL_i) / ρL_i
        HR = (ER_i + pR_i) / ρR_i
        H_avg = (sqrt_ρL * HL + sqrt_ρR * HR) / (sqrt_ρL + sqrt_ρR)
        a_avg = sqrt((γ - 1) * (H_avg - 0.5 * u_avg^2))

        ρ_avg = max(ρ_avg, 1e-6)
        a_avg = max(a_avg, 1e-6)


        # 计算特征值
        λ1 = u_avg - a_avg
        λ2 = u_avg
        λ3 = u_avg + a_avg

        # Entropy fix
        δ = 0.1 * a_avg
        if abs(λ1) < δ
            λ1 = 0.5 * (λ1^2 / δ + δ)
        end
        if abs(λ3) < δ
            λ3 = 0.5 * (λ3^2 / δ + δ)
        end

        # 计算特征向量
        K1 = [1, u_avg - a_avg, H_avg - u_avg * a_avg]
        K2 = [1, u_avg, 0.5 * u_avg^2]
        K3 = [1, u_avg + a_avg, H_avg + u_avg * a_avg]

        # 计算变量增量
        ΔU = UR_i - UL_i

        # 计算波强度（α 系数）
        Δp = pR_i - pL_i
        Δu = uR_i - uL_i
        Δρ = ρR_i - ρL_i

        α2 = ( (γ - 1) * (Δρ * (H_avg - u_avg^2) + u_avg * Δu * (2 * u_avg) - Δp) ) / (a_avg^2)
        α1 = (Δρ * (u_avg + a_avg) - Δu * ρ_avg - a_avg * α2) / (2 * a_avg)
        α3 = Δρ - α1 - α2

        # 计算数值通量
        flux[:, i] = 0.5 * (FL_i + FR_i) - 0.5 * (α1 * abs(λ1) * K1 + α2 * abs(λ2) * K2 + α3 * abs(λ3) * K3)
    end

    return flux
end

end # module Solvers
