module TimeIntegration

export TimeIntegrator, EulerForward, RK3, advance_in_time

using ..ProblemDefinition
using ..Solvers
using ..BoundaryConditions
using ..Reconstruction

# 抽象时间推进器类型
abstract type TimeIntegrator end

# 具体时间推进器类型
struct EulerForward <: TimeIntegrator end
struct RK3 <: TimeIntegrator end

# 时间推进函数接口
function advance_in_time(
    integrator::TimeIntegrator,
    state,
    dt::Float64,
    dx::Float64,
    problem::ShockTubeProblem,
    solver::Solver,
    bc::BoundaryCondition,
    recon::ReconstructionMethod
)
    if integrator isa EulerForward
        advance_euler_forward!(state, dt, dx, problem, solver, bc, recon)
    elseif integrator isa RK3
        advance_rk3!(state, dt, dx, problem, solver, bc, recon)
    else
        error("Unknown time integrator type")
    end
end

# 欧拉显式方法
function advance_euler_forward!(state, dt::Float64, dx::Float64, problem::ShockTubeProblem, solver::Solver, bc::BoundaryCondition, recon::ReconstructionMethod)
    nx = length(state.ρ) - 2  # 去掉幽灵单元
    γ = problem.γ

    # 应用边界条件
    apply_boundary_conditions!(state, bc)

    # 计算通量
    flux = compute_flux(solver, state, problem, recon)

    # 提取守恒变量
    U = zeros(3, nx + 2)  # 包含幽灵单元
    U[1, :] = state.ρ
    U[2, :] = state.ρ .* state.u
    U[3, :] = state.p ./ (γ - 1) + 0.5 * state.ρ .* state.u.^2

    # 更新守恒变量
    for i in 2:nx+1  # 跳过幽灵单元
        U[:, i] = U[:, i] - dt / dx * (flux[:, i] - flux[:, i - 1])
    end

    # 更新状态变量
    state.ρ .= U[1, :]
    state.u .= U[2, :] ./ state.ρ
    state.p .= (γ - 1) .* (U[3, :] .- 0.5 .* state.ρ .* state.u.^2)
end

# 三阶 TVD Runge-Kutta 方法
function advance_rk3!(state, dt::Float64, dx::Float64, problem::ShockTubeProblem, solver::Solver, bc::BoundaryCondition, recon::ReconstructionMethod)
    nx = length(state.ρ) - 2  # 去掉幽灵单元
    γ = problem.γ

    # Stage 1
    apply_boundary_conditions!(state, bc)
    flux = compute_flux(solver, state, problem, recon)

    U = zeros(3, nx + 2)  # 包含幽灵单元
    U[1, :] = state.ρ
    U[2, :] = state.ρ .* state.u
    U[3, :] = state.p ./ (γ - 1) + 0.5 * state.ρ .* state.u.^2

    U1 = similar(U)
    for i in 2:nx+1
        U1[:, i] = U[:, i] - dt / dx * (flux[:, i] - flux[:, i - 1])
    end

    # 更新状态变量
    state1 = deepcopy(state)
    state1.ρ .= U1[1, :]
    state1.u .= U1[2, :] ./ state1.ρ
    state1.p .= (γ - 1) .* (U1[3, :] .- 0.5 .* state1.ρ .* state1.u.^2)

    # Stage 2
    apply_boundary_conditions!(state1, bc)
    flux = compute_flux(solver, state1, problem, recon)

    U2 = similar(U)
    for i in 2:nx+1
        U2[:, i] = 0.75 * U[:, i] + 0.25 * (U1[:, i] - dt / dx * (flux[:, i] - flux[:, i - 1]))
    end

    # 更新状态变量
    state2 = deepcopy(state)
    state2.ρ .= U2[1, :]
    state2.u .= U2[2, :] ./ state2.ρ
    state2.p .= (γ - 1) .* (U2[3, :] .- 0.5 .* state2.ρ .* state2.u.^2)

    # Stage 3
    apply_boundary_conditions!(state2, bc)
    flux = compute_flux(solver, state2, problem, recon)

    U3 = similar(U)
    for i in 2:nx+1
        U3[:, i] = (1.0/3.0) * U[:, i] + (2.0/3.0) * (U2[:, i] - dt / dx * (flux[:, i] - flux[:, i - 1]))
    end

    # 更新状态变量
    state.ρ .= U3[1, :]
    state.u .= U3[2, :] ./ state.ρ
    state.p .= (γ - 1) .* (U3[3, :] .- 0.5 .* state.ρ .* state.u.^2)
end

end # module TimeIntegration
