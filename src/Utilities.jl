module Utilities

export simulate, compare_solutions, get_analytical_solution

using ..ProblemDefinition
using ..Solvers
using ..BoundaryConditions
using ..Reconstruction
using ..TimeIntegration
using ..AnalyticalShockTube

# 导入必要的包用于绘图（如果需要）
using Plots

# 主模拟函数
"""
    simulate(
        problem::ShockTubeProblem,
        solver::Solver,
        bc::BoundaryCondition,
        recon::ReconstructionMethod,
        integrator::TimeIntegrator
    ) -> (x, ρ, u, p)

运行数值模拟，返回位置和状态变量。

- `problem`: 问题定义，包括几何、初始条件、物理参数等
- `solver`: 求解器类型实例（如 `LaxFriedrichs()`）
- `bc`: 边界条件类型实例（如 `ReflectiveBoundary()`）
- `recon`: 重构方法实例（如 `LinearReconstruction()`）
- `integrator`: 时间推进方案实例（如 `EulerForward()`）

返回值：

- `x`: 网格点位置数组
- `ρ`, `u`, `p`: 最终时刻的状态变量数组
"""
function simulate(
    problem::ShockTubeProblem,
    solver::Solver,
    bc::BoundaryCondition,
    recon::ReconstructionMethod,
    integrator::TimeIntegrator
)
    # 初始化
    nx = problem.nx
    CFL = problem.CFL
    t_end = problem.t_end
    γ = problem.γ
    x_left, x_right, x_interface = problem.geometry
    dx = (x_right - x_left) / (nx - 1)
    x = range(x_left, stop=x_right, length=nx)
    
    # 扩展网格以包含幽灵单元
    x_ext = [x_left - dx; x; x_right + dx]
    nx_ext = length(x_ext)
    
    # 初始化状态变量，包括幽灵单元
    ρ = zeros(nx_ext)
    u = zeros(nx_ext)
    p = zeros(nx_ext)
    
    # 初始化条件
    for i in 2:nx+1  # 跳过第一个幽灵单元
        xi = x[i - 1]  # 对应的物理位置
        if xi < x_interface
            ρ[i] = problem.left_state.ρ
            u[i] = problem.left_state.u
            p[i] = problem.left_state.p
        else
            ρ[i] = problem.right_state.ρ
            u[i] = problem.right_state.u
            p[i] = problem.right_state.p
        end
        # E[i] = p[i] / (γ - 1) + 0.5 * ρ[i] * u[i]^2
    end
    
    # 创建状态结构体
    state = (ρ=ρ, u=u, p=p)
    
    # 时间推进
    t = 0.0
    while t < t_end
        # 应用边界条件
        apply_boundary_conditions!(state, bc)
        
        # 计算最大波速以确定时间步长
        a = sqrt.(γ * state.p ./ state.ρ)
        max_speed = maximum(abs.(state.u) + a)
        dt = CFL * dx / max_speed
        if t + dt > t_end
            dt = t_end - t
        end
        
        # 调用时间推进函数
        advance_in_time(
            integrator,
            state,
            dt,
            dx,
            problem,
            solver,
            bc,
            recon
        )

        # 更新状态变量后
        state.ρ .= max.(state.ρ, 1e-6)
        state.p .= max.(state.p, 1e-6)

        # 更新时间
        t += dt
        println("Time: ", t)
    end
    
    # 返回去除幽灵单元的状态变量
    return x, state.ρ[2:end-1], state.u[2:end-1], state.p[2:end-1]
end

# 函数：获取解析解
function get_analytical_solution(problem::ShockTubeProblem, x_arr)
    ana_problem = AnaShockTubeProblem(
        geometry=problem.geometry,
        left_state=problem.left_state,
        right_state=problem.right_state,
        t=problem.t_end,
        γ=problem.γ
    )
    positions, regions, values = anasolve(ana_problem, x_arr)
    return values
end

# 函数：比较数值解和解析解
function compare_solutions(x, ρ, u, p, analytical_values)
    # 绘制并保存密度比较图
    density_plot = plot(x, ρ, label="Numerical Density", xlabel="x", ylabel="Density")
    plot!(density_plot, analytical_values.x, analytical_values.ρ, label="Analytical Density", linestyle=:dash)
    savefig(density_plot, "density_comparison.png")
    
    # 绘制并保存速度比较图
    velocity_plot = plot(x, u, label="Numerical Velocity", xlabel="x", ylabel="Velocity")
    plot!(velocity_plot, analytical_values.x, analytical_values.u, label="Analytical Velocity", linestyle=:dash)
    savefig(velocity_plot, "velocity_comparison.png")
    
    # 绘制并保存压力比较图
    pressure_plot = plot(x, p, label="Numerical Pressure", xlabel="x", ylabel="Pressure")
    plot!(pressure_plot, analytical_values.x, analytical_values.p, label="Analytical Pressure", linestyle=:dash)
    savefig(pressure_plot, "pressure_comparison.png")
end

end # module Utilities
