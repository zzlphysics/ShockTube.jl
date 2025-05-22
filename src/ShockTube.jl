module ShockTube

# Include other modules
include("ProblemDefinition.jl")
include("BoundaryConditions.jl")
include("Reconstruction.jl")
include("Solvers.jl")
include("TimeIntegration.jl")
include("AnalyticalShockTube.jl")
include("Utilities.jl")


# Using statements for submodules
using .ProblemDefinition
using .Solvers
using .BoundaryConditions
using .Reconstruction
using .TimeIntegration
using .Utilities
using .AnalyticalShockTube

# Exported symbols
export ShockTubeProblem, Solver, BoundaryCondition, ReconstructionMethod, TimeIntegrator
export LaxFriedrichs, HLL, HLLC, Roe, ReflectiveBoundary, PeriodicBoundary, NeumannBoundary, DirichletBoundary
export WENOReconstruction, LinearReconstruction, EulerForward, RK3
export simulate, get_analytical_solution, compare_solutions

end # module ShockTube



# module ShockTube

# using Plots
# include("AnalyticalShockTube.jl")  # Ensure this line is included to access AnalyticalShockTube
# using .AnalyticalShockTube

# export simulate, plot_results, ShockTubeProblem, Solver, LaxFriedrichs, HLL, initialize_conditions, boundary_conditions!, compare_solutions, get_analytical_solution

# # Abstract Solver type for flexibility
# abstract type Solver end

# # Lax-Friedrichs Solver
# struct LaxFriedrichs <: Solver end

# struct HLL <: Solver end


# # Shock Tube Problem structure
# Base.@kwdef struct ShockTubeProblem
#     geometry::Tuple{Float64, Float64, Float64} = (0.0, 1.0, 0.5) # (x_left, x_right, x_interface)
#     left_state::NamedTuple{(:ρ, :u, :p), Tuple{Float64, Float64, Float64}} = (ρ=1.0, u=0.0, p=1.0)
#     right_state::NamedTuple{(:ρ, :u, :p), Tuple{Float64, Float64, Float64}} = (ρ=0.125, u=0.0, p=0.1)
#     t_end::Float64 = 0.2
#     γ::Float64 = 1.4
#     nx::Int = 100
#     CFL::Float64 = 0.7
# end

# # Boundary Conditions
# function boundary_conditions!(ρ, u, p)
#     ρ[1] = ρ[2]
#     ρ[end] = ρ[end-1]
#     u[1] = u[2]
#     u[end] = u[end-1]
#     p[1] = p[2]
#     p[end] = p[end-1]
# end

# # Initial conditions setup
# function initialize_conditions(x, left_state, right_state)
#     if x < 0.5
#         return left_state.ρ, left_state.u, left_state.p
#     else
#         return right_state.ρ, right_state.u, right_state.p
#     end
# end

# function lax_friedrichs!(ρ, u, p, dx, γ, CFL)
#     nx = length(ρ)
#     # Extend arrays to include ghost cells
#     ρ_ext = [ρ[1]; ρ; ρ[end]]
#     u_ext = [u[1]; u; u[end]]
#     p_ext = [p[1]; p; p[end]]
#     E_ext = p_ext / (γ - 1) + 0.5 * ρ_ext .* u_ext.^2

#     # Calculate speed of sound
#     a_ext = sqrt.(γ * p_ext ./ ρ_ext)

#     # Check for negative density or pressure before proceeding
#     if any(ρ_ext .<= 0) || any(p_ext .<= 0)
#         error("Negative density or pressure encountered in extended arrays.")
#     end

#     # Determine the maximum time step based on CFL condition
#     max_speed = maximum(abs.(u_ext) + a_ext)
#     dt = CFL * dx / max_speed  # CFL condition

#     # Calculate fluxes at interfaces
#     F = zeros(3, nx + 1)
#     for i in 1:nx + 1
#         # Compute average states
#         ρ_avg = 0.5 * (ρ_ext[i] + ρ_ext[i+1])
#         u_avg = 0.5 * (u_ext[i] + u_ext[i+1])
#         E_avg = 0.5 * (E_ext[i] + E_ext[i+1])
#         p_avg = (γ - 1) * (E_avg - 0.5 * ρ_avg * u_avg^2)

#         # Ensure positive density and pressure
#         if ρ_avg <= 0 || p_avg <= 0
#             error("Negative density or pressure encountered at interface $i.")
#         end

#         # Compute fluxes
#         F[1, i] = ρ_avg * u_avg
#         F[2, i] = ρ_avg * u_avg^2 + p_avg
#         F[3, i] = (E_avg + p_avg) * u_avg
#     end

#     # Update conserved quantities
#     ρ_new = similar(ρ)
#     u_new = similar(u)
#     p_new = similar(p)

#     for i in 1:nx
#         # Original conserved quantities
#         ρ_old = ρ[i]
#         m_old = ρ[i] * u[i]
#         E_old = p[i] / (γ - 1) + 0.5 * ρ[i] * u[i]^2

#         # Update conserved quantities
#         ρ_new[i] = ρ_old - dt / dx * (F[1, i+1] - F[1, i])
#         m_new = m_old - dt / dx * (F[2, i+1] - F[2, i])
#         E_new = E_old - dt / dx * (F[3, i+1] - F[3, i])

#         # Ensure positive density
#         if ρ_new[i] <= 0
#             error("Negative density encountered at index $i.")
#         end

#         # Update primitive variables
#         u_new[i] = m_new / ρ_new[i]
#         p_new[i] = (γ - 1) * (E_new - 0.5 * ρ_new[i] * u_new[i]^2)

#         # Ensure positive pressure
#         if p_new[i] <= 0
#             error("Negative pressure encountered at index $i.")
#         end
#     end

#     # Apply boundary conditions
#     boundary_conditions!(ρ_new, u_new, p_new)

#     return ρ_new, u_new, p_new, dt  # Returning dt for tracking
# end


# # Corrected HLL Solver Implementation
# function hll!(ρ, u, p, dx, γ, CFL)
#     nx = length(ρ)
#     # Extend arrays to include ghost cells
#     ρ_ext = [ρ[1]; ρ; ρ[end]]
#     u_ext = [u[1]; u; u[end]]
#     p_ext = [p[1]; p; p[end]]
#     E_ext = p_ext / (γ - 1) + 0.5 * ρ_ext .* u_ext.^2

#     # Calculate speed of sound
#     a_ext = sqrt.(γ * p_ext ./ ρ_ext)

#     # Determine the maximum time step based on CFL condition
#     max_speed = maximum(abs.(u_ext) + a_ext)
#     dt = CFL * dx / max_speed  # CFL condition

#     # Initialize flux arrays
#     F = zeros(3, nx + 1)

#     for i in 1:nx + 1
#         # Left and right states at interface
#         UL = [ρ_ext[i], ρ_ext[i] * u_ext[i], E_ext[i]]
#         UR = [ρ_ext[i+1], ρ_ext[i+1] * u_ext[i+1], E_ext[i+1]]

#         # Compute wave speeds
#         aL = sqrt(γ * p_ext[i] / ρ_ext[i])
#         aR = sqrt(γ * p_ext[i+1] / ρ_ext[i+1])
#         sL = min(u_ext[i] - aL, u_ext[i+1] - aR)
#         sR = max(u_ext[i] + aL, u_ext[i+1] + aR)

#         # Compute fluxes for left and right states
#         FL = [ρ_ext[i] * u_ext[i],
#               ρ_ext[i] * u_ext[i]^2 + p_ext[i],
#               (E_ext[i] + p_ext[i]) * u_ext[i]]
#         FR = [ρ_ext[i+1] * u_ext[i+1],
#               ρ_ext[i+1] * u_ext[i+1]^2 + p_ext[i+1],
#               (E_ext[i+1] + p_ext[i+1]) * u_ext[i+1]]

#         # Compute HLL flux at interface
#         if sL >= 0
#             F[:, i] = FL
#         elseif sR <= 0
#             F[:, i] = FR
#         else
#             F[:, i] = (sR * FL - sL * FR + sL * sR * (UR - UL)) / (sR - sL)
#         end
#     end

#     # Update conserved quantities
#     for i in 1:nx
#         U = [ρ[i], ρ[i] * u[i], p[i] / (γ - 1) + 0.5 * ρ[i] * u[i]^2]
#         U_new = U - dt / dx * (F[:, i+1] - F[:, i])

#         # Update variables
#         ρ[i] = U_new[1]
#         u[i] = U_new[2] / U_new[1]
#         E = U_new[3]
#         p[i] = (γ - 1) * (E - 0.5 * ρ[i] * u[i]^2)
#     end

#     # Apply boundary conditions
#     boundary_conditions!(ρ, u, p)

#     return ρ, u, p, dt  # Returning dt for tracking
# end

# # Adjusted simulate function
# function simulate(problem::ShockTubeProblem, solver::Solver=LaxFriedrichs())
#     nx = problem.nx
#     CFL = problem.CFL
#     t_end = problem.t_end
#     γ = problem.γ
#     dx = (problem.geometry[2] - problem.geometry[1]) / (nx - 1)
#     x = range(problem.geometry[1], stop=problem.geometry[2], length=nx)
#     ρ = zeros(nx)
#     u = zeros(nx)
#     p = zeros(nx)

#     # Initialize conditions
#     for i in 1:nx
#         ρ[i], u[i], p[i] = initialize_conditions(x[i], problem.left_state, problem.right_state)
#     end

#     # Time-stepping loop
#     t = 0.0
#     nstep = 0
#     while t < t_end
#         if solver isa LaxFriedrichs
#             ρ, u, p, dt = lax_friedrichs!(ρ, u, p, dx, γ, CFL)
#         elseif solver isa HLL
#             ρ, u, p, dt = hll!(ρ, u, p, dx, γ, CFL)
#         else
#             error("Unknown solver type")
#         end

#         if t + dt > t_end
#             dt = t_end - t  # Adjust dt to not exceed t_end
#         end
        
#         nstep +=1
#         t += dt

#         println(nstep, " ", dt)
#     end

#     return x, ρ, u, p
# end


# # Function to call the analytical solver from AnalyticalShockTube module
# function get_analytical_solution(problem::ShockTubeProblem, x_arr)
#     # Use AnaShockTubeProblem structure to pass parameters to the analytical solver
#     ana_problem = AnaShockTubeProblem(
#         geometry=problem.geometry,
#         left_state=problem.left_state,
#         right_state=problem.right_state,
#         t=problem.t_end,
#         γ=problem.γ
#     )
#     positions, regions, values = anasolve(ana_problem, x_arr)
#     return values
# end

# # Plot results and compare with analytical solution, and save each plot as a file
# function compare_solutions(x, ρ, u, p, analytical_values)
#     # Plot and save Density comparison
#     density_plot = plot(x, ρ, label="Numerical Density", xlabel="x", ylabel="Density")
#     plot!(density_plot, analytical_values.x, analytical_values.ρ, label="Analytical Density", linestyle=:dash)
#     savefig(density_plot, "density_comparison.png")

#     # Plot and save Velocity comparison
#     velocity_plot = plot(x, u, label="Numerical Velocity", xlabel="x", ylabel="Velocity")
#     plot!(velocity_plot, analytical_values.x, analytical_values.u, label="Analytical Velocity", linestyle=:dash)
#     savefig(velocity_plot, "velocity_comparison.png")

#     # Plot and save Pressure comparison
#     pressure_plot = plot(x, p, label="Numerical Pressure", xlabel="x", ylabel="Pressure")
#     plot!(pressure_plot, analytical_values.x, analytical_values.p, label="Analytical Pressure", linestyle=:dash)
#     savefig(pressure_plot, "pressure_comparison.png")
# end



# end # module ShockTube



