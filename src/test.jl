# Load the ShockTube module
include("ShockTube.jl")
using .ShockTube

# Define the problem setup
left_state = (ρ = 1.0, u = 0.0, p = 1.0)    # Left initial state: density, velocity, pressure
right_state = (ρ = 0.125, u = 0.0, p = 0.1) # Right initial state: density, velocity, pressure
geometry = (0.0, 1.0, 0.5)                 # Geometry: (x_left, x_right, x_interface)

# Create the ShockTubeProblem instance with the defined parameters
problem = ShockTubeProblem(
    geometry = geometry,
    left_state = left_state,
    right_state = right_state,
    t_end = 0.1,    # Simulation end time
    γ = 1.4,        # Specific heat ratio
    nx = 500,       # Number of grid points
    CFL = 0.1      # Time step size
)

# 选择组件
solver = HLLC() #LaxFriedrichs HLL HLLC Roe
bc = NeumannBoundary()
recon = LinearReconstruction()  # 可以选择 LinearReconstruction() WENOReconstruction()
integrator = RK3()  # 或者 RK3()  EulerForward

# 运行模拟
x, ρ, u, p = simulate(problem, solver, bc, recon, integrator)

# Generate the x array for plotting
x_arr = range(problem.geometry[1], stop=problem.geometry[2], length=problem.nx)

# Run the numerical simulation using Lax-Friedrichs method
# x, ρ_num, u_num, p_num = simulate(problem, LaxFriedrichs())

# Obtain the analytical solution from AnalyticalShockTube
analytical_values = get_analytical_solution(problem, x_arr)

# Compare the numerical and analytical solutions
compare_solutions(x, ρ, u, p, analytical_values)

