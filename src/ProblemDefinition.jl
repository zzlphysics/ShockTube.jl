module ProblemDefinition

export ShockTubeProblem, initialize_conditions

# Define the ShockTubeProblem structure
Base.@kwdef struct ShockTubeProblem
    geometry::Tuple{Float64, Float64, Float64} = (0.0, 1.0, 0.5) # (x_left, x_right, x_interface)
    left_state::NamedTuple{(:ρ, :u, :p), Tuple{Float64, Float64, Float64}} = (ρ=1.0, u=0.0, p=1.0)
    right_state::NamedTuple{(:ρ, :u, :p), Tuple{Float64, Float64, Float64}} = (ρ=0.125, u=0.0, p=0.1)
    t_end::Float64 = 0.2
    γ::Float64 = 1.4
    nx::Int = 100
    CFL::Float64 = 0.7
end

# Function to initialize conditions
function initialize_conditions(x, left_state, right_state)
    if x < left_state.x_interface
        return left_state.ρ, left_state.u, left_state.p
    else
        return right_state.ρ, right_state.u, right_state.p
    end
end

end # module ProblemDefinition
