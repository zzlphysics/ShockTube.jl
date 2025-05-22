module BoundaryConditions

export BoundaryCondition, ReflectiveBoundary, PeriodicBoundary, NeumannBoundary, DirichletBoundary
export apply_boundary_conditions!

# Abstract BoundaryCondition type
abstract type BoundaryCondition end

# Concrete BoundaryCondition types

"""
    ReflectiveBoundary()

反射边界条件。在边界处，速度分量取反，其余变量保持不变。
"""
struct ReflectiveBoundary <: BoundaryCondition end

"""
    PeriodicBoundary()

周期边界条件。边界处的变量与对侧边界的变量相等。
"""
struct PeriodicBoundary <: BoundaryCondition end

"""
    NeumannBoundary()

Neumann（自然）边界条件。在边界处，变量的空间导数为零，即变量在边界处的值与相邻内部点的值相同。
"""
struct NeumannBoundary <: BoundaryCondition end

"""
    DirichletBoundary(ρ::Float64, u::Float64, p::Float64)

Dirichlet（固定值）边界条件。在边界处，变量取指定的固定值。
"""
struct DirichletBoundary <: BoundaryCondition
    ρ::Float64
    u::Float64
    p::Float64
end

# Function to apply boundary conditions

"""
    apply_boundary_conditions!(state::State, bc::BoundaryCondition)

根据指定的边界条件类型，更新状态变量 `state` 的边界值。
"""
function apply_boundary_conditions!(state, bc::BoundaryCondition)
    if bc isa ReflectiveBoundary
        apply_reflective_boundary!(state)
    elseif bc isa PeriodicBoundary
        apply_periodic_boundary!(state)
    elseif bc isa NeumannBoundary
        apply_neumann_boundary!(state)
    elseif bc isa DirichletBoundary
        apply_dirichlet_boundary!(state, bc)
    else
        error("Unknown boundary condition type")
    end
end

# Implementations of boundary condition applications

# Reflective Boundary Condition
function apply_reflective_boundary!(state)
    # Assuming ghost cells are at index 1 and end
    # Left boundary
    state.ρ[1] = state.ρ[2]
    state.u[1] = -state.u[2]     # Reflect velocity
    state.p[1] = state.p[2]
    # Right boundary
    state.ρ[end] = state.ρ[end-1]
    state.u[end] = -state.u[end-1]  # Reflect velocity
    state.p[end] = state.p[end-1]
end

# Periodic Boundary Condition
function apply_periodic_boundary!(state)
    # Left boundary
    state.ρ[1] = state.ρ[end-1]
    state.u[1] = state.u[end-1]
    state.p[1] = state.p[end-1]
    # Right boundary
    state.ρ[end] = state.ρ[2]
    state.u[end] = state.u[2]
    state.p[end] = state.p[2]
end

# Neumann Boundary Condition (Zero Gradient)
function apply_neumann_boundary!(state)
    # Left boundary
    state.ρ[1] = state.ρ[2]
    state.u[1] = state.u[2]
    state.p[1] = state.p[2]
    # Right boundary
    state.ρ[end] = state.ρ[end-1]
    state.u[end] = state.u[end-1]
    state.p[end] = state.p[end-1]
end

# Dirichlet Boundary Condition (Fixed Values)
function apply_dirichlet_boundary!(state, bc::DirichletBoundary)
    # Left boundary
    state.ρ[1] = bc.ρ
    state.u[1] = bc.u
    state.p[1] = bc.p
    # Right boundary
    state.ρ[end] = bc.ρ
    state.u[end] = bc.u
    state.p[end] = bc.p
end

end # module BoundaryConditions
