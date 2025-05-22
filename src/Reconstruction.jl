module Reconstruction

export ReconstructionMethod, LinearReconstruction, WENOReconstruction, reconstruct

# 抽象重构方法类型
abstract type ReconstructionMethod end

# 具体重构方法类型
struct LinearReconstruction <: ReconstructionMethod end
struct WENOReconstruction <: ReconstructionMethod end

# 重构函数接口
"""
    reconstruct(recon::ReconstructionMethod, q_ext::Vector{Float64})

对变量数组 `q_ext` 进行重构，返回左、右侧的变量值。

- `recon`: 重构方法的实例（如 `LinearReconstruction()`）
- `q_ext`: 扩展后的变量数组，包括幽灵单元

返回值：

- `(qL, qR)`: 左、右侧的变量数组，大小为 `nx + 1`
"""
function reconstruct(recon::ReconstructionMethod, q_ext::Vector{Float64})
    if recon isa LinearReconstruction
        return linear_reconstruction(q_ext)
    elseif recon isa WENOReconstruction
        return weno5_reconstruction(q_ext)
    else
        error("Unknown reconstruction method")
    end
end

# 实现具体的重构方法

# 线性重构（一阶精度）
function linear_reconstruction(q_ext::Vector{Float64})
    nx = length(q_ext) - 2  # 去掉幽灵单元
    qL = zeros(nx + 1)
    qR = zeros(nx + 1)

    for i in 1:nx + 1
        qL[i] = q_ext[i]
        qR[i] = q_ext[i + 1]
    end

    return qL, qR
end

function weno_reconstruction(q_ext::Vector{Float64})
    nx = length(q_ext) - 2  # 去掉幽灵单元
    qL = zeros(nx + 1)
    qR = zeros(nx + 1)

    # WENO5 系数
    ε = 1e-6  # 防止除零的微小数
    γ = [1.0/10.0, 6.0/10.0, 3.0/10.0]  # 线性权重

    # 计算右界面值 qR
    for i in 3:nx
        # 构建模板
        s0 = q_ext[i - 2:i]
        s1 = q_ext[i - 1:i + 1]
        s2 = q_ext[i:i + 2]

        # 计算平滑指标 β
        β = zeros(3)
        β[1] = (13/12)*(s0[1] - 2*s0[2] + s0[3])^2 + (1/4)*(s0[1] - 4*s0[2] + 3*s0[3])^2
        β[2] = (13/12)*(s1[1] - 2*s1[2] + s1[3])^2 + (1/4)*(s1[1] - s1[3])^2
        β[3] = (13/12)*(s2[1] - 2*s2[2] + s2[3])^2 + (1/4)*(3*s2[1] - 4*s2[2] + s2[3])^2

        # 计算 α 权重
        α = γ ./ ((ε .+ β).^2)
        ω = α / sum(α)

        # 多项式重构
        p0 = (2*s0[1] - 7*s0[2] + 11*s0[3]) / 6
        p1 = (-s1[1] + 5*s1[2] + 2*s1[3]) / 6
        p2 = (2*s2[1] + 5*s2[2] - s2[3]) / 6

        # 重构值
        qR[i] = ω[1]*p0 + ω[2]*p1 + ω[3]*p2
    end

    # 边界处采用三阶外推方法
    qR[1] = 3 * q_ext[2] - 3 * q_ext[3] + q_ext[4]  # 三阶外推
    qR[2] = 3 * q_ext[3] - 3 * q_ext[4] + q_ext[5]  # 三阶外推
    qR[nx + 1] = 3 * q_ext[nx + 1] - 3 * q_ext[nx] + q_ext[nx - 1]  # 三阶外推

    # 计算左界面值 qL
    for i in 3:nx
        # 构建模板（反向）
        s0 = q_ext[i + 2:-1:i]
        s1 = q_ext[i + 1:-1:i - 1]
        s2 = q_ext[i:-1:i - 2]

        # 计算平滑指标 β
        β = zeros(3)
        β[1] = (13/12)*(s0[1] - 2*s0[2] + s0[3])^2 + (1/4)*(s0[1] - 4*s0[2] + 3*s0[3])^2
        β[2] = (13/12)*(s1[1] - 2*s1[2] + s1[3])^2 + (1/4)*(s1[1] - s1[3])^2
        β[3] = (13/12)*(s2[1] - 2*s2[2] + s2[3])^2 + (1/4)*(3*s2[1] - 4*s2[2] + s2[3])^2

        # 计算 α 权重
        α = γ ./ ((ε .+ β).^2)
        ω = α / sum(α)

        # 多项式重构
        p0 = (2*s0[1] - 7*s0[2] + 11*s0[3]) / 6
        p1 = (-s1[1] + 5*s1[2] + 2*s1[3]) / 6
        p2 = (2*s2[1] + 5*s2[2] - s2[3]) / 6

        # 重构值
        qL[i] = ω[1]*p0 + ω[2]*p1 + ω[3]*p2
    end

    # 边界处采用三阶外推方法
    qL[1] = 3 * q_ext[2] - 3 * q_ext[3] + q_ext[4]  # 三阶外推
    qL[2] = 3 * q_ext[3] - 3 * q_ext[4] + q_ext[5]  # 三阶外推
    qL[nx + 1] = 3 * q_ext[nx + 1] - 3 * q_ext[nx] + q_ext[nx - 1]  # 三阶外推

    return qL, qR
end


function weno5_reconstruction(q_ext::Vector{Float64})
    nx = length(q_ext) - 2  # q_ext 包含 nx 个内部单元和 2 个幽灵单元
    qL = zeros(nx + 1)
    qR = zeros(nx + 1)

    # 扩展 q_ext，使其左右各有 3 个幽灵单元
    q_ext_extended = zeros(nx + 6)
    
    # 将原始数据复制到扩展数组中
    q_ext_extended[3:nx + 4] = q_ext

    # 左边界外插（使用已有的幽灵单元）
    q_ext_extended[2] = 2 * q_ext_extended[3] - q_ext_extended[4]
    q_ext_extended[1] = 2 * q_ext_extended[2] - q_ext_extended[3]

    # 右边界外插
    q_ext_extended[nx + 5] = 2 * q_ext_extended[nx + 4] - q_ext_extended[nx + 3]
    q_ext_extended[nx + 6] = 2 * q_ext_extended[nx + 5] - q_ext_extended[nx + 4]

    # 循环遍历每个界面
    for i in 1:nx + 1
        j = i + 2  # 调整索引以匹配 q_ext_extended

        # 左侧重构（qL）
        beta0 = (13/12)*(q_ext_extended[j - 2] - 2*q_ext_extended[j - 1] + q_ext_extended[j])^2 +
                (1/4)*(q_ext_extended[j - 2] - 4*q_ext_extended[j - 1] + 3*q_ext_extended[j])^2
        beta1 = (13/12)*(q_ext_extended[j - 1] - 2*q_ext_extended[j] + q_ext_extended[j + 1])^2 +
                (1/4)*(q_ext_extended[j - 1] - q_ext_extended[j + 1])^2
        beta2 = (13/12)*(q_ext_extended[j] - 2*q_ext_extended[j + 1] + q_ext_extended[j + 2])^2 +
                (1/4)*(3*q_ext_extended[j] - 4*q_ext_extended[j + 1] + q_ext_extended[j + 2])^2

        gamma0, gamma1, gamma2 = 1/10, 6/10, 3/10
        epsilon = 1e-6

        alpha0 = gamma0 / (epsilon + beta0)^2
        alpha1 = gamma1 / (epsilon + beta1)^2
        alpha2 = gamma2 / (epsilon + beta2)^2
        alpha_sum = alpha0 + alpha1 + alpha2

        w0 = alpha0 / alpha_sum
        w1 = alpha1 / alpha_sum
        w2 = alpha2 / alpha_sum

        p0 = (1/3)*q_ext_extended[j - 2] - (7/6)*q_ext_extended[j - 1] + (11/6)*q_ext_extended[j]
        p1 = (-1/6)*q_ext_extended[j - 1] + (5/6)*q_ext_extended[j] + (1/3)*q_ext_extended[j + 1]
        p2 = (1/3)*q_ext_extended[j] + (5/6)*q_ext_extended[j + 1] - (1/6)*q_ext_extended[j + 2]

        qL[i] = w0*p0 + w1*p1 + w2*p2

        # 右侧重构（qR）
        beta0 = (13/12)*(q_ext_extended[j + 3] - 2*q_ext_extended[j + 2] + q_ext_extended[j + 1])^2 +
                (1/4)*(q_ext_extended[j + 3] - 4*q_ext_extended[j + 2] + 3*q_ext_extended[j + 1])^2
        beta1 = (13/12)*(q_ext_extended[j + 2] - 2*q_ext_extended[j + 1] + q_ext_extended[j])^2 +
                (1/4)*(q_ext_extended[j + 2] - q_ext_extended[j])^2
        beta2 = (13/12)*(q_ext_extended[j + 1] - 2*q_ext_extended[j] + q_ext_extended[j - 1])^2 +
                (1/4)*(3*q_ext_extended[j + 1] - 4*q_ext_extended[j] + q_ext_extended[j - 1])^2

        alpha0 = gamma0 / (epsilon + beta0)^2
        alpha1 = gamma1 / (epsilon + beta1)^2
        alpha2 = gamma2 / (epsilon + beta2)^2
        alpha_sum = alpha0 + alpha1 + alpha2

        w0 = alpha0 / alpha_sum
        w1 = alpha1 / alpha_sum
        w2 = alpha2 / alpha_sum

        p0 = (1/3)*q_ext_extended[j + 3] - (7/6)*q_ext_extended[j + 2] + (11/6)*q_ext_extended[j + 1]
        p1 = (-1/6)*q_ext_extended[j + 2] + (5/6)*q_ext_extended[j + 1] + (1/3)*q_ext_extended[j]
        p2 = (1/3)*q_ext_extended[j + 1] + (5/6)*q_ext_extended[j] - (1/6)*q_ext_extended[j - 1]

        qR[i] = w0*p0 + w1*p1 + w2*p2
    end

    return qL, qR
end






end # module Reconstruction
