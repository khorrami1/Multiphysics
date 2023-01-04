
using Tensors
using Ferrite
using Rotations
using SIAMFANLEquations


struct Phenomenological_power_low{T}
    C_11 :: T # Elastic modulii
    C_12 :: T # Elastic modulii
    C_44 :: T # Elastic modulii
    γ_dot_0 :: T # Reference shear rate
    m :: T # Rate sensitivity
    τ_hat_0 :: T # Initial slip hardness
    τ_hat_sat :: T # Saturation slip hardness
    h_0 :: T # Reference hardening rate
    a :: T # Hardening exponent 
    q_lat :: T # Latent Hardening
    num_slip :: Int
end


mat_Al = Phenomenological_power_low{Float64}(106.75*1e3, 60.41*1e3, 28.34*1e3,
    0.001, 20, 31.0, 63.0, 75.0, 2.25, 1.4, 12)

# F and F_e are given, then we gonna solve F_p:
# gamma_dot_0

# E = 200e9
# ν = 0.3
# dim = 3
# λ = E*ν / ((1 + ν) * (1 - 2ν))
# μ = E / (2(1 + ν))
# δ(i,j) = i == j ? 1.0 : 0.0

# function get_func_elasticity(λ, μ)
#     (i,j,k,l) -> λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
# end

# func_elas = get_func_elasticity(λ::Number, μ::Number)

# units: Stress: [MPa], Length: [μm], Force: [μN]
C_11 = mat_Al.C_11 # C_22 = C_33
C_12 = mat_Al.C_12
C_44 = mat_Al.C_44
get_C_qubic(C_11::Number, C_12::Number, C_44::Number) = [C_11 C_12 C_12 0 0 0; C_12 C_11 C_12 0 0 0 ; C_12 C_12 C_11 0 0 0; 0 0 0 C_44 0 0;  0 0 0 0 C_44 0;  0 0 0 0 0 C_44]

get_C_qubic(C_11, C_12, C_44)

function get_C_qubic_from_ijkl(C_11::Number, C_12::Number, C_44::Number)
    return function (i,j,k,l)
        if i==j
            if k==l
                if k==i
                    C_11
                else
                    C_12
                end
            else
                0.0
            end
        else
            if i==k && j==l
                C_44
            else
                0.0
            end
        end
    end
end



func_elas = get_C_qubic_from_ijkl(C_11, C_12, C_44)

C = SymmetricTensor{4, 3}(func_elas)

q = rand(QuatRotation)
axis_angle = AngleAxis(q)

axis_angle_vec3d = Vec((axis_angle.axis_x, axis_angle.axis_y, axis_angle.axis_z ))

C_elas = Tensors.rotate(C, axis_angle_vec3d, axis_angle.theta )

# E_e = zero(Tensor{2, 3})
# E_e[1,2] = 1e-1

# E_e = Tensor{2, 3}((i,j) -> i == 1 && j == 2 ? 0.01 : 0.0)

# S_e = C_rotated ⊡ E_e

# gamma_dot_0 = 1e-3
# n = 20
# gamma_dot_p = 0.0
# dt = 1
# gamma_p = 1e-3
# h0 = 50.0*1000
# xi_0 = 120.0

# struct MaterialState_isoJ2{T, S<:SecondOrderTensor{3, T}}
#     F_p :: S
#     F_e :: S
#     S_e :: S
#     L_p :: S
#     gamma_p :: T
#     xi :: T
# end

# function MaterialState_isoJ2()
#     return MaterialState_isoJ2(
#         one(Tensor{2, 3}),
#         one(Tensor{2, 3}),
#         zero(Tensor{2, 3}),
#         one(Tensor{2, 3}), 0.0, 10.0
#     )
# end

# num_q = 4
# num_cell = 2

# matstates_old = [MaterialState_isoJ2() for _ in 1:num_q, _ in 1:num_cell]
# matstates_new = [MaterialState_isoJ2() for _ in 1:num_q, _ in 1:num_cell]

# matstate = @view matstates_new[:,2]

# matstate[1] = MaterialState_isoJ2()

# matstates_old = matstates_new

# F = Matrix{Float64}(I, 3, 3)
# F_p = Matrix{Float64}(I, 3, 3)
# F_p[1,2] = 0.5
# I_3_3 = Matrix{Float64}(I, 3, 3)

# F = one(Tensor{2,3})
# F_p = Tensor{2, 3}((i,j) -> (i == 1 && j == 2) ? 0.005 : 0.0) + F

# res_F_p_vec = zeros(9)
# F_p_vec = collect(F_p.data)

# function Res_Fp!(res_F_p_vec,F_p_vec, gamma_dot_p, gamma_p, F, h0, xi_0, n, F_p)
#     F_p_next = Tensor{2,3}(F_p_vec)
#     gamma_p_next = gamma_dot_p*dt + gamma_p
#     xi_next = xi_0 + h0*gamma_p_next
#     F_e = F ⋅ inv(F_p_next)
#     E_e = 0.5*(tdot(F_e) - one(Tensor{2,3}))
#     S_e = C_rotated ⊡ E_e
#     S_e_dev = dev(S_e)
#     norm_S_e_dev = norm(S_e_dev)
#     gamma_dot_p = gamma_dot_0 * (norm_S_e_dev/xi_next)^n
#     L_p = gamma_dot_p * S_e_dev/norm_S_e_dev
#     F_dot_p = L_p ⋅ F_p_next
#     res_F_p = F_p_next .- F_dot_p * dt .- F_p
#     res_F_p_vec .= res_F_p[:]
#     return  res_F_p_vec
# end

#Res_Fp!(res_F_p_vec,F_p_vec, gamma_dot_p, gamma_p, F, I_3_3, h0, xi_0, F_p)

# Res_Fp_nsoli!(res_F_p_vec, F_p_vec) = Res_Fp!(res_F_p_vec,F_p_vec, gamma_dot_p, gamma_p, F, h0, xi_0, n, F_p)

#Res_Fp_nsoli!(F_p_vec, res_F_p_vec)

# krylov_dims = 4
# sol = nsoli(Res_Fp_nsoli!, F_p_vec, F_p_vec, zeros(length(F_p_vec), krylov_dims))
# println()
# @show sol.history



# Plastic Flow and Hardening at a Material Point

function get_FCC_slip_system()
    s_1 = Vec((1, -1, 0))/sqrt(2)
    n_1 = Vec((1, 1, 1))/sqrt(3)
    s_2 = Vec((1, 0, -1))/sqrt(2)
    n_2 = Vec((1, 1, 1))/sqrt(3)
    s_3 = Vec((0, 1, -1))/sqrt(2)
    n_3 = Vec((1, 1, 1))/sqrt(3)
    s_4 = Vec((1, 1, 0))/sqrt(2)
    n_4 = Vec((1, -1, -1))/sqrt(3)
    s_5 = Vec((1, 0, 1))/sqrt(2)
    n_5 = Vec((1, -1, -1))/sqrt(3)
    s_6 = Vec((0, 1, -1))/sqrt(2)
    n_6 = Vec((1, -1, -1))/sqrt(3)
    s_7 = Vec((1, 1, 0))/sqrt(2)
    n_7 = Vec((1, -1, 1))/sqrt(3)
    s_8 = Vec((1, 0, -1))/sqrt(2)
    n_8 = Vec((1, -1, 1))/sqrt(3)
    s_9 = Vec((0, 1, 1))/sqrt(2)
    n_9 = Vec((1, -1, 1))/sqrt(3)
    s_10 = Vec((1, -1, 0))/sqrt(2)
    n_10 = Vec((-1, -1, 1))/sqrt(3)
    s_11 = Vec((1, 0, 1))/sqrt(2)
    n_11 = Vec((-1, -1, 1))/sqrt(3)
    s_12 = Vec((0, 1, 1))/sqrt(2)
    n_12 = Vec((-1, -1, 1))/sqrt(3)
    return [s_1, s_2, s_3, s_4, s_5, s_6, s_7, s_8, s_9, s_10, s_11, s_12],
     [n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8, n_9, n_10, n_11, n_12]
end

function get_FCC_slip_system_DAMASK()
    s_1 = Vec((0, 1, -1))/sqrt(2)
    n_1 = Vec((1, 1, 1))/sqrt(3)
    s_2 = Vec((-1, 0, 1))/sqrt(2)
    n_2 = Vec((1, 1, 1))/sqrt(3)
    s_3 = Vec((1, -1, 0))/sqrt(2)
    n_3 = Vec((1, 1, 1))/sqrt(3)
    s_4 = Vec((0, -1, -1))/sqrt(2)
    n_4 = Vec((-1, -1, 1))/sqrt(3)
    s_5 = Vec((1, 0, 1))/sqrt(2)
    n_5 = Vec((-1, -1, 1))/sqrt(3)
    s_6 = Vec((-1, 1, 0))/sqrt(2)
    n_6 = Vec((-1, -1, 1))/sqrt(3)
    s_7 = Vec((0, -1, 1))/sqrt(2)
    n_7 = Vec((1, -1, -1))/sqrt(3)
    s_8 = Vec((-1, 0, -1))/sqrt(2)
    n_8 = Vec((1, -1, -1))/sqrt(3)
    s_9 = Vec((1, 1, 0))/sqrt(2)
    n_9 = Vec((1, -1, -1))/sqrt(3)
    s_10 = Vec((0, 1, 1))/sqrt(2)
    n_10 = Vec((-1, 1, -1))/sqrt(3)
    s_11 = Vec((1, 0, -1))/sqrt(2)
    n_11 = Vec((-1, 1, -1))/sqrt(3)
    s_12 = Vec((-1, -1, 0))/sqrt(2)
    n_12 = Vec((-1, 1, -1))/sqrt(3)
    return [s_1, s_2, s_3, s_4, s_5, s_6, s_7, s_8, s_9, s_10, s_11, s_12],
     [n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8, n_9, n_10, n_11, n_12]
end

s_slip, n_slip = get_FCC_slip_system_DAMASK()

Schmid_tensors = [s_slip[a] ⊗ n_slip[a] for a in 1:mat_Al.num_slip]
symmetric_Schmid = symmetric.(Schmid_tensors)

s_a = zero(Tensor{1,3})
n_a = zero(Tensor{1,3})
S = zero(Tensor{2,3})

F_e = one(Tensor{2,3})

C_e = F_e' ⋅ F_e
# τ = [C_e ⊡ (S ⋅ symmetric_Schmid[a]) for a in 1:num_slip]

# τ_hat = [mat_Al.τ_hat_0 for a in 1:num_slip]

# h = [mat_Al.h_0*(1-τ_hat[a]/mat_Al.τ_hat_sat)^mat_Al.a for a in 1:num_slip]

function get_hardening_matrix(q::Number)
    A = Matrix([1.0 0 0 ; 0 1.0 0; 0 0 1.0])
    qA = Matrix([1.0 0 0 ; 0 1.0 0; 0 0 1.0])*q
    return [A qA qA qA ; qA A qA qA ; qA qA A qA; qA qA qA A]
end

H_hardening = get_hardening_matrix(mat_Al.q_lat) 

# γ_dot = [mat_Al.γ_dot_0 for _ in 1:num_slip]

# τ_hat_dot = [sum([H_hardening[a,b]*norm(γ_dot[b]) for b in 1:num_slip]) for a in 1:num_slip]

# Δt = 1e-3

# τ_hat_new = τ_hat .+ τ_hat_dot .* Δt

# γ_dot = [γ_dot_0 * (norm(τ[a]/τ_hat_new[a]))^mat_Al.m * sign(τ[a]) for a in 1:num_slip]

# L_p = sum([γ_dot[a] * s_slip[a] ⊗ n_slip[a] for a in 1:num_slip])

# L_p

# F_p = one(Tensor{2,3}) + rand(Tensor{2,3})
# F_p = F_p / det(F_p)^(1/3) 
# det(F_p)


# F = one(Tensor{2,3})
# # F_p = Tensor{2, 3}((i,j) -> (i == 1 && j == 2) ? 0.005 : 0.0) + F
# F_p = one(Tensor{2,3})

# res_F_p_vec = zeros(9)
# F_p_vec = collect(F_p.data)

struct material_state_pheno_power{T, S<:SecondOrderTensor{3, T}}
    F_e :: S
    F_p :: S
    τ_hat :: Vector{<:T}
    γ_dot :: Vector{<:T}
end

function material_state_pheno_power(τ_hat, γ_dot)
    return material_state_pheno_power(one(Tensor{2, 3}), one(Tensor{2, 3}), τ_hat, γ_dot)
end

τ_hat_0 = [mat_Al.τ_hat_0 for a in 1:mat_Al.num_slip]
γ_dot_0 = [mat_Al.γ_dot_0 for _ in 1:mat_Al.num_slip]

mat_state_old = material_state_pheno_power(τ_hat_0, γ_dot_0)
mat_state_new = material_state_pheno_power(τ_hat_0, γ_dot_0)

function Residual_Lp!(residual_vec, F_p_vec, F_p_0, F, F_e, τ_hat_old, τ_hat, γ_dot, C_elas, symmetric_Schmid, mat, Δt)

    F_p = Tensor{2,3}(F_p_vec)
    # F_p = F_p / (det(F_p))^(1/3) # what if det(F_p)<0, you can't ^(1/3)
    inv_F_p = inv(F_p)
    F_e = F ⋅ inv_F_p
    C_e = F_e' ⋅ F_e
    S = 0.5*C_elas ⊡ (C_e - one(Tensor{2,3}))
    τ = [C_e ⊡ (S ⋅ symmetric_Schmid[a]) for a in 1:mat.num_slip]
    h = [mat.h_0*(1-τ_hat[a]/mat.τ_hat_sat)^mat.a for a in 1:mat.num_slip]
    τ_hat_dot = [sum([H_hardening[a,b]*h[b]*norm(γ_dot[b]) for b in 1:mat.num_slip]) for a in 1:mat.num_slip]
    τ_hat .= τ_hat_old .+ τ_hat_dot .* Δt

    γ_dot .= [mat.γ_dot_0 * (abs(τ[a]/τ_hat[a]))^mat_Al.m * sign(τ[a]) for a in 1:mat.num_slip]

    L_p = sum([γ_dot[a] * Schmid_tensors[a] for a in 1:mat.num_slip])
    residual = F_p - F_p_0 - L_p ⋅ F_p *Δt
    # residual = L_p - (one(Tensor{2,3}) - F_p ⋅ inv_F_p)/Δt
    residual_vec .= collect(residual)[:]

end

F = Tensor{2, 3}((i,j) -> (i == 1 && j == 2) ? 0.005 : 0.0)
F = F + one(Tensor{2,3})
F_e = one(Tensor{2,3})
F_p_0 = one(Tensor{2,3})
# F_p_0 = Tensor{2,3}(RotMatrix(q))
F_p_vec = collect(F_p_0.data)
Δt = 0.001
residual = zeros(9)

Residual_Lp!(residual, F_p_vec, F, F_e, F_p_0, τ_hat_0, τ_hat, γ_dot_0, C_elas, symmetric_Schmid, mat_Al, Δt)

Res_Lp_nsoli!(residual, F_p_vec) = Residual_Lp!(residual, F_p_vec, F, F_e, F_p_0, τ_hat_0, τ_hat_0, γ_dot_0, C_elas, symmetric_Schmid, mat_Al, Δt)

Res_Lp_nsoli!(residual, F_p_vec)

krylov_dims = 10
sol = nsoli(Res_Lp_nsoli!, F_p_vec, F_p_vec, zeros(length(F_p_vec), krylov_dims))
@show sol.history

a = [F_e]

b = [F_p_0]

b .= a
