

# F and F_e are given, then we gonna solve F_p:
# gamma_dot_0

using LinearAlgebra

function calc_deviatoric(M::Matrix{<:Number})
    return M .- LinearAlgebra.tr(M)/3
end

mutable struct J2Plasticity_linearHardeing
    n :: Number
    gamma_dot_0 :: Number
    h0 :: Number
    gamma_p :: Number
    xi :: Number
end


gamma_dot_0 = 1e-3
n = 20
gamma_dot_p = 0.0
dt = 1
gamma_p = 1e-3
h0 = 50
xi_0 = 10

F = Matrix{Float64}(I, 3, 3)
F_p = Matrix{Float64}(I, 3, 3)
F_p[1,2] = 0.5
I_3_3 = Matrix{Float64}(I, 3, 3)

res_F_p_vec = zeros(9)
F_p_vec = F_p[:]

function Res_Fp!(res_F_p_vec,F_p_vec, gamma_dot_p, gamma_p, F, I_3_3, h0, xi_0, F_p)
    F_p_next = reshape(F_p_vec, 3,3)
    gamma_p_next = gamma_dot_p*dt + gamma_p
    xi_next = xi_0 + h0*gamma_p_next
    F_e = F*inv(F_p_next)
    E_e = 0.5*(transpose(F_e)*F_e .- I_3_3)
    S_e = 25*E_e
    S_e_dev = calc_deviatoric(S_e)
    norm_S_e_dev = norm(S_e_dev)
    gamma_dot_p = gamma_dot_0 * (norm_S_e_dev/xi_next)^n
    L_p = gamma_dot_p * S_e_dev/norm_S_e_dev
    F_dot_p = L_p * F_p_next
    res_F_p = F_p_next .- F_dot_p * dt .- F_p
    res_F_p_vec .= res_F_p[:]
    return  res_F_p_vec
end

#Res_Fp!(res_F_p_vec,F_p_vec, gamma_dot_p, gamma_p, F, I_3_3, h0, xi_0, F_p)

Res_Fp_nsoli!(res_F_p_vec, F_p_vec) = Res_Fp!(res_F_p_vec,F_p_vec, gamma_dot_p, gamma_p, F, I_3_3, h0, xi_0, F_p)

#Res_Fp_nsoli!(F_p_vec, res_F_p_vec)

krylov_dims = 4
sol = nsoli(Res_Fp_nsoli!, F_p_vec, F_p_vec, zeros(length(F_p_vec), krylov_dims))
println()
@show sol.history

