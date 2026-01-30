# ==============================================================================
# FD_Evolution_Core.jl
# Time-Domain Evolution Solver 
#
# Reference:
#  - [QNM families: classification and competition](https://arxiv.org/pdf/2510.02033.)
# ==============================================================================

using DelimitedFiles
using MultiFloats
using Printf
using Dates
using Plots

MultiFloats.use_bigfloat_transcendentals()

# ==============================================================================
# Configuration & Constants
# ==============================================================================

# Fixed Precision: Float64x4 (Quad-Double, ~63 digits)
const MF = Float64x4
const BF_BITS = 256

# Physical Constants
const M_BH   = BigFloat("1.0")
const MLO    = BigFloat("1.908396946564885496183206106870229007633587786259541984732824427") # 1/(M-l0/2)
const ALPHA  = BigFloat("7.0")
const L_NUM  = 6
const S_SPIN = 1

# Grid Constants
const N_STEPS   = 280000
const J_OFFSET  = 30800
const DT        = BigFloat("0.002")
const DR        = BigFloat("0.004")

# Gaussian Wave Packet Constants
const P_A = BigFloat("78.0") # Center position
const P_B = BigFloat("2.0")  # Wave Packet Width
const P_C = BigFloat("1.0")  # Amplitude

# Paths
const ROOT_DIR = @__DIR__
const G_FILE   = joinpath(ROOT_DIR, "r(rs)_grid_M1_L6_a700_lo952.csv") # Precomputed r(r*) grid via numerical integration

# ==============================================================================
# Kernels
# ==============================================================================

mf(x) = MF(x)

function load_g(fp::AbstractString)
    if !isfile(fp); error("Grid file not found"); end
    setprecision(BigFloat, BF_BITS)
    raw = readdlm(fp, ',', String; header=false)
    v = vec(raw)
    g = Vector{BigFloat}(undef, length(v))
    for i in eachindex(v)
        g[i] = parse(BigFloat, v[i])
    end
    println("Grid loaded: $(length(g)) pts")
    return g
end

function gen_V(rg::Vector{BigFloat})
    N = length(rg)
    setprecision(BigFloat, BF_BITS)
    V_b = Vector{BigFloat}(undef, N)
    l_v = BigFloat(L_NUM)
    s_v = BigFloat(S_SPIN)
    d_s2 = (S_SPIN == 2) ? BigFloat(1.0) : BigFloat(0.0)
    
    @inbounds for i in 1:N
        r = rg[i]
        if abs(r) < 1e-100; V_b[i] = 0; continue; end
        
        e_mr = exp(MLO * r)
        e_n2mr = exp(-2 * MLO * r)
        
        t1 = e_n2mr * (e_mr * (r - 2*M_BH) + r * ALPHA)
        t2 = e_mr * (l_v*(l_v+1)*r - 2*M_BH*(s_v-1) - 4*M_BH*d_s2) + MLO * r^2 * ALPHA * ((s_v-1) + MLO*r*d_s2)
        V_b[i] = t1 * t2 / (r^4)
    end
    return [mf(v) for v in V_b]
end

function init_psi(N::Int)
    setprecision(BigFloat, BF_BITS)
    p_b = Vector{BigFloat}(undef, N)
    i_bf = BigFloat(N_STEPS)
    j_bf = BigFloat(J_OFFSET)
    
    @inbounds for j in 1:N
        c_idx = (-(i_bf - j_bf) + (BigFloat(j) - 1))
        arg = DR * c_idx - P_A
        p_b[j] = P_C * exp(-(arg^2) / (2 * P_B^2))
    end
    
    curr = [mf(x) for x in p_b]
    prev = zeros(MF, N)
    next = zeros(MF, N)
    return prev, curr, next
end

function core_step!(nex::Vector{MF}, cur::Vector{MF}, pre::Vector{MF}, V::Vector{MF}, l2::MF, d2::MF)
    N = length(cur)
    TWO = MF(2.0)
    # Single-threaded loop
    for j in 2:(N-1)
        @inbounds begin
            nex[j] = -pre[j] + TWO*cur[j] - (V[j] * d2 * cur[j]) + l2 * (cur[j-1] - TWO*cur[j] + cur[j+1])
        end
    end
end

# ==============================================================================
# Driver
# ==============================================================================

function exec_run()
    println("--- Solver Start (Single Thread) ---")
    println("Prec: $(MF)")
    
    cf = DT / DR
    if cf >= 1.0; println("WARN: CFL >= 1.0"); end

    rg = load_g(G_FILE)
    V = gen_V(rg)
    p_pre, p_cur, p_nex = init_psi(length(rg))
    
    oidx = N_STEPS
    if oidx > length(rg); oidx = div(length(rg), 2); end
    
    wform = Vector{MF}(undef, N_STEPS)
    wform[1] = p_cur[oidx]

    dt_mf = mf(DT)
    dr_mf = mf(DR)
    dt2 = dt_mf^2
    lam2 = dt2 / (dr_mf^2)
    
    b_pre = p_pre
    b_cur = p_cur
    b_nex = p_nex
    
    t0 = now()
    pi = div(N_STEPS, 100)
    
    println("Evolving...")
    for t in 1:(N_STEPS-1)
        core_step!(b_nex, b_cur, b_pre, V, lam2, dt2)
        
        wform[t+1] = b_nex[oidx]
        
        tmp = b_pre; b_pre = b_cur; b_cur = b_nex; b_nex = tmp
        
        if t % pi == 0
            el = now() - t0
            @printf("\r%.1f%% | T: %s | L: %.4f", t/N_STEPS*100, el, log10(abs(Float64(wform[t+1])) + 1e-100))
        end
    end
    println("\nDone. T: ", now() - t0)
    return wform
end

res = exec_run()

function viz(w::Vector{MF})
    println("Plotting...")
    y = [log10(abs(Float64(x)) + 1e-100) for x in w]
    t = [Float64(DT) * (i-1) for i in 1:length(y)]
    
    s_i = Int(round(0.07 * length(y)))
    e_i = Int(round(0.99 * length(y)))
    
    p = plot(t[s_i:e_i], y[s_i:e_i],
             xlabel="t", ylabel="log10|psi|", legend=false,
             lw=1.5, color=:blue, grid=true)
    display(p)
    return p
end

p = viz(res)

f_out = joinpath(ROOT_DIR, "TD_M1_L6_a700_lo952.pdf")
savefig(f_out)
println("Fig saved: $f_out")

d_out = joinpath(ROOT_DIR, "TD_data_M1L6S1_Echo.csv")
writedlm(d_out, string.(res))
println("Data saved: $d_out")