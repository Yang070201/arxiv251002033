# ==============================================================================
# MP_QNFextraction_ky.jl
# QNM & Amplitude Extractor 
#
# Reference:
#  - [QNM families: classification and competition](https://arxiv.org/pdf/2510.02033.)
#  - [Mining information from binary black hole mergers: a comparison of estimation methods for complex exponentials in noise](https://arxiv.org/pdf/gr-qc/0701086)
# ==============================================================================

using DelimitedFiles, GenericLinearAlgebra, LinearAlgebra, Plots, Printf
plotlyjs()

# --- Config ---
setprecision(BigFloat, 256); const BF = BigFloat

# --- Core Functions ---
function ld_f(fp::AbstractString)
    !isfile(fp) && error("File missing: $fp")
    raw = readdlm(fp, ',', String; header=false)
    d = Vector{BF}(undef, length(raw))
    try
        for (i, s) in enumerate(vec(raw)); d[i] = parse(BF, strip(s)); end
    catch; error("Parse err"); end
    println("Data loaded. Len: $(length(d))")
    return d
end

function proc_A(v::AbstractVector, k::Integer)
    println(">> P1: SVD...")
    s_d = v[1:k:end]; N = length(s_d); L = floor(Int, N/2)
    R, C = N-L, L
    Y = Matrix{eltype(v)}(undef, R, C)
    for j in 1:C, i in 1:R; Y[i,j] = s_d[i+j-1]; end
    F = svd(Y)
    println("Max SV: $(Float64(F.S[1]))")
    p = scatter(Float64.(F.S), yaxis=:log, label="SVD", title="Cliff Plot", marker=:circle, size=(800,500))
    return s_d, F, L, p
end

function proc_B(tup::Tuple, M::Integer, dt::Number, k::Integer)
    println(">> P2: Modes (p=$M)...")
    sd, F, L, _ = tup; U, S, V = F.U, F.S, F.V
    N = length(sd); (L <= M) && @warn "Unstable: L<=M"
    rows, cols = N-L, L
    Y1 = Matrix{eltype(sd)}(undef, rows, cols)
    for j in 1:cols, i in 1:rows; Y1[i,j] = sd[i+j]; end
    iS = Diagonal([1/s for s in S[1:M]])
    Z = iS * U[:,1:M]' * Y1 * V[:,1:M]
    z_vals = eigen(ComplexF64.(Complex{BF}.(Z))).values
    dt_eff = dt * k; im_u = Complex(zero(eltype(sd)), one(eltype(sd)))
    w = (im_u / dt_eff) .* log.(z_vals)
    Vm = Matrix{Complex{eltype(dt)}}(undef, N, M)
    for m in 1:M
        zk=z_vals[m]; c=one(zk)
        for i in 1:N; Vm[i,m]=c; (i<N) && (c*=zk); end
    end
    h = Vm \ sd
    println("Done. $(length(w)) modes.")
    return collect(zip(w, h))
end

function proc_C(res::Vector, T::Number)
    println(">> P3: Analysis...")
    phy = filter(x->real(x[1])>0, res)
    isempty(phy) && return []
    srt = sort(phy, by=x->abs(imag(x[1])))
    w = [x[1] for x in srt]; h = [x[2] for x in srt]
    im_w = abs.(imag.(w)); re_w = abs.(real.(w)); h2 = abs2.(h)
    cf = 1 .- exp.(-2 .* im_w .* T)
    ec = (h2 .* (re_w.^2) ./ (im_w .+ 1e-100)) .* cf
    ef = ec ./ sum(ec)
    println("-"^60)
    @printf("%-4s | %-26s | %-11s | %-9s\n", "n", "Omega", "|h|", "Energy Fraction")
    println("-"^60)
    for i in 1:length(srt)
        @printf("%-4d | %.5e %s %.5ei | %.5e | %.5e\n", i-1, Float64(real(w[i])), imag(w[i])<0 ? "-" : "+", Float64(abs(imag(w[i]))), Float64(abs(h[i])), Float64(ef[i]))
    end
    return [(n=i-1, w=w[i], h=h[i], ef=ef[i], cf=cf[i]) for i in 1:length(srt)]
end

function exp_html(res::Vector)
    isempty(res) && return
    fs(v) = (s=@sprintf("%.45e", v); contains(s,"e") ? replace(s, "e"=>" * 10^(") * ")" : s)
    h_str = """<!DOCTYPE html><html><head><style>body{font-family:sans-serif;padding:20px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #ddd;padding:8px;text-align:left}th{background:#f2f2f2}.n{font-family:monospace}</style></head><body><h3>Results</h3><table><thead><tr><th>n</th><th>Re(w)</th><th>Im(w)</th><th>|h|</th><th>Energy Fraction</th><th>Correction Factor</th></tr></thead><tbody>"""
    for r in res
        h_str *= "<tr><td>$(r.n)</td><td class='n'>$(fs(real(r.w)))</td><td class='n'>$(fs(imag(r.w)))</td><td class='n'>$(fs(abs(r.h)))</td><td class='n'>$(fs(r.ef))</td><td class='n'>$(@sprintf("%.20f",Float64(r.cf)))</td></tr>"
    end
    h_str *= "</tbody></table></body></html>"
    p = joinpath(tempdir(), "res_$(round(Int,time())).html")
    open(f->write(f, h_str), p, "w")
    try; Sys.iswindows() ? run(`cmd /c start $p`) : run(`xdg-open $p`); catch; println("Saved: $p"); end
end

# ================= EXECUTION =================
# Path setup
const D_DIR = @__DIR__
const FP = joinpath(D_DIR, "TD_data_M1L6S1_Echo_myver.csv") 
const W = ld_f(FP) # Precomputed time-domain data from "FD_Evolution_kys.jl"

# Params
const IMAX = 280000
const DT = BF(2)/1000

# ---Sampling Window (Adjustable) ---
tt1 = ceil(Int, 0.12 * IMAX)
tt2 = ceil(Int, 1.00 * IMAX)
WIN_T = BF(tt2 - tt1) * DT

println("Plotting preview...")
# Plotting logic kept raw for adjustment
pt = (0:length(W)-1) .* Float64(DT)
py = log10.(abs.(Float64.(W)) .+ 1e-100)
# View ratio (Adjust here)
r_s = 0.07
r_e = 1.00
idx_a = max(1, round(Int, r_s*length(py))); idx_b = min(length(py), round(Int, r_e*length(py)))
plt = plot(pt[idx_a:idx_b], py[idx_a:idx_b], lab="W", lw=0.8, c=:blue, size=(800,400), grid=true)
vline!(plt, [tt1*Float64(DT), tt2*Float64(DT)], l=(:red,:dash), lab="Win")
display(plt)

# --- MP Process ---
sub_W = W[tt1:tt2]
const K_STEP = 958  # Sampling step size

# SVD
out_1 = proc_A(sub_W, K_STEP);
display(out_1[4])

# Modes
const P_ORD = 34    # Mode order
res_raw = proc_B(out_1, P_ORD, DT, K_STEP);
res_fin = proc_C(res_raw, WIN_T);

# Export
exp_html(res_fin)

# The results can be checked in Table.I in https://arxiv.org/pdf/2510.02033