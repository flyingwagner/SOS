using JLD2
using Printf
using Statistics

const DATA_FILE = "sos_data.jld2"
const OUT_DIR   = "fit_results"
mkpath(OUT_DIR)

data       = load(DATA_FILE, "results")
meta       = data["_meta"]
L_list_all = sort(meta["L_list"])   # 升序：[4,8,32,64,128,256,512]
sigma_list = meta["sigma_list"]

# L_min 起点序列，依次忽略小尺寸
L_min_list = [4, 8, 16, 32, 64]

get_W2(L, s)   = data["L$(L)_s$(s)"]["W2_mean"],   data["L$(L)_s$(s)"]["W2_err"]
get_chi(L, s)  = data["L$(L)_s$(s)"]["chi_mean"],  data["L$(L)_s$(s)"]["chi_err"]
get_chip(L, s) = data["L$(L)_s$(s)"]["chip_mean"], data["L$(L)_s$(s)"]["chip_err"]

# -------------------------------------------------------
# 加权线性回归：y = a*x + b，权重 w = 1/err^2
# 返回 (a, b, a_err, b_err, R2, chi2_red)
# -------------------------------------------------------
function weighted_linreg(x, y, yerr)
    w   = 1.0 ./ yerr.^2
    W   = sum(w)
    Wx  = sum(w .* x)
    Wy  = sum(w .* y)
    Wxx = sum(w .* x.^2)
    Wxy = sum(w .* x .* y)
    D   = W * Wxx - Wx^2

    a = (W * Wxy - Wx * Wy) / D
    b = (Wxx * Wy - Wx * Wxy) / D

    a_err = sqrt(W  / D)
    b_err = sqrt(Wxx / D)

    y_pred = a .* x .+ b
    ss_res = sum(((y .- y_pred) ./ yerr).^2)
    ss_tot = sum(w .* (y .- (Wy/W)).^2)
    R2     = 1.0 - sum((y .- y_pred).^2) / sum((y .- mean(y)).^2)
    n      = length(x)
    chi2_red = (n > 2) ? ss_res / (n - 2) : NaN

    return a, b, a_err, b_err, R2, chi2_red
end

# -------------------------------------------------------
# 写表头
# -------------------------------------------------------
function write_header_log(io)
    @printf(io, "%-8s  %-6s  %-5s  %10s  %10s  %10s  %10s  %8s  %10s\n",
        "sigma", "L_min", "N_pts",
        "a(slope)", "a_err", "b(intercept)", "b_err",
        "R2", "chi2_red")
    println(io, repeat("-", 90))
end

function write_header_pow(io)
    @printf(io, "%-8s  %-6s  %-5s  %10s  %10s  %10s  %10s  %8s  %10s\n",
        "sigma", "L_min", "N_pts",
        "alpha", "alpha_err", "log_A", "log_A_err",
        "R2", "chi2_red")
    println(io, repeat("-", 90))
end

# -------------------------------------------------------
# 对三个量分别做拟合，每个 sigma 一个文件
# -------------------------------------------------------
for sigma in sigma_list
    fname = joinpath(OUT_DIR, "fit_sigma$(sigma).dat")
    open(fname, "w") do io
        println(io, "# Fitting results for sigma = $sigma")
        println(io, "# All fits go from L_min to L=512")
        println(io)

        # ---------- M² ~ a·ln(L) + b ----------
        println(io, "## M2 ~ a*ln(L) + b")
        write_header_log(io)
        for L_min in L_min_list
            Ls = filter(L -> L >= L_min, L_list_all)
            length(Ls) < 3 && continue
            W2v  = [get_W2(L, sigma)[1] for L in Ls]
            W2e  = [get_W2(L, sigma)[2] for L in Ls]
            a, b, ae, be, R2, chi2 = weighted_linreg(log.(Float64.(Ls)), W2v, W2e)
            @printf(io, "%-8.1f  %-6d  %-5d  %10.5f  %10.5f  %10.5f  %10.5f  %8.5f  %10.5f\n",
                sigma, L_min, length(Ls), a, ae, b, be, R2, chi2)
        end
        println(io)

        # ---------- χ ~ A·L^α ----------
        println(io, "## chi ~ A * L^alpha  (log-log fit)")
        write_header_pow(io)
        for L_min in L_min_list
            Ls = filter(L -> L >= L_min, L_list_all)
            length(Ls) < 3 && continue
            cv  = [get_chi(L, sigma)[1] for L in Ls]
            ce  = [get_chi(L, sigma)[2] for L in Ls]
            # 误差传播：d(ln y) ≈ dy/y
            loge = ce ./ abs.(cv)
            a, b, ae, be, R2, chi2 = weighted_linreg(log.(Float64.(Ls)), log.(cv), loge)
            @printf(io, "%-8.1f  %-6d  %-5d  %10.5f  %10.5f  %10.5f  %10.5f  %8.5f  %10.5f\n",
                sigma, L_min, length(Ls), a, ae, b, be, R2, chi2)
        end
        println(io)

        # ---------- χ' ~ A·L^α ----------
        println(io, "## chi_prime ~ A * L^alpha  (log-log fit)")
        write_header_pow(io)
        for L_min in L_min_list
            Ls = filter(L -> L >= L_min, L_list_all)
            length(Ls) < 3 && continue
            cpv = [get_chip(L, sigma)[1] for L in Ls]
            cpe = [get_chip(L, sigma)[2] for L in Ls]
            loge = cpe ./ abs.(cpv)
            a, b, ae, be, R2, chi2 = weighted_linreg(log.(Float64.(Ls)), log.(cpv), loge)
            @printf(io, "%-8.1f  %-6d  %-5d  %10.5f  %10.5f  %10.5f  %10.5f  %8.5f  %10.5f\n",
                sigma, L_min, length(Ls), a, ae, b, be, R2, chi2)
        end
        println(io)
    end
    println("已保存 $fname")
end

println("\n全部完成，结果在 $OUT_DIR/")
