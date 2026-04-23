using JLD2
using PyCall
using PyPlot

plt = pyimport("matplotlib.pyplot")
pyimport("scienceplots")
plt.style.use(["science"])

const SIGMA    = 2.0
const L_EXTRA  = 512
const L_CG     = 1024
const FIG_DIR  = joinpath("figures", "benchmark")
mkpath(FIG_DIR)

fig, ax = subplots(figsize=(5, 4))
ax.set_title("\$\\sigma = $SIGMA\$")
ax.set_xlabel("\$|\\Delta h|\$")
ax.set_ylabel("\$P(|\\Delta h|)\$")
ax.set_yscale("log")

# --- extra 数据，L=512 ---
fname_extra = "data_extra/L$(L_EXTRA)_s$(SIGMA)_extra.jld2"
if isfile(fname_extra)
    d  = load(fname_extra)
    bc = d["bin_centers"]
    hc = d["hist_counts"]
    bw = d["BIN_WIDTH"]
    prob = hc ./ (sum(hc) * bw)
    mask = prob .> 1e-12
    ax.plot(bc[mask], prob[mask], "-", lw=1.2, label="run\\_measure\\_extra (\$L=$L_EXTRA\$)")
else
    println("警告: $fname_extra 不存在")
end

# --- cg 数据，L=1024，b=1（level 1）---
fname_cg = "data_cg/L$(L_CG)_s$(SIGMA)_cg.jld2"
if isfile(fname_cg)
    d   = load(fname_cg)
    bc  = d["bin_centers"]
    hc  = d["hist"][1]      # level 1 = b=1
    bw  = d["BIN_WIDTH"]
    prob = hc ./ (sum(hc) * bw)
    mask = prob .> 1e-12
    ax.plot(bc[mask], prob[mask], "--", lw=1.2, label="run\\_coarse\\_grain \$b=1\$ (\$L=$L_CG\$)")
else
    println("警告: $fname_cg 不存在")
end

ax.legend()
fig.tight_layout()
outpath = joinpath(FIG_DIR, "histogram_bench.pdf")
fig.savefig(outpath)
plt.close(fig)
println("已保存 $outpath")
