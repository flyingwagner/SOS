using JLD2
using Glob

const DATA_DIR  = "data"
const OUT_FILE  = "sos_data.jld2"

# 扫描 data/ 下所有 L*_s*.jld2 文件
files = glob("L*_s*.jld2", DATA_DIR)
if isempty(files)
    error("data/ 下没有找到任何数据文件")
end

results = Dict{String, Any}()
L_set     = Set{Int}()
sigma_set = Set{Float64}()

for f in files
    d = load(f)
    L     = d["L"]
    sigma = d["sigma"]
    key   = "L$(L)_s$(sigma)"
    results[key] = d
    push!(L_set, L)
    push!(sigma_set, sigma)
    println("已加载 $f  →  $key")
end

L_list     = sort(collect(L_set), rev=true)
sigma_list = sort(collect(sigma_set))

# 从任意一条记录取公共参数
sample = first(values(results))
results["_meta"] = Dict(
    "L_list"     => L_list,
    "sigma_list" => sigma_list,
    "K"          => sample["K"],
)

jldsave(OUT_FILE; results)
println("\n合并完成：$(length(results)-1) 条记录 → $OUT_FILE")
println("  L_list     = $L_list")
println("  sigma_list = $sigma_list")
