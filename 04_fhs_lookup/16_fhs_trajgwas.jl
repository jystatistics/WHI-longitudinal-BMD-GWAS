# 2026-11-25. FHS look-up: TrajGWAS score test (SPA) + Wald test for the look-up SNPs, same model as the rebuilt WHI analysis.
# Usage: julia 16_fhs_trajgwas.jl <label> <outcome column> <phenotype csv> <look-up plink prefix> <eigenvec> <outdir>
#   mean:  y ~ 1 + time + AGE + BMI + PC1..PC10   (FHS file has baseline BMI, not height and weight: BMI replaces HEIGHTX + WEIGHTX;
#          no cohort term, one cohort)
#   WSV:   ~ 1 + AGE + BMI          random effects: 1 + time        AGE = AGE20, BMI = BMI20; covariates z-scored
#   Same fixed optimizer cascade and acceptance rule as 07_trajgwas.jl.
using CSV, DataFrames, SnpArrays, TrajGWAS, Statistics, Printf
const SM = first(m for (k, m) in Base.loaded_modules if k.name == "StatsModels")
term, FormulaTerm = SM.term, SM.FormulaTerm
W  = first(m for (k, m) in Base.loaded_modules if k.name == "WiSER")
IP = first(m for (k, m) in Base.loaded_modules if k.name == "Ipopt")
lab, outcome, phenofile, prefix, pcfile, outdir = ARGS
y = Symbol(outcome); mkpath(outdir)
tofloat(v) = [x === nothing ? missing : x for x in tryparse.(Float64, string.(v))]
df = CSV.read(phenofile, DataFrame); df.IID = string.(df.IID)
# TrajGWAS re-reads the .fam and may parse numeric IDs as numbers (leading zeros lost), which broke the ID check.
# Give every ID a letter prefix in the phenotype, the PCs and a relabelled copy of the look-up .fam.
println("IIDs all digits: $(all(id -> all(isdigit, id), df.IID)); with a leading zero: $(count(id -> startswith(id, "0"), unique(df.IID)))")
df.IID = "S" .* df.IID
n0 = nrow(df); df = unique(df, [:IID, :time]); println("duplicate (IID, time) rows dropped: $(n0 - nrow(df))")
rename!(df, :AGE20 => :AGE, :BMI20 => :BMI)
select!(df, Not(intersect(names(df), ["PC$i" for i in 1:20])))   # use the new PCs, not any old PC columns in the file
pcs = CSV.read(pcfile, DataFrame; delim = ' ', ignorerepeated = true, header = vcat(["FID", "IID"], ["PC$i" for i in 1:10]))
pcs.IID = "S" .* string.(pcs.IID); allunique(pcs.IID) || error("duplicated IIDs in $pcfile")
df = innerjoin(df, select(pcs, Not(:FID)), on = :IID)
pcsyms = [Symbol("PC$i") for i in 1:10]
vars = vcat([y, :time, :AGE, :BMI], pcsyms)
for v in vars; df[!, v] = tofloat(df[!, v]); end
df = df[completecases(df, vars), :]; disallowmissing!(df, vars)
# FHS "time" is the exam number (20, 22, 24). Shift so baseline = 0: an uncentred time column (mean ~22) next to the
# intercept made the model numerically singular (Infs/NaNs). SNP tests do not depend on this shift.
t0 = minimum(df.time); df.time = df.time .- t0; println("time = exam number - $t0")
src = prefix; prefix = joinpath(outdir, basename(src) * "_sid")
for ext in (".bed", ".bim"); isfile(prefix * ext) || symlink(src * ext, prefix * ext); end
open(prefix * ".fam", "w") do io
    for l in eachline(src * ".fam"); f = split(l); f[1] = "S" * f[2]; f[2] = f[1]; println(io, join(f, " ")); end
end
famids = [split(l)[2] for l in eachline(prefix * ".fam")]
df = filter(r -> r.IID in Set(famids), df); sort!(df, [:IID, :time])
subj = unique(df.IID); rowinds = [findfirst(==(s), famids) for s in subj]
nv = [count(==(s), df.IID) for s in subj]
println("$lab: $(length(subj)) unrelated women, $(nrow(df)) visits; visits per woman: ", join(["$k=$(count(==(k), nv))" for k in sort(unique(nv))], ", "))
for v in vars
    x = df[!, v]; @printf("  %-22s mean %9.4f  sd %8.4f  min %9.4f  max %9.4f\n", v, mean(x), std(x), minimum(x), maximum(x))
end
for v in vcat([:AGE, :BMI], pcsyms); x = df[!, v]; df[!, v] = (x .- mean(x)) ./ std(x); end
meanf = FormulaTerm(term(y), term(1) + term(:time) + term(:AGE) + term(:BMI) + sum(term.(pcsyms)))
wsvf  = FormulaTerm(term(y), term(1) + term(:AGE) + term(:BMI))
ref   = FormulaTerm(term(y), term(1) + term(:time))
statusof(txt) = (ms = collect(eachmatch(r"status = ([A-Z_]+)", txt)); isempty(ms) ? "UNKNOWN" : ms[end].captures[1])
finiteparams(m) = try all(isfinite, m.β) && all(isfinite, m.τ) catch; true end
plausible(m) = try -20 < m.τ[1] < 0 && all(abs.(m.τ[2:end]) .<= 5) catch; true end
settings = [("WiSER default, LS start", nothing, :ls), ("WiSER default, MoM start", nothing, :mom),
            ("Ipopt max_iter 100, LS start", ["print_level" => 0, "max_iter" => 100], :ls),
            ("Ipopt watchdog, LS start", ["print_level" => 0, "watchdog_shortened_iter_trigger" => 3, "max_iter" => 100], :ls),
            ("Ipopt max_iter 1000, LS start", ["print_level" => 0, "max_iter" => 1000], :ls),
            ("Ipopt max_iter 100, MoM start", ["print_level" => 0, "max_iter" => 100], :mom),
            ("Ipopt watchdog, MoM start", ["print_level" => 0, "watchdog_shortened_iter_trigger" => 3, "max_iter" => 100], :mom),
            ("Ipopt max_iter 1000, MoM start", ["print_level" => 0, "max_iter" => 1000], :mom),
            ("Ipopt mehrotra+warm start 1000, MoM start", ["print_level" => 0, "mehrotra_algorithm" => "yes", "warm_start_init_point" => "yes", "max_iter" => 1000], :mom)]
function tryfit(re, opts, ini)
    m = W.WSVarLmmModel(meanf, re, wsvf, :IID, df)
    path, io = mktemp()
    redirect_stdout(io) do
        initm = ini == :mom ? W.init_mom!(m) : W.init_ls!(m)
        opts === nothing ? W.fit!(m; init = initm) : W.fit!(m, IP.Optimizer(); init = initm, solver_config = Dict{String,Any}(opts))
    end
    close(io); txt = read(path, String); rm(path); return m, statusof(txt)
end
nm = nothing; nmlab = ""; almost = nothing; reused = ref
# random intercept + slope first (the WHI model); only if no setting converges, random intercept (logged as MODEL USED)
for (re, sl, opts, ini) in vcat([(ref, s...) for s in settings], [(FormulaTerm(term(y), term(1)), s...) for s in settings])
    nm === nothing || break
    re !== ref && almost !== nothing && break
    try
        m, st = tryfit(re, opts, ini); ok = finiteparams(m) && plausible(m); sl = sl * (re === ref ? "" : " [RE intercept only]")
        println("  -> $sl: final status $st", ok ? "" : " (rejected: non-finite/implausible)", "  tau = ", try round.(m.τ, digits = 3) catch; "?" end)
        if ok && st in ("LOCALLY_SOLVED", "OPTIMAL"); global nm = m; global nmlab = "$sl ($st)"; global reused = re
        elseif ok && st == "ALMOST_LOCALLY_SOLVED" && almost === nothing; global almost = (m, "$sl ($st)", re); end
    catch e; println("  -> $sl: fit failed: ", first(sprint(showerror, e), 200)); end
end
if nm === nothing && almost !== nothing
    global nm = almost[1]; global nmlab = almost[2]; global reused = almost[3]
end
nm === nothing && error("FHS null model did not converge under any optimizer setting")
println("MODEL USED: mean = $meanf | RE = $reused | WSV = $wsvf | optimizer: $nmlab")
open(joinpath(outdir, "$lab.null.txt"), "w") do io; println(io, "optimizer: ", nmlab); show(io, nm); end
nsnp = countlines(prefix * ".bim"); println("look-up variants in $prefix.bim: $nsnp; fam rows $(length(famids)); geneticrowinds $(length(rowinds))")
try
    trajgwas(nm, prefix; geneticrowinds = rowinds, pvalfile = joinpath(outdir, "$lab.score.txt"))
    trajgwas(nm, prefix; geneticrowinds = rowinds, pvalfile = joinpath(outdir, "$lab.wald.txt"), test = :wald, snpinds = 1:nsnp)
catch e
    println("trajgwas failed: ", first(sprint(showerror, e), 600)); rethrow()
end
println("done: $nsnp look-up variants -> $outdir/$lab.{score,wald}.txt")
