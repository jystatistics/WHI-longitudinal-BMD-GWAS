# 2026-11-24. Unified TrajGWAS for one group x outcome x chromosome.
# Usage: julia 07_trajgwas.jl <group> <outcome> <chr> <phenotype csv>
#   group   : label used in output names (black, hispanic, ai_an, ...)
#   outcome : column in the phenotype file (hipbmd, spinebmd, or any other biomarker)
# Same model for every group and outcome:
#   TIME = "time" column = BMD exam year since baseline (0, 1, 3, 6, 9 ...)
#   mean:                     y ~ 1 + TIME + AGE + HEIGHTX + WEIGHTX + cohort + PC1..PC10
#   within-subject variance:  ~ 1 + AGE + BMI (default; env WSV), decided 2026-11-25. PC1-2 were dropped because no
#     optimizer setting converged with them for Black hip or AI/AN hip (DIAG 7); TIME was dropped because with it the
#     Hispanic spine and AI/AN null models were degenerate/not converged and their tau tests invalid. PCs and TIME stay in
#     the mean model, and TIME stays in the random effects.
#   random effects: 1 + TIME (random slope) in every set (env RE=intercept for a sensitivity run only).
#   Only optimizer settings / starting values vary between sets (fixed, logged order); the model never changes.
#   AGE, HEIGHTX, WEIGHTX and the PCs are z-scored (2026-11-25): unscaled covariates overflow the log-linear
#   within-subject variance model (Ipopt INVALID_MODEL). SNP tests are unaffected; time stays in years.
# PCs come from WHI_PCA_v5_r08.eigenvec (the rebuilt merge, duplicates removed); any PC columns in the phenotype file are ignored.
using CSV, DataFrames, SnpArrays, TrajGWAS, Statistics, Printf, LinearAlgebra
# StatsModels is not in our Julia environment (only loaded as a TrajGWAS dependency); take it from
# the loaded modules instead of Pkg.add, so the working package versions are not changed.
const SM = first(m for (k, m) in Base.loaded_modules if k.name == "StatsModels")
term, FormulaTerm = SM.term, SM.FormulaTerm

group, outcome, chr, phenofile = ARGS[1], ARGS[2], parse(Int, ARGS[3]), ARGS[4]
const B = joinpath(ENV["WHI_DIR"], "10_Rebuild_v2")
timevar = Symbol(get(ENV, "TIMEVAR", "time"))
y = Symbol(outcome)
tofloat(v) = [x === nothing ? missing : x for x in tryparse.(Float64, string.(v))]

df = CSV.read(phenofile, DataFrame)
select!(df, Not(intersect(names(df), ["PC$i" for i in 1:20])))
df.IID = string.(df.IID)
pcs = CSV.read(joinpath(B, "WHI_PCA_v5_r08.eigenvec"), DataFrame; delim = ' ', ignorerepeated = true,
               header = vcat(["FID", "IID"], ["PC$i" for i in 1:20]))
pcs.IID = string.(pcs.IID)
coh = CSV.read(joinpath(B, "cohort_by_iid.txt"), DataFrame; delim = ' ', header = ["IID", "cohort"])
coh.IID = string.(coh.IID)
allunique(pcs.IID) && allunique(coh.IID) || error("duplicated IIDs in the PCA or cohort file")
# women genotyped twice under two IDs (04e): phenotype rows filed under the dropped ID go to the kept ID
aliasfile = joinpath(B, "alias_all.txt")
if isfile(aliasfile)
    alias = Dict(String(split(l)[1]) => String(split(l)[2]) for l in eachline(aliasfile) if length(split(l)) == 2)
    function resolve(id)
        n = 0
        while haskey(alias, id) && n < 5
            id = alias[id]; n += 1
        end
        return id
    end
    nmapped = count(id -> haskey(alias, id), df.IID)
    df.IID = [resolve(id) for id in df.IID]
    n0 = nrow(df); df = unique(df, [:IID, timevar])
    println("ID aliases: $nmapped phenotype rows re-labelled to the kept ID; $(n0 - nrow(df)) duplicate (IID, time) rows dropped")
end
nvisits = nrow(df)
df = innerjoin(df, select(pcs, Not(:FID)), on = :IID)
df = innerjoin(df, coh, on = :IID)
println("visits in phenotype file: $nvisits; with genotype PCs and cohort: $(nrow(df))")

pcsyms = [Symbol("PC$i") for i in 1:10]
vars = vcat([y, timevar, :AGE, :HEIGHTX, :WEIGHTX], pcsyms)
for v in vars; df[!, v] = tofloat(df[!, v]); end
df = df[completecases(df, vars), :]
disallowmissing!(df, vars)

# genotype rows in the order subjects first appear in the (sorted) data
# GENOPREFIX=<plink prefix>: use another genotype file (e.g. a single-SNP look-up); outputs then carry its name as a tag
prefix = get(ENV, "GENOPREFIX", joinpath(B, "geno", "chr$chr"))
gtag = haskey(ENV, "GENOPREFIX") ? "." * basename(ENV["GENOPREFIX"]) : ""
famids = [split(l)[2] for l in eachline(prefix * ".fam")]
allunique(famids) || error("duplicated IIDs in $prefix.fam")
df = filter(r -> r.IID in Set(famids), df)
sort!(df, [:IID, timevar])
subj = unique(df.IID)
rowinds = [findfirst(==(s), famids) for s in subj]
for v in vars
    x = df[!, v]; @printf("  %-8s mean %10.4f  sd %9.4f  min %10.4f  max %10.4f\n", v, mean(x), std(x), minimum(x), maximum(x))
end
nv = [count(==(s), df.IID) for s in subj]
println("  visits per woman: ", join(["$k=$(count(==(k), nv))" for k in sort(unique(nv))], ", "))
df.BMI = df.WEIGHTX ./ (df.HEIGHTX ./ 100) .^ 2      # kg/m2, from the same baseline height and weight
for v in vcat([:AGE, :HEIGHTX, :WEIGHTX, :BMI], pcsyms)
    x = df[!, v]; df[!, v] = (x .- mean(x)) ./ std(x)
end
percoh = unique(df[:, [:IID, :cohort]])
println("$group $outcome chr$chr: $(length(subj)) women, $(nrow(df)) visits; women per cohort: ",
        join(["$c=$(count(==(c), percoh.cohort))" for c in sort(unique(percoh.cohort))], ", "))

rhs = term(1) + term(timevar) + term(:AGE) + term(:HEIGHTX) + term(:WEIGHTX) + sum(term.(pcsyms))
length(unique(df.cohort)) > 1 && (rhs = rhs + term(:cohort))
meanf = FormulaTerm(term(y), rhs)
# Within-subject variance (WSV) model: a reduced covariate set, because WiSER does not converge when the WSV model
# carries the full mean-model list (2026-11-25 diagnostics). Default time + AGE + BMI + PC1 + PC2; override with
# WSV="time,AGE" etc. ("time" means the TIMEVAR column).
wsvars = [v == "time" ? timevar : Symbol(v) for v in split(get(ENV, "WSV", "AGE,BMI"), ",")]
wsvf = FormulaTerm(term(y), term(1) + sum(term.(wsvars)))
println("mean model: ", meanf); println("WSV model:  ", wsvf)

# DIAG=7: grid of solver settings x starting values x rescaling for the null model, compact status matrix, then stop.
if get(ENV, "DIAG", "0") == "7"
    W   = first(m for (k, m) in Base.loaded_modules if k.name == "WiSER")
    IP  = first(m for (k, m) in Base.loaded_modules if k.name == "Ipopt")
    MOI = first(m for (k, m) in Base.loaded_modules if k.name == "MathOptInterface")
    for (k, m) in Base.loaded_modules
        k.name in ("WiSER", "TrajGWAS", "Ipopt", "MathOptInterface") && println("version ", k.name, " ", try pkgversion(m) catch; "?" end)
    end
    kws = try reduce(vcat, [Base.kwarg_decl(mm) for mm in methods(W.fit!)]) catch; Symbol[] end
    println("fit! keywords: ", kws); foreach(println, methods(W.fit!))
    function mksolver(opts)
        if :solver_config in kws
            return (IP.Optimizer(), Dict{String,Any}(opts))
        end
        o = IP.Optimizer()
        A = isdefined(MOI, :RawOptimizerAttribute) ? MOI.RawOptimizerAttribute : MOI.RawParameter
        for (k, v) in opts; MOI.set(o, A(k), v); end
        return (o, nothing)
    end
    function cap(fn)
        path, io = mktemp()
        redirect_stdout(fn, io); close(io); txt = read(path, String); rm(path); return txt
    end
    code(st) = st in ("LOCALLY_SOLVED", "OPTIMAL") ? "S" : st == "ALMOST_LOCALLY_SOLVED" ? "a" :
               st == "NUMERICAL_ERROR" ? "N" : st == "ITERATION_LIMIT" ? "I" : "?"
    solvers = [("dflt", nothing),
               ("p100", ["print_level" => 0, "max_iter" => 100]),
               ("p1k",  ["print_level" => 0, "max_iter" => 1000]),
               ("wdog", ["print_level" => 0, "watchdog_shortened_iter_trigger" => 3, "max_iter" => 100]),
               ("mw1k", ["print_level" => 0, "mehrotra_algorithm" => "yes", "warm_start_init_point" => "yes", "max_iter" => 1000])]
    # rescaled copies: time in SD units (time_s) and/or the outcome z-scored. SNP tests are invariant to these.
    d0 = copy(df); d0.time_s = d0[!, timevar] ./ std(d0[!, timevar])
    ysd = std(d0[!, y]); d0.y_s = (d0[!, y] .- mean(d0[!, y])) ./ ysd
    println("outcome sd $(round(ysd, sigdigits=4)); time sd $(round(std(d0[!, timevar]), sigdigits=4))")
    gf(lhs, vs) = FormulaTerm(term(lhs), term(1) + sum(term.(vs)))
    hasc = length(unique(df.cohort)) > 1
    X = hcat(ones(nrow(df)), Matrix(df[:, vcat([timevar, :AGE, :HEIGHTX, :WEIGHTX], pcsyms)]),
             hasc ? hcat([Float64.(df.cohort .== c) for c in sort(unique(df.cohort))[2:end]]...) : zeros(nrow(df), 0))
    sv = svdvals(X)
    println("mean design: $(size(X, 2)) columns, rank $(rank(X)), condition number $(round(sv[1] / sv[end], sigdigits=3))")
    println("legend: S solved, a almost solved, N NUMERICAL_ERROR, I ITERATION_LIMIT, F exception, x non-finite estimates")
    println("columns: ", join(["$(s[1])/$(i)" for s in solvers for i in ("ls", "mom")], " "))
    for (re, wsvlab, wsvv) in (("RE 1+time", "WSV t+AGE+BMI", [:AGE, :BMI]), ("RE 1+time", "WSV t+AGE+BMI+PC1-2", [:AGE, :BMI, :PC1, :PC2]),
                               ("RE 1",      "WSV t+AGE+BMI", [:AGE, :BMI]), ("RE 1",      "WSV t+AGE+BMI+PC1-2", [:AGE, :BMI, :PC1, :PC2]))
        for (sc, tv, yv) in (("raw", timevar, y), ("time_s", :time_s, y), ("y_s", timevar, :y_s), ("both", :time_s, :y_s))
            mvars = vcat([tv, :AGE, :HEIGHTX, :WEIGHTX], pcsyms); hasc && push!(mvars, :cohort)
            mf = gf(yv, mvars); wf = gf(yv, vcat([tv], wsvv))
            rf = re == "RE 1" ? FormulaTerm(term(yv), term(1)) : gf(yv, [tv])
            cells = String[]
            for (sl, opts) in solvers, ini in (:ls, :mom)
                c = try
                    m = W.WSVarLmmModel(mf, rf, wf, :IID, d0)
                    txt = cap() do
                        initm = ini == :mom ? W.init_mom!(m) : W.init_ls!(m)
                        if opts === nothing
                            W.fit!(m; init = initm)
                        else
                            o, cfg = mksolver(opts)
                            cfg === nothing ? W.fit!(m, o; init = initm) : W.fit!(m, o; init = initm, solver_config = cfg)
                        end
                    end
                    ms = collect(eachmatch(r"status = ([A-Z_]+)", txt))
                    st = isempty(ms) ? "?" : code(ms[end].captures[1])
                    (st == "S" || st == "a") && !(all(isfinite, m.β) && all(isfinite, m.τ)) ? "x" : st
                catch e
                    global lasterr = first(sprint(showerror, e), 150); "F"
                end
                push!(cells, rpad(c, 9))
            end
            println(rpad("$group $outcome | $re | $wsvlab | $sc", 66), " ", join(cells))
        end
    end
    isdefined(Main, :lasterr) && println("last exception: ", lasterr)
    exit(0)
end
# DIAG=1: fit only the null model under several specifications to find what breaks the fit, then stop
if get(ENV, "DIAG", "0") == "1"
    W = first(m for (k, m) in Base.loaded_modules if k.name == "WiSER")
    f(vs) = FormulaTerm(term(y), term(1) + sum(term.(vs)))
    base = [:AGE, :HEIGHTX, :WEIGHTX]
    pc2, pc5, pc10 = pcsyms[1:2], pcsyms[1:5], pcsyms
    T = [timevar]; C = [:cohort]
    re1, reslope = FormulaTerm(term(y), term(1)), FormulaTerm(term(y), term(1) + term(timevar))
    f0 = FormulaTerm(term(y), term(1))
    keepv(k) = df[in.(df.IID, Ref(Set(subj[nv .>= k]))), :]
    # DIAG round 6: confirm the candidate within-subject variance (WSV) models in every group/outcome.
    g(lhs, vs) = isempty(vs) ? FormulaTerm(term(lhs), term(1)) : FormulaTerm(term(lhs), term(1) + sum(term.(vs)))
    A3 = [:AGE, :HEIGHTX, :WEIGHTX]; full = vcat(T, A3, pc10)
    length(unique(df.cohort)) > 1 && (full = vcat(full, C))
    TA = vcat(T, [:AGE])
    rs = FormulaTerm(term(y), term(1) + term(timevar)); ri = FormulaTerm(term(y), term(1))
    specs = [
      ("1 WSV time+AGE+BMI+PC1-2, RE 1+time", 1, g(y, full), g(y, vcat(TA, [:BMI], pc2)), rs, false),
      ("2 WSV time+AGE+BMI, RE 1+time",       1, g(y, full), g(y, vcat(TA, [:BMI])),      rs, false),
      ("3 WSV time+AGE+PC1-2, RE 1+time",     1, g(y, full), g(y, vcat(TA, pc2)),         rs, false),
      ("4 WSV time+AGE, RE 1+time",           1, g(y, full), g(y, TA),                    rs, false),
      ("5 WSV time+AGE+BMI+PC1-2, RE 1",      1, g(y, full), g(y, vcat(TA, [:BMI], pc2)), ri, false)]
    for (lab, k, mf, wf, rf, mom) in specs
        d = keepv(k)
        try
            m = W.WSVarLmmModel(mf, rf, wf, :IID, d)
            mom ? W.fit!(m; init = W.init_mom!(m)) : W.fit!(m)
            println("OK   ", lab, "  [", length(unique(d.IID)), " women, ", nrow(d), " visits]")
        catch e
            println("FAIL ", lab, "  [", length(unique(d.IID)), " women, ", nrow(d), " visits]   ", first(sprint(showerror, e), 160))
        end
    end
    exit(0)
end

# WiSER.fit! does not throw when Ipopt fails; it prints "status = ..." for each run. Capture that text and accept
# the model only if the final run converged and all parameters are finite.
const GOOD = ("LOCALLY_SOLVED", "OPTIMAL", "ALMOST_LOCALLY_SOLVED")
function fitlog(m)
    W = first(x for (k, x) in Base.loaded_modules if k.name == "WiSER")
    path, io = mktemp()
    redirect_stdout(io) do
        W.fit!(m)
    end
    close(io); txt = read(path, String); rm(path); return txt
end
statusof(txt) = (ms = collect(eachmatch(r"status = ([A-Z_]+)", txt)); isempty(ms) ? "UNKNOWN" : ms[end].captures[1])
finiteparams(m) = try all(isfinite, m.β) && all(isfinite, m.τ) catch; true end
# A "solved" fit can still be degenerate (2026-11-25: AI/AN tau for time = 29.6 or -14175 with SE ~ 0). Within-subject
# log-variance intercept must be in a sensible range and the other WSV effects (z-scored covariates; time per year) modest.
plausible(m) = try -20 < m.τ[1] < 0 && all(abs.(m.τ[2:end]) .<= 5) catch; true end
# a non-default WSV (env WSV) writes to its own folder, e.g. black_hipbmd_wsv_AGE_BMI, so test runs never overwrite the main ones
outdir = joinpath(B, "trajgwas", "$(group)_$(outcome)" * (haskey(ENV, "WSV") ? "_wsv_" * replace(ENV["WSV"], "," => "_") : "")); mkpath(outdir)
# One model for every set (2026-11-25): mean model above, random intercept + slope, WSV ~ 1 + AGE + BMI.
# Only the optimizer settings may change between sets (they do not change the model): try them in this fixed order and
# keep the first fit whose last run is LOCALLY_SOLVED/OPTIMAL with finite estimates; if none, the first ALMOST_LOCALLY_SOLVED
# (Ipopt "solved to acceptable level"). The setting used is printed and written to the null-model file.
W  = first(m for (k, m) in Base.loaded_modules if k.name == "WiSER")
IP = first(m for (k, m) in Base.loaded_modules if k.name == "Ipopt")
println("IIDs all digits: $(all(id -> all(isdigit, id), subj)); with a leading zero: $(count(id -> startswith(id, "0"), subj))")
reslope = FormulaTerm(term(y), term(1) + term(timevar))
reformula = get(ENV, "RE", "slope") == "intercept" ? FormulaTerm(term(y), term(1)) : reslope
# SPA=0: score tests without the saddlepoint approximation (sensitivity check; written to chrN.nospa.pval.txt)
nospa = get(ENV, "SPA", "1") == "0"
nullfile = joinpath(outdir, "chr$chr$gtag.null.txt")
pvalfile = joinpath(outdir, nospa ? "chr$chr$gtag.nospa.pval.txt" : "chr$chr$gtag.pval.txt")
# WALD=<file with one snpid per line>: Wald test (effect sizes + SEs) for the listed SNPs on this chromosome only,
# written to chrN.wald.txt; the genome-wide score-test files are not touched.
waldfile = get(ENV, "WALD", "")
waldinds = Int[]
if !isempty(waldfile)
    want = Set(strip.(filter(!isempty, readlines(waldfile))))
    bimids = [split(l)[2] for l in eachline(prefix * ".bim")]
    waldinds = findall(in(want), bimids)
    println("Wald mode: $(length(waldinds)) of $(length(want)) listed SNPs are on chr$chr")
    isempty(waldinds) && exit(0)
    pvalfile = joinpath(outdir, "chr$chr$gtag.wald.txt")
end
println("trajgwas keywords: ", try unique(reduce(vcat, [Base.kwarg_decl(mm) for mm in methods(trajgwas)])) catch; "?" end)
settings = [("WiSER default, LS start",                      nothing, :ls),
            ("WiSER default, MoM start",                     nothing, :mom),
            ("Ipopt max_iter 100, LS start",                 ["print_level" => 0, "max_iter" => 100], :ls),
            ("Ipopt watchdog, LS start",                     ["print_level" => 0, "watchdog_shortened_iter_trigger" => 3, "max_iter" => 100], :ls),
            ("Ipopt max_iter 1000, LS start",                ["print_level" => 0, "max_iter" => 1000], :ls),
            ("Ipopt max_iter 100, MoM start",                ["print_level" => 0, "max_iter" => 100], :mom),
            ("Ipopt watchdog, MoM start",                    ["print_level" => 0, "watchdog_shortened_iter_trigger" => 3, "max_iter" => 100], :mom),
            ("Ipopt max_iter 1000, MoM start",               ["print_level" => 0, "max_iter" => 1000], :mom),
            ("Ipopt mehrotra+warm start 1000, MoM start",    ["print_level" => 0, "mehrotra_algorithm" => "yes", "warm_start_init_point" => "yes", "max_iter" => 1000], :mom)]
function tryfit(opts, ini)
    m = W.WSVarLmmModel(meanf, reformula, wsvf, :IID, df)
    path, io = mktemp()
    redirect_stdout(io) do
        initm = ini == :mom ? W.init_mom!(m) : W.init_ls!(m)
        opts === nothing ? W.fit!(m; init = initm) :
                           W.fit!(m, IP.Optimizer(); init = initm, solver_config = Dict{String,Any}(opts))
    end
    close(io); txt = read(path, String); rm(path)
    return m, txt, statusof(txt)
end
nm = nothing; nmlab = ""; almost = nothing
for (lab, opts, ini) in settings
    try
        m, txt, st = tryfit(opts, ini)
        print(txt)
        fin = finiteparams(m); pl = fin && plausible(m); ok = fin && pl
        println("  -> $lab: final status $st", !fin ? " (non-finite estimates)" : !pl ? " (implausible WSV estimates, rejected)" : "",
                "  tau = ", try round.(m.τ, digits = 3) catch; "?" end)
        if ok && st in ("LOCALLY_SOLVED", "OPTIMAL")
            global nm = m; global nmlab = "$lab (final status $st)"; break
        elseif ok && st == "ALMOST_LOCALLY_SOLVED" && almost === nothing
            global almost = (m, "$lab (final status $st)")
        end
    catch e
        println("  -> $lab: fit failed: ", first(sprint(showerror, e), 200))
    end
end
if nm === nothing && almost !== nothing
    global nm = almost[1]; global nmlab = almost[2]
end
nm === nothing && error("null model did not converge under any optimizer setting")
println("MODEL USED: mean = $meanf | RE = $reformula | WSV = $wsvf | optimizer: $nmlab")
open(nullfile, "w") do io; println(io, "optimizer: ", nmlab); show(io, nm); end
println("null model id order matches famids[rowinds]: ", try string.(nm.ids) == famids[rowinds] catch err; "n/a ($(typeof(err)))" end)
try
    if !isempty(waldinds)
        trajgwas(nm, prefix; geneticrowinds = rowinds, pvalfile = pvalfile, test = :wald, snpinds = waldinds)
    elseif nospa
        trajgwas(nm, prefix; geneticrowinds = rowinds, pvalfile = pvalfile, usespa = false)
    else
        trajgwas(nm, prefix; geneticrowinds = rowinds, pvalfile = pvalfile)
    end
catch e
    println("trajgwas scan failed: ", first(sprint(showerror, e), 300))
    println("available trajgwas methods:"); foreach(println, methods(trajgwas))
    rethrow()
end
println("pval file: $(countlines(pvalfile) - 1) variants  ($pvalfile)")
