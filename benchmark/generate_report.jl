using Dates, Printf

const SCENARIOS = [
    (name="small",      nsites=120,  ntaxa=80,  total=5_000,  repeats=5, inner=100),
    (name="default",    nsites=400,  ntaxa=200, total=10_000, repeats=5, inner=20),
    (name="many_sites", nsites=800,  ntaxa=250, total=10_000, repeats=3, inner=3),
    (name="wide",       nsites=300,  ntaxa=1_000, total=20_000, repeats=3, inner=3),
    (name="large",      nsites=1_200, ntaxa=300, total=20_000, repeats=2, inner=1),
]

const BASE_TASKS = [
    "richness",
    "shannon_entropy",
    "alpha_diversity",
    "bray_curtis_distance_matrix",
    "jaccard_distance_matrix",
    "hellinger_distance_matrix",
]

const TASK_DISPLAY = Dict(
    "richness"                  => "richness",
    "shannon_entropy"           => "Shannon entropy",
    "alpha_diversity"           => "alpha diversity",
    "bray_curtis_distance_matrix" => "Bray-Curtis matrix",
    "jaccard_distance_matrix"   => "Jaccard matrix",
    "hellinger_distance_matrix" => "Hellinger matrix",
)

function benchmark_python()
    configured = get(ENV, "DIVERSITY_BENCH_PYTHON", "")
    if !isempty(configured)
        return configured
    end
    venv_python = joinpath(homedir(), ".venvs", "diversity-bench", "bin", "python")
    return isfile(venv_python) ? venv_python : "python3"
end

function run_capture(command, scenario)
    env = [
        "DIVERSITY_BENCH_NSITES"   => string(scenario.nsites),
        "DIVERSITY_BENCH_NTAXA"    => string(scenario.ntaxa),
        "DIVERSITY_BENCH_TOTAL"    => string(scenario.total),
        "DIVERSITY_BENCH_REPEATS"  => string(scenario.repeats),
        "DIVERSITY_BENCH_INNER"    => string(scenario.inner),
    ]
    return readchomp(setenv(command, env))
end

function parse_rows(output, scenario)
    rows = NamedTuple[]
    for line in split(output, '\n')
        isempty(strip(line)) && continue
        startswith(line, "language,") && continue
        parts = split(strip(line), ',')
        length(parts) in (7, 8) || continue
        inner = length(parts) == 8 ? parse(Int, parts[7]) : 1
        best_seconds = length(parts) == 8 ? String(parts[8]) : String(parts[7])
        push!(rows, (
            scenario=scenario.name,
            language=String(parts[1]),
            package=String(parts[2]),
            task=String(parts[3]),
            nsites=parse(Int, parts[4]),
            ntaxa=parse(Int, parts[5]),
            repeats=parse(Int, parts[6]),
            inner=inner,
            best_seconds=best_seconds,
        ))
    end
    return rows
end

function run_all()
    python = benchmark_python()
    julia = Base.julia_cmd()
    commands = [
        `$julia --project=. benchmark/julia_benchmark.jl`,
        `$python benchmark/python_benchmark.py`,
        `Rscript benchmark/r_vegan_benchmark.R`,
    ]

    rows = NamedTuple[]
    failures = String[]
    for scenario in SCENARIOS
        for command in commands
            try
                append!(rows, parse_rows(run_capture(command, scenario), scenario))
            catch err
                push!(failures, "$(scenario.name): $(command) failed with $(typeof(err)): $(err)")
            end
        end
    end
    return rows, failures, python
end

function benchmark_environment(python)
    cpu_name = isdefined(Sys, :CPU_NAME) ? getfield(Sys, :CPU_NAME) : "unknown"
    cpu_model = cpu_model_name()
    return (
        julia=string(VERSION),
        julia_threads=Threads.nthreads(),
        cpu=cpu_name,
        cpu_model=cpu_model,
        machine=Sys.MACHINE,
        os=string(Sys.KERNEL),
        python=python,
    )
end

function cpu_model_name()
    if Sys.islinux() && isfile("/proc/cpuinfo")
        for line in eachline("/proc/cpuinfo")
            if startswith(line, "model name")
                parts = split(line, ':'; limit=2)
                length(parts) == 2 && return strip(parts[2])
            end
        end
    end
    return "unknown"
end

function write_raw_csv(path, rows)
    open(path, "w") do io
        println(io, "scenario,language,package,task,nsites,ntaxa,repeats,inner,best_seconds")
        for row in rows
            println(io, join((
                row.scenario, row.language, row.package, row.task,
                row.nsites, row.ntaxa, row.repeats, row.inner, row.best_seconds,
            ), ","))
        end
    end
end

# ── display helpers ────────────────────────────────────────────────────────────

function is_prevalidated(row)
    return endswith(row.task, "_prevalidated")
end

function row_label(row)
    row.language == "Julia" || return row.language
    return is_prevalidated(row) ? "Julia (pre)" : "Julia (safe)"
end

function row_color(row)
    row.language == "Julia" && is_prevalidated(row) && return "#8fa8f4"
    row.language == "Julia"   && return "#4063d8"
    row.language == "Python"  && return "#389826"
    row.language == "R"       && return "#9558b2"
    return "#666666"
end

function xml_escape(value)
    text = string(value)
    text = replace(text, "&"  => "&amp;")
    text = replace(text, "<"  => "&lt;")
    text = replace(text, ">"  => "&gt;")
    text = replace(text, "\"" => "&quot;")
    return replace(text, "'"  => "&apos;")
end

function seconds_value(row)
    value = tryparse(Float64, row.best_seconds)
    return value === nothing ? NaN : value
end

function format_time(seconds::Float64)
    isnan(seconds)      && return "—"
    seconds == 0.0      && return "0"
    seconds < 1e-3      && return string(round(seconds * 1e6, sigdigits=3)) * " µs"
    seconds < 1.0       && return string(round(seconds * 1e3, sigdigits=3)) * " ms"
    return string(round(seconds, sigdigits=3)) * " s"
end

function scenario_heading(scenario)
    label = titlecase(replace(scenario.name, "_" => " "))
    return "$(label) ($(scenario.nsites) sites × $(scenario.ntaxa) taxa)"
end

# ── figure generation ─────────────────────────────────────────────────────────

function ordered_task_rows(rows, task)
    ordered = NamedTuple[]
    pre_task = task * "_prevalidated"
    for scenario in SCENARIOS
        srows = [r for r in rows if r.scenario == scenario.name]
        # Julia safe then Julia prevalidated, then Python, then R
        append!(ordered, [r for r in srows if r.language == "Julia" && r.task == task])
        append!(ordered, [r for r in srows if r.language == "Julia" && r.task == pre_task])
        for language in ("Python", "R")
            append!(ordered, [r for r in srows if r.language == language && r.task == task])
        end
    end
    return ordered
end

function write_task_figure(path, rows, task, title)
    task_rows = ordered_task_rows(rows, task)
    values = [seconds_value(row) for row in task_rows]
    positive_values = [v for v in values if isfinite(v) && v > 0]
    isempty(positive_values) && return false

    min_positive = minimum(positive_values)
    max_value    = maximum(positive_values)
    plot_min  = min_positive / 2
    plot_max  = max_value
    log_min   = log10(plot_min)
    log_max   = log10(plot_max)
    log_span  = max(log_max - log_min, 1e-9)

    width       = 1060
    left        = 270
    right       = 185
    top         = 72
    row_height  = 24
    gap         = 6
    chart_width = width - left - right
    height      = top + length(task_rows) * row_height + 70
    axis_y      = height - 48

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(width)" height="$(height)" viewBox="0 0 $(width) $(height)">""")
        println(io, """<rect width="100%" height="100%" fill="#ffffff"/>""")
        println(io, """<text x="24" y="32" font-family="sans-serif" font-size="20" font-weight="700" fill="#222222">$(xml_escape(title))</text>""")
        println(io, """<text x="24" y="54" font-family="sans-serif" font-size="12" fill="#555555">Log-scaled best elapsed time in seconds. Zero values reflect timer resolution and are drawn at half the smallest positive time.</text>""")

        for power in floor(Int, log_min):ceil(Int, log_max)
            tick = 10.0^power
            (tick < plot_min || tick > plot_max) && continue
            x = left + (log10(tick) - log_min) / log_span * chart_width
            println(io, """<line x1="$(x)" y1="$(top - 8)" x2="$(x)" y2="$(axis_y)" stroke="#eeeeee" stroke-width="1"/>""")
            println(io, """<text x="$(x)" y="$(axis_y + 18)" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#666666">1e$(power)</text>""")
        end

        for (index, row) in enumerate(task_rows)
            y           = top + (index - 1) * row_height
            value       = seconds_value(row)
            finite_val  = isfinite(value) ? value : 0.0
            plot_value  = finite_val > 0 ? finite_val : plot_min
            bar_width   = (log10(plot_value) - log_min) / log_span * chart_width
            label       = "$(row.scenario) / $(row_label(row))"
            shown       = finite_val == 0 ? "0" : row.best_seconds

            println(io, """<text x="$(left - 12)" y="$(y + 15)" text-anchor="end" font-family="sans-serif" font-size="12" fill="#333333">$(xml_escape(label))</text>""")
            println(io, """<rect x="$(left)" y="$(y + gap / 2)" width="$(max(bar_width, 1.0))" height="$(row_height - gap)" fill="$(row_color(row))" rx="2"/>""")
            println(io, """<text x="$(left + max(bar_width, 1.0) + 8)" y="$(y + 15)" font-family="sans-serif" font-size="11" fill="#333333">$(xml_escape(shown))</text>""")
        end

        legend_x = width - right + 25
        legend_items = [
            ("Julia (safe)", "#4063d8"),
            ("Julia (pre)",  "#8fa8f4"),
            ("Python",       "#389826"),
            ("R",            "#9558b2"),
        ]
        for (index, (label, color)) in enumerate(legend_items)
            y = top + (index - 1) * 22
            println(io, """<rect x="$(legend_x)" y="$(y)" width="12" height="12" fill="$(color)" rx="2"/>""")
            println(io, """<text x="$(legend_x + 18)" y="$(y + 11)" font-family="sans-serif" font-size="12" fill="#333333">$(label)</text>""")
        end

        println(io, """<line x1="$(left)" y1="$(axis_y)" x2="$(left + chart_width)" y2="$(axis_y)" stroke="#bbbbbb" stroke-width="1"/>""")
        println(io, """</svg>""")
    end
    return true
end

function write_figures(results_dir, rows)
    figures_dir = joinpath(results_dir, "figures")
    mkpath(figures_dir)
    specs = [
        ("richness",                  "Richness Across Scenarios",
         "Row-wise observed richness."),
        ("shannon_entropy",           "Shannon Entropy Across Scenarios",
         "Row-wise Shannon entropy."),
        ("alpha_diversity",           "Alpha-Diversity Summary Across Scenarios",
         "Compact exploratory alpha-diversity summary. R/vegan is not included because the R benchmark uses separate vegan calls."),
        ("bray_curtis_distance_matrix", "Bray-Curtis Distance Matrix Across Scenarios",
         "Dense pairwise Bray-Curtis dissimilarity matrix."),
        ("jaccard_distance_matrix",   "Jaccard Distance Matrix Across Scenarios",
         "Dense pairwise incidence Jaccard distance matrix."),
        ("hellinger_distance_matrix", "Hellinger Distance Matrix Across Scenarios",
         "Dense pairwise Hellinger distance matrix. R/vegan is not included for this direct helper comparison."),
    ]
    figures = NamedTuple[]
    for (task, title, caption) in specs
        filename = "$(task).svg"
        path = joinpath(figures_dir, filename)
        if write_task_figure(path, rows, task, title)
            push!(figures, (title=title, path=joinpath("figures", filename), caption=caption))
        end
    end
    return figures
end

# ── markdown report ────────────────────────────────────────────────────────────

function lookup_seconds(rows, scenario_name, language, task)
    for row in rows
        if row.scenario == scenario_name && row.language == language && row.task == task
            v = tryparse(Float64, row.best_seconds)
            return v === nothing ? NaN : v
        end
    end
    return NaN
end

function lookup_python_package(rows, scenario_name, task)
    for row in rows
        row.scenario == scenario_name && row.language == "Python" && row.task == task || continue
        return row.package
    end
    return ""
end

function write_scenario_table(io, rows, scenario)
    julia_safe_col  = String[]
    julia_pre_col   = String[]
    python_col      = String[]
    r_col           = String[]
    python_packages = Set{String}()

    for task in BASE_TASKS
        push!(julia_safe_col, format_time(lookup_seconds(rows, scenario.name, "Julia", task)))
        push!(julia_pre_col,  format_time(lookup_seconds(rows, scenario.name, "Julia", task * "_prevalidated")))
        push!(python_col,     format_time(lookup_seconds(rows, scenario.name, "Python", task)))
        push!(r_col,          format_time(lookup_seconds(rows, scenario.name, "R", task)))
        pkg = lookup_python_package(rows, scenario.name, task)
        isempty(pkg) || push!(python_packages, pkg)
    end

    pkg_note = isempty(python_packages) ? "" : " (" * join(sort(collect(python_packages)), "/") * ")"
    println(io, "| Task | Julia (safe) | Julia (pre) | Python$(pkg_note) | R/vegan |")
    println(io, "|---|---:|---:|---:|---:|")
    for (i, task) in enumerate(BASE_TASKS)
        display = get(TASK_DISPLAY, task, task)
        println(io, "| $(display) | $(julia_safe_col[i]) | $(julia_pre_col[i]) | $(python_col[i]) | $(r_col[i]) |")
    end
end

function write_markdown(path, rows, failures, python, figures, environment)
    generated = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    open(path, "w") do io
        println(io, "# DiversityAndDissimilarity Benchmark Report")
        println(io)
        println(io, "Generated: $(generated)")
        println(io)
        println(io, "Python interpreter: `$(python)`")
        println(io)
        println(io, "## Environment")
        println(io)
        println(io, "| Item | Value |")
        println(io, "|---|---|")
        println(io, "| CPU | `$(environment.cpu)` |")
        println(io, "| CPU model | `$(environment.cpu_model)` |")
        println(io, "| Machine | `$(environment.machine)` |")
        println(io, "| Operating system | `$(environment.os)` |")
        println(io, "| Julia | `$(environment.julia)` |")
        println(io, "| Julia threads | `$(environment.julia_threads)` |")
        println(io, "| Python interpreter | `$(environment.python)` |")
        println(io)
        println(io, "The benchmark matrix is generated deterministically at runtime for each scenario. Timings are best per-call elapsed time over the listed repeats and inner repetitions.")
        println(io)
        println(io, "## Scenarios")
        println(io)
        println(io, "| Scenario | Sites | Taxa | Total abundance | Repeats | Inner repetitions |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for scenario in SCENARIOS
            println(io, "| $(scenario.name) | $(scenario.nsites) | $(scenario.ntaxa) | $(scenario.total) | $(scenario.repeats) | $(scenario.inner) |")
        end

        if !isempty(figures)
            println(io)
            println(io, "## Figures")
            println(io)
            println(io, "Each figure shows Julia (safe), Julia (pre-validated), Python, and R bars per scenario on a shared log-scaled time axis. \"Julia (safe)\" validates on every call; \"Julia (pre)\" validates once before the loop.")
            for figure in figures
                println(io)
                println(io, "### $(figure.title)")
                println(io)
                println(io, "![$(figure.title)]($(figure.path))")
                println(io)
                println(io, figure.caption)
            end
        end

        println(io)
        println(io, "## Results")
        println(io)
        println(io, "\"Julia (safe)\" validates on every call (default user experience). \"Julia (pre)\" calls `validate` once before the loop, matching the Python and R calling convention and providing a fair computational comparison.")
        for scenario in SCENARIOS
            println(io)
            println(io, "### $(scenario_heading(scenario))")
            println(io)
            write_scenario_table(io, rows, scenario)
        end

        if !isempty(failures)
            println(io)
            println(io, "## Failures")
            println(io)
            for failure in failures
                println(io, "- `$(failure)`")
            end
        end

        println(io)
        println(io, "## Notes")
        println(io)
        println(io, "- Dense pairwise distance matrices scale quadratically in the number of sites.")
        println(io, "- Reported times are per call. Small scenarios use larger inner repetition counts to avoid coarse timer-resolution artifacts.")
        println(io, "- The Julia benchmark runs every timed function once before measurement so JIT compilation work is not included in the reported timings.")
        println(io, "- Julia and R timings call package-level APIs. Python uses NumPy/SciPy reference workflows in `benchmark/python_benchmark.py`.")
        println(io, "- The benchmark is intended to compare practical workflows, not isolated kernel implementations.")
    end
end

# ── LaTeX paper artifacts ──────────────────────────────────────────────────────

# Tasks included in the paper table: (scenario, task_key, display_name)
const LATEX_TABLE_ROWS = [
    ("default", "richness",                    "richness"),
    ("default", "shannon_entropy",             "Shannon entropy"),
    ("default", "bray_curtis_distance_matrix", "Bray--Curtis matrix"),
    ("default", "jaccard_distance_matrix",     "Jaccard matrix"),
    ("default", "hellinger_distance_matrix",   "Hellinger matrix"),
    ("large",   "bray_curtis_distance_matrix", "Bray--Curtis matrix"),
    ("large",   "jaccard_distance_matrix",     "Jaccard matrix"),
]

function format_latex_seconds(s::Float64)
    isnan(s) && return "---"
    s == 0.0 && return "0"
    exp_val = floor(Int, log10(abs(s)))
    mantissa = s / 10.0^exp_val
    m = round(mantissa, sigdigits=3)
    m_str = if m == round(m)
        string(round(Int, m))
    else
        raw = @sprintf("%.2f", m)
        raw = rstrip(raw, '0')
        endswith(raw, '.') ? raw[1:end-1] : raw
    end
    return "\$$(m_str) \\times 10^{$(exp_val)}\$"
end

function write_latex_table(path, rows, environment)
    open(path, "w") do io
        println(io, "\\begin{table}[t]")
        println(io, "\\centering")
        println(io, "\\caption{Selected warmed, per-call benchmark results on $(environment.cpu_model) with Julia $(environment.julia) using $(environment.julia_threads) Julia thread(s). Times in seconds. \\emph{Julia (safe)} validates on every call (the default user experience); \\emph{Julia (pre)} calls \\texttt{validate} once before the timing loop, matching the Python and R calling convention. Dashes indicate no direct equivalent in that stack.}")
        println(io, "\\label{tab:benchmarks}")
        println(io, "\\small")
        println(io, "\\resizebox{\\textwidth}{!}{%")
        println(io, "\\begin{tabular}{llrrrr}")
        println(io, "\\toprule")
        println(io, "Scenario & Task & Julia (safe) & Julia (pre) & Python & R/\\texttt{vegan}\\\\")
        println(io, "\\midrule")
        last_scen = ""
        for (scen, task, display) in LATEX_TABLE_ROWS
            scen != last_scen && last_scen != "" && println(io, "\\midrule")
            last_scen = scen
            js = format_latex_seconds(lookup_seconds(rows, scen, "Julia",   task))
            jp = format_latex_seconds(lookup_seconds(rows, scen, "Julia",   task * "_prevalidated"))
            py = format_latex_seconds(lookup_seconds(rows, scen, "Python",  task))
            rv = format_latex_seconds(lookup_seconds(rows, scen, "R",       task))
            println(io, "$(scen) & $(display) & $(js) & $(jp) & $(py) & $(rv)\\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "}")
        println(io, "\\end{table}")
    end
end

function write_latex_figure(path, rows)
    open(path, "w") do io
        println(io, "\\begin{figure}[t]")
        println(io, "\\centering")
        println(io, "\\definecolor{juliasafe}{HTML}{4063d8}")
        println(io, "\\definecolor{juliapre}{HTML}{8fa8f4}")
        println(io, "\\definecolor{pythoncol}{HTML}{389826}")
        println(io, "\\definecolor{rcol}{HTML}{9558b2}")
        println(io, "\\begin{tikzpicture}")
        println(io, "\\begin{axis}[")
        println(io, "  ybar,")
        println(io, "  bar width=8pt,")
        println(io, "  ymode=log,")
        println(io, "  ylabel={seconds},")
        println(io, "  symbolic x coords={default,large},")
        println(io, "  xtick=data,")
        println(io, "  legend style={at={(0.5,1.05)},anchor=south,legend columns=4},")
        println(io, "  width=0.88\\linewidth,")
        println(io, "  height=6cm,")
        println(io, "  nodes near coords,")
        println(io, "  nodes near coords align={vertical},")
        println(io, "  every node near coord/.append style={font=\\tiny},")
        println(io, "]")
        for (label, color, task) in [
            ("Julia (safe)", "juliasafe", "jaccard_distance_matrix"),
            ("Julia (pre)",  "juliapre",  "jaccard_distance_matrix_prevalidated"),
            ("SciPy",        "pythoncol", "jaccard_distance_matrix"),
            ("R/\\texttt{vegan}", "rcol", "jaccard_distance_matrix"),
        ]
            lang = startswith(label, "Julia") ? "Julia" :
                   startswith(label, "SciPy") ? "Python" : "R"
            d = lookup_seconds(rows, "default", lang, task)
            l = lookup_seconds(rows, "large",   lang, task)
            if !isnan(d) && !isnan(l)
                println(io, "\\addplot[fill=$(color)] coordinates {(default,$(d)) (large,$(l))};")
            end
        end
        println(io, "\\legend{Julia (safe),Julia (pre),SciPy,R/\\texttt{vegan}}")
        println(io, "\\end{axis}")
        println(io, "\\end{tikzpicture}")
        println(io, "\\caption{Jaccard distance matrix timings after bitset incidence encoding. \\emph{Julia (safe)} includes per-call validation; \\emph{Julia (pre)} validates once before the loop, matching the SciPy and R/\\texttt{vegan} calling convention. The benchmark is workflow-level and hardware-dependent, but illustrates the algorithmic benefit of packed incidence representation.}")
        println(io, "\\label{fig:jaccard}")
        println(io, "\\end{figure}")
    end
end

# ── entry point ────────────────────────────────────────────────────────────────

function main()
    rows, failures, python = run_all()
    results_dir = joinpath("benchmark", "results")
    mkpath(results_dir)
    write_raw_csv(joinpath(results_dir, "benchmark-results.csv"), rows)
    figures = write_figures(results_dir, rows)
    environment = benchmark_environment(python)
    write_markdown(joinpath(results_dir, "benchmark-report.md"), rows, failures, python, figures, environment)
    write_latex_table(joinpath("notes", "benchmark-results.tex"), rows, environment)
    write_latex_figure(joinpath("notes", "Figs", "benchmark-jaccard.tex"), rows)
    println("Wrote $(length(rows)) benchmark rows to $(joinpath(results_dir, "benchmark-report.md"))")
    for failure in failures
        println("FAILED: ", failure)
    end
    isempty(failures) || exit(1)
end

main()
