using Dates

const SCENARIOS = [
    (name="small", nsites=120, ntaxa=80, total=5_000, repeats=5, inner=100),
    (name="default", nsites=400, ntaxa=200, total=10_000, repeats=5, inner=20),
    (name="many_sites", nsites=800, ntaxa=250, total=10_000, repeats=3, inner=3),
    (name="wide", nsites=300, ntaxa=1_000, total=20_000, repeats=3, inner=3),
    (name="large", nsites=1_200, ntaxa=300, total=20_000, repeats=2, inner=1),
]

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
        "DIVERSITY_BENCH_NSITES" => string(scenario.nsites),
        "DIVERSITY_BENCH_NTAXA" => string(scenario.ntaxa),
        "DIVERSITY_BENCH_TOTAL" => string(scenario.total),
        "DIVERSITY_BENCH_REPEATS" => string(scenario.repeats),
        "DIVERSITY_BENCH_INNER" => string(scenario.inner),
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
                row.scenario,
                row.language,
                row.package,
                row.task,
                row.nsites,
                row.ntaxa,
                row.repeats,
                row.inner,
                row.best_seconds,
            ), ","))
        end
    end
end

function xml_escape(value)
    text = string(value)
    text = replace(text, "&" => "&amp;")
    text = replace(text, "<" => "&lt;")
    text = replace(text, ">" => "&gt;")
    text = replace(text, "\"" => "&quot;")
    return replace(text, "'" => "&apos;")
end

function seconds_value(row)
    value = tryparse(Float64, row.best_seconds)
    return value === nothing ? NaN : value
end

function language_color(language)
    language == "Julia" && return "#4063d8"
    language == "Python" && return "#389826"
    language == "R" && return "#9558b2"
    return "#666666"
end

function ordered_task_rows(rows, task)
    ordered = NamedTuple[]
    for scenario in SCENARIOS
        scenario_rows = [row for row in rows if row.scenario == scenario.name && row.task == task]
        for language in ("Julia", "Python", "R")
            append!(ordered, [row for row in scenario_rows if row.language == language])
        end
    end
    return ordered
end

function write_task_figure(path, rows, task, title)
    task_rows = ordered_task_rows(rows, task)
    values = [seconds_value(row) for row in task_rows]
    positive_values = [value for value in values if isfinite(value) && value > 0]
    isempty(positive_values) && return false

    min_positive = minimum(positive_values)
    max_value = maximum(positive_values)
    plot_min = min_positive / 2
    plot_max = max_value
    log_min = log10(plot_min)
    log_max = log10(plot_max)
    log_span = max(log_max - log_min, 1e-9)

    width = 1040
    left = 270
    right = 170
    top = 68
    row_height = 24
    gap = 6
    chart_width = width - left - right
    height = top + length(task_rows) * row_height + 70
    axis_y = height - 48

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(width)" height="$(height)" viewBox="0 0 $(width) $(height)">""")
        println(io, """<rect width="100%" height="100%" fill="#ffffff"/>""")
        println(io, """<text x="24" y="32" font-family="sans-serif" font-size="20" font-weight="700" fill="#222222">$(xml_escape(title))</text>""")
        println(io, """<text x="24" y="54" font-family="sans-serif" font-size="12" fill="#555555">Log-scaled best elapsed time in seconds. Zero values reflect timer resolution and are drawn at half the smallest positive time.</text>""")

        for power in floor(Int, log_min):ceil(Int, log_max)
            tick = 10.0^power
            if tick < plot_min || tick > plot_max
                continue
            end
            x = left + (log10(tick) - log_min) / log_span * chart_width
            println(io, """<line x1="$(x)" y1="$(top - 8)" x2="$(x)" y2="$(axis_y)" stroke="#eeeeee" stroke-width="1"/>""")
            println(io, """<text x="$(x)" y="$(axis_y + 18)" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#666666">1e$(power)</text>""")
        end

        for (index, row) in enumerate(task_rows)
            y = top + (index - 1) * row_height
            value = seconds_value(row)
            finite_value = isfinite(value) ? value : 0.0
            plot_value = finite_value > 0 ? finite_value : plot_min
            bar_width = (log10(plot_value) - log_min) / log_span * chart_width
            label = "$(row.scenario) / $(row.language)"
            shown = finite_value == 0 ? "0" : row.best_seconds

            println(io, """<text x="$(left - 12)" y="$(y + 15)" text-anchor="end" font-family="sans-serif" font-size="12" fill="#333333">$(xml_escape(label))</text>""")
            println(io, """<rect x="$(left)" y="$(y + gap / 2)" width="$(max(bar_width, 1.0))" height="$(row_height - gap)" fill="$(language_color(row.language))" rx="2"/>""")
            println(io, """<text x="$(left + max(bar_width, 1.0) + 8)" y="$(y + 15)" font-family="sans-serif" font-size="11" fill="#333333">$(xml_escape(shown))</text>""")
        end

        legend_x = width - right + 25
        for (index, language) in enumerate(("Julia", "Python", "R"))
            y = top + (index - 1) * 22
            println(io, """<rect x="$(legend_x)" y="$(y)" width="12" height="12" fill="$(language_color(language))" rx="2"/>""")
            println(io, """<text x="$(legend_x + 18)" y="$(y + 11)" font-family="sans-serif" font-size="12" fill="#333333">$(language)</text>""")
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
        ("richness", "Richness Across Scenarios", "Row-wise observed richness."),
        ("shannon_entropy", "Shannon Entropy Across Scenarios", "Row-wise Shannon entropy."),
        ("alpha_diversity", "Alpha-Diversity Summary Across Scenarios", "Compact exploratory alpha-diversity summary. R/vegan is not included because the R benchmark uses separate vegan calls."),
        ("bray_curtis_distance_matrix", "Bray-Curtis Distance Matrix Across Scenarios", "Dense pairwise Bray-Curtis dissimilarity matrix."),
        ("jaccard_distance_matrix", "Jaccard Distance Matrix Across Scenarios", "Dense pairwise incidence Jaccard distance matrix."),
        ("hellinger_distance_matrix", "Hellinger Distance Matrix Across Scenarios", "Dense pairwise Hellinger distance matrix. R/vegan is not included for this direct helper comparison."),
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
            println(io, "The figures use a log-scaled time axis so small row-wise summaries and larger pairwise matrix calculations can be read on the same page. The raw values remain in the results table below.")
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
        println(io, "| Scenario | Language | Package | Task | Sites | Taxa | Repeats | Inner | Best seconds |")
        println(io, "|---|---|---|---|---:|---:|---:|---:|---:|")
        for row in rows
            println(io, "| $(row.scenario) | $(row.language) | $(row.package) | $(row.task) | $(row.nsites) | $(row.ntaxa) | $(row.repeats) | $(row.inner) | $(row.best_seconds) |")
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

function main()
    rows, failures, python = run_all()
    results_dir = joinpath("benchmark", "results")
    mkpath(results_dir)
    write_raw_csv(joinpath(results_dir, "benchmark-results.csv"), rows)
    figures = write_figures(results_dir, rows)
    environment = benchmark_environment(python)
    write_markdown(joinpath(results_dir, "benchmark-report.md"), rows, failures, python, figures, environment)
    println("Wrote $(length(rows)) benchmark rows to $(joinpath(results_dir, "benchmark-report.md"))")
    for failure in failures
        println("FAILED: ", failure)
    end
    isempty(failures) || exit(1)
end

main()
