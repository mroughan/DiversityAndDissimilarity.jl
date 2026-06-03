using DiversityAndDissimilarity

const OUTPUT_DIR = joinpath(@__DIR__, "src", "assets")
const DOT_PATH = joinpath(OUTPUT_DIR, "type-tree.dot")
const SVG_PATH = joinpath(OUTPUT_DIR, "type-tree.svg")
const PDF_PATH = joinpath(OUTPUT_DIR, "type-tree.pdf")

function package_types(mod::Module)
    values = Any[]
    for name in names(mod; all=false, imported=false)
        binding = getproperty(mod, name)
        if binding isa DataType || binding isa UnionAll
            type = Base.unwrap_unionall(binding)
            parentmodule(type) === mod && push!(values, type)
        end
    end
    return sort(unique(values); by=string)
end

function dot_id(type::DataType)
    return replace(string(type), r"[^A-Za-z0-9_]" => "_")
end

function dot_label(type::DataType)
    return string(type)
end

function dot_shape(type::DataType)
    return isabstracttype(type) ? "box" : "ellipse"
end

function write_type_tree_dot(path::AbstractString, types)
    type_set = Set(types)
    roots = [type for type in types if !(supertype(type) in type_set)]
    open(path, "w") do io
        println(io, "digraph DiversityAndDissimilarityTypes {")
        println(io, "  graph [rankdir=TB, bgcolor=\"transparent\", pad=\"0.2\", nodesep=\"0.35\", ranksep=\"0.55\"];")
        println(io, "  node [style=\"filled\", fillcolor=\"#f8fafc\", color=\"#334155\", fontname=\"Helvetica\", fontsize=11];")
        println(io, "  edge [color=\"#64748b\", arrowsize=0.7];")
        println(io)
        for root in roots
            println(io, "  subgraph cluster_", dot_id(root), " {")
            println(io, "    label=\"", dot_label(root), " tree\";")
            println(io, "    color=\"#cbd5e1\";")
            println(io, "    fontname=\"Helvetica\";")
            for type in types
                if type == root || _has_ancestor(type, root)
                    println(io, "    ", dot_id(type), " [label=\"", dot_label(type), "\", shape=", dot_shape(type), "];")
                end
            end
            println(io, "  }")
            println(io)
        end
        for type in types
            parent = supertype(type)
            parent in type_set || continue
            println(io, "  ", dot_id(parent), " -> ", dot_id(type), ";")
        end
        println(io, "}")
    end
    return path
end

function _has_ancestor(type::DataType, ancestor::DataType)
    parent = supertype(type)
    while parent !== Any
        parent == ancestor && return true
        parent = supertype(parent)
    end
    return false
end

function run_graphviz(dot_path::AbstractString, svg_path::AbstractString, pdf_path::AbstractString)
    dot = Sys.which("dot")
    if dot === nothing
        @warn "Graphviz 'dot' was not found. Wrote DOT file only." dot_path
        return false
    end
    run(`$dot -Tsvg -o $svg_path $dot_path`)
    run(`$dot -Tpdf -o $pdf_path $dot_path`)
    return true
end

function main()
    mkpath(OUTPUT_DIR)
    types = package_types(DiversityAndDissimilarity)
    write_type_tree_dot(DOT_PATH, types)
    run_graphviz(DOT_PATH, SVG_PATH, PDF_PATH)
    return DOT_PATH
end

main()
