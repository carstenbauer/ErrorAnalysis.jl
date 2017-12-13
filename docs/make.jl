using Documenter, ErrorAnalysis

makedocs(
    # options
    
)

makedocs(
    modules = [ErrorAnalysis],
    format = :html,
    sitename = "ErrorAnalysis.jl",
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md"
        # "Subsection" => [
        #     ...
        # ]
    ]
)

deploydocs(
    repo   = "github.com/crstnbr/ErrorAnalysis.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    julia  = "release",
    osname = "linux"
)