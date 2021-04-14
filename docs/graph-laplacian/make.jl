using Documenter, RLinearAlgebra

makedocs(sitename="Randomized Graph Laplacian Algorithms",
    pages=[
        "Overview" => "index.md",
        "Graphs and Graph Laplacians" => "graphs.md",
        "Cholesky Decompositions" => "cholesky.md"
    ])

# Graphs are generated using tikz and pdflatex, see figures/
# Once graphs are generated, output is saved to figures/png
