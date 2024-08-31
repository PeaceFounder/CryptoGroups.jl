dirname(@__DIR__) in LOAD_PATH || Base.push!(LOAD_PATH, dirname(@__DIR__))

using Documenter
using CryptoGroups

makedocs(
    sitename = "CryptoGroups.jl",
    repo = Documenter.Remotes.GitHub("PeaceFounder", "CryptoGroups.jl"),
    format = Documenter.HTML(),
    modules = [CryptoGroups],
    warnonly = true
)


deploydocs(repo = "github.com/PeaceFounder/CryptoGroups.jl.git")

