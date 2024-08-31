using Documenter
using CryptoGroups

makedocs(
    sitename = "CryptoGroups",
    repo = Documenter.Remotes.GitHub("PeaceFounder", "CryptoGroups.jl"),
    format = Documenter.HTML(),
    modules = [CryptoGroups],
    warnonly = true
)


deploydocs(repo = "github.com/PeaceFounder/CryptoGroups.jl.git")

