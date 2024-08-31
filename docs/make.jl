using Documenter
using CryptoGroups

makedocs(
    sitename = "CryptoGroups",
    format = Documenter.HTML(),
    modules = [CryptoGroups]
)


deploydocs(repo = "github.com/PeaceFounder/CryptoGroups.jl.git")

