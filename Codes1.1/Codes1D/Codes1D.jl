module Codes1D
using LinearAlgebra
using SpecialFunctions
using Revise
using Test
struct QuadratureFormula
    points::Array{Float64}
    weights::Array{Float64}
end


include("Jacobi.jl")

end
