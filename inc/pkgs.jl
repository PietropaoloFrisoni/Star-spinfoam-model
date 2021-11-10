using Distributions    
using Random
using HalfIntegers    
using ParallelDataTransfer
using Tullio
using LoopVectorization    
using BenchmarkTools 
using JLD2   
using LinearAlgebra
using KahanSummation    
using DelimitedFiles    
using ElasticArrays 
using CSV
using DataFrames
using Dates


# useful command to install all required pkgs in a single loop

#=
vec = ["Distributions", "Random", "HalfIntegers", "ParallelDataTransfer", "LoopVectorization", "BenchmarkTools", "JLD2", "LinearAlgebra", "KahanSummation", "DelimitedFiles", "ElasticArrays", "CSV", "DataFrames", "Dates", "Tullio"]

import Pkg

for p in vec
Pkg.add("$p")
end
=#
