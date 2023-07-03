import Pkg
Pkg.activate(dirname(Base.current_project()))

using PyCall
# Update SSL certificate to avoid issue in GitHub Actions CI
certifi = pyimport("certifi")
ENV["SSL_CERT_FILE"] = certifi.where()
# println("Current SSL certificate: ", ENV["SSL_CERT_FILE"])

using Revise
using Test
using JLD2
using Infiltrator

include(joinpath("..", "src/Sleipnir.jl"))

include(joinpath(Sleipnir.root_dir, "test/params_construction.jl"))
include(joinpath(Sleipnir.root_dir, "test/glacier_construction.jl"))

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Parameters constructors" params_constructor_specified(save_refs=true)

@testset "Parameters constructors" params_constructor_default(save_refs=true)

@testset "Glaciers and Climate constructors" glacier_construction(save_refs=true)
