# Defining dumb inputs is easier than creating a whole `Simulation` object for testing `generate_input`
module DumbTestInputs
    using Sleipnir: AbstractInput
    import Sleipnir: default_name, get_input

    struct A <: AbstractInput end
    default_name(::A) = :a
    get_input(::A, simulation, glacier_idx, t) = ones(3)

    struct B <: AbstractInput end
    default_name(::B) = :b
    get_input(::B, simulation, glacier_idx, t) = range(1, 3)

    struct C <: AbstractInput end
    default_name(::C) = :c
    get_input(::C, simulation, glacier_idx, t) = t
end

using Sleipnir: _normalize_law_inputs, generate_inputs, Law, apply_law!, init_cache, cache_type

generate_inputs_testset() = @testset "generate_inputs" begin
    (;A, B, C) = DumbTestInputs

    @test generate_inputs((a = A(), b = B(), c = C()), nothing, nothing, 3.) == (a = ones(3), b = range(1, 3), c = 3.)    
    @test generate_inputs((x = A(), y = B(), z = C()), nothing, nothing, 2.) == (x = ones(3), y = range(1, 3), z = 2.)

    # type stability
    @inferred generate_inputs((a = A(), b = B(), c = C()), nothing, nothing, 3.)
    @inferred generate_inputs((x = A(), y = B(), z = C()), nothing, nothing, 2.)
end

normalize_law_inputs_testset() = @testset "_normalize_law_inputs" begin
    (;A, B, C) = DumbTestInputs

    @test _normalize_law_inputs((A(), B(), C())) == (a = A(), b = B(), c = C())
    @test _normalize_law_inputs((a = A(), b = B(), c = C())) == (a = A(), b = B(), c = C())
    @test _normalize_law_inputs((x = A(), y = B(), z = C())) == (x = A(), y = B(), z = C())
end

apply_law_testset() = @testset "apply_law" begin 
    (;A, B, C) = DumbTestInputs

    @testset "ConstantLaw" begin
        # fake simulation
        simulation = (;
            glacier = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),                
            ]
        )

        law = ConstantLaw{Matrix{Float64}}(
            function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glacier[glacier_idx]
                fill(6., nx, ny)
            end,
        )

        glacier_idx = 2
        t = 2.
        θ = (;a = 3.)

        @test cache_type(law) == Matrix{Float64}
        cache = @inferred init_cache(law, simulation, glacier_idx, θ)


        apply_law!(law, cache, simulation, glacier_idx, t, θ)

        @test cache == fill(6., 2, 3)
    end

    @testset "Law without inputs" begin
        # fake simulation
        simulation = (;
            glacier = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),                
            ]
        )

        law = Law{Matrix{Float64}}(;
            f! = function (cache, simulation, glacier_idx, t, θ)
                @. cache = θ.a * t
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glacier[glacier_idx]
                zeros(nx, ny)
            end,
        )

        glacier_idx = 2
        t = 2.
        θ = (;a = 3.)

        @test cache_type(law) == Matrix{Float64}
        cache = @inferred init_cache(law, simulation, glacier_idx, θ)


        apply_law!(law, cache, simulation, glacier_idx, t, θ)

        @test cache == fill(6., 2, 3)
    end

    @testset "Law with inputs" begin
        # fake simulation
        simulation = (;
            glacier = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),                
            ]
        )

        law = Law{Matrix{Float64}}(;
            inputs = (A(), B(), C()),
            f! = function (cache, inputs, θ)
                @. cache = θ.a * inputs.c
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glacier[glacier_idx]
                zeros(nx, ny)
            end,
        )

        glacier_idx = 2
        t = 2.
        θ = (;a = 3.)

        @test cache_type(law) == Matrix{Float64}
        cache = @inferred init_cache(law, simulation, glacier_idx, θ)

        apply_law!(law, cache, simulation, glacier_idx, t, θ)

        @test cache == fill(6., 2, 3)


        # fake simulation
        simulation = (;
            glacier = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),                
            ]
        )

        law = Law{Matrix{Float64}}(;
            inputs = (x = A(), y = B(), z = C()),
            f! = function (cache, inputs, θ)
                @. cache = θ.a * inputs.z
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glacier[glacier_idx]
                zeros(nx, ny)
            end,
        )

        glacier_idx = 2
        t = 2.
        θ = (;a = 3.)

        @test cache_type(law) == Matrix{Float64}
        cache = @inferred init_cache(law, simulation, glacier_idx, θ)

        apply_law!(law, cache, simulation, glacier_idx, t, θ)

        @test cache == fill(6., 2, 3)
    end
end
