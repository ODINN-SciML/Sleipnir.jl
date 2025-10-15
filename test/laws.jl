# Using mock inputs is simpler than building a full `Simulation` object for testing `generate_input`
module MockTestInputs
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

module MockVjpsPrepLaw
    import Sleipnir: AbstractPrepVJP
    struct FakeVjpsPrepLaw <: AbstractPrepVJP end
end

using Sleipnir: _normalize_law_inputs, generate_inputs

generate_inputs_testset() = @testset "generate_inputs" begin
    (;A, B, C) = MockTestInputs

    @test generate_inputs((a = A(), b = B(), c = C()), nothing, nothing, 3.) == (a = ones(3), b = range(1, 3), c = 3.)
    @test generate_inputs((x = A(), y = B(), z = C()), nothing, nothing, 2.) == (x = ones(3), y = range(1, 3), z = 2.)

    # type stability
    JET.@test_opt generate_inputs((a = A(), b = B(), c = C()), nothing, nothing, 3.)
    JET.@test_opt generate_inputs((x = A(), y = B(), z = C()), nothing, nothing, 2.)
end

normalize_law_inputs_testset() = @testset "_normalize_law_inputs" begin
    (;A, B, C) = MockTestInputs

    @test _normalize_law_inputs((A(), B(), C())) == (a = A(), b = B(), c = C())
    @test _normalize_law_inputs((a = A(), b = B(), c = C())) == (a = A(), b = B(), c = C())
    @test _normalize_law_inputs((x = A(), y = B(), z = C())) == (x = A(), y = B(), z = C())
end

function test_law(;
    law,
    simulation,
    glacier_idx,
    θ,
    t,
    expected_cache_type,
    expected_cache,
    inputs_defined,
    expected_inputs,
    expected_apply_law_in_model,
)
    @test cache_type(law) == expected_cache_type
    JET.@test_opt init_cache(law, simulation, glacier_idx, θ)

    # test apply_law!
    let
        cache = init_cache(law, simulation, glacier_idx, θ)

        apply_law!(law, cache, simulation, glacier_idx, t, θ)
        @test cache == expected_cache
    end

    # test apply_law! stability
    let
        cache = init_cache(law, simulation, glacier_idx, θ)

        JET.@test_opt apply_law!(law, cache, simulation, glacier_idx, t, θ)
    end

    # test build_affect
    let
        # fake integrator
        integrator = (;
            p = simulation,
            t = t,
        )

        cache = init_cache(law, simulation, glacier_idx, θ)
        affect! = build_affect(law, cache, glacier_idx, θ)

        affect!(integrator)
        @test cache == expected_cache
    end

    # test build_affect stability
    let
        # fake integrator
        integrator = (;
            p = simulation,
            t = t,
        )

        cache = init_cache(law, simulation, glacier_idx, θ)
        affect! = build_affect(law, cache, glacier_idx, θ)

        JET.@test_opt affect!(integrator)
    end

    # test inputs
    if inputs_defined
        @test inputs(law) == expected_inputs
        JET.@test_opt inputs(law)
    else
        # Check that it fails
        @test_throws "Inputs are not defined." inputs(law)
    end

    # Test whether to apply the law inside the model
    @test apply_law_in_model(law) == expected_apply_law_in_model
end

apply_law_testset() = @testset "Law" begin
    (;A, B, C) = MockTestInputs

    @testset "ConstantLaw" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )

        law = ConstantLaw{Matrix{Float64}}(
            function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                fill(6., nx, ny)
            end,
        )

        test_law(;
            law,
            simulation,
            glacier_idx = 2,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = fill(6., 2, 3),
            expected_cache_type = Matrix{Float64},
            inputs_defined = true,
            expected_inputs = (;),
            expected_apply_law_in_model = false,
        )

        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4, n=4.),
                (; nx=2, ny=3, n=5.),
            ]
        )

        law = ConstantLaw{Float64}(
            function (simulation, glacier_idx, θ)
                simulation.glaciers[glacier_idx].n
            end,
        )

        test_law(;
            law,
            simulation,
            glacier_idx = 2,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = 5.,
            expected_cache_type = Float64,
            inputs_defined = true,
            expected_inputs = (;),
            expected_apply_law_in_model = false,
        )
    end

    @testset "Law without inputs" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )

        law = Law{Matrix{Float64}}(;
            f! = function (cache, simulation, glacier_idx, t, θ)
                @. cache = θ.a * t
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                zeros(nx, ny)
            end,
        )

        test_law(;
            law,
            simulation,
            glacier_idx = 2,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = fill(6., 2, 3),
            expected_cache_type = Matrix{Float64},
            inputs_defined = false,
            expected_inputs = (;),
            expected_apply_law_in_model = true,
        )
    end

    @testset "Law with inputs" begin
        # fake simulation
        simulation = (;
            glaciers = [
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
                (; nx, ny) = simulation.glaciers[glacier_idx]
                zeros(nx, ny)
            end,
        )

        test_law(;
            law,
            simulation,
            glacier_idx = 2,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = fill(6., 2, 3),
            expected_cache_type = Matrix{Float64},
            inputs_defined = true,
            expected_inputs = (a=A(), b=B(), c=C()),
            expected_apply_law_in_model = true,
        )

        # fake simulation
        simulation = (;
            glaciers = [
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
                (; nx, ny) = simulation.glaciers[glacier_idx]
                zeros(nx, ny)
            end,
        )

        test_law(;
            law,
            simulation,
            glacier_idx = 2,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = fill(6., 2, 3),
            expected_cache_type = Matrix{Float64},
            inputs_defined = true,
            expected_inputs = (x = A(), y = B(), z = C()),
            expected_apply_law_in_model = true,
        )
    end
end


function test_law_vjp(;
    law,
    simulation,
    glacier_idx,
    θ,
    t,
    expected_cache_type,
    expected_cache,
    expected_inputs,
    precomputed_vjp,
)
    @test cache_type(law) == expected_cache_type
    cache = init_cache(law, simulation, glacier_idx, θ)

    @test inputs(law) == expected_inputs
    JET.@test_opt inputs(law)

    apply_law!(law, cache, simulation, glacier_idx, t, θ)
    @test cache == expected_cache

    @test is_precomputable_law_VJP(law) == precomputed_vjp
    if precomputed_vjp
        (; FakeVjpsPrepLaw) = MockVjpsPrepLaw
        vjpsPrepLaw = FakeVjpsPrepLaw()
        precompute_law_VJP(law, cache, vjpsPrepLaw, simulation, glacier_idx, t, θ)
    end
    law_VJP_input(law, cache, simulation, glacier_idx, t, θ)
    law_VJP_θ(law, cache, simulation, glacier_idx, t, θ)
end

apply_vjp_law_testset() = @testset "Law VJPs" begin
    (;A, B, C) = MockTestInputs

    @testset "Law without VJPs" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )
        glacier_idx = 2
        θ = (;a = 3.0)
        t = 2.0

        law = Law{Matrix{Float64}}(;
            inputs = (A(), B(), C()),
            f! = function (cache, inputs, θ)
                @. cache = θ.a * inputs.c
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                zeros(nx, ny)
            end,
        )
        cache = init_cache(law, simulation, glacier_idx, θ)
        @test !is_precomputable_law_VJP(law)
        @test law_VJP_input(law, cache, simulation, glacier_idx, t, θ) === nothing
        @test law_VJP_θ(law, cache, simulation, glacier_idx, t, θ) === nothing
    end

    @testset "Law with custom VJP" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )

        law = Law{MatrixCache}(;
            inputs = (A(), B(), C()),
            f! = function (cache, inputs, θ)
                @. cache.value = θ.a * inputs.c
            end,
            f_VJP_input! = function (cache, inputs, θ)
                @. cache.vjp_inp = θ.a
            end,
            f_VJP_θ! = function (cache, inputs, θ)
                @. cache.vjp_θ = sum(inputs.c)
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                MatrixCache(zeros(nx, ny), zeros(nx, ny), zeros(1))
            end,
        )

        glacier_idx = 2
        (; nx, ny) = simulation.glaciers[glacier_idx]
        test_law_vjp(;
            law,
            simulation,
            glacier_idx = glacier_idx,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = MatrixCache(ones(nx, ny)*6, zeros(nx, ny), zeros(1)),
            expected_cache_type = MatrixCache,
            expected_inputs = (a=A(), b=B(), c=C()),
            precomputed_vjp = false,
        )
    end

    @testset "Law with VJP precomputation" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )

        law = Law{MatrixCache}(;
            inputs = (A(), B(), C()),
            f! = function (cache, inputs, θ)
                @. cache.value = θ.a * inputs.c
            end,
            f_VJP_input! = function (cache, inputs, θ)
            end,
            f_VJP_θ! = function (cache, inputs, θ)
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                MatrixCache(zeros(nx, ny), zeros(nx, ny), zeros(1))
            end,
            p_VJP! = function (cache, vjpsPrepLaw, inputs, θ)
                @. cache.vjp_inp = θ.a
                @. cache.vjp_θ = sum(inputs.c)
            end,
        )

        glacier_idx = 2
        (; nx, ny) = simulation.glaciers[glacier_idx]
        test_law_vjp(;
            law,
            simulation,
            glacier_idx = glacier_idx,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = MatrixCache(ones(nx, ny)*6, zeros(nx, ny), zeros(1)),
            expected_cache_type = MatrixCache,
            expected_inputs = (a=A(), b=B(), c=C()),
            precomputed_vjp = true,
        )
    end

    @testset "Law with VJP precomputation w/o parameters" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )

        law = Law{MatrixCache}(;
            inputs = (A(), B(), C()),
            f! = function (cache, inputs, θ)
                @. cache.value = inputs.c
            end,
            f_VJP_input! = function (cache, inputs, θ)
            end,
            f_VJP_θ! = function (cache, inputs, θ)
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                MatrixCache(zeros(nx, ny), zeros(nx, ny), zeros(0))
            end,
            p_VJP! = function (cache, vjpsPrepLaw, inputs, θ)
                @. cache.vjp_inp = 1.0
            end,
        )

        glacier_idx = 2
        (; nx, ny) = simulation.glaciers[glacier_idx]
        test_law_vjp(;
            law,
            simulation,
            glacier_idx = glacier_idx,
            θ = (;),
            t = 2.0,
            expected_cache = MatrixCache(ones(nx, ny)*2, zeros(nx, ny), zeros(0)),
            expected_cache_type = MatrixCache,
            expected_inputs = (a=A(), b=B(), c=C()),
            precomputed_vjp = true,
        )
    end

    @testset "Law with VJP precomputation w/o f_VJP_*!" begin
        # fake simulation
        simulation = (;
            glaciers = [
                (; nx=5, ny=4),
                (; nx=2, ny=3),
            ]
        )

        law = Law{MatrixCache}(;
            inputs = (A(), B(), C()),
            f! = function (cache, inputs, θ)
                @. cache.value = θ.a * inputs.c
            end,
            init_cache = function (simulation, glacier_idx, θ)
                (; nx, ny) = simulation.glaciers[glacier_idx]
                MatrixCache(zeros(nx, ny), zeros(nx, ny), zeros(1))
            end,
            p_VJP! = function (cache, vjpsPrepLaw, inputs, θ)
                @. cache.vjp_inp = θ.a
                @. cache.vjp_θ = sum(inputs.c)
            end,
        )

        glacier_idx = 2
        (; nx, ny) = simulation.glaciers[glacier_idx]
        test_law_vjp(;
            law,
            simulation,
            glacier_idx = glacier_idx,
            θ = (;a = 3.0),
            t = 2.0,
            expected_cache = MatrixCache(ones(nx, ny)*6, zeros(nx, ny), zeros(1)),
            expected_cache_type = MatrixCache,
            expected_inputs = (a=A(), b=B(), c=C()),
            precomputed_vjp = true,
        )
    end
end
