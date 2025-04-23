# Defining dumb inputs is easier than creating a whole `Simulation` object for testing `generate_input`
module DumbTestInputs
    using Sleipnir
    import Sleipnir: default_name, get_input

    struct A <: Sleipnir.AbstractInput end
    default_name(::A) = :a
    get_input(::A, simulation, glacier_idx, t) = ones(3)

    struct B <: Sleipnir.AbstractInput end
    default_name(::B) = :b
    get_input(::B, simulation, glacier_idx, t) = range(1, 3)

    struct C <: Sleipnir.AbstractInput end
    default_name(::C) = :c
    get_input(::C, simulation, glacier_idx, t) = t
end

using Sleipnir: generate_inputs, apply_law

generate_inputs_testset() = @testset "generate_inputs" begin
    (;A, B, C) = DumbTestInputs

    @test generate_inputs((a = A(), b = B(), c = C()), nothing, nothing, 3.) == (a = ones(3), b = range(1, 3), c = 3.)    
    @test generate_inputs((x = A(), y = B(), z = C()), nothing, nothing, 2.) == (x = ones(3), y = range(1, 3), z = 2.)

    # type stability
    @inferred generate_inputs((a = A(), b = B(), c = C()), nothing, nothing, 3.)
    @inferred generate_inputs((x = A(), y = B(), z = C()), nothing, nothing, 2.)
end

apply_law_testset() = @testset "apply_law" begin 
    (;A, B, C) = DumbTestInputs

    law = Sleipnir.Law(;
        inputs = (A(), B(), C()),
        f = (inputs, θ) -> inputs.a .* inputs.b .+ inputs.c .+ θ.a,
    )

    theta = (;a = 3.)

    @test apply_law(law, nothing, nothing, 2., theta) == [6., 7., 8.]
    @inferred apply_law(law, nothing, nothing, 2., theta)

    law = Sleipnir.Law(;
        inputs = (x = A(), y = B(), z = C()),
        f = (inputs, θ) -> inputs.x .* inputs.y .+ inputs.z .+ θ.a,
    )

    theta = (;a = 3.)

    @test apply_law(law, nothing, nothing, 2., theta) == [6., 7., 8.]
    @inferred apply_law(law, nothing, nothing, 2., theta)
end
