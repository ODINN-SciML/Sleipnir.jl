export prepare_vjp_law
export CustomVJP, DIVJP
export apply_all_non_callback_laws!, apply_all_callback_laws!, precompute_all_VJPs_laws!

"""
    emptyVJP(cache, simulation, glacier_idx, t, θ)
    emptyVJPWithInputs(cache, inputs, θ)

Function that defines an empty `!`-style function for the VJPs.
The two methods define the two possible signatures for the function that updates the cache in-place.
It is for internal use only and it isn't exposed to the user.
"""
emptyVJP(cache, simulation, glacier_idx, t, θ) = nothing
emptyVJPWithInputs(cache, inputs, θ) = nothing

"""
    emptyPrepVJP(cache, vjpPrep, simulation, glacier_idx, t, θ)
    emptyPrepVJPWithInputs(cache, vjpPrep, inputs, θ)

Function that defines an empty `!`-style function for the preparation of the VJPs.
The two methods define the two possible signatures for the function that updates the cache in-place.
Trying to apply this function will yield an error.
It is for internal use only and it isn't exposed to the user.
"""
emptyPrepVJP(cache, vjpPrep, simulation, glacier_idx, t, θ) = throw("This VJP preparation has not been defined.")
emptyPrepVJPWithInputs(cache, vjpPrep, inputs, θ) = throw("This VJP preparation has not been defined.")

"""
    ∂law∂θ!()

This function serves as a placeholder and should be replaced by other implementations in ODINN.
This implementation throws an error.
It is for internal use only and is not exposed to the user.
"""
∂law∂θ!() = throw("Function ∂law∂θ! not implemented. Please provide an implementation.")

"""
    ∂law∂inp!()

This function serves as a placeholder and should be replaced by other implementations in ODINN.
This implementation throws an error.
It is for internal use only and is not exposed to the user.
"""
∂law∂inp!() = throw("Function ∂law∂inp! not implemented. Please provide an implementation.")

"""
    prepare_vjp_law(simulation, law::AbstractLaw, law_cache, θ, glacier_idx)

Function used to prepare the VJPs at the initialization of the model cache.
It is used for example to compile VJPs of the laws to be differentiated using DifferentiationInterface.jl.
"""
prepare_vjp_law(simulation, law::AbstractLaw, law_cache, θ, glacier_idx) = throw("This prepare function has not been defined.")

"""
    VJPType

Abstract type representing the mode of Vector-Jacobian Product (VJP) computation used in a law.
Subtypes of `VJPType` define how the VJP is evaluated (e.g., custom or using DifferentiationInterface.jl).
"""
abstract type VJPType end

"""
    CustomVJP <: VJPType

Indicates that a law uses a custom-defined function for VJP computation.
This is used when the VJP is provided manually rather than computed automatically.
"""
struct CustomVJP <: VJPType end

"""
    DIVJP <: VJPType

Indicates that a law uses the default VJP computation provided by DifferentiationInterface.jl.
This is used when no custom VJP function is provided, and the Vector-Jacobian Product (VJP) is computed automatically using DifferentiationInterface.jl.
"""
struct DIVJP <: VJPType end



"""
    apply_all_non_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ)

This function is a placeholder and must be implemented for your custom model type.

It is intended to apply all non-callback laws in the simulation to update the model cache for a given glacier at time `t` with parameters `θ`. By default, calling this function will throw an error to indicate that the user should provide their own implementation tailored to their model.

# Arguments
- `model::AbstractModel`: The model instance.
- `cache`: The cache object storing state variables.
- `simulation`: The simulation context.
- `glacier_idx`: Index identifying the glacier.
- `t`: The current simulation time.
- `θ`: The parameter vector.

# Throws
- Always throws an error: `"This function should not be called. Implement apply_all_non_callback_laws! for your own model."`
"""
apply_all_non_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ) = throw("This function should not be called. Implement apply_all_non_callback_laws! for your own model.")

"""
    precompute_all_VJPs_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ)

This function is a placeholder and must be implemented for your custom model type.

It is intended to precompute the VJPs for all the laws that are used in a model. By default, calling this function will throw an error to indicate that the user should provide their own implementation tailored to their model.

# Arguments
- `model::AbstractModel`: The model instance.
- `cache`: The cache object storing state variables.
- `simulation`: The simulation context.
- `glacier_idx`: Index identifying the glacier.
- `t`: The current simulation time.
- `θ`: The parameter vector.

# Throws
- Always throws an error: `"This function should not be called. Implement precompute_all_VJPs_laws! for your own model."`
"""
precompute_all_VJPs_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ) = throw("This function should not be called. Implement precompute_all_VJPs_laws! for your own model.")

"""
    apply_all_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ)

This function is a placeholder and must be implemented for your custom model type.

It is intended to apply all callback laws in the simulation to update the model cache for a given glacier at time `t` with parameters `θ`. By default, calling this function will throw an error to indicate that the user should provide their own implementation tailored to their model.

# Arguments
- `model::AbstractModel`: The model instance.
- `cache`: The cache object storing state variables.
- `simulation`: The simulation context.
- `glacier_idx`: Index identifying the glacier.
- `t`: The current simulation time.
- `θ`: The parameter vector.

# Throws
- Always throws an error: `"This function should not be called. Implement apply_all_callback_laws! for your own model."`
"""
apply_all_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ) = throw("This function should not be called. Implement apply_all_callback_laws! for your own model.")
