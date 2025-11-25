export AbstractLaw
export apply_law!, init_cache, law_VJP_input, law_VJP_θ, precompute_law_VJP, cache_type,
       is_callback_law, is_precomputable_law_VJP, callback_freq, inputs, inputs_defined,
       apply_law_in_model
export Container
export build_affect

"""
    AbstractLaw

Abstract type representing a synthetic law.
Currently it's only used for testing by making easier to create dumb laws, but in the future it may be cleaner to use different concrete type of laws
(for example `CallbackLaw`, `ContinuousLaw`, or `LearnableLaw`)

Concrete subtypes must implement:

  - `apply_law!(::ConcreteLaw, state, simulation, glacier_idx, t, θ)`
  - `init_cache(::ConcreteLaw, glacier, glacier_idx)`
  - `law_VJP_input(::ConcreteLaw, cache, simulation, glacier_idx, t, θ)`
  - `law_VJP_θ(::ConcreteLaw, cache, simulation, glacier_idx, t, θ)`
  - `precompute_law_VJP(::ConcreteLaw, cache, vjpsPrepLaw, simulation, glacier_idx, t, θ)`
  - `cache_type(::ConcreteLaw)`
  - `is_callback_law(::ConcreteLaw)`
  - `is_precomputable_law_VJP(::ConcreteLaw)`
  - `callback_freq(::ConcreteLaw)`
  - `inputs(::ConcreteLaw)`
  - `inputs_defined(::ConcreteLaw)`
  - `apply_law_in_model(::ConcreteLaw)`
"""
abstract type AbstractLaw end
function apply_law!(law::AbstractLaw, cache, simulation, glacier_idx, t, θ)
    throw(error("Concrete subtypes of AbstractLaw must implement apply_law!. Please provide an implementation for $(typeof(law))."))
end
function init_cache(law::AbstractLaw, simulation, glacier_idx, θ)
    throw(error("Concrete subtypes of AbstractLaw must implement init_cache. Please provide an implementation for $(typeof(law))."))
end
function law_VJP_input(law::AbstractLaw, cache, simulation, glacier_idx, t, θ)
    throw(error("Concrete subtypes of AbstractLaw must implement law_VJP_input. Please provide an implementation for $(typeof(law))."))
end
function law_VJP_θ(law::AbstractLaw, cache, simulation, glacier_idx, t, θ)
    throw(error("Concrete subtypes of AbstractLaw must implement law_VJP_θ. Please provide an implementation for $(typeof(law))."))
end
function precompute_law_VJP(
        law::AbstractLaw, cache, vjpsPrepLaw, simulation, glacier_idx, t, θ)
    throw(error("Concrete subtypes of AbstractLaw must implement precompute_law_VJP. Please provide an implementation for $(typeof(law))."))
end
function cache_type(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement cache_type. Please provide an implementation for $(typeof(law))."))
end
function is_callback_law(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement is_callback_law. Please provide an implementation for $(typeof(law))."))
end
function is_precomputable_law_VJP(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement is_precomputable_law_VJP. Please provide an implementation for $(typeof(law))."))
end
function callback_freq(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement callback_freq. Please provide an implementation for $(typeof(law))."))
end
function inputs(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement inputs. Please provide an implementation for $(typeof(law))."))
end
function inputs_defined(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement inputs_defined. Please provide an implementation for $(typeof(law))."))
end
function apply_law_in_model(law::AbstractLaw)
    throw(error("Concrete subtypes of AbstractLaw must implement apply_law_in_model. Please provide an implementation for $(typeof(law))."))
end

"""
    Container

Abstract type that defines a container to be used in the PDE solver.
It is useful to retrieve the `simulation` object when applying callback laws.
"""
abstract type Container end

"""
    retrieve_simulation(p)
    retrieve_simulation(p::Container)

Function that retrieves the `simulation` object from `integrator.p` when called from a callback.
If `p` is a subtype of `Container`, then `p.simulation` is returned, otherwise it returns `p`.
It is for internal use only and it isn't exposed to the user.
"""
retrieve_simulation(p) = p
retrieve_simulation(p::Container) = p.simulation

"""
    build_affect(law::AbstractLaw, cache, glacier_idx, θ)

Return a `!`-style function suitable for use in a callback, which applies the given `law`
to update the `cache` for a specific glacier and parameters `θ`, using the simulation time.
"""
function build_affect(law::AbstractLaw, cache, glacier_idx, θ)
    # The let block make sure that every variable are type stable
    return let law = law, cache = cache, glacier_idx = glacier_idx, θ = θ
        function affect!(integrator)
            simulation = retrieve_simulation(integrator.p)
            t = integrator.t

            apply_law!(law, cache, simulation, glacier_idx, t, θ)
        end
    end
end
