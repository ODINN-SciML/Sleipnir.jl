export safe_approx
export check_concrete_types, check_field_types

# Function to override "≈" to handle nothing values and dict.
# `rtol` is forwarded to `isapprox` for fields whose values vary slightly across
# platforms (e.g. projection-derived coordinates); when `nothing`, the default `≈` is used.
function safe_approx(a, b; rtol = nothing)
    if isnothing(a) && isnothing(b)
        return true
    elseif isnothing(a) || isnothing(b)
        return false
    elseif (a isa Dict) && (b isa Dict)
        keys(a) == keys(b) || return false
        return all(_approx(a[k], b[k], rtol) for k in keys(a))
    elseif (a isa Vector) && (b isa Vector)
        return all(_approx(ai, bi, rtol) for (ai, bi) in zip(a, b))
    else
        return _approx(a, b, rtol)
    end
end
function _approx(x, y, rtol)
    if typeof(x)==DateTime && typeof(y)==DateTime
        return x==y
    else
        return isnothing(rtol) ? (x ≈ y) : isapprox(x, y; rtol = rtol)
    end
end

"""
    check_concrete_types(x, indent=0, show=true)

Check recursively that the fields of a nested struct are all concrete.
Print information for each field and return false if any of the fields is not concrete.
"""
function check_concrete_types(x, indent = 0; show = true)
    T = typeof(x)
    prefix = " " ^ indent
    all_concrete = true
    if !isconcretetype(T)
        if show
            println("$(prefix)⚠ NON-CONCRETE: $(T)")
        end
        all_concrete = false
    else
        if show
            println("$(prefix)✓ $(T)")
        end
    end
    for field in fieldnames(T)
        val = getfield(x, field)
        if isstructtype(typeof(val)) && !isprimitivetype(typeof(val))
            if show
                println("$(prefix)  .$(field):")
            end
            all_concrete &= check_concrete_types(val, indent + 4; show = show)
        else
            FT = typeof(val)
            if isconcretetype(FT)
                if show
                    println("$(prefix)  .$(field): ✓ $(FT)")
                end
            else
                if show
                    println("$(prefix)  .$(field): ⚠ NON-CONCRETE: $(FT)")
                end
                all_concrete = false
            end
        end
    end
    return all_concrete
end

"""
    check_field_types(x, indent=0; show=true)

Check recursively that the fields of type are all concrete.
Print information for each field and return false if any of the fields is not concrete.
"""
function check_field_types(T::Type, indent = 0; show = true)
    prefix = " " ^ indent
    all_concrete = true
    for (name, FT) in zip(fieldnames(T), fieldtypes(T))
        if isconcretetype(FT)
            if show
                println("$(prefix)✓ .$(name)::$(FT)")
            end
        else
            if show
                println("$(prefix)⚠ .$(name)::$(FT)")
            end
            all_concrete = false
        end
        if isstructtype(FT) && !isprimitivetype(FT)
            all_concrete &= check_field_types(FT, indent + 4; show = show)
        end
    end
    return all_concrete
end
