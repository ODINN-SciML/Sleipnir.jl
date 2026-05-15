export safe_approx

# Function to override "≈" to handle nothing values and dict
function safe_approx(a, b)
    if isnothing(a) && isnothing(b)
        return true
    elseif isnothing(a) || isnothing(b)
        return false
    elseif (a isa Dict) && (b isa Dict)
        keys(a) == keys(b) || return false
        return all(a[k] ≈ b[k] for k in keys(a))
    else
        return a ≈ b
    end
end
