export safe_approx

# Function to override "≈" to handle nothing values
function safe_approx(a, b)
    if isnothing(a) && isnothing(b)
        return true
    elseif isnothing(a) || isnothing(b)
        return false
    else
        return a ≈ b
    end
end

function get_property_or_zero_PyObject(obj, prop_name::Symbol)
    if PyCall.hasproperty(obj, prop_name)
        return PyCall.getproperty(obj, prop_name)
    else
        return 0.0
    end
end
    
