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

