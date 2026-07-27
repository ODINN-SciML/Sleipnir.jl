export velocityProduct

# Abstract type as a parent type for all losses
abstract type GeneralAbstractLoss end

# Basic losses whose purpose is to be used within other losses
abstract type AbstractSimpleLoss <: GeneralAbstractLoss end

# More advanced losses that are used in the code
abstract type AbstractLoss <: GeneralAbstractLoss end

velocityProduct(lossType::AbstractLoss) = nothing
