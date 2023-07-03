
"""
    stop_condition_tstops(u,t,integrator, tstops)  

Function that iterates through the tstops, with a closure including `tstops`
"""
function stop_condition_tstops(u,t,integrator, tstops) 
    t in tstops
end

