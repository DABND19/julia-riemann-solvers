module FluxSolver

using StaticArrays

abstract type State end

function get_wave_velocities(solution::State)::Tuple{Float64, Float64} 
  error("You must implement get_wave_velocities method for $(typeof(solution)).")
end

function get_flux(solution::State)::SVector{3, Float64} 
  error("You must implement get_flux method for $(typeof(solution)).")
end

end
