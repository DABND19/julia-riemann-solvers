module HLL

using StaticArrays

import ..RiemannSolvers.FluxSolver
import ..RiemannSolvers.GasFlow

struct State <: FluxSolver.State
  left::GasFlow.Params
  right::GasFlow.Params
end

@inline function FluxSolver.get_wave_velocities(state::State)::Tuple{Float64,Float64}
  u_l = state.left.velocity
  c_l = GasFlow.sound_speed(state.left)

  u_r = state.right.velocity
  c_r = GasFlow.sound_speed(state.right)

  return (
    min(u_l, u_r) - max(c_l, c_r), max(u_l, u_r) + max(c_l, c_r)
  )
end

function FluxSolver.get_flux(state::State)::SVector{3,Float64}
  s::Float64 = 0.0
  S_l, S_r = FluxSolver.get_wave_velocities(state)

  u_l = GasFlow.to_conservative(state.left)
  F_l = GasFlow.to_flux(state.left)

  u_r = GasFlow.to_conservative(state.right)
  F_r = GasFlow.to_flux(state.right)

  if s <= S_l
    return F_l
  end

  if s >= S_r
    return F_r
  end

  return (S_r * F_l - S_l * F_r + S_l * S_r * (u_r - u_l)) / (S_r - S_l)
end

end
