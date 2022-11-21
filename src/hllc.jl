module HLLC

using StaticArrays

import Main.FluxSolver
import Main.GasFlow

struct State <: FluxSolver.State
  left::GasFlow.Params
  right::GasFlow.Params
end

function U_c_k(flow::GasFlow.Params, S_k::Float64, S_c::Float64)::SVector{3,Float64}
  r_k = flow.density
  p_k = flow.pressure
  u_k = flow.velocity
  E_k = GasFlow.total_energy(flow)

  c = r_k * (S_k - u_k) / (S_k - S_c)
  return c * @SVector [1.0, S_c, E_k + (S_c - u_k) * (S_c + p_k / (r_k * (S_k - u_k)))]
end

function get_contact_velocity(left::GasFlow.Params, right::GasFlow.Params, S_l::Float64, S_r::Float64)::Float64
  p_l = left.pressure
  r_l = left.density
  u_l = left.velocity

  p_r = right.pressure
  r_r = right.density
  u_r = right.velocity

  m_l = r_l * (S_l - u_l)
  m_r = r_r * (S_r - u_r)

  return (p_r - p_l + u_l * m_l - u_r * m_r) / (m_l - m_r)
end

function FluxSolver.get_wave_velocities(state::State)::Tuple{Float64,Float64}
  u_l = state.left.velocity
  c_l = GasFlow.sound_speed(state.left)

  u_r = state.right.velocity
  c_r = GasFlow.sound_speed(state.right)

  return (min(u_l, u_r) - max(c_l, c_r), max(u_l, u_r) + max(c_l, c_r))
end

function FluxSolver.get_flux(state::State)::SVector{3,Float64}
  s::Float64 = 0.0

  F_l = GasFlow.to_flux(state.left)
  F_r = GasFlow.to_flux(state.right)

  S_l, S_r = FluxSolver.get_wave_velocities(state)
  S_c = get_contact_velocity(state.left, state.right, S_l, S_r)

  if s <= S_l
    return F_l
  end

  if s >= S_r
    return F_r
  end

  U_l = GasFlow.to_conservative(state.left)
  U_r = GasFlow.to_conservative(state.right)
  U_c_l = U_c_k(state.left, S_l, S_c)
  U_c_r = U_c_k(state.right, S_r, S_c)

  if s < S_c
    return F_l + S_l * (U_c_l - U_l)
  else
    return F_r + S_r * (U_c_r - U_r)
  end
end

end
