module GasFlow

import JSON

const GAMMA::Float64 = 5.0 / 3.0

struct Params
  pressure::Float64
  density::Float64
  velocity::Float64
end

function sound_speed(flow::Params)::Float64
  if flow.pressure == 0.0 || flow.density == 0.0
    return 0.0
  end
  return sqrt(GAMMA * flow.pressure / flow.density)
end

function internal_energy(flow::Params)::Float64
  return flow.pressure / ((GAMMA - 1.0) * flow.density)
end

function total_energy(flow::Params)::Float64
  return internal_energy(flow) + 0.5 * flow.velocity^2
end

function mach_number(flow::Params)::Float64
  return abs(flow.velocity) / sound_speed(flow)
end

function from_conservative(u::Vector{Float64})::Params
  @assert length(u) == 3
  u_1, u_2, u_3 = u

  density::Float64 = u_1
  velocity::Float64 = u_2 / u_1
  pressure::Float64 = (GAMMA - 1.0) * (u_3 - 0.5 * density * velocity^2)

  @assert density >= 0.0
  @assert pressure >= 0.0

  return Params(pressure, density, velocity)
end

function to_conservative(flow::Params)::Vector{Float64}
  return [
    flow.density,
    flow.density * flow.velocity,
    flow.density * total_energy(flow)
  ]
end

function to_flux(flow::Params)::Vector{Float64}
  return [
    flow.density * flow.velocity,
    flow.pressure + flow.density * flow.velocity^2,
    flow.density * total_energy(flow) * flow.velocity + flow.pressure * flow.velocity
  ]
end

function JSON.lower(params::Params)
  return Dict(
    "pressure" => params.pressure,
    "density" => params.density,
    "velocity" => params.velocity,
    "mach_number" => mach_number(params),
  )
end

end
