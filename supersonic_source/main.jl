include("../riemann_solvers/RiemannSolvers.jl")

using JSON
using StaticArrays

using .RiemannSolvers

const R_0::Float64 = 0.05
const R_INF::Float64 = 1.0

const M_0::Float64 = 10.0

const SOURCE_FLOW::GasFlow.Params = GasFlow.Params(
  1.0 / (GasFlow.GAMMA * M_0^2),
  1.0,
  1.0
)

function exact_solution(r::Float64)::GasFlow.Params
  return GasFlow.Params(
    SOURCE_FLOW.pressure * (R_0 / r)^(2 * GasFlow.GAMMA),
    SOURCE_FLOW.density * (R_0 / r)^2,
    SOURCE_FLOW.velocity
  )
end

function source_boundary_condition(state::Godunov.State)::SVector{3,Float64}
  return GasFlow.to_flux(SOURCE_FLOW)
end

function main()
  if length(ARGS) != 2
    print("You must provide time and cells count.")
    println()
    return
  end
  t_end::Float64 = parse(Float64, ARGS[1])
  cells_count = parse(Int, ARGS[2])
  r = collect(LinRange(R_0, R_INF, cells_count + 1))
  initial_flow = [exact_solution(0.5 * (r[i+1] + r[i])) for i in range(1, cells_count)]

  state = Godunov.State{HLLC.State}(r, initial_flow)

  Godunov.run!(
    state,
    t_end,
    Godunov.spherical_difference_schema,
    source_boundary_condition,
    Godunov.right_soft_boundary_condition
  )

  print(json(state))
  println()
end

main()
