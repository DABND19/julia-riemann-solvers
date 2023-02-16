include("../riemann_solvers/RiemannSolvers.jl")

using JSON
using .RiemannSolvers

const X_LEFT::Float64 = -1.0
const X_RIGHT::Float64 = 1.0
const X_DIAPH::Float64 = 0.0

function main()
  if length(ARGS) != 8
    print(stderr, "You must provide 8 positional arguments:\n")
    print(stderr, "1) Left pressure.\n")
    print(stderr, "2) Left density.\n")
    print(stderr, "3) Left velocity.\n")
    print(stderr, "4) Right pressure.\n")
    print(stderr, "5) Right density.\n")
    print(stderr, "6) Right velocity.\n")
    print(stderr, "7) Time.\n")
    print(stderr, "8) Cells count.\n")
    return
  end
  left = GasFlow.Params(
    parse(Float64, ARGS[1]),
    parse(Float64, ARGS[2]),
    parse(Float64, ARGS[3])
  )
  right = GasFlow.Params(
    parse(Float64, ARGS[4]),
    parse(Float64, ARGS[5]),
    parse(Float64, ARGS[6])
  )
  t_end = parse(Float64, ARGS[7])
  cells_count::Int = parse(Int, ARGS[8])

  x = collect(LinRange(X_LEFT, X_RIGHT, cells_count + 1))
  initial_flow = [0.5 * (x[i] + x[i+1]) < X_DIAPH ? left : right for i in range(1, cells_count)]
  state = Godunov.State(x, initial_flow)

  Godunov.run!(
    state,
    t_end,
    Godunov.plain_difference_schema,
    Godunov.left_soft_boundary_condition,
    Godunov.right_soft_boundary_condition
  )

  print(JSON.json(state))
  println()
end

main()
