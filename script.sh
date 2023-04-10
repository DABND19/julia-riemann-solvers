julia -t 8 stellar_wind/main.jl --t-end=100 --Nh=1.429e-3 > stellar_wind/data/test_nh_small.json
julia -t 8 stellar_wind/main.jl --t-end=100 --Nh=1.429e-2 > stellar_wind/data/test_nh_medium.json
julia -t 8 stellar_wind/main.jl --t-end=100 --Nh=5.7e-2 > stellar_wind/data/test_nh_large.json
julia -t 8 stellar_wind/main.jl --t-end=100 --Nh=1.14e-1 > stellar_wind/data/test_nh_xlarge.json
julia -t 8 stellar_wind/main.jl --t-end=100 --Nh=0.5 > stellar_wind/data/test_nh_xxlarge.json
