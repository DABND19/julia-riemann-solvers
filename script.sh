time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=0.0 \
  --A=0.1 \
  --frames-count=1000 \
  --input-file=stellar_wind/data/test_nh_zero.json \
  > data/a_10_nh_zero.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=0.0 \
  --A=0.33 \
  --frames-count=1000 \
  --input-file=stellar_wind/data/test_nh_zero.json \
  > data/a_33_nh_zero.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=0.0 \
  --A=0.5 \
  --frames-count=1000 \
  --input-file=stellar_wind/data/test_nh_zero.json \
  > data/a_50_nh_zero.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=0.0 \
  --A=0.66 \
  --frames-count=1000 \
  --input-file=stellar_wind/data/test_nh_zero.json \
  > data/a_66_nh_zero.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=0.0 \
  --A=0.8 \
  --frames-count=1000 \
  --input-file=stellar_wind/data/test_nh_zero.json \
  > data/a_80_nh_zero.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=0.0 \
  --A=0.95 \
  --frames-count=1000 \
  --input-file=stellar_wind/data/test_nh_zero.json \
  > data/a_95_nh_zero.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=1.429e-2 \
  --A=0.1 \
  --frames-count=1000 \
  --input-file=inputs/test_nh_medium.json \
  > data/a_10_nh_medium.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=1.429e-2 \
  --A=0.33 \
  --frames-count=1000 \
  --input-file=inputs/test_nh_medium.json \
  > data/a_33_nh_medium.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=1.429e-2 \
  --A=0.5 \
  --frames-count=1000 \
  --input-file=inputs/test_nh_medium.json \
  > data/a_50_nh_medium.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=1.429e-2 \
  --A=0.66 \
  --frames-count=1000 \
  --input-file=inputs/test_nh_medium.json \
  > data/a_66_nh_medium.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=1.429e-2 \
  --A=0.8 \
  --frames-count=1000 \
  --input-file=inputs/test_nh_medium.json \
  > data/a_80_nh_medium.jsonl
time julia -t 8 stellar_wind/main.jl \
  --t-end=100 \
  --Nh=1.429e-2 \
  --A=0.95 \
  --frames-count=1000 \
  --input-file=inputs/test_nh_medium.json \
  > data/a_95_nh_medium.jsonl
