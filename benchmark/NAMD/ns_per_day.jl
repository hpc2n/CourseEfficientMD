println(typeof(ARGS))
println(ARGS[1])
timing_str="TIMING:"
time_step_str="TIMESTEP"

tlist=[]
global dt = 1.0
for line in eachline(ARGS[1])
      word = split(line)
      s = length(word)
      if s>1 && word[2]==time_step_str 
          global dt = parse(Float64,word[3])
      end
      if s>0 && word[1]==timing_str
          t = parse(Float64,split(word[8],'/')[1])
          push!(tlist,t)
      end
end

# calculate the mean
n = length(tlist)
global t_avg = 0.0
for t in tlist
    global t_avg +=t
end
if n > 0
    t_avg /= n
end

# calculate variance
t_var = 0.0
for t in tlist
   global t_var += (t - t_avg)^2
end
if n > 0
    t_var /= n
end

# calculate standard deviation
t_std = sqrt(t_var)

# calculate nanoseconds per day
if t_avg > 0.0 
   ns_per_day = (dt / t_avg) * (60.0 * 60.0 * 24.0 * 1e-6)
end

println("Nanoseconds per day ",ns_per_day)
println("Mean time per step ",t_avg)
println("Standard deviation ",t_std)
