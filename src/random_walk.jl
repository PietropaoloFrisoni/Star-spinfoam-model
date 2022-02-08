function random_walk_function(j::Float64, D::Int64, d::Int64, A::Array{Float64, 5}, N::Int64, b::Int64, σ::Float64, draws_folder::String, ampls_folder::String, chain_id::Int64, verbosity::Int64, number_of_existing_draws::Int64=0)

  # containers for contraction (some allocations)
  v1 = zeros(Float64, D)
  v2 = zeros(Float64, D)
  v3 = zeros(Float64, D)
  v4 = zeros(Float64, D)
  v5 = zeros(Float64, D)
    
  draws = ElasticArray{Int64}(undef, 21, 0)     # some allocations    
  ampls = ElasticArray{Float64}(undef, 0)      
  
  Truncated_Coefficients = zeros(Float64, D)    # some allocations   
  Normal_distribution = Normal(0, σ)
 
  for i = 0:d
  r = 0.0
    for n = (-i):1:(d-i)            
    r += (cdf(Normal_distribution,n+0.5)-cdf(Normal_distribution,n-0.5))  
    end 
  Truncated_Coefficients[i+1] = r 
  end     
    
  if (myid() == 1) 
  if(verbosity > 1)  
  println("Truncated coefficients for j=$(j) are $(Truncated_Coefficients)\n")         
  end 
  end  

  #Initial draw and gaussian 
  draw = Array{Int64}(undef, 21)           # some allocations --- # 1 final slot for molteplicity
  gaussian_draw = Array{Int64}(undef, 20)  # some allocations 
    
  amp = 0.0    
  while(amp == 0) 
    for i = 1:20
    @inbounds draw[i] = rand((1:D))     # some allocations        
    end 
  amp = star_amplitude(D, A, v1, v2, v3, v4, v5, draw)  
  end
    
  draw[21] = 1 # Initial molteplicity      
  
  if (myid() == 1) 
  if(verbosity > 1)  
  println("Initial draw is $(draw[1:20]) with amp $(amp)\n")       
  end 
  end      
      
  #Proposed draw  
  prop_draw = Array{Int64}(undef, 20)      # some allocations
 
  acceptance_ratio = 0      
  molteplicity = 1                              # initial molteplicity counter
  draw_float_sample = Array{Float64}(undef,1)   # some allocations --- Distribution pkg does not allow it to be scalar without allocating memory     
    
  RW_monitor = true  # to test if the RW actually moved
    
  for n=1:1:N 
       
      RW_monitor = true  
        
      # Sampling proposed draw   
      Cx=Cx_prop=1.0 
      for i=1:1:20 
           while true
           rand!(Normal_distribution, draw_float_sample)    
           @inbounds gaussian_draw[i] = round(Int64, draw_float_sample[1])   
           @inbounds prop_draw[i] = draw[i] + gaussian_draw[i] 
           @inbounds !(1 <= prop_draw[i] && prop_draw[i] <= (D)) || break
           end
        if (gaussian_draw[i] != 0) RW_monitor = false; end    
           @inbounds Cx*=Truncated_Coefficients[draw[i]]    
           @inbounds Cx_prop*=Truncated_Coefficients[prop_draw[i]]  
      end 
        
      if (RW_monitor == true) 
      
          if (myid() == 1) 
          if(verbosity > 1)  
          println("Iteration $(n)---------------------------------------------------------------\n")          
          println("The prop_draw below:\n$(prop_draw[1:20])\nturns out to be equal to the current draw:\n$(draw[1:20])\nso that the molteplicity of the current draw is raised to $(molteplicity + 1)\n")
          end 
          end    
            
        acceptance_ratio += 1   
        molteplicity += 1
        continue    
            
      else      
        
          if (myid() == 1) 
          if(verbosity > 1)  
          println("Iteration $(n)---------------------------------------------------------------\n")        
          println("draw is:\n$(draw[1:20])\nprop_draw is:\n$(prop_draw[1:20])\namp is $(amp)\n")      
          end 
          end    
        

        Prop_amp = star_amplitude(D, A, v1, v2, v3, v4, v5, prop_draw) 
        
        p=min(1.0,(((Prop_amp^2)/(amp^2)))*(Cx/Cx_prop))       
      
        if (isnan(p))
        error("got NaN while computing densities ratio: prop_draw = $(prop_draw), amp = $(amp)")           
        end        
    
        r = rand()    
          
        if(r<p)
        
                if (myid() == 1) 
                if(verbosity > 1)  
                println("prop_draw $(prop_draw) was accepted, since p=$(p) and r=$(r)\n")
                end 
                end   

                    if(n>b)
                    # TODO avoid resizing at every iteration (pre-allocating assuming >~ 30% accept. rate and re-allocated only at the end) 
                    # and eventually replace matrices with 1-d arrays (more efficient?)
                    resize!(draws, 21, size(draws)[2]+1)     # some allocations
                    resize!(ampls, size(ampls)[1]+1)         # some allocations
                    draw[21] = molteplicity
                    draws[:, end] = draw[:]
                    ampls[end] = amp

                        if (myid() == 1) 
                        if(verbosity > 1)  
                        println("The old draw $(draw[1:20]) has been stored with molteplicity $(draw[21])\nthe corresponding amplitude $(amp) has been stored as well")  
                        end 
                        end 
                                
                    end
                    
                molteplicity = 1
                
             @turbo for i = 1:20
                    draw[i] = prop_draw[i]
                    end
                    
                amp = Prop_amp
                acceptance_ratio += 1
                
                if (myid() == 1) 
                if(verbosity > 1) 
                println("Now the new draw is $(draw[1:20])\nthe new amp is $(amp)\n")
                end 
                end 
        else      
                molteplicity += 1    
                if (myid() == 1) 
                if(verbosity > 1)  
                println("prop_draw $(prop_draw) was rejected, since p=$(p) and r=$(r)\nThe current draw $(draw[1:20]) remains the same and its molteplicity is $(molteplicity)\namp remains $(amp)\n")
                end 
                end
        end # if condition  (r<p)
            
        end # if condition  RW_monitor == true  
        
      # final storage  
      if (n == N)  
      resize!(draws, 21, size(draws)[2]+1)     # some allocations 
      resize!(ampls, size(ampls)[1]+1)         # some allocations 
      draw[21] = molteplicity
      draws[:, end] .= draw[:]
      ampls[end] = amp      
        if (myid() == 1) 
        if(verbosity > 1) 
        println("The last draw $(draw[1:20]) has been stored with molteplicity $(draw[21])\n")  
        end 
        end                     
      end  
        
  end # n cycle      
   
  #This is such that draws has structure: [N,21], so the operators will be computed faster  
  draws = transpose(draws)  # some allocations  
   
  if (chain_id == 1)  
  println("Done! $(acceptance_ratio*100/N)% of proposed draws have been accepted in master chain")  
  end               
        
  draws_number = size(draws)[1] 
  
  if (chain_id == 1) 
  if(verbosity > 0) 
  println("$(draws_number) draws stored in master chain") 
  end 
  end        
    
  @save "$(draws_folder)/draws_chain=$(chain_id+number_of_existing_draws).jld2" draws    # some allocations 
  @save "$(ampls_folder)/ampls_chain=$(chain_id+number_of_existing_draws).jld2" ampls    # some allocations 
 
end
