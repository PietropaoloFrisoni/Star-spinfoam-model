function compute_volumes_function!(j, D, vertex, draws, number_of_draws, ampls, volumes, volumes_sq, volumes_matrix_values, N, b, volumes_folder, volumes_sq_folder, chain_id, total_volumes_already_stored=0)
 
# Here draws are trasposed to allow a faster broadcast 

  #containers for contraction
  v1 = zeros(Float64, D)
  v2 = zeros(Float64, D)
  v3 = zeros(Float64, D)
  v4 = zeros(Float64, D)
  v5 = zeros(Float64, D)
 
  indices = zeros(Int64,20)
 
  all_volumes = zeros(Float64, number_of_draws, 1)   
  all_squared_volumes = zeros(Float64, number_of_draws, 1)  
    
     for n = 1:number_of_draws   
     
     outer_step_vol = 0.0 
     outer_step_sq_vol = 0.0 
     indices[:] .= draws[1:20, n]
     
       for i_n = 1:D
       
       indices[1] = i_n
       amp = star_amplitude(D, vertex, v1, v2, v3, v4, v5, indices)         
       @inbounds outer_step_vol += amp*volumes_matrix_values[1][j][i_n, draws[1, n]]
       @inbounds outer_step_sq_vol += amp*volumes_matrix_values[2][j][i_n, draws[1, n]]
                  
       end
       
     @inbounds all_volumes[1] += outer_step_vol*draws[21,n]/(ampls[n]) 
     @inbounds all_squared_volumes[1] += outer_step_sq_vol*draws[21,n]/(ampls[n]) 
                 
     end    
     
  volumes[1] = sum_kbn(all_volumes[:])/(N - b)  
  volumes_sq[1] = sum_kbn(all_squared_volumes[:])/(N - b) 
  
  @save "$(volumes_folder)/volumes_chain=$(chain_id + total_volumes_already_stored).jld2" volumes 
  @save "$(volumes_sq_folder)/volumes_sq_chain=$(chain_id + total_volumes_already_stored).jld2" volumes_sq
  
end







function compute_volumes_pseudo_correlations_function!(j, D, draws, number_of_draws, ampls, vertex, volumes_pseudo_correlations, volumes_matrix_values, N, b, volumes_pseudo_correlations_folder, chain_id, node_1, node_2, total_volumes_pseudo_correlations_already_stored=0)

# Here draws are trasposed to allow a faster broadcast 

  #containers for contraction
  v1 = zeros(Float64, D)
  v2 = zeros(Float64, D)
  v3 = zeros(Float64, D)
  v4 = zeros(Float64, D)
  v5 = zeros(Float64, D)

  index_1 = node_1 
  index_2 = node_2 

  all_volumes_pseudo_correlations_for_two_specific_nodes = zeros(Float64, number_of_draws, 1)   

  indices = Array{Int64}(undef, 20)  

    for n = 1:number_of_draws

    @inbounds indices[:] .= draws[1:20, n]
    outer_step_vol_1 = 0.0

        for i_n = 1:D 

        indices[index_1] = i_n
        outer_step_vol_2 = 0.0

            for i_m = 1:D

            indices[index_2] = i_m
            amp = star_amplitude(D, vertex, v1, v2, v3, v4, v5, indices)   
            @inbounds outer_step_vol_2 += amp*volumes_matrix_values[1][j][i_m, draws[index_2, n]]

            end   

        @inbounds outer_step_vol_1 += outer_step_vol_2*volumes_matrix_values[1][j][i_n, draws[index_1, n]]
        
        end

    @inbounds all_volumes_pseudo_correlations_for_two_specific_nodes[1] += outer_step_vol_1*draws[21,n]/(ampls[n])     
 
    end

  volumes_pseudo_correlations[1] = sum_kbn(all_volumes_pseudo_correlations_for_two_specific_nodes[:])/(N - b)  
  
  @save "$(volumes_pseudo_correlations_folder)/volumes_pseudo_correlations_node_1=$(node_1)_node_2=$(node_2)_chain=$(chain_id + total_volumes_pseudo_correlations_already_stored).jld2" volumes_pseudo_correlations 
  
end








# this should run only on master process
function volumes_assemble(conf::Configuration, chains_to_assemble::Int64)

  volumes_all_chains = zeros(Float64, 1, chains_to_assemble) 
  volumes_sq_all_chains = zeros(Float64, 1, chains_to_assemble)    
  
  volumes_spread = zeros(Float64, 1)  
    
  for id_chain=1:chains_to_assemble
    
    @load "$(conf.volumes_folder)/volumes_chain=$(id_chain).jld2" volumes 
    @load "$(conf.volumes_sq_folder)/volumes_sq_chain=$(id_chain).jld2" volumes_sq
      
    volumes_all_chains[:, id_chain] = volumes[:]
    volumes_sq_all_chains[:, id_chain] = volumes_sq[:]
      
  end # cycle in id_chain
      
  volumes_all_chains = sum(volumes_all_chains, dims = 2)
  volumes_sq_all_chains = sum(volumes_sq_all_chains, dims = 2)
  volumes_all_chains[:] ./= chains_to_assemble
  volumes_sq_all_chains[:] ./= chains_to_assemble
  
  volumes_all_chains[1] = round(volumes_all_chains[1], digits = 5)
  volumes_sq_all_chains[1] = round(volumes_sq_all_chains[1], digits = 5)
  
  
  volumes_all_chains = vec(volumes_all_chains) # otherwise dataframe it's impossible to create
  volumes_sq_all_chains = vec(volumes_sq_all_chains)  
  
  volumes_dataframe = DataFrame(to_rename = volumes_all_chains)
  volumes_sq_dataframe = DataFrame(to_rename = volumes_sq_all_chains)
  column_name = "j=$(conf.j)"
  rename!(volumes_dataframe, :to_rename => column_name) # julia is weird
  rename!(volumes_sq_dataframe, :to_rename => column_name)   
 
  volumes_table_name = "/volumes_$(chains_to_assemble)_chains_combined.csv"
  volumes_sq_table_name = "/volumes_sq_$(chains_to_assemble)_chains_combined.csv"
      
  volumes_table_full_path = conf.tables_folder*volumes_table_name
  volumes_sq_table_full_path = conf.tables_folder*volumes_sq_table_name
      
  CSV.write(volumes_table_full_path, volumes_dataframe)
  CSV.write(volumes_sq_table_full_path, volumes_sq_dataframe)
  
  # spread is here
  # I compute spread by combining squared and volumes which were previously combined between multiple chains 
  # for j=0.5 sometimes the difference is extremely close to zero but negative (degen. case), this is why the abs() is necessary
  volumes_spread[1] = sqrt(abs(volumes_sq_all_chains[1] - volumes_all_chains[1]^2)) 
  
  @save "$(conf.volumes_spread_folder)/volumes_spread_$(chains_to_assemble)_chains_combined.jld2" volumes_spread
  
  volumes_spread_dataframe = DataFrame(to_rename = volumes_spread)
  rename!(volumes_spread_dataframe, :to_rename => column_name)  
  
  volumes_spread_table_name = "/volumes_spread_$(chains_to_assemble)_chains_combined.csv"
  volumes_spread_table_full_path = conf.tables_folder*volumes_spread_table_name
  
  CSV.write(volumes_spread_table_full_path, volumes_spread_dataframe)
  
  return volumes_dataframe,volumes_spread_dataframe

end





function volumes_correlations_assemble(conf::Configuration, chains_to_assemble::Int64)

  volumes_table_name = "/volumes_$(chains_to_assemble)_chains_combined.csv" 
  volumes_table_full_path = conf.tables_folder*volumes_table_name
  volumes_all_chains = vec(Matrix(DataFrame(CSV.File(volumes_table_full_path))))

  volumes_pseudo_correlations_all_chains = zeros(Float64, 1, chains_to_assemble)   
  
  for id_chain=1:chains_to_assemble
  
    @load "$(conf.volumes_pseudo_correlations_folder)/volumes_pseudo_correlations_node_1=$(conf.volumes_correlations_node_1)_node_2=$(conf.volumes_correlations_node_2)_chain=$(id_chain).jld2" volumes_pseudo_correlations 
  
    volumes_pseudo_correlations_all_chains[:, id_chain] = volumes_pseudo_correlations[:]
    
  end 
  
  volumes_pseudo_correlations_all_chains = sum(volumes_pseudo_correlations_all_chains, dims = 2)
  volumes_pseudo_correlations_all_chains[:] ./= chains_to_assemble
  
  volumes_pseudo_correlations_all_chains = vec(volumes_pseudo_correlations_all_chains)
    
  volumes_pseudo_correlations_all_chains[1] = round(volumes_pseudo_correlations_all_chains[1], digits = 5)
   
  @load "$(conf.volumes_spread_folder)/volumes_spread_$(chains_to_assemble)_chains_combined.jld2" volumes_spread
    
  # now we have everything to compute full correlations
  
  volumes_correlations = zeros(Float64, 1)  
  
  # this makes sense only because all volumes and spreads are equal
  volumes_correlations[1] = (volumes_pseudo_correlations_all_chains[1] - volumes_all_chains[1]^2)/(volumes_spread[1]^2)
  
  @save "$(conf.volumes_correlations_folder)/volumes_correlations_$(chains_to_assemble)_chains_combined.jld2" volumes_correlations
  
  volumes_correlations_dataframe = DataFrame(to_rename = volumes_correlations)
  column_name = "j=$(conf.j)"
  rename!(volumes_correlations_dataframe, :to_rename => column_name)  
  
  volumes_correlations_table_name = "/volumes_correlations_node_1=$(conf.volumes_correlations_node_1)_node_2=$(conf.volumes_correlations_node_2)_$(chains_to_assemble)_chains_combined.csv"
  volumes_correlations_table_full_path = conf.tables_folder*volumes_correlations_table_name
  
  CSV.write(volumes_correlations_table_full_path, volumes_correlations_dataframe)
  
  return volumes_correlations_dataframe

end


