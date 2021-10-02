function check_on_preliminary_parameters(data_folder_path::String, verbosity_random_walk::Int64, verbosity_flux::Int64)

  # data folder check
  if (isdir(data_folder_path) == false) error("data folder path does not exists") end    


  # check on printing but not assembling
  if (print_final_data_on_terminal == true && assemble_and_save_final_data == false)
  printstyled("warning: you chose not to combine chains and therefore final data will not be assembled. The request to print the latter on terminal will be ignored\n"; color = :red) 
  end 

end


# check input configuration
function check_configuration!(user_conf::Vector{Any}, data_folder_path::String, verbosity_random_walk::Int64, verbosity_flux::Int64, chain_id::Int64)
 
  j = user_conf[1]
  M = user_conf[2]
  N = user_conf[3]
  b = user_conf[4]
  σ = user_conf[5]
  
  random_walk = user_conf[6][1]
  add_chains = user_conf[6][2]
  
  compute_angles = user_conf[7][1] 
  compute_angles_correlations = user_conf[7][2] 
  
  compute_volumes = user_conf[8][1]   
  compute_volumes_correlations = user_conf[8][2]
  volumes_correlations_node_1 = user_conf[8][3]
  volumes_correlations_node_2 = user_conf[8][4]
  
  # set spin to float ending in .0 if it is integer
  if (typeof(j) == Int64) 
  user_conf[1] = convert(Float64, user_conf[1]) 
  j = convert(Float64, j) 
  end 
  
  # check on presence of vertex amplitude
  if (isfile("vertex_ampls/$(M)/vertex_j=$(j).jld2") == false) error("The vertex amplitude for $(M) model and j=$(j) does not exists in the current directory. Unless the latter was removed, try to run the code from inside the main folder") end
     
  # (j, M, N, b, σ)  check
  if (typeof(j) != Float64 || j < 0) error("Please assign spin as positive float in user_conf $(user_conf)\n") end
  check_model = (M == "BF" || M == "EPRL") 
  if (check_model == false) error("I don't understand the model you want to compute in user_conf $(user_conf)\n") end 
  if (typeof(N) != Int64 || N <= 0) error("I don't understand number of iterations in user_conf $(user_conf)\nPlease use a positive integer number") end   
  if (b > N) error("burnin is greater than number of iterations in user_conf $(user_conf)") end   
  if (typeof(σ) != Float64 || σ <= 0) error("Please assign standard deviation as positive float in user_conf $(user_conf)\n") end  
  
  # random walk flags check
  if (typeof(add_chains) != Bool) error("add_chains must have a bool value in user_conf $(user_conf)") end 
  if (add_chains == true && random_walk == false) error("You can't pretend to add chains if you don't do random walk in user_conf $(user_conf)") end   
  
  # angles flags check
  if (typeof(compute_angles) != Bool) error("compute_angles must have a bool value in user_conf $(user_conf)") end    
  if (typeof(compute_angles_correlations) != Bool) error("compute_angles_correlations must have a bool value in user_conf $(user_conf)") end   
  
  # volumes flags check  
  if (typeof(volumes_correlations_node_1) != Int64 || typeof(volumes_correlations_node_2) != Int64) error("In user_conf $(user_conf) the nodes for the volumes correlations must be integers") end     
  if (volumes_correlations_node_1 < 1 || volumes_correlations_node_1 > 20) error("In user_conf $(user_conf), node 1 for the volumes correlations must be between 1 and 20") end
  if (volumes_correlations_node_2 < 1 || volumes_correlations_node_2 > 20) error("In user_conf $(user_conf), node 2 for the volumes correlations must be between 1 and 20") end
    
end # end of function check_configuration!





function check_stored_operators!(user_conf::Vector{Any}, conf::Configuration, data_folder_path::String, verbosity_flux::Int64, chain_id::Int64, number_of_chains::Int64)

  j = user_conf[1]
  M = user_conf[2]
  volumes_correlations_node_1 = user_conf[8][3]
  volumes_correlations_node_2 = user_conf[8][4]  
   
  # DRAWS CHECK
   
  push!(user_conf[6], isfile("$(conf.draws_folder)/draws_chain=$(chain_id).jld2")) # this is user_conf[6][3] 
      
  # check if there's an equal number of ampls and draws  
  number_of_existing_draws = file_count(conf.draws_folder)
  number_of_existing_ampls = file_count(conf.ampls_folder) 
  if (number_of_existing_draws != number_of_existing_ampls) error("There are different numbers of stored draws and ampls of user_conf $(user_conf)") end  
  push!(user_conf, number_of_existing_draws) # this is user_conf[9]           
  
  if (verbosity_flux > 1 && chain_id == 1) println("Found $(number_of_existing_draws) draws previously stored for the $(M) config with j=$(j)") end
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  # ANGLES CHECK
  
  push!(user_conf[7], isfile("$(conf.angles_folder)/angles_chain=$(chain_id).jld2")) # this is user_conf[7][3] 
      
  number_of_existing_angles = file_count(conf.angles_folder)
  number_of_existing_angles_sq = file_count(conf.angles_sq_folder)
  number_of_existing_angles_pseudo_correlations = file_count(conf.angles_pseudo_correlations_folder)
  if (number_of_existing_angles != number_of_existing_angles_sq) error("There are different numbers of stored angles and angles sq of user_conf $(user_conf)") end                                 
  push!(user_conf, number_of_existing_angles) # this is user_conf[10] 
  
  push!(user_conf[7], isfile("$(conf.angles_pseudo_correlations_folder)/angles_pseudo_correlations_chain=$(chain_id).jld2")) # this is user_conf[7][4] 
  push!(user_conf, number_of_existing_angles_pseudo_correlations) # this is user_conf[11]
  
  if (conf.compute_angles_correlations == true && (number_of_existing_angles < number_of_chains) && conf.compute_angles == false) 
  error("for user_conf $(user_conf), you chose to compute angles correlations and not angles, but there are $(number_of_existing_angles) angles previously stored and you are using $(number_of_chains) chains. Since angles correlations require also the computation of angles, you should compute angles as well") 
  end
  
  if (conf.compute_angles_correlations == true && (isfile("$(conf.tables_folder)/angles_$(number_of_chains)_chains_combined.csv") == false) && conf.compute_angles == false && assemble_and_save_final_data == true) 
  error("for user_conf $(user_conf), you chose to compute angles correlations and not angles, then assemble the final data. The problem is that there are no previously combined angles with $(number_of_chains) chains, therefore you should compute angles as well and assemble final data for such a number of chains") 
  end  
  
  if (conf.random_walk == false && conf.compute_angles == true && number_of_existing_angles < number_of_chains && number_of_existing_angles != 0 && chain_id == 1) 
  println("warning: for user_conf $(user_conf) you chose to don't do RW and compute angles. Since there are $(number_of_existing_angles) angles previously stored and $(number_of_chains) chains in this run, $(number_of_chains - number_of_existing_angles) angles will be computed") 
  end
  
  if (conf.random_walk == false && conf.compute_angles_correlations == true && number_of_existing_angles_pseudo_correlations < number_of_chains && number_of_existing_angles_pseudo_correlations != 0 && chain_id == 1) 
  println("warning: for user_conf $(user_conf) you chose to don't do RW and compute angles_pseudo_correlations. Since there are $(number_of_existing_angles_pseudo_correlations) operators <A_n,A_m> previously stored and $(number_of_chains) chains in this run, $(number_of_chains - number_of_existing_angles_pseudo_correlations) operators <A_n,A_m> will be computed") 
  end

  if (verbosity_flux > 1 && chain_id == 1) println("Found $(number_of_existing_angles) angles previously stored for the $(M) config with j=$(j)") end
  if (verbosity_flux > 1 && chain_id == 1) println("Found $(number_of_existing_angles_pseudo_correlations) operators <A_n,A_m> previously stored for the $(M) config with j=$(j)") end
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------                    

  # VOLUMES CHECK
  
  push!(user_conf[8], isfile("$(conf.volumes_folder)/volumes_chain=$(chain_id).jld2")) # this is user_conf[8][5]    
      
  number_of_existing_volumes = file_count(conf.volumes_folder)
  number_of_existing_volumes_sq = file_count(conf.volumes_sq_folder) 
  number_of_existing_volumes_pseudo_correlations = file_count(conf.volumes_pseudo_correlations_folder)
  if (number_of_existing_volumes != number_of_existing_volumes_sq) error("There are different numbers of stored volumes and volumes sq of user_conf $(user_conf)") end   
  push!(user_conf, number_of_existing_volumes) # this is user_conf[12]  
  
  push!(user_conf[8], isfile("$(conf.volumes_pseudo_correlations_folder)/volumes_pseudo_correlations_node_1=$(conf.volumes_correlations_node_1)_node_2=$(conf.volumes_correlations_node_2)_chain=$(chain_id).jld2")) # this is user_conf[8][6] 
  push!(user_conf, number_of_existing_volumes_pseudo_correlations) # this is user_conf[13]       
  
  if (conf.compute_volumes_correlations == true && (number_of_existing_volumes < number_of_chains) && conf.compute_volumes == false) 
  error("for user_conf $(user_conf), you chose to compute volumes correlations and not volumes, but there are $(number_of_existing_volumes) volumes previously stored and you are using $(number_of_chains) chains. Since volumes correlations require also the computation of volumes, you should compute volumes as well") 
  end  
  
  if (conf.compute_volumes_correlations == true && (isfile("$(conf.tables_folder)/volumes_$(number_of_chains)_chains_combined.csv") == false) && conf.compute_volumes == false && assemble_and_save_final_data == true) 
  error("for user_conf $(user_conf), you chose to compute volumes correlations and not volumes, then assemble the final data. The problem is that there are no previously combined volumes with $(number_of_chains) chains, therefore you should compute volumes as well and assemble final data for such a number of chains") 
  end    
  
  if (conf.random_walk == false && conf.compute_volumes == true && number_of_existing_volumes < number_of_chains && number_of_existing_volumes != 0 && chain_id == 1) 
  println("warning: for user_conf $(user_conf) you chose to don't do RW and compute volumes. Since there are $(number_of_existing_volumes) volumes previously stored and $(number_of_chains) chains in this run, $(number_of_chains - number_of_existing_volumes) volumes will be computed") 
  end  
  
  if (conf.random_walk == false && conf.compute_volumes_correlations == true && number_of_existing_volumes_pseudo_correlations < number_of_chains && number_of_existing_volumes_pseudo_correlations != 0 && chain_id == 1) 
  println("warning: for user_conf $(user_conf) you chose to don't do RW and compute volumes_pseudo_correlations. Since there are $(number_of_existing_volumes_pseudo_correlations) operators <V_$(volumes_correlations_node_1),V_$(volumes_correlations_node_2)> previously stored and $(number_of_chains) chains in this run, $(number_of_chains - number_of_existing_volumes_pseudo_correlations) operators <V_$(volumes_correlations_node_1),V_$(volumes_correlations_node_2)> will be computed") 
  end   
  
  if (verbosity_flux > 1 && chain_id == 1) println("Found $(number_of_existing_volumes) volumes previously stored for the $(M) config with j=$(j)") end
  if (verbosity_flux > 1 && chain_id == 1) println("Found $(number_of_existing_volumes_pseudo_correlations) operators <V_$(volumes_correlations_node_1),V_$(volumes_correlations_node_2)> previously stored for the $(M) config with j=$(j)") end 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
  
  # if user doesn't want to do RW
  
  if (conf.random_walk == false)  
  if (number_of_existing_draws < number_of_chains) error("In user_conf $(user_conf) you said that you don't want to do random walk, but there are $(conf.total_draws_already_stored) draws already stored and you are using $(number_of_chains) chains") end
  end # if check on RW
  
              
end # end of function check_stored_operators!
