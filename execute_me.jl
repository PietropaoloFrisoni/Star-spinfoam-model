current_folder = pwd() 

printstyled("\nM-H algorithm for star model with $(nprocs()) Markov chain(s)\n\n"; bold = true, color = :blue) 

using Distributed

if (nprocs() > length(Sys.cpu_info())) printstyled("WARNING: you are using more chains than available cores on this system. Performances will be affected\n"; bold = true, color = :red) end
sleep(2)

println("Precompiling files and packages...") 

@everywhere begin
include("configs_to_compute.jl")
include("inc/pkgs.jl")
include("src/init.jl")
include("src/random_walk.jl")
include("src/star_amplitude.jl")
include("src/utilities.jl")
include("src/check.jl")
include("src/operators_def.jl")
include("src/angles.jl")
include("src/volumes.jl")
end

println("done\n") 
sleep(1)
println("checking configurations to compute...") 

  @everywhere begin  
  
    chain_id = myid()
    number_of_chains = nprocs()
    number_conf = size(Configurations)[1] 
    
    for user_conf in Configurations  
    check_configuration!(user_conf, data_folder_path, verbosity_random_walk, verbosity_flux, chain_id)
    end # cycle on configurations    
      
  end # end everywhere
  
println("done\n") 
sleep(1)  

println("creating folders...")

    for user_conf in Configurations    
    conf = init_config(user_conf, data_folder_path)        
    make_folders(conf)
    end    
    
println("done\n") 

println("checking stored draws and operators...")   

  @everywhere begin  
    
    for user_conf in Configurations     
    conf = init_config(user_conf, data_folder_path)        
    check_stored_operators!(user_conf, conf, data_folder_path, verbosity_flux, chain_id, number_of_chains)  
    end # cycle on configurations    
      
  end # end everywhere

println("done\n") 
sleep(1)
println("-------------------------------------------------------------------------\n") 

@everywhere begin

  for user_conf in Configurations
  
    conf = init_config(user_conf, data_folder_path, user_conf[6][3], user_conf[9], user_conf[7][3], user_conf[7][4], user_conf[10], user_conf[11], user_conf[8][5], user_conf[8][6], user_conf[12], user_conf[13]) 
    
    if (chain_id == 1) 
    printstyled("Starting with configuration:\nM=$(conf.M), j=$(conf.j), N=$(conf.N), b=$(conf.b), σ=$(conf.σ)\n\n"; bold = true, color = :bold) 
    sleep(1)
    end
    
    @load "vertex_ampls/$(conf.M)/vertex_j=$(conf.j).jld2" vertex
    
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
  
    if (conf.random_walk == true) 
    
        if (chain_id == 1)
        printstyled("Starting the random walk...\n"; color = :blue) 
        end
    
        if (conf.add_chains == true)
            if (chain_id == 1) 
            println("There are $(conf.total_draws_already_stored) chains already stored for this configuration, and $(number_of_chains) will be added\n") 
            @time random_walk_function(conf.j, conf.D, conf.d, vertex, conf.N, conf.b, conf.σ, conf.draws_folder, conf.ampls_folder, chain_id, verbosity_random_walk, conf.total_draws_already_stored)  
            else
            random_walk_function(conf.j, conf.D, conf.d, vertex, conf.N, conf.b, conf.σ, conf.draws_folder, conf.ampls_folder, chain_id, verbosity_random_walk, conf.total_draws_already_stored)
            end   
        else  
            if (chain_id == 1) 
            @time random_walk_function(conf.j, conf.D, conf.d, vertex, conf.N, conf.b, conf.σ, conf.draws_folder, conf.ampls_folder, chain_id, verbosity_random_walk)  
            else
            random_walk_function(conf.j, conf.D, conf.d, vertex, conf.N, conf.b, conf.σ, conf.draws_folder, conf.ampls_folder, chain_id, verbosity_random_walk)
            end               
        end  # check on add chains
        
    end  # if on random walk
    
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
    if (conf.compute_angles == true)     
    
      if (chain_id == 1)
      printstyled("\nStarting computing angles <A_n> for n=1...20\n...\n"; color = :green) 
      end
    
        if (conf.random_walk == false && conf.angles_already_stored_for_this_chain == true)       
        # skip angles computation for this chain            
        else
       
          angles = Array{Float64}(undef, 20)  
          angles_sq = Array{Float64}(undef, 20) 
          angles_vector_values = angle_vector(conf.j, conf.D, conf.d)     
   
          @load "$(conf.draws_folder)/draws_chain=$(chain_id).jld2" draws 
          number_of_draws = size(draws)[1]   
        
          if (conf.add_chains == true)
              if (chain_id == 1) 
              println("There are $(conf.total_angles_already_stored) angles chains already stored for this configuration, and $(number_of_chains) will be added\n") 
              @time compute_angles_function!(conf.j_half_int, conf.D, draws, number_of_draws, angles, angles_sq, angles_vector_values, conf.N, conf.b, conf.angles_folder, conf.angles_sq_folder, chain_id, conf.total_angles_already_stored)
              else
              compute_angles_function!(conf.j_half_int, conf.D, draws, number_of_draws, angles, angles_sq, angles_vector_values, conf.N, conf.b, conf.angles_folder, conf.angles_sq_folder, chain_id, conf.total_angles_already_stored)
              end   
          else    
              if (chain_id == 1) 
              @time compute_angles_function!(conf.j_half_int, conf.D, draws, number_of_draws, angles, angles_sq, angles_vector_values, conf.N, conf.b, conf.angles_folder, conf.angles_sq_folder, chain_id)
              else
              compute_angles_function!(conf.j_half_int, conf.D, draws, number_of_draws, angles, angles_sq, angles_vector_values, conf.N, conf.b, conf.angles_folder, conf.angles_sq_folder, chain_id)
              end                  
          end  # check on add chains     
        
        end # check on previously computed angles for this chain and random walk false
        
      if (chain_id == 1)
      println("Done!\n")
      end
      
    end # if on angles computation

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
    
    if (conf.compute_angles_correlations == true) # only if user wants to compute angles correlations
    
        if (chain_id == 1)
        printstyled("\nStarting computing operators <A_n,A_m> for all possible combinations\n...\n"; color = :yellow) 
        end    
    
          if (conf.random_walk == false && conf.angles_pseudo_correlations_already_stored_for_this_chain == true)   
          # skip angles pseudo corr computation for this chain            
          else
    
            angles_pseudo_correlations = zeros(20,20)   
            angles_vector_values = angle_vector(conf.j, conf.D, conf.d)  
    
            @load "$(conf.draws_folder)/draws_chain=$(chain_id).jld2" draws 
            number_of_draws = size(draws)[1]   
            
            if (conf.add_chains == true)
                if (chain_id == 1) 
                println("There are $(conf.total_angles_pseudo_correlations_already_stored) operators <A_n,A_m> chains already stored for this configuration, and $(number_of_chains) will be added\n") 
                @time compute_angles_pseudo_correlations_function!(conf.j_half_int, conf.d, draws, number_of_draws, angles_pseudo_correlations, angles_vector_values, conf.N, conf.b, conf.angles_pseudo_correlations_folder, chain_id, conf.total_angles_pseudo_correlations_already_stored)
                else
                compute_angles_pseudo_correlations_function!(conf.j_half_int, conf.d, draws, number_of_draws, angles_pseudo_correlations, angles_vector_values, conf.N, conf.b, conf.angles_pseudo_correlations_folder, chain_id, conf.total_angles_pseudo_correlations_already_stored)
                end     
            else    
                if (chain_id == 1) 
                @time compute_angles_pseudo_correlations_function!(conf.j_half_int, conf.d, draws, number_of_draws, angles_pseudo_correlations, angles_vector_values, conf.N, conf.b, conf.angles_pseudo_correlations_folder, chain_id)
                else
                compute_angles_pseudo_correlations_function!(conf.j_half_int, conf.d, draws, number_of_draws, angles_pseudo_correlations, angles_vector_values, conf.N, conf.b, conf.angles_pseudo_correlations_folder, chain_id)
                end                   
            end  # check on add chains      
               
          end # check on previously computed angles pseudo correlations and random walk false    
            
          
        if (chain_id == 1)
        println("Done!\n")
        end            
  
    end # if on angles pseudo correlations
    
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
    
    if (conf.compute_volumes == true)   
    
        if (chain_id == 1)
        printstyled("\nStarting computing volumes <V_1>\n...\n"; color = :magenta) 
        end      
    
          if (conf.random_walk == false && conf.volumes_already_stored_for_this_chain == true)     
          # skip volumes computation for this chain
          else    
        
            volumes = Array{Float64}(undef, 1)    # julia is weird and it doesn't allow single scalar values to be modified unless they are unidimensional arrays 
            volumes_sq = Array{Float64}(undef, 1) 
            volumes_matrix_values = volume_matrix(conf.j_half_int, conf.D, conf.d)           
         
            @load "$(conf.draws_folder)/draws_chain=$(chain_id).jld2" draws 
       
            #This is such that draws has structure: [21,N], so the volumes will be computed faster  
            draws = transpose(draws)  
            number_of_draws = size(draws)[2]           
        
            @load "$(conf.ampls_folder)/ampls_chain=$(chain_id).jld2" ampls     
            @load "vertex_ampls/$(conf.M)/vertex_j=$(conf.j).jld2" vertex
        
            if (conf.add_chains == true)       
                if (chain_id == 1) 
                println("There are $(conf.total_volumes_already_stored) volumes chains already stored for this configuration, and $(number_of_chains) will be added\n") 
                @time compute_volumes_function!(conf.j_half_int, conf.D, vertex, draws, number_of_draws, ampls, volumes, volumes_sq, volumes_matrix_values, conf.N, conf.b, conf.volumes_folder, conf.volumes_sq_folder, chain_id, conf.total_volumes_already_stored)  
                else
                compute_volumes_function!(conf.j_half_int, conf.D, vertex, draws, number_of_draws, ampls, volumes, volumes_sq, volumes_matrix_values, conf.N, conf.b, conf.volumes_folder, conf.volumes_sq_folder, chain_id, conf.total_volumes_already_stored) 
                end   
            else      
                if (chain_id == 1) 
                @time compute_volumes_function!(conf.j_half_int, conf.D, vertex, draws, number_of_draws, ampls, volumes, volumes_sq, volumes_matrix_values, conf.N, conf.b, conf.volumes_folder, conf.volumes_sq_folder, chain_id)  
                else
                compute_volumes_function!(conf.j_half_int, conf.D, vertex, draws, number_of_draws, ampls, volumes, volumes_sq, volumes_matrix_values, conf.N, conf.b, conf.volumes_folder, conf.volumes_sq_folder, chain_id) 
                end               
            end  # check on add chains                         
        
          end # check on previously computed volumes and random walk false
          
        if (chain_id == 1)
        println("Done!\n")
        end              
    
    end # if on volumes computation
    
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
    
    if (conf.compute_volumes_correlations == true) # only if user wants to compute volumes correlations
    
        if (chain_id == 1)
        printstyled("\nStarting computing operators <V_$(conf.volumes_correlations_node_1),V_$(conf.volumes_correlations_node_2)>\n...\n"; color = :light_red) 
        end    
    
          if (conf.random_walk == false && conf.volumes_pseudo_correlations_already_stored_for_this_chain == true)      
          # skip volumes pseudo corr computation for this chain        
          else    
  
            volumes_pseudo_correlations = Array{Float64}(undef, 1)  
            volumes_matrix_values = volume_matrix(conf.j_half_int, conf.D, conf.d)  # I'm using F. old functions   
    
            @load "$(conf.draws_folder)/draws_chain=$(chain_id).jld2" draws 
       
            #This is such that draws has structure: [21,N], so the volumes correlations will be computed faster  
            draws = transpose(draws)  
            number_of_draws = size(draws)[2]           
        
            @load "$(conf.ampls_folder)/ampls_chain=$(chain_id).jld2" ampls 
            @load "vertex_ampls/$(conf.M)/vertex_j=$(conf.j).jld2" vertex
  
            if (conf.add_chains == true)
                if (chain_id == 1) 
                println("There are $(conf.total_volumes_pseudo_correlations_already_stored) operators <V_$(conf.volumes_correlations_node_1),V_$(conf.volumes_correlations_node_2)> chains already stored for this configuration, and $(number_of_chains) will be added\n")
                @time compute_volumes_pseudo_correlations_function!(conf.j_half_int, conf.D, draws, number_of_draws, ampls, vertex, volumes_pseudo_correlations, volumes_matrix_values, conf.N, conf.b, conf.volumes_pseudo_correlations_folder, chain_id, conf.volumes_correlations_node_1, conf.volumes_correlations_node_2, conf.total_volumes_pseudo_correlations_already_stored)
                else
                compute_volumes_pseudo_correlations_function!(conf.j_half_int, conf.D, draws, number_of_draws, ampls, vertex, volumes_pseudo_correlations, volumes_matrix_values, conf.N, conf.b, conf.volumes_pseudo_correlations_folder, chain_id, conf.volumes_correlations_node_1, conf.volumes_correlations_node_2, conf.total_volumes_pseudo_correlations_already_stored)
                end   
            else
                if (chain_id == 1) 
                @time compute_volumes_pseudo_correlations_function!(conf.j_half_int, conf.D, draws, number_of_draws, ampls, vertex, volumes_pseudo_correlations, volumes_matrix_values, conf.N, conf.b, conf.volumes_pseudo_correlations_folder, chain_id, conf.volumes_correlations_node_1, conf.volumes_correlations_node_2)
                else
                compute_volumes_pseudo_correlations_function!(conf.j_half_int, conf.D, draws, number_of_draws, ampls, vertex, volumes_pseudo_correlations, volumes_matrix_values, conf.N, conf.b, conf.volumes_pseudo_correlations_folder, chain_id, conf.volumes_correlations_node_1, conf.volumes_correlations_node_2)
                end              
            end  # check on add chains 
                        
          end # check on previously computed and random walk false      
  
        if (chain_id == 1)
        println("Done!\n")
        end      
  
    end # if on volumes correlations        
    
  if (chain_id == 1)
  println("\n-------------------------------------------------------------------------\n") 
  end      
     
  end # end conf cycle
  
end # end everywhere --- PARALLELIZATION ENDS HERE

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

  printstyled("All computations completed successfully!\n\n\n"; bold = true, color = :blue)
    
  user_BF_angles_table = DataFrame()
  user_BF_angles_spread_table = DataFrame()
  user_BF_volumes_table = DataFrame()
  user_BF_volumes_spread_table = DataFrame()  
  
  user_EPRL_angles_table = DataFrame()
  user_EPRL_angles_spread_table = DataFrame()  
  user_EPRL_volumes_table = DataFrame()
  user_EPRL_volumes_spread_table = DataFrame()
  
  user_BF_volumes_correlations_table = DataFrame()
  user_EPRL_volumes_correlations_table = DataFrame() 
    
  array = String[];
  for i=1:20, j=1:20 
  string = "C($(i),$(j))"
  push!(array, string)        
  end
  user_BF_angles_correlations_table = DataFrame(nodes = array)
  user_EPRL_angles_correlations_table = DataFrame(nodes = array)
 
  
  for user_conf in Configurations
  
    conf = init_config(user_conf, data_folder_path) 
    
    printstyled("\nStart assembling $(number_of_chains) chain(s) for:\nM=$(conf.M), j=$(conf.j), N=$(conf.N), b=$(conf.b), σ=$(conf.σ)\n"; color = :bold) 
    
      if (conf.add_chains == true)
      println("\nEven if you chose to add chains, only $(number_of_chains) chains are assembled for every operator. This will change in future updates!\n")        
      end        
  
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------             
            
            if (conf.compute_angles == true) 
            # printstyled("\nAssembling $(number_of_chains) angles <A_n> chains\n...\n"; color = :green) 
            angles_dataframes = angles_assemble(conf, number_of_chains)      
                if (conf.M == "BF")
            global user_BF_angles_table = hcat(user_BF_angles_table, angles_dataframes[1]) 
            global user_BF_angles_spread_table = hcat(user_BF_angles_spread_table, angles_dataframes[2]) 
                else
            global user_EPRL_angles_table = hcat(user_EPRL_angles_table, angles_dataframes[1]) 
            global user_EPRL_angles_spread_table = hcat(user_EPRL_angles_spread_table, angles_dataframes[2])          
                end 
            end # if on angle computation
    
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------             
    
            if (conf.compute_angles_correlations == true)
            # printstyled("\nAssembling $(number_of_chains) operators <A_n,A_m> chains\n...\n"; color = :yellow) 
            angles_correlations_dataframe = angles_correlations_assemble(conf, number_of_chains)         
                if (conf.M == "BF")  
            global user_BF_angles_correlations_table = hcat(user_BF_angles_correlations_table, angles_correlations_dataframe)     
                else
            global user_EPRL_angles_correlations_table = hcat(user_EPRL_angles_correlations_table, angles_correlations_dataframe)   
                end 
            end  # if on angle correlation computation
      
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------         

            if (conf.compute_volumes == true) 
            # printstyled("\nAssembling $(number_of_chains) volumes <V_1> chains\n...\n"; color = :magenta) 
            volumes_dataframes = volumes_assemble(conf, number_of_chains)      
                if (conf.M == "BF")
            global user_BF_volumes_table = hcat(user_BF_volumes_table, volumes_dataframes[1]) 
            global user_BF_volumes_spread_table = hcat(user_BF_volumes_spread_table, volumes_dataframes[2])     
                else 
            global user_EPRL_volumes_table = hcat(user_EPRL_volumes_table, volumes_dataframes[1]) 
            global user_EPRL_volumes_spread_table = hcat(user_EPRL_volumes_spread_table, volumes_dataframes[2])     
                end 
            end # if on angle computation      
      
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------   

            if (conf.compute_volumes_correlations == true)
            # printstyled("\nAssembling $(number_of_chains) operators <V_$(conf.volumes_correlations_node_1),V_$(conf.volumes_correlations_node_2)> chains\n...\n"; color = :light_red) 
            volumes_correlations_dataframe = volumes_correlations_assemble(conf, number_of_chains)
                if (conf.M == "BF")  
            global user_BF_volumes_correlations_table = hcat(user_BF_volumes_correlations_table, volumes_correlations_dataframe)
                else        
            global user_EPRL_volumes_correlations_table = hcat(user_EPRL_volumes_correlations_table, volumes_correlations_dataframe)          
                end 
            end  # if on volumes correlation computation      
                 
    println("Done!\n")  
    println("\n-------------------------------------------------------------------------")      
      
  end # end conf cycle

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

  current_date = now()
  store_final_data_path = "$(current_folder)/final_data/run_on:$(current_date)"
  mkpath(store_final_data_path)

  if (isempty(user_BF_angles_table) == false)
 
   CSV.write("$(store_final_data_path)/user_BF_angles_table.csv", user_BF_angles_table)   
   CSV.write("$(store_final_data_path)/user_BF_angles_spread_table.csv", user_BF_angles_spread_table)
 
   if (print_final_data_on_terminal == true)
     println("These are the values of the BF angles that you asked for:\n", user_BF_angles_table, "\n")
     println("...and these are the values of corresponding spread:\n", user_BF_angles_spread_table, "\n\n\n")
   end
   
  end
 
  if (isempty(user_EPRL_angles_table) == false)
 
   CSV.write("$(store_final_data_path)/user_EPRL_angles_table.csv", user_EPRL_angles_table)   
   CSV.write("$(store_final_data_path)/user_EPRL_angles_spread_table.csv", user_EPRL_angles_spread_table)
 
   if (print_final_data_on_terminal == true)
     println("These are the values of the EPRL angles that you asked for:\n", user_EPRL_angles_table, "\n")
     println("...and these are the values of corresponding spread:\n", user_EPRL_angles_spread_table, "\n\n\n")
   end
   
 end
 
  if (isempty(user_BF_volumes_table) == false)
 
   CSV.write("$(store_final_data_path)/user_BF_volumes_table.csv", user_BF_volumes_table)   
   CSV.write("$(store_final_data_path)/user_BF_volumes_spread_table.csv", user_BF_volumes_spread_table)
 
   if (print_final_data_on_terminal == true)
     println("These are the values of the BF volumes that you asked for:\n", user_BF_volumes_table, "\n")
     println("...and these are the values of corresponding spread:\n", user_BF_volumes_spread_table, "\n\n\n")
   end
  
 end
 
  if (isempty(user_EPRL_volumes_table) == false)
 
   CSV.write("$(store_final_data_path)/user_EPRL_volumes_table.csv", user_EPRL_volumes_table)   
   CSV.write("$(store_final_data_path)/user_EPRL_volumes_spread_table.csv", user_EPRL_volumes_spread_table)
 
   if (print_final_data_on_terminal == true)
     println("These are the values of the EPRL volumes that you asked for:\n", user_EPRL_volumes_table, "\n")
     println("...and these are the values of corresponding spread:\n", user_EPRL_volumes_spread_table, "\n\n\n")
   end
   
 end 
 
 if (size(user_BF_angles_correlations_table)[2] != 1)
 
   CSV.write("$(store_final_data_path)/user_BF_angles_correlations_table.csv", user_BF_angles_correlations_table)
 
   if (print_final_data_on_terminal == true)
     println("These are the values of the BF correlations that you asked for:\n", user_BF_angles_correlations_table, "\n")
   end
  
 end
 
  if (size(user_EPRL_angles_correlations_table)[2] != 1)
 
   CSV.write("$(store_final_data_path)/user_EPRL_angles_correlations_table.csv", user_EPRL_angles_correlations_table)
 
   if (print_final_data_on_terminal == true)
     println("These are the values of the EPRL correlations that you asked for:\n", user_EPRL_angles_correlations_table, "\n")
   end
 
 end
 
  if (isempty(user_BF_volumes_correlations_table) == false)
 
   CSV.write("$(store_final_data_path)/user_BF_volumes_correlations_table.csv", user_BF_volumes_correlations_table)   
   
   if (print_final_data_on_terminal == true)
     println("These are the values of the BF volumes correlations that you asked for:\n", user_BF_volumes_correlations_table, "\n")
   end
  
 end  
 
  if (isempty(user_EPRL_volumes_correlations_table) == false)
 
   CSV.write("$(store_final_data_path)/user_EPRL_volumes_correlations_table.csv", user_EPRL_volumes_correlations_table)   
  
   if (print_final_data_on_terminal == true)
     println("These are the values of the EPRL volumes correlations that you asked for:\n", user_EPRL_volumes_correlations_table, "\n")
   end
 
 end   
 
 # saving Configurations chosen by user to final_data
 configs_to_compute_file = current_folder*"/configs_to_compute.jl"
 Configurations_computed_in_this_run_file = store_final_data_path*"/Configurations_computed_in_this_run.txt"
 cp(configs_to_compute_file, Configurations_computed_in_this_run_file)

printstyled("The execution terminated successfully!\n\n"; color = :bold)

printstyled("The data have been saved in the 'final_data' folder\n\n"; bold = true, color = :blue) 
 
