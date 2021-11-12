function compute_density_matrix_function!(j, D, vertex, draws, number_of_draws, ampls, density_matrix, N, b, density_matrices_folder, chain_id, subsystem, number_of_nodes_in_subsystem, density_matrix_linear_dim, total_density_matrices_already_stored=0)

  # containers for contraction
  v1 = zeros(Float64, D)
  v2 = zeros(Float64, D)
  v3 = zeros(Float64, D)
  v4 = zeros(Float64, D)
  v5 = zeros(Float64, D)
 
  indices = zeros(Int64,20) 
  vector_unflatted_indices = zeros(Int, number_of_nodes_in_subsystem)
  comparison_vector = zeros(Int, number_of_nodes_in_subsystem)
 
  for flat_index_row = 1:density_matrix_linear_dim
  
      for n = 1:number_of_draws 
      
          from_flat_index_to_vector_indices!(D, vector_unflatted_indices, number_of_nodes_in_subsystem, flat_index_row)
          
          for k=1:number_of_nodes_in_subsystem
          @inbounds comparison_vector[k] = abs(vector_unflatted_indices[k] - draws[subsystem[k], n])
          end
          
          if maximum(comparison_vector) > 0
          continue
          end
          
          for i=1:20
          @inbounds indices[i] = draws[i, n]
          end

              for flat_index_column = flat_index_row:density_matrix_linear_dim

                  from_flat_index_to_vector_indices!(D, vector_unflatted_indices, number_of_nodes_in_subsystem, flat_index_column)

                  for k=1:number_of_nodes_in_subsystem
                  @inbounds indices[subsystem[k]] = vector_unflatted_indices[k] 
                  end
              
                    if (flat_index_column == flat_index_row) 
                    @inbounds density_matrix[flat_index_row, flat_index_column] += 1*draws[21,n] 
                    continue
                    end   
              
                  @inbounds density_matrix[flat_index_row, flat_index_column] += star_amplitude(D, vertex, v1, v2, v3, v4, v5, indices)*draws[21,n]/ampls[n]
              
              end
               
      end # draw cycle
      
  end # flat_index_row cycle
 
  density_matrix = (density_matrix + transpose(density_matrix) - Diagonal(density_matrix))/sum(Diagonal(density_matrix))
  
  @save "$(density_matrices_folder)/density_matrix_chain=$(chain_id + total_density_matrices_already_stored).jld2" density_matrix 

end







function compute_density_matrix_function_multithread!(j, D, vertex, draws, number_of_draws, ampls, density_matrix, N, b, density_matrices_folder, chain_id, subsystem, number_of_nodes_in_subsystem, density_matrix_linear_dim, total_density_matrices_already_stored=0)

  number_of_threads = Threads.nthreads()

  # containers for contraction
  v1 = [zeros(Float64, D) for i in 1:number_of_threads]
  v2 = [zeros(Float64, D) for i in 1:number_of_threads]
  v3 = [zeros(Float64, D) for i in 1:number_of_threads]
  v4 = [zeros(Float64, D) for i in 1:number_of_threads]
  v5 = [zeros(Float64, D) for i in 1:number_of_threads]
 
  indices = [zeros(Int64, 20)  for i in 1:number_of_threads] 
  vector_unflatted_indices = [zeros(Int64, number_of_nodes_in_subsystem)  for i in 1:number_of_threads] 
  comparison_vector = [zeros(Int64, number_of_nodes_in_subsystem)  for i in 1:number_of_threads] 

  Threads.@threads for flat_index_row = 1:density_matrix_linear_dim

  thread_id = Threads.threadid()
   
      for n = 1:number_of_draws 
      
          from_flat_index_to_vector_indices!(D, vector_unflatted_indices[thread_id], number_of_nodes_in_subsystem, flat_index_row)
          
          for k=1:number_of_nodes_in_subsystem
          @inbounds comparison_vector[thread_id][k] = abs(vector_unflatted_indices[thread_id][k] - draws[subsystem[k], n])
          end
          
          if maximum(comparison_vector[thread_id]) > 0
          continue
          end
          
          for i=1:20
          @inbounds indices[thread_id][i] = draws[i, n]
          end

              for flat_index_column = flat_index_row:density_matrix_linear_dim

                  from_flat_index_to_vector_indices!(D, vector_unflatted_indices[thread_id], number_of_nodes_in_subsystem, flat_index_column)

                  for k=1:number_of_nodes_in_subsystem
                  @inbounds indices[thread_id][subsystem[k]] = vector_unflatted_indices[thread_id][k] 
                  end
              
                    if (flat_index_column == flat_index_row) 
                    @inbounds density_matrix[flat_index_row, flat_index_column] += 1*draws[21,n] 
                    continue
                    end   
              
                  @inbounds density_matrix[flat_index_row, flat_index_column] += star_amplitude(D, vertex, v1[thread_id], v2[thread_id], v3[thread_id], v4[thread_id], v5[thread_id], indices[thread_id])*draws[21,n]/ampls[n]
              
              end
               
      end # draw cycle
      
  end # flat_index_row cycle
 
  density_matrix = (density_matrix + transpose(density_matrix) - Diagonal(density_matrix))/sum(Diagonal(density_matrix))
  
  @save "$(density_matrices_folder)/density_matrix_chain=$(chain_id + total_density_matrices_already_stored).jld2" density_matrix 

end







function entropy_assemble(conf::Configuration, chains_to_assemble::Int64, density_matrix_linear_dim)

  entropy_numerical_fluctuations = zeros(Float64, 3) 
   
  vec_of_entropy = zeros(Float64, chains_to_assemble)

  for id_chain=1:chains_to_assemble
    
    @load "$(conf.density_matrices_folder)/density_matrix_chain=$(id_chain).jld2" density_matrix 

    # it would be much better to diagonalize during the parallel phase, but for small matrices the difference is negligible 
    density_matrix_eigenvalues = eigvals(Symmetric(density_matrix))
    
    entropy = 0.0
  
      for i=1:density_matrix_linear_dim
      entropy -= density_matrix_eigenvalues[i]*log(Complex(density_matrix_eigenvalues[i]))
      end 
  
    # discard small imaginary values due to numerical errors
    entropy = real(entropy)
    
    vec_of_entropy[id_chain] = entropy

  end # cycle in id_chain
      
  # average_entropy 
  entropy_numerical_fluctuations[1] = mean(vec_of_entropy[:])  
  
  # std_dev_entropy
  entropy_numerical_fluctuations[2] = std(vec_of_entropy[:])    
  
  # number of combined chains
  entropy_numerical_fluctuations[3] = chains_to_assemble       

  entropy_dataframe = DataFrame(to_rename = entropy_numerical_fluctuations[1])
  entropy_numerical_fluctuations_dataframe = DataFrame(to_rename = entropy_numerical_fluctuations)
  column_name = "j=$(conf.j)"
  rename!(entropy_dataframe, :to_rename => column_name) # julia is weird   
  rename!(entropy_numerical_fluctuations_dataframe, :to_rename => column_name)  

  entropy_table_name = "/entropy_$(chains_to_assemble)_chains_combined.csv"
  entropy_numerical_fluctuations_table_name = "/entropy_numerical_fluctuations_$(chains_to_assemble)_chains_combined.csv"
     
  entropy_table_full_path = conf.tables_folder*entropy_table_name    
  entropy_numerical_fluctuations_full_path = conf.tables_folder*entropy_numerical_fluctuations_table_name
      
  CSV.write(entropy_table_full_path, entropy_dataframe)    
  CSV.write(entropy_numerical_fluctuations_full_path, entropy_numerical_fluctuations_dataframe)  
  
  return entropy_dataframe, entropy_numerical_fluctuations_dataframe
  
end





