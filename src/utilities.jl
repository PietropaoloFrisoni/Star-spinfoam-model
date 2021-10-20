# dimension of spin representations
dimension(x::Int64) = 2x 
dimension(x::HalfInteger) = twice(x) 
dimension(x::Float64) = round(Int64, 2x)



Dimension(x::Int64) = 2x + 1
Dimension(x::HalfInteger) = twice(x) + 1
Dimension(x::Float64) = round(Int64, 2x + 1)



# returns number of files in a folder
function file_count(folder::String)

  files_and_dirs = readdir(folder)     
  number_of_files = size(files_and_dirs)[1]
  return number_of_files
  
end



# maps flat to vector indices in the density matrix
function from_flat_index_to_vector_indices!(D, vector_of_unflatted_indices, number_of_nodes_in_subsystem, flat_index)

  flat_index -= 1
      
    for k in 1:number_of_nodes_in_subsystem
    vector_of_unflatted_indices[k] = flat_index % D + 1
    flat_index = div(flat_index, D)
    end

end         





















