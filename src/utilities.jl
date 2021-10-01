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























