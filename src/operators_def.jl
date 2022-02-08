#"Precompute the angles cosine operatore in the intw basis for given j"
function angle_vector(j::Float64, D::Int64, d::Int64)
   
    A = Dict{HalfInteger, Vector{Float64}}()   # Dict{K,V}() constructs a hash table with keys of type K and values of type V
    A2 = Dict{HalfInteger, Vector{Float64}}() 
          
    Angles = Vector{Float64}(undef, D)  
    Angles2 = Vector{Float64}(undef, D)     
        
    for i = 0:d
    Angles[i+1] = (i*(i+1) - 2j*(j+1)) / (2j*(j+1))
    Angles2[i+1] = Angles[i+1]^2        
    end
        
    A[j] = Angles[:]
    A2[j] = Angles2[:]       
    
 return A,A2
   
end




#"Precompute the matrices Vj and V2j (vol squared) in the intw basis for all js for the volume operator"
function volume_matrix(js, D::Int64, d::Int64)
   
    V = Dict{HalfInteger, Matrix{Float64}}()   # Dict{K,V}() constructs a hash table with keys of type K and values of type V
    V2 = Dict{HalfInteger, Matrix{Float64}}()
   
    for j in js
      
        Vj = Matrix{Complex{Float64}}(undef, D, D)       # Matrix{T}(undef, m, n) Construct an uninitialized Matrix{T} of size m×n. See undef.
        V2j = Matrix{Complex{Float64}}(undef, D, D)

        # build A matrix
        ud = zeros(Complex{Float64}, d)  
        for k in 1:d
            ud[k] =  -1im * (1/4) * (k^2 * (D^2 - k^2) / sqrt(4*(k^2) - 1))
        end
        ld = .-ud                            
        zd = zeros(Complex{Float64}, D)
        A = convert(Array, Tridiagonal(ld, zd, ud)) 
        
        # eigvals and vectors
        decomp = eigen(Hermitian(A))
        
        qs = decomp.values
        Qs = decomp.vectors 
       
        for k in 1:D
            for kp in 1:D

                s = zero(Complex{Float64})
                s2 = zero(Complex{Float64})

                for i in 1:D
                    v = sqrt(abs(qs[i]))
                    v2 = abs(qs[i])
                    s = s + v * Qs[k,i] * conj(Qs[kp, i])
                    s2 = s2 + v2 * Qs[k,i] * conj(Qs[kp, i])
                end

                Vj[k, kp] = s
                V2j[k, kp] = s2

            end
        end
       
        V[j] = copy(Float64.(real(Vj))) # discard small imaginary values due to numerical errors
        V2[j] = copy(Float64.(real(V2j))) # this is equal to Vj squared 
    end
   
    V 
    V2
    
    return V,V2
   
end
