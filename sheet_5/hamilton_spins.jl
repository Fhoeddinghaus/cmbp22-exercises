using SparseArrays
using LinearAlgebra

function reset_Js(N) 
    Js = [sparse(zeros(N,N)) for α in 1:3]; # 3 N×N matrices (α = x,y,z)
    return Js
end

# convert [n]₁₀ to [n]₂ (in reversed order of notation, that is also used for |n⟩)
n_dual(n) = digits(n, base=2, pad=N)

# Function to calculate |n⟩ ↦ |l⟩ for the action of SᵢSⱼ for α = x,y
l(n,i,j) = n + (1 - 2 * n_dual(n)[i]) * 2^(i-1) + (1 - 2 * n_dual(n)[j]) * 2^(j-1)

# prefactor λˣ for the x-component of the interaction of spins i and j of |n⟩ 
#    -JˣᵢⱼSˣᵢSˣⱼ|n⟩ ↦ λʸ|l⟩
function x_link_prefactor(Js, n, i, j)
    α = 1 # x =̂ 1
    return -1//4 * Js[α][i,j]
end

# prefactor λʸ for the y-component of the interaction of spins i and j of |n⟩ 
#    -JʸᵢⱼSʸᵢSʸⱼ|n⟩ ↦ λʸ|l⟩
function y_link_prefactor(Js, n, i, j)
    α = 2 # y =̂ 2
    # calculate the sign, +1 for nᵢ != nⱼ, -1 for nᵢ == nⱼ
    s = - (2 * n_dual(n)[i] - 1) * (2 * n_dual(n)[j] - 1)
    return - 1//4 * s * Js[α][i,j]
end

# prefactor λᶻ for the z-component of the interaction of spins i and j of |n⟩ 
#    -JᶻᵢⱼSᶻᵢSᶻⱼ|n⟩ ↦ λᶻ|l⟩
function z_link_prefactor(Js, n, i, j)
    α = 3 # z =̂ 3
    # calculate the sign, +1 for nᵢ == nⱼ, -1 for nᵢ != nⱼ
    s = (2 * n_dual(n)[i] - 1) * (2 * n_dual(n)[j] - 1)
    return -1//4 * s * Js[α][i,j]
end

function calculate_hamilton_matrix(Js, N)
    H̄ = sparse(zeros(2^N, 2^N)) # empty 2ᴺ×2ᴺ matrix
    for n in 0:(2^N-1)
        # for each |n⟩ set the non-zero matrix elements
        
        # 1. loop over all components α = x,y,z 
        for α in [1,2,3]
            # 2. loop over only non-zero links in Js
            # iterate over non-zero elements by using findnz()
            for (i, j, J) in zip(findnz(Js[α])...)
                
                if α == 1 # α =̂ x
                    m = l(n, i, j) # calculate |l⟩ from |n⟩ for i,j
                    # add to the matrix element ⟨l|H|n⟩ the x-link-term
                    H̄[m+1,n+1] += x_link_prefactor(Js, n, i, j)
                elseif α == 2 # α =̂ y
                    m = l(n, i, j) # calculate |l⟩ from |n⟩ for i,j
                    # add to the matrix element ⟨l|H|n⟩ the y-link-term
                    H̄[m+1,n+1] += y_link_prefactor(Js, n, i, j)
                else
                    # α =̂ z, diagonal element
                    H̄[n+1,n+1] += z_link_prefactor(Js, n, i, j)
                end
            end
        end
    end
    return H̄
end