using Pkg
Pkg.activate(".")

using LinearAlgebra, Random

# Construct random weighted adjacency matrix (returns a dense matrix)
function generate_adjacency(
    nv::Int, #Number of nodes in graph
    ne_max::Int, #Maximum number of "directed" edgesfrom a node
    we_max::Float64, #Maximum weight of an edge
)
    # Empty Adjacency Matrix
    A = zeros(nv,nv)

    for i = 1:nv
        # Generate a random permutation of the remaining nodes, and
        # select the first between 1 and ne_max (random value) nodes as
        # neighbors
        neighbors = vcat(1:i-1,i+1:nv)[randperm(nv-1)[1:rand(1:ne_max)]]

        # Calculate the number of edges originating from current node
        degree = length(neighbors)

        # Generate the weights of the edges
        weights = we_max*rand(degree)

        #Add the neighbors and weights to the adjacency graph
        for j = 1:degree
            A[neighbors[j],i] = A[i,neighbors[j]] = weights[j]
        end
    end

    return A
end

# One step cholesky decomposition (outer product formulation)
function one_step_cholesky(A, method_args=nothing)
    n = size(A,1)
    L11 = sum(A[:,1])
    B = -A[2:end,1]/L11
    A1 = zeros(n-1,n-1)
    for i = 1:(n-1)
        for j = i+1:(n-1)
            A1[i,j] = A1[j,i] = A[i+1,j+1] + A[i+1,1]*A[j+1,1]/L11
        end
    end

    return B, L11, A1
end

# One step of sampled cholesky decomposition
function one_step_sampled_cholesky(A, bern_p=0.5)
    n = size(A,1)
    L11 = sum(A[:,1])

    # If row sum is zero
    if abs(L11) < 1e-15
        B = zeros(n-1)
        A1 = copy(A[2:end,2:end])

        return B, 1.0, A1
    end

    # Else if row sum is not zero
    B = -A[2:end,1]/L11
    A1 = zeros(n-1,n-1)
    for i = 1:(n-1)
        for j = i+1:(n-1)
            if abs(A[i+1,j+1]) > 1e-15
                A1[i,j] = A1[j,i] = A[i+1,j+1] + A[i+1,1]*A[j+1,1]/L11
            elseif rand() <= bern_p
                A1[i,j] = A1[j,i] = A[i+1,1]*A[j+1,1]/(L11*bern_p)
            end
        end
    end

    return B, L11, A1
end

# Repeated cholesky decomposition
function chol(A; method=one_step_cholesky, method_args = nothing)
    n = size(A,1)
    C = zeros(n,n) + I
    D = zeros(size(A,1))
    A1 = copy(A)
    for i = 1:n
        B, L11, A1 = method(A1, method_args)
        C[i+1:end,i] = B
        D[i] = L11
    end
    return C, D
end

# Computes excess number of entries that are filled in relative to A in
# Cholesky decomposition term C
function fillin(A, C)
    count = 0
    for i = 1:size(A,1)
        for j = 1:i-1
            if A[i,j] < 1e-15 &&  abs(C[i,j]) > 1e-15
                count += 1
            end
        end
    end
    return count
end


# Comparisons
nv = 100      # Number of Nodes
ne_max = 5    # Max number of edges "originating" from a node
we_max = 2.0  # Maximum weight of an edge (uniformly selected from 0 to we_max)

# Adjacency Matrix
A = generate_adjacency(nv, ne_max, we_max)
sum(A,dims=2)           # Degree of nodes
sum( A .> 0.0, dims=2)  # Number of edges per node
nnz = sum( A .> 0.0)/2  # Number edges

# Standard and Sampled Cholesky Decomposition
# method_args is the probability of filling in a nonzero entry in A when
# computing the sampled Cholesky decomposition
Cd, Dd = chol(A);
Cs, Ds = chol(A, method=one_step_sampled_cholesky, method_args=0.01);

# Dd and Ds are the diagonal matrices of the Cholesky Decompositions
# Dd is from the standard cholesky decomposition
# Ds is from the sampled cholesky decomposition
hcat(Dd, Ds)

# Cd and Cs are the unit lower triangular matrices; Cd for
# Cd is from the standard cholesky decomposition
# Cs is from the sampled cholesky decomposition
round.(Cd, sigdigits=4)
round.(Cs, sigdigits=4)

# Calculates number of excess nonzeros in Cholesky vs in original adjacency
# Divided by the original number of edges
fillin(A,Cd)/nnz |> println
fillin(A,Cs)/nnz |> println

precond = Cs * diagm(sqrt.(Ds))

L = precond \ Cd
sort(Dd)
eigvals( diagm( sum(A,dims=2)[:,1]) - A)
eigvals(L * diagm(Dd) * L')
