
function firstNByGroup(vect::AbstractVector, N::Integer)
    """
    firstNByGroupIndices(vec::AbstractVector, N::Integer)

    Return the indices of the first `N` occurrences of each unique element in `vec`.

    This is useful when you want to subset another array based on limited occurrences per group.

    # Arguments
    - `vec`: A vector of values to group by (e.g., transcription factors).
    - `N`: The maximum number of elements to select per unique group.

    # Returns
    - A vector of indices corresponding to the first `N` entries per group.
    """
    seen = Dict{eltype(vect), Int}()
    idxs = Int[]
    for (i, v) in enumerate(vect)
        seen[v] = get(seen, v, 0) + 1
        if seen[v] <= N
            push!(idxs, i)
        end
    end
    return idxs
end



# function firstNByGroup(vect::AbstractVector, N::Integer)
#     """
#     firstNByGroupMask(vect::AbstractVector, N::Integer)

#     Return a boolean mask selecting the first `N` occurrences of each unique element in `vec`.

#     This is useful for logical indexing when subsetting arrays with the same length as `vec`.

#     # Arguments
#     - `vect`: A vector of values to group by.
#     - `N`: The maximum number of elements to select per group.

#     # Returns
#     - A boolean mask vector with `true` at positions of the first `N` occurrences per group.
#     """
#     seen = Dict{eltype(vect), Int}()
#     mask = falses(length(vect))
#     for i in eachindex(vect)
#         seen[vect[i]] = get(seen, vect[i], 0) + 1
#         mask[i] = seen[vect[i]] <= N
#     end
#     return mask
# end



