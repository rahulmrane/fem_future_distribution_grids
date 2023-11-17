module FastSparse

    using StaticArrays

    export FastSparseMatrix, add!, update!

    """
    # FastSparseMatrix

    A struct that contains the data for a sparse matrix, so that it can be constructed efficiently.

    Attributes:
    - i_row: a vector of row indices
    - i_col: a vector of column indices
    - value: the values corresponding to the row and column index
    """
    mutable struct FastSparseMatrix
        i_row::Vector{Int64}
        i_col::Vector{Int64}
        value::Vector{Float64}
        FastSparseMatrix(nelements::Int64) = new(Vector{Int64}(undef, 9*nelements), Vector{Int64}(undef, 9*nelements), Vector{Float64}(undef, 9*nelements))
    end
    
    function add!(fsp::FastSparseMatrix, element_id::Int64, nodes::Vector{Int64}, matrix::StaticArraysCore.SMatrix{3, 3, Float64, 9})
        for (num, (index, element)) in enumerate(pairs(matrix))
            (i, j) = Tuple(index)
            fsp.i_row[(element_id-1)*9 + num] = nodes[i]
            fsp.i_col[(element_id-1)*9 + num] = nodes[j]
            fsp.value[(element_id-1)*9 + num] = element
        end
    end

end # module FastSparse