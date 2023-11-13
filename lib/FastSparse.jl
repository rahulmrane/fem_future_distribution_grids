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
        i_row::Vector{Int}
        i_col::Vector{Int}
        value::Vector{Float64}
        FastSparseMatrix(nelements::Int) = new(Vector{Int}(undef, 9*nelements), Vector{Int}(undef, 9*nelements), Vector{Complex{Float64}}(undef, 9*nelements))
    end

    # function add!(fsp::FastSparseMatrix, id::Int, nodes::Vector{Int}, matrix::StaticArraysCore.SMatrix{3, 3, Float64, 9})
    function add!(fsp::FastSparseMatrix, id::Int, nodes::Vector{Int}, matrix::Matrix{Float64})
        for (num, (index, element)) in enumerate(pairs(matrix))
            (i, j) = Tuple(index)
            fsp.i_row[(id-1)*9 + num] = nodes[i]
            fsp.i_col[(id-1)*9 + num] = nodes[j]
            fsp.value[(id-1)*9 + num] = element
        end
        return nothing
    end

    function update!(fsp::FastSparseMatrix, id::Int, matrix::Matrix{Float64})
        for (num, (index, element)) in enumerate(pairs(matrix))
            fsp.value[(id-1)*9 + num] = element
        end
        return nothing
    end

end # module FastSparse