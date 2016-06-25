

type ColumnIterator
   A::Matrix
end
columns(A::Matrix) = ColumnIterator(A)
Base.start(::ColumnIterator) = 0::Int
Base.done(iter::ColumnIterator, state::Int) = (state == size(iter.A, 2))
Base.next(iter::ColumnIterator, state::Int) = iter.A[:, state+1], state+1
