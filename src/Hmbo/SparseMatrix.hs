module Hmbo.SparseMatrix 
    (zeroMatrix
    ,numEntries
    ) where

data Entry = Entry Int Int Double

newtype SparseMatrix = SparseMatrix [Entry]

zeroMatrix :: SparseMatrix
zeroMatrix = SparseMatrix []

numEntries :: SparseMatrix -> Int
numEntries (SparseMatrix entries) = length entries
