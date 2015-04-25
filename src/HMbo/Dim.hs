{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module HMbo.Dim (
    Dim(unPositive),
    toDim,
    fromDim
    ) where

newtype Dim = Dim { unPositive :: Int}
  deriving(Eq, Show, Ord, Num)

toDim :: Int -> Maybe Dim
toDim d = if d < 1 then Nothing else Just (Dim d)

fromDim :: Dim -> Int
fromDim (Dim d) = d

