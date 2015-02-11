module HMbo.LinearOp(
    zero,
    getDim,
    Dim(unPositive),
    toDim,
    fromDim
    ) where

newtype Amplitude = Amplitude Double
  deriving(Show)
newtype MatrixElement = MatrixElement Double
  deriving(Show)
newtype Dim = Dim { unPositive :: Int}
  deriving(Eq, Show, Ord)

data LinearOp = Kron Dim LinearOp LinearOp
              | Plus Dim LinearOp LinearOp
              | KetBra Dim MatrixElement
              | ScaledId Dim Amplitude
              | Null Dim
  deriving(Show)

zero :: Dim -> LinearOp
zero d = Null d

getDim :: LinearOp -> Dim
getDim (Kron d _ _) = d
getDim (Plus d _ _) = d
getDim (KetBra d _) = d
getDim (ScaledId d _) = d
getDim (Null d) = d

toDim :: Int -> Maybe Dim
toDim d = if (d < 1) then Nothing else Just (Dim d)

fromDim :: Dim -> Int
fromDim (Dim d) = d
