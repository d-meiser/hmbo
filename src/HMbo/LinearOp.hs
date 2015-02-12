{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module HMbo.LinearOp(
    zero,
    identity,
    getDim,
    Dim(unPositive),
    toDim,
    fromDim,
    mul,
    add,
    scale,
    kron
    ) where

type Amplitude = Double
newtype MatrixElement = MatrixElement Double
  deriving(Show)
newtype Dim = Dim { unPositive :: Int}
  deriving(Eq, Show, Ord, Num)

data SparseMatrixEntry = SparseMatrixEntry Int Int Amplitude
  deriving(Show, Eq)

data LinearOp = Kron Dim LinearOp LinearOp
              | Plus Dim LinearOp LinearOp
              | KetBra Dim SparseMatrixEntry
              | ScaledId Dim Amplitude
              | Null Dim
  deriving(Show, Eq)

zero :: Dim -> LinearOp
zero d = Null d
identity :: Dim -> LinearOp
identity d = ScaledId d 1.0

getDim :: LinearOp -> Dim
getDim (Kron d _ _) = d
getDim (Plus d _ _) = d
getDim (KetBra d _) = d
getDim (ScaledId d _) = d
getDim (Null d) = d

mul :: LinearOp -> LinearOp -> Maybe LinearOp
mul op1 op2 | d1 == d2 = Just $ mul' op1 op2
            | otherwise = Nothing
            where
              d1 = getDim op1
              d2 = getDim op2
mul' :: LinearOp -> LinearOp -> LinearOp
mul' op@(Null _) _ = op
mul' _ op@(Null _) = op
mul' (ScaledId _ a) op = scale a op
mul' op (ScaledId _ a) = scale a op
mul' (Plus d op1 op2) op = Plus d (mul' op1 op) (mul' op2 op)
mul' op (Plus d op1 op2) = Plus d (mul' op op1) (mul' op op2)
mul' (KetBra d m1) (KetBra _ m2) =
  case theProduct of
    Just m3 -> KetBra d m3
    Nothing -> Null d
    where
      mProd :: SparseMatrixEntry -> SparseMatrixEntry -> Maybe SparseMatrixEntry
      mProd (SparseMatrixEntry r1 c1 a1) (SparseMatrixEntry r2 c2 a2)
        | c1 == r2 = Just $ SparseMatrixEntry r1 c2 (a1 * a2)
        | otherwise = Nothing
      theProduct = mProd m1 m2



add :: LinearOp -> LinearOp -> Maybe LinearOp
add op1 op2 | d1 == d2 = Just $ Plus d1 op1 op2
            | otherwise = Nothing
            where
              d1 = getDim op1
              d2 = getDim op2

scale :: Amplitude -> LinearOp -> LinearOp
scale a (Kron d op1 op2) = Kron d (scale a op1) op2
scale a (Plus d op1 op2) = Plus d (scale a op1) (scale a op2)
scale a (KetBra d (SparseMatrixEntry r c me)) = KetBra d $
  SparseMatrixEntry r c (a * me)
scale a (ScaledId d me) = ScaledId d (a * me)
scale _ n@(Null _) = n

kron :: LinearOp -> LinearOp -> LinearOp
kron (Null d) op = Null (d * (getDim op))
kron op (Null d) = Null (d * (getDim op))
kron (ScaledId d1 a1) (ScaledId d2 a2) = ScaledId (d1 * d2) (a1 * a2)
kron op1 op2 = Kron ((getDim op1) * (getDim op2)) op1 op2

toDim :: Int -> Maybe Dim
toDim d = if (d < 1) then Nothing else Just (Dim d)

fromDim :: Dim -> Int
fromDim (Dim d) = d
