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
    kron,
    apply,
    identityMatrix,
    isZero,
    transpose,
    Ket
    ) where


import Data.Maybe
import qualified Data.Vector.Unboxed as VU
type Ket = VU.Vector Amplitude

type Amplitude = Double
newtype MatrixElement = MatrixElement Double
  deriving(Show)
newtype Dim = Dim { unPositive :: Int}
  deriving(Eq, Show, Ord, Num)

data SparseMatrixEntry = SparseMatrixEntry Int Int Amplitude
  deriving(Show, Eq)

data LinearOp = Kron Dim LinearOp LinearOp
              | Plus Dim [LinearOp]
              | KetBra Dim SparseMatrixEntry
              | ScaledId Dim Amplitude
  deriving(Show, Eq)

zero :: Dim -> LinearOp
zero d = Plus d []
identity :: Dim -> LinearOp
identity d = ScaledId d 1.0

getDim :: LinearOp -> Dim
getDim (Kron d _ _) = d
getDim (Plus d _) = d
getDim (KetBra d _) = d
getDim (ScaledId d _) = d

mul :: LinearOp -> LinearOp -> Maybe LinearOp
mul op1 op2 | d1 == d2 = Just $ mul' op1 op2
            | otherwise = Nothing
            where
              d1 = getDim op1
              d2 = getDim op2
mul' :: LinearOp -> LinearOp -> LinearOp
mul' (ScaledId _ a) op = scale a op
mul' op (ScaledId _ a) = scale a op
mul' (Plus d ops) op = Plus d (map ($ op) (fmap mul' ops))
mul' op (Plus d ops) = Plus d (map (mul' op) ops)
mul' (KetBra d m1) (KetBra _ m2) =
  case theProduct of
    Just m3 -> KetBra d m3
    Nothing -> Plus d []
    where
      mProd :: SparseMatrixEntry -> SparseMatrixEntry -> Maybe SparseMatrixEntry
      mProd (SparseMatrixEntry r1 c1 a1) (SparseMatrixEntry r2 c2 a2)
        | c1 == r2 = Just $ SparseMatrixEntry r1 c2 (a1 * a2)
        | otherwise = Nothing
      theProduct = mProd m1 m2
mul' (Kron dim a b) (Kron _ c d) 
        | dimsCompatible = Kron  dim (mul' a c) (mul' b d)
        | otherwise = Plus dim []
  where
    dimsCompatible = getDim a == getDim c && getDim b == getDim d

mul' op1@(KetBra dim _) op2@(Kron _ _ _) = Plus dim products
  where
    products = filter isNonZero $
      [mul' op1 op2' | op2' <- 
        [KetBra dim entry | entry <- (toSparseEntries op2)]]
    isNonZero = not . isZero

mul' op1@(Kron dim _ _ ) op2@(KetBra _ _) = Plus dim products
  where
    products = filter isNonZero $
      [mul' op1' op2 | op1' <- 
        [KetBra dim entry | entry <- (toSparseEntries op1)]]
    isNonZero = not . isZero

toSparseEntries :: LinearOp -> [SparseMatrixEntry]
toSparseEntries (Kron _ op1 op2) = [combine se1 se2 |
                                     se1 <- (toSparseEntries op1),
                                     se2 <- (toSparseEntries op2)]
  where
    combine :: SparseMatrixEntry -> SparseMatrixEntry -> SparseMatrixEntry
    combine
      (SparseMatrixEntry r1 c1 a1)
      (SparseMatrixEntry r2 c2 a2) = SparseMatrixEntry (r1 * d2 + r2) (c1 * d2 + c2) (a1 * a2)
      where
        d2 = fromDim $ getDim op2
toSparseEntries (Plus _ ops) = concat $ map toSparseEntries ops
toSparseEntries (KetBra _ op) = [op]
toSparseEntries (ScaledId d a) = [SparseMatrixEntry i i a | i <- [0..((fromDim d) - 1)]]

identityMatrix :: Dim -> LinearOp
identityMatrix d = ScaledId d 1.0

isZero ::LinearOp -> Bool
isZero (Plus _ []) = True
isZero _ = False

add :: LinearOp -> LinearOp -> Maybe LinearOp
add op1 op2 | d1 == d2 = Just $ Plus d1 [op1, op2]
            | otherwise = Nothing
            where
              d1 = getDim op1
              d2 = getDim op2

scale :: Amplitude -> LinearOp -> LinearOp
scale a (Kron d op1 op2) = Kron d (scale a op1) op2
scale a (Plus d ops) = Plus d (map (scale a) ops)
scale a (KetBra d (SparseMatrixEntry r c me)) = KetBra d $
  SparseMatrixEntry r c (a * me)
scale a (ScaledId d me) = ScaledId d (a * me)

kron :: LinearOp -> LinearOp -> LinearOp
kron (ScaledId d1 a1) (ScaledId d2 a2) = ScaledId (d1 * d2) (a1 * a2)
kron op1 op2 | isZero op1 || isZero op2 = Plus ((getDim op1) * (getDim op2)) []
             | otherwise = Kron ((getDim op1) * (getDim op2)) op1 op2

toDim :: Int -> Maybe Dim
toDim d = if (d < 1) then Nothing else Just (Dim d)

fromDim :: Dim -> Int
fromDim (Dim d) = d

apply :: LinearOp -> Ket -> Maybe Ket
apply (ScaledId d a) x | fromDim d == VU.length x = Just $ VU.map (a *) x
                       | otherwise = Nothing
apply (KetBra d (SparseMatrixEntry i j a)) x
  | fromDim d == VU.length x =
      Just $  VU.generate (fromDim d) yfunc
  | otherwise = Nothing
    where
      yfunc :: Int -> Amplitude
      yfunc i' | i' == i = a * ((VU.!) x j)
               | otherwise = 0
apply (Plus d ops) x | fromDim d == VU.length x =
                        Just $ foldl addVecs zeroVec (map (\op -> fromJust $ apply op x) ops)
                     | otherwise = Nothing
     where
      addVecs :: Ket -> Ket -> Ket
      addVecs = VU.zipWith (+)
      zeroVec = VU.replicate (fromDim d) 0

apply (Kron _ op1 op2) x = Just $
    transpose d2 $
    applyOne op1 d1 $
    transpose d1 $
    applyOne op2 d2 x
  where
    d1 = fromDim $ getDim op1
    d2 = fromDim $ getDim op2
    subVectors :: Int -> Ket -> [Ket]
    subVectors d v = [VU.slice (d * j) d v | j <- [0..(numVecs - 1)]]
      where
        numVecs = (VU.length v) `div` d
    applyOne op d v = VU.concat $
      map (\x' -> fromJust $ apply op x') (subVectors d v)

transpose :: Int -> Ket -> Ket
transpose n v = VU.fromList $
                  [(VU.!) v (i * n + j) | j <- [0..(n - 1)], i <- [0..(m - 1)]]
  where
    m = (VU.length v) `div` n

