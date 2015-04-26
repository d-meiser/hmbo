\section{Linear Operators}
\label{sec:LinearOperators}

\begin{code}
module HMbo.LinearOp(
    LinearOp,
    zero,
    eye,
    sigmaX,
    sigmaY,
    sigmaZ,
    sigmaPlus,
    sigmaMinus,
    ketBra,
    getDim,
    toDim,
    fromDim,
    mul,
    add,
    scale,
    kron,
    apply,
    isZero,
    SimpleOperator,
    simplify,
    showSO,
    sApply,
    checkSimpleOp
    ) where

import HMbo.Dim
import HMbo.Amplitude
import HMbo.Ket
import qualified Data.Vector.Unboxed as VU
import Data.Complex
import Data.List (intercalate,foldl')
import Control.Monad (liftM2)
import Data.Maybe ()

data SparseMatrixEntry = SparseMatrixEntry Int Int Amplitude
  deriving(Show, Eq)

data LinearOp = Kron Dim LinearOp LinearOp
              | Plus Dim [LinearOp]
              | KetBra Dim SparseMatrixEntry
              | ScaledId Dim Amplitude
  deriving(Show, Eq)

zero :: Dim -> LinearOp
zero d = Plus d []
eye :: Dim -> LinearOp
eye d = ScaledId d 1.0
ketBra:: Dim -> Int -> Int -> Amplitude -> Maybe LinearOp
ketBra d i j a | i >= 0 && i < d' && j >= 0 && j < d' = Just $
                  KetBra d (SparseMatrixEntry i j a)
               | otherwise = Nothing
  where
    d' = fromDim d

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

mul' op1@(KetBra dim _) op2@(Kron{}) = Plus dim products
  where
    products = filter isNonZero
      [mul' op1 op2' | op2' <-
        [KetBra dim entry | entry <- toSparseEntries op2]]
    isNonZero = not . isZero

mul' op1@(Kron dim _ _ ) op2@(KetBra _ _) = Plus dim products
  where
    products = filter isNonZero
      [mul' op1' op2 | op1' <-
        [KetBra dim entry | entry <- toSparseEntries op1]]
    isNonZero = not . isZero

toSparseEntries :: LinearOp -> [SparseMatrixEntry]
toSparseEntries (Kron _ op1 op2) = [combine se1 se2 |
                                     se1 <- toSparseEntries op1,
                                     se2 <- toSparseEntries op2]
  where
    combine :: SparseMatrixEntry -> SparseMatrixEntry -> SparseMatrixEntry
    combine
      (SparseMatrixEntry r1 c1 a1)
      (SparseMatrixEntry r2 c2 a2) = SparseMatrixEntry (r1 * d2 + r2) (c1 * d2 + c2) (a1 * a2)
      where
        d2 = fromDim $ getDim op2
toSparseEntries (Plus _ ops) = concatMap toSparseEntries ops
toSparseEntries (KetBra _ op) = [op]
toSparseEntries (ScaledId d a) = [SparseMatrixEntry i i a | i <- [0..(fromDim d - 1)]]

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
kron op1 op2 | isZero op1 || isZero op2 = Plus (getDim op1 * getDim op2) []
             | otherwise = Kron (getDim op1 * getDim op2) op1 op2

apply :: LinearOp -> Ket -> Maybe Ket
apply (ScaledId d a) x | fromDim d == VU.length x = Just $ VU.map (a *) x
                       | otherwise = Nothing
apply (KetBra d (SparseMatrixEntry i j a)) x
  | fromDim d == VU.length x =
      Just $  VU.generate (fromDim d) yfunc
  | otherwise = Nothing
    where
      yfunc :: Int -> Amplitude
      yfunc i' | i' == i = a * (VU.!) x j
               | otherwise = 0
apply (Plus d ops) x | fromDim d == VU.length x =
                        foldl' addVecs
                          (Just zeroVec)
                          (map (`apply` x) ops)
                     | otherwise = Nothing
     where
      addVecs :: Maybe Ket -> Maybe Ket -> Maybe Ket
      addVecs _ Nothing = Nothing
      addVecs Nothing _ = Nothing
      addVecs (Just a) (Just b) = Just $ VU.zipWith (+) a b
      zeroVec = VU.replicate (fromDim d) 0

apply (Kron _ op1 op2) x = do
    xp <- applyOne op2 d2 x
    xpp <- applyOne op1 d1 (transpose d1 xp)
    return (transpose d2 xpp)
  where
    d1 = fromDim $ getDim op1
    d2 = fromDim $ getDim op2
    subVectors :: Int -> Ket -> [Ket]
    subVectors d v = [VU.slice (d * j) d v | j <- [0..(numVecs - 1)]]
      where
        numVecs = VU.length v `div` d
    applyOne :: LinearOp -> Int -> Ket -> Maybe Ket
    applyOne op d v = VU.concat `fmap` mapM (apply op) (subVectors d v)


sigmaX :: LinearOp
Just (Just sigmaX) = liftM2 add (ketBra d 0 1 1.0) (ketBra d 1 0 1.0)
  where
    Just d = toDim 2

sigmaY :: LinearOp
Just (Just sigmaY) = liftM2 add (ketBra d 0 1 (0.0 :+ 1.0))
                                (ketBra d 1 0 (0.0 :+ (-1.0)))
  where
    Just d = toDim 2

sigmaZ :: LinearOp
Just (Just sigmaZ) = liftM2 add (ketBra d 1 1 1.0) (ketBra d 0 0 (-1.0))
  where
    Just d = toDim 2

sigmaPlus :: LinearOp
Just sigmaPlus = ketBra d 1 0 1.0
  where
    Just d = toDim 2

sigmaMinus :: LinearOp
Just sigmaMinus = ketBra d 0 1 1.0
  where
    Just d = toDim 2


data SimpleOperator = SimpleOperator Amplitude [Slice]
  deriving (Show)

data Slice = IdentityMatrix Dim
           | NonZeroLocation Dim Int Int
  deriving (Show)

simplify :: LinearOp -> [SimpleOperator]
simplify (ScaledId d a) = [SimpleOperator a [IdentityMatrix d]]
simplify (KetBra d (SparseMatrixEntry i j a)) =
  [SimpleOperator a [NonZeroLocation d i j]]
simplify (Plus _ ops) = concatMap simplify ops
simplify (Kron _ op1 op2) = map mergePairs
    [simpleKron sop1 sop2 | sop1 <- simplify op1 , sop2 <- simplify op2]
  where
    simpleKron (SimpleOperator a1 s1) (SimpleOperator a2 s2) =
      SimpleOperator (a1 * a2) (s1 ++ s2)
    mergePairs :: SimpleOperator -> SimpleOperator
    mergePairs (SimpleOperator a slices) = SimpleOperator a $
      mergePairs' slices
    mergePairs' ((IdentityMatrix d1):(IdentityMatrix d2):rest) =
      mergePairs' ((IdentityMatrix (d1 * d2)):rest)
    mergePairs' ((NonZeroLocation d1 i1 j1):(NonZeroLocation d2 i2 j2):rest) =
      mergePairs'
        ((NonZeroLocation (d1 * d2) (i1 * d2' + i2) (j1 * d2' + j2)):rest)
      where
        d2' = fromDim d2
    mergePairs' [] = []
    mergePairs' (s:ss) = s:(mergePairs' ss)


sApply :: [SimpleOperator] -> Ket -> Maybe Ket
sApply ops k
  | dimsOk = Just $ (foldl' vAdd vZero) (map ((flip sApply') k) ops)
  | otherwise = Nothing
  where
    vAdd :: Ket -> Ket -> Ket
    vAdd = VU.zipWith (+)
    vZero = VU.replicate dim 0
    dim = VU.length k
    dimsOk :: Bool
    dimsOk = and $ map ((VU.length k ==) . totalDim) ops
    totalDim :: SimpleOperator -> Int
    totalDim (SimpleOperator _ ss) = product (map sDim ss)

-- | Checks whether a simple operator can be applied to a Ket
checkSimpleOp :: SimpleOperator -> Ket -> Bool
checkSimpleOp (SimpleOperator _ []) k
  | VU.length k == 1 = True
  | otherwise = False
checkSimpleOp (SimpleOperator a (s:ss)) k
  | totalDim == VU.length k =
    case s of
      IdentityMatrix _ ->
        checkSimpleOp (SimpleOperator a ss) (VU.take blockSize k)
      NonZeroLocation _ i j ->
        i >= 0 && i < d && j >= 0 && j < d &&
          checkSimpleOp (SimpleOperator a ss) (VU.take blockSize k)
  | otherwise = False
  where
    blockSize = product (map sDim ss)
    d = sDim s
    totalDim = d * blockSize

sDim :: Slice -> Int
sDim (IdentityMatrix d') = fromDim d'
sDim (NonZeroLocation d' _ _) = fromDim d'
sApply' :: SimpleOperator -> Ket -> Ket
sApply' (SimpleOperator a []) k = VU.map (a *) k
sApply' (SimpleOperator a [s]) k =
  case s of
    IdentityMatrix _ -> VU.map (a *) k
    NonZeroLocation _ i j ->
      VU.concat
        [VU.replicate i 0
        ,VU.singleton (a * ((VU.!) k j))
        ,VU.replicate (d - (i + 1)) 0
        ]
  where
    d = sDim s
sApply' (SimpleOperator a (s:ss)) k =
  case s of
    IdentityMatrix _ ->
      VU.concat $
        map (sApply' (SimpleOperator a ss)) $
        [VU.slice (i * blockSize) blockSize k | i <- [0..(d - 1)]]
    NonZeroLocation _ i j ->
        VU.concat $
        [VU.replicate (i * blockSize) 0
        ,sApply'
           (SimpleOperator a ss)
           (VU.slice (j * blockSize) blockSize k)
        ,VU.replicate (totalDim - (i + 1) * blockSize) 0
        ]
  where
    blockSize = product (map sDim ss)
    d = sDim s
    totalDim = d * blockSize

showSO :: SimpleOperator -> String
showSO (SimpleOperator a slices) = show a ++ " (" ++
  intercalate " X " (map showSlice slices) ++ ")"
showSlice :: Slice -> String
showSlice (IdentityMatrix d) = "I(" ++ show (fromDim d) ++ ")"
showSlice (NonZeroLocation d i j) =
  "|" ++ show i ++ "><" ++ show j ++ "|_(" ++ show (fromDim d) ++ ")"

\end{code}
