\section{Many body operators}
\label{sec:ManyBodyOperators}

\ignore{
\begin{code}
module HMbo.ManyBodyOperator(
    ManyBodyOperator,
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

import qualified Data.Vector.Unboxed as VU
import Data.Complex
import Data.List (intercalate,foldl')
import Control.Monad (liftM2)
import Data.Maybe ()

import HMbo.Dim
import HMbo.Amplitude
import HMbo.Ket
\end{code}
}

Computations in quantum mechanical many body calculations
can be divided into two distinct classes.  On the one hand there are
algebraic computations at the operator level similar to what one does
with pencil and paper.  On the other hand there are numerical
computations in which the quantum mechanical operator is applied to a
state vector.  The latter quickly become much too time consuming for
humans to carry out.  They must be implemented with efficient computer
code.  Typically, the computational cost of the algebraic operations is
of polynomial order in the number $N$ of constituent elementary quantum
systems, e.g. $\mathcal{O}(N)$ or $\mathcal{O}(N^2)$.  The elementary
quantum systems might be individual spins or multi-level atoms.
However application of the operators to state vectors is of exponential
complexity, $\mathcal{O}(e^N)$, owing to the exponential scaling of the
dimension of the Hilbert space of the quantum system.

Importantly, the needs of these two classes of computation are very
different from one another.  For algebraic operator manipulations we
require an expressive set of operations for composing simpler operators
into more complex ones.  Computational efficiency is of secondary
importance for these manipulations.  For the second type of computations
we only need one operation: applying a many particle operator to a state
vector.  But computational efficiency is of utmost importance for this
latter type of computation.

In order to accommodate the needs of these two disparate types of
computations we introduce two separate concepts, namely that of an
algebraic operator and that of a numerical operator.  A {\tt
ManyBodyOperator} has both these representations so it can offer a
rich set of algebraic operations without sacrificing performance when
the {\tt ManyBodyOperator} is applied to a {\tt Ket}:
\begin{code}
data ManyBodyOperator = ManyBodyOperator { aOp :: AOperator
                                         , nOp :: NOperator
                                         }
  deriving (Show, Eq)
\end{code}
The type {\tt AOperator} is tailored towards algebraic operations {\tt
NOperator} is tailored towards numerical application of the operator to
state vectors.  {\tt NOperator}s are created from {\tt AOperator}s by
means of a compilation process.  The compilation can be performed with
computational cost comparable to other algebraic operations, e.g.
$\mathcal{O}(N)$.

The {\tt AOperator} is structured in such a way that it supports the
typical algebraic operations such as operator addition or outer
products:
\begin{code}
data AOperator = Kron Dim AOperator AOperator
               | Plus Dim [AOperator]
               | KetBra Dim SparseMatrixEntry
               | ScaledId Dim Amplitude
  deriving (Show, Eq)

data SparseMatrixEntry = SparseMatrixEntry Int Int Amplitude
  deriving(Show, Eq)
\end{code}

\begin{code}
data NOperator = NOperator
  deriving (Show, Eq)
\end{code}

\begin{code}

zero :: Dim -> ManyBodyOperator
zero d = ManyBodyOperator (Plus d []) NOperator

eye :: Dim -> ManyBodyOperator
eye d = ManyBodyOperator (ScaledId d 1.0) NOperator

ketBra:: Dim -> Int -> Int -> Amplitude -> Maybe ManyBodyOperator
ketBra d i j a | i >= 0 && i < d' && j >= 0 && j < d' = Just $
                  ManyBodyOperator
                    (KetBra d (SparseMatrixEntry i j a)) NOperator
               | otherwise = Nothing
  where
    d' = fromDim d

getDim :: ManyBodyOperator -> Dim
getDim (ManyBodyOperator aop _) = aGetDim aop

aGetDim :: AOperator -> Dim
aGetDim (Kron d _ _) = d
aGetDim (Plus d _) = d
aGetDim (KetBra d _) = d
aGetDim (ScaledId d _) = d

mul :: ManyBodyOperator -> ManyBodyOperator -> Maybe ManyBodyOperator
mul op1 op2 | d1 == d2 = Just $
              ManyBodyOperator (aMul (aOp op1) (aOp op2)) NOperator
            | otherwise = Nothing
            where
              d1 = getDim op1
              d2 = getDim op2
aMul :: AOperator -> AOperator -> AOperator
aMul (ScaledId _ a) op = aScale a op
aMul op (ScaledId _ a) = aScale a op
aMul (Plus d ops) op = Plus d (map ($ op) (fmap aMul ops))
aMul op (Plus d ops) = Plus d (map (aMul op) ops)
aMul (KetBra d m1) (KetBra _ m2) =
  case theProduct of
    Just m3 -> KetBra d m3
    Nothing -> Plus d []
    where
      mProd :: SparseMatrixEntry -> SparseMatrixEntry -> Maybe SparseMatrixEntry
      mProd (SparseMatrixEntry r1 c1 a1) (SparseMatrixEntry r2 c2 a2)
        | c1 == r2 = Just $ SparseMatrixEntry r1 c2 (a1 * a2)
        | otherwise = Nothing
      theProduct = mProd m1 m2
aMul (Kron dim a b) (Kron _ c d)
        | dimsCompatible = Kron  dim (aMul a c) (aMul b d)
        | otherwise = Plus dim []
  where
    dimsCompatible = aGetDim a == aGetDim c && aGetDim b == aGetDim d

aMul op1@(KetBra dim _) op2@(Kron{}) = Plus dim products
  where
    products = filter isNonZero
      [aMul op1 op2' | op2' <-
        [KetBra dim entry | entry <- toSparseEntries op2]]
    isNonZero = not . isZero

aMul op1@(Kron dim _ _ ) op2@(KetBra _ _) = Plus dim products
  where
    products = filter isNonZero
      [aMul op1' op2 | op1' <-
        [KetBra dim entry | entry <- toSparseEntries op1]]
    isNonZero = not . isZero

toSparseEntries :: AOperator -> [SparseMatrixEntry]
toSparseEntries (Kron _ op1 op2) = [combine se1 se2 |
                                     se1 <- toSparseEntries op1,
                                     se2 <- toSparseEntries op2]
  where
    combine :: SparseMatrixEntry -> SparseMatrixEntry -> SparseMatrixEntry
    combine
      (SparseMatrixEntry r1 c1 a1)
      (SparseMatrixEntry r2 c2 a2) = SparseMatrixEntry (r1 * d2 + r2) (c1 * d2 + c2) (a1 * a2)
      where
        d2 = fromDim $ aGetDim op2
toSparseEntries (Plus _ ops) = concatMap toSparseEntries ops
toSparseEntries (KetBra _ op) = [op]
toSparseEntries (ScaledId d a) = [SparseMatrixEntry i i a | i <- [0..(fromDim d - 1)]]

isZero :: AOperator -> Bool
isZero (Plus _ []) = True
isZero _ = False

add :: ManyBodyOperator -> ManyBodyOperator -> Maybe ManyBodyOperator
add op1 op2 = case (aAdd (aOp op1) (aOp op2)) of
  Just result -> Just $ ManyBodyOperator result NOperator
  Nothing -> Nothing
aAdd :: AOperator -> AOperator -> Maybe AOperator
aAdd op1 op2 | d1 == d2 = Just $ Plus d1 [op1, op2]
            | otherwise = Nothing
            where
              d1 = aGetDim op1
              d2 = aGetDim op2

scale :: Amplitude -> ManyBodyOperator -> ManyBodyOperator
scale a op = ManyBodyOperator (aScale a (aOp op)) NOperator

aScale :: Amplitude -> AOperator -> AOperator
aScale a (Kron d op1 op2) = Kron d (aScale a op1) op2
aScale a (Plus d ops) = Plus d (map (aScale a) ops)
aScale a (KetBra d (SparseMatrixEntry r c me)) = KetBra d $
  SparseMatrixEntry r c (a * me)
aScale a (ScaledId d me) = ScaledId d (a * me)

kron :: ManyBodyOperator -> ManyBodyOperator -> ManyBodyOperator
kron op1 op2 = ManyBodyOperator (aKron (aOp op1) (aOp op2)) NOperator
aKron :: AOperator -> AOperator -> AOperator
aKron (ScaledId d1 a1) (ScaledId d2 a2) = ScaledId (d1 * d2) (a1 * a2)
aKron op1 op2 | isZero op1 || isZero op2 = Plus (aGetDim op1 * aGetDim op2) []
              | otherwise = Kron (aGetDim op1 * aGetDim op2) op1 op2

apply :: ManyBodyOperator -> Ket -> Maybe Ket
apply op k = aApply (aOp op) k

aApply :: AOperator -> Ket -> Maybe Ket
aApply (ScaledId d a) x | fromDim d == VU.length x = Just $ VU.map (a *) x
                        | otherwise = Nothing
aApply (KetBra d (SparseMatrixEntry i j a)) x
  | fromDim d == VU.length x =
      Just $  VU.generate (fromDim d) yfunc
  | otherwise = Nothing
    where
      yfunc :: Int -> Amplitude
      yfunc i' | i' == i = a * (VU.!) x j
               | otherwise = 0
aApply (Plus d ops) x | fromDim d == VU.length x =
                        foldl' addVecs
                          (Just zeroVec)
                          (map (`aApply` x) ops)
                      | otherwise = Nothing
     where
      addVecs :: Maybe Ket -> Maybe Ket -> Maybe Ket
      addVecs _ Nothing = Nothing
      addVecs Nothing _ = Nothing
      addVecs (Just a) (Just b) = Just $ VU.zipWith (+) a b
      zeroVec = VU.replicate (fromDim d) 0

aApply (Kron _ op1 op2) x = do
    xp <- applyOne op2 d2 x
    xpp <- applyOne op1 d1 (transpose d1 xp)
    return (transpose d2 xpp)
  where
    d1 = fromDim $ aGetDim op1
    d2 = fromDim $ aGetDim op2
    subVectors :: Int -> Ket -> [Ket]
    subVectors d v = [VU.slice (d * j) d v | j <- [0..(numVecs - 1)]]
      where
        numVecs = VU.length v `div` d
    applyOne :: AOperator -> Int -> Ket -> Maybe Ket
    applyOne op d v = VU.concat `fmap` mapM (aApply op) (subVectors d v)


sigmaX :: ManyBodyOperator
Just (Just sigmaX) = liftM2 add (ketBra d 0 1 1.0) (ketBra d 1 0 1.0)
  where
    Just d = toDim 2

sigmaY :: ManyBodyOperator
Just (Just sigmaY) = liftM2 add (ketBra d 0 1 (0.0 :+ 1.0))
                                (ketBra d 1 0 (0.0 :+ (-1.0)))
  where
    Just d = toDim 2

sigmaZ :: ManyBodyOperator
Just (Just sigmaZ) = liftM2 add (ketBra d 1 1 1.0) (ketBra d 0 0 (-1.0))
  where
    Just d = toDim 2

sigmaPlus :: ManyBodyOperator
Just sigmaPlus = ketBra d 1 0 1.0
  where
    Just d = toDim 2

sigmaMinus :: ManyBodyOperator
Just sigmaMinus = ketBra d 0 1 1.0
  where
    Just d = toDim 2


data SimpleOperator = SimpleOperator Amplitude [Slice]
  deriving (Show)

data Slice = IdentityMatrix Dim
           | NonZeroLocation Dim Int Int
  deriving (Show)

simplify :: AOperator -> [SimpleOperator]
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
