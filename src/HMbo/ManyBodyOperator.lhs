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
    numberOperator,
    annihilationOperator,
    creationOperator,
    ketBra,
    getDim,
    toDim,
    fromDim,
    mul,
    add,
    scale,
    kron,
    apply,
    isZero
    ) where

import qualified Data.Vector.Unboxed as VU
import Data.Complex
import Data.List (foldl')
import Control.Monad (liftM2,foldM)
import Data.Maybe (fromJust)

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
The type {\tt AOperator} is tailored towards algebraic operations and
{\tt NOperator} is tailored towards numerical application of the
operator to state vectors.  {\tt NOperator}s are created from {\tt
AOperator}s by means of a compilation process.  In the compilation
analogy the {\tt AOperator} corresponds to a high level language while
the {\tt NOperator} is a form of assembly language.  Thus, we can
construct a full {\tt ManyBodyOperator} from an {\tt AOperator} as
follows:
\begin{code}
fromAOperator :: AOperator -> ManyBodyOperator
fromAOperator op = ManyBodyOperator op (compile op)
\end{code}
The {\tt compile} function will be discussed in
section~\ref{ssec:NOperator}.  The compilation can be performed with
computational cost comparable to other algebraic operations, e.g.
$\mathcal{O}(N)$.


\subsection{Construction of many body operators}
\label{ssec:ManyBodyOperatorConstructors}

Any {\tt ManyBodyOperator} can be constructed from the following
primitive operators along with the algebraic operations discussed in
section \ref{ssec:ManyBodyOperatorAlgebra}.  The details of the {\tt
AOperator} constructions is discussed in section \ref{ssec:AOperator}.

Perhaps the simplest operator is the null operator,
\begin{code}
-- | Construct a null operator
zero :: Dim -> ManyBodyOperator
zero d = fromAOperator (Plus d [])
\end{code}

Similarly, the function {\tt eye} constructs an identity operator:
\begin{code}
-- | Construct an identity operator
eye :: Dim -> ManyBodyOperator
eye d = fromAOperator (ScaledId d 1.0)
\end{code}

The function {\tt ketBra} for constructing an operator with a single
non-zero entry is slightly more involved:
\begin{code}
-- | Construct an operator with a single non-zero matrix element.
--   This function returns Nothing when i or j are out of bounds.
ketBra :: Dim       -- ^ Dimension of the operator
       -> Int       -- ^ Row index of non-zero entry
       -> Int       -- ^ Column index of non-zero entry
       -> Amplitude -- ^ Non-zero matrix element
       -> Maybe ManyBodyOperator
ketBra d i j a | i >= 0 && i < d' && j >= 0 && j < d' = Just $
                 fromAOperator (KetBra d (SparseMatrixEntry i j a))
               | otherwise = Nothing
  where
    d' = fromDim d
\end{code}
The main reason for the more complicated code is the fact that this
function can fail.  To ensure that ill-formed {\tt ManyBodyOperator}
instances can never exist we return {\tt Nothing} in that case.


\subsection{Algebra of many body operator}
\label{ssec:ManyBodyOperatorAlgebra}

The basic algebraic operations on many body operators are multiplication
by a scalar ({\tt scale}), operator addition ({\tt add}), operator
multiplicaiton ({\tt mul}), and Kronecker products of operators ({\tt
kron}).  As discussed in the introduction to this section, the {\tt
AOperator} is specifically designed for these operations.  Therefore it
is no surprise that the actual work involved in these algebraic
operations is carried out by functions on the {\tt AOperator} field.  We
prefix these operators with ``a'', e.g. {\tt aMul} for {\tt mul}.  Once
the {\tt AOperator} of the result is computed we find the {\tt
NOperator} simply by means of compilation\footnote{It may be worth
pointing out that this compilation is performed lazily.  As a
consequence, the {\tt NOperator} is not actually constructed for
intermediate results that are never applied to a {\tt Ket}.}.

\begin{code}
-- | Multiply an operator by a scalar.
scale :: Amplitude -> ManyBodyOperator -> ManyBodyOperator
scale a op = fromAOperator (aScale a (aOp op))

-- | Add two operators.
--   This function returns Nothing when operators with different
--   dimensions are supplied.
add :: ManyBodyOperator -> ManyBodyOperator -> Maybe ManyBodyOperator
add op1 op2 = case (aAdd (aOp op1) (aOp op2)) of
  Just result -> Just $ fromAOperator result
  Nothing -> Nothing

-- | Product of two operators (matrix multiplication).
--   This function returns Nothing when operators with different
--   dimensions are supplied.
mul :: ManyBodyOperator -> ManyBodyOperator -> Maybe ManyBodyOperator
mul op1 op2 | d1 == d2 = Just $
              fromAOperator (aMul (aOp op1) (aOp op2))
            | otherwise = Nothing
            where
              d1 = getDim op1
              d2 = getDim op2

-- | Compute the Kronecker product of two operators.
kron :: ManyBodyOperator -> ManyBodyOperator -> ManyBodyOperator
kron op1 op2 = fromAOperator (aKron (aOp op1) (aOp op2))

-- | Compute the dimension of an operator.
getDim :: ManyBodyOperator -> Dim
getDim (ManyBodyOperator aop _) = aGetDim aop
\end{code}

The functions {\tt scale} and {\tt kron} trivially forward to {\tt
aScale} and {\tt aKron}, respectively.  The functions {\tt add} and {\tt
mul} on the other hand are a little bit trickier because they can fail.
The sum or product of two operators of different dimension is
meaningless.  In the event of incompatible operators, {\tt add} and {\tt
mul} return {\tt Nothing}.  In separate efforts we have attempted to
encode this failure mode in the type system.  But we found that the
result was unduly cumbersome in the common case where the operators are
in fact compatible.  In our current design, there are three {\tt
ManyBodyOperator} related functions that may fail: {\tt ketBra}, {\tt
add}, and {\tt mul}.

Finally, to apply {\tt ManyBodyOperators} to {\tt Ket}s we use the {\tt
NOperator} part. \textbf{Properly implement this}:
\begin{code}
-- | Apply an operator to a vector
--   This function fails if the dimension of the operator does not match
--   the dimension of the vector.
apply :: ManyBodyOperator -> Ket -> Maybe Ket
apply op k = nApply (nOp op) k
\end{code}


\subsection{Algebraic operators}
\label{ssec:AOperator}

The {\tt AOperator} is designed to make it easy to implement the
typical algebraic operations such as addition or outer products:
\begin{code}
data AOperator = Kron Dim AOperator AOperator
               | Plus Dim [AOperator]
               | KetBra Dim SparseMatrixEntry
               | ScaledId Dim Amplitude
  deriving (Show, Eq)
\end{code}
The data constructor {\tt Kron} represents the outer product of two
operators. {\tt Plus} is used to represent sums of operators.  In
contrast to Kron we use a list to represent the summands.  This is to
allow the representation of the null operator.  The {\tt KetBra}
constructor is used to represent elementary sparse matrix entries.  Its
name derives from the notation customary in quantum mechanics,
$|i\rangle\langle j|$, where $|i\rangle$ and $|j\rangle$ are basis
vectors.  Such operators are represented naively by the following data
type:
\begin{code}
data SparseMatrixEntry = SparseMatrixEntry Int Int Amplitude
  deriving(Show, Eq)
\end{code}
Finally, {\tt ScaledId} is the identity operator scaled by a
complex number, $\alpha \mathds{1}$.  The scale factor is needed so that
any {\tt AOperator} can be scaled by a complex number.  All forms of
AOperator keep track of the dimension of the Hilbert space on which they
act.  The Hilbert space is used to decide if two operators are
compatible with one another for purposes of addition and multiplication.


\subsubsection{Implementation of many body operator algebra}

Next we turn to a discussion of how the operator algebra is implemented.
We start with the multiplication of an operator by a scalar:
\begin{code}
aScale :: Amplitude -> AOperator -> AOperator
aScale a (Kron d op1 op2) = Kron d (aScale a op1) op2
aScale a (Plus d ops) = Plus d (map (aScale a) ops)
aScale a (KetBra d (SparseMatrixEntry r c me)) = KetBra d $
  SparseMatrixEntry r c (a * me)
aScale a (ScaledId d me) = ScaledId d (a * me)
\end{code}
For the elementary operators {\tt ScaledId}s and {\tt KetBra}s we simply
multiply the operator's scale factor by {\tt a}.  For a sum of operators
we scale each summand and for an outer product we scale the first
factor.

Addition of operators is rather simpler because one of the data
constructors is exactly the sum of operators.  Nearly all the code in
{\tt aAdd} is error checking to make sure that we never end up with the
sum of two incompatible operators:
\begin{code}
aAdd :: AOperator -> AOperator -> Maybe AOperator
aAdd op1 op2 | d1 == d2 = Just $ Plus d1 [op1, op2]
            | otherwise = Nothing
            where
              d1 = aGetDim op1
              d2 = aGetDim op2

aGetDim :: AOperator -> Dim
aGetDim (Kron d _ _) = d
aGetDim (Plus d _) = d
aGetDim (KetBra d _) = d
aGetDim (ScaledId d _) = d
\end{code}

Multiplicatino of operators is quite complicated and we have to
distinguish many different cases.  Multiplication by an identity amounts
to scaling an operator:
\begin{code}
aMul :: AOperator -> AOperator -> AOperator
aMul (ScaledId _ a) op = aScale a op
aMul op (ScaledId _ a) = aScale a op
\end{code}
When one of the factors is a sum we use the distributive law of operator
multiplication:
\begin{code}
aMul (Plus d ops) op = Plus d (map ($ op) (fmap aMul ops))
aMul op (Plus d ops) = Plus d (map (aMul op) ops)
\end{code}
When both factors are a {\tt KetBra} we have to check if the column
index of the first matches the row index of the second.  If they match
we have a non-zero entry.  If they don't match the product is zero which
we represent by {\tt Plus d []}:
\begin{code}
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
\end{code}
When multiplying two outer products we have to distinguish two main
sub-cases.  If the dimensions of the Kronecker factors in the two matrix
multiplication factors match we can simply recursively multiply the
Kronecker factors.  If the dimensions do not match we balk. \textbf{Fix
this}.
\begin{code}
aMul (Kron dim a b) (Kron _ c d)
        | dimsCompatible = Kron  dim (aMul a c) (aMul b d)
        | otherwise =
          error "Can't handle products of Krons when dimensions don't match"
  where
    dimsCompatible = aGetDim a == aGetDim c && aGetDim b == aGetDim d
\end{code}
When a {\tt KetBra} is multiplied with a Kronecker product we proceed by
turning the Kronecker Product into a list of sparse entries.  We can
then use the distributive law to compute the operator product by
multiplying the {\tt KetBra} with each sparse entry.
\begin{code}
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
\end{code}
Turning the different types of algebraic operators into lists of sparse
entries is straight forward but tedious:
\begin{code}
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
\end{code}
Note that many of the elementary products are zero.  We filter those
with the following predicate:
\begin{code}
isZero :: AOperator -> Bool
isZero (Plus _ []) = True
isZero _ = False
\end{code}

Thanks to the design of the {\tt AOperator} it is trivial to take the
outer product of operators.  As an optimization, we merge adjacent
identity operators into a larger identity operator:
\begin{code}
aKron :: AOperator -> AOperator -> AOperator
aKron (ScaledId d1 a1) (ScaledId d2 a2) = ScaledId (d1 * d2) (a1 * a2)
aKron op1 op2 | isZero op1 || isZero op2 = Plus (aGetDim op1 * aGetDim op2) []
              | otherwise = Kron (aGetDim op1 * aGetDim op2) op1 op2
\end{code}


\subsection{Numerical operators}
\label{ssec:NOperator}

In order to efficiently apply operators to states we need a simple
representation of operators, an ``assembly language'' for many body
operators.  We choose the following representation:
\begin{equation}
\hat A = \sum _j \alpha_j a_j \otimes b_j \otimes \ldots\;.
\label{eqn:NOperator}
\end{equation}
An operator is a sum of simple operators, each of which is the tensor
product of elementary operators $a_j$, $b_j$, etc.  The elementary
operators are of one of two forms.  They are either an identity operator
or they are a matrix with a single non-zero entry.  The matrix element
of the non-zero entry is one.

These ideas are straight forward to translate into Haskell.  An {\tt
NOperator} is a list of simple operators corresponding to the sum in
Eq.~(\ref{eqn:NOperator}:
\begin{code}
type NOperator = [SimpleOperator]
\end{code}

A {\tt SimpleOperator} has an amplitude, $\alpha_j$ in
Eq.~(\ref{eqn:NOperator}), and a list of {\tt Slice}s representing the
tensor product:
\begin{code}
data SimpleOperator = SimpleOperator Amplitude [Slice]
  deriving (Show, Eq)
\end{code}

A {\tt Slice} can be either an identity matrix or a matrix with a single
non-zero matrix entry:
\begin{code}
data Slice = IdentityMatrix Dim
           | NonZeroLocation Dim Int Int
  deriving (Show, Eq)
\end{code}

The compilation process is rather straight forward. We flatten Kronecker
products of sums into sums of Kronecker products using the
distributivity of the Kronecker product and collect scale factors in
each term.  As an optimization, we combine adjacent identity matrices
into larger identity matrices and adjacent {\tt NonZeroLocation}s into a
single {\tt NonZeroLocation} using the {\tt mergePairs} function.
\begin{code}
compile :: AOperator -> NOperator
compile (ScaledId d a) = [SimpleOperator a [IdentityMatrix d]]
compile (KetBra d (SparseMatrixEntry i j a)) =
  [SimpleOperator a [NonZeroLocation d i j]]
compile (Plus _ ops) = concatMap compile ops
compile (Kron _ op1 op2) = map mergePairs
    [simpleKron sop1 sop2 | sop1 <- compile op1 , sop2 <- compile op2]
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
\end{code}

Now that we know how to compile {\tt AOperator}s into {\tt NOperator}s
all that's left to do is to implement application of {\tt NOperator}s to
{\tt Ket}s.  Operator application is achieved by applying each simple
operator to the {\tt Ket} and by adding the results together:
\begin{code}
nApply :: NOperator -> Ket -> Maybe Ket
nApply ops k
  | dimsOk = Just $ (foldl' vAdd vZero) (map ((flip nApply') k) ops)
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
\end{code}

The {\tt sDim} function simply extracts the dimension of a {\tt Slice}:
\begin{code}
sDim :: Slice -> Int
sDim (IdentityMatrix d') = fromDim d'
sDim (NonZeroLocation d' _ _) = fromDim d'
\end{code}

We are left with computing the application of the individual {\tt
SimpleOperator}s to a {\tt Ket} which we do by recursively traversing
the list of slices:
\begin{code}
nApply' :: SimpleOperator -> Ket -> Ket
nApply' (SimpleOperator a []) k = VU.map (a *) k
nApply' (SimpleOperator a [s]) k =
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
nApply' (SimpleOperator a (s:ss)) k =
  case s of
    IdentityMatrix _ ->
      VU.concat $
        map (nApply' (SimpleOperator a ss)) $
        [VU.slice (i * blockSize) blockSize k | i <- [0..(d - 1)]]
    NonZeroLocation _ i j ->
        VU.concat $
        [VU.replicate (i * blockSize) 0
        ,nApply'
           (SimpleOperator a ss)
           (VU.slice (j * blockSize) blockSize k)
        ,VU.replicate (totalDim - (i + 1) * blockSize) 0
        ]
  where
    blockSize = product (map sDim ss)
    d = sDim s
    totalDim = d * blockSize
\end{code}
Note that we have unrolled the case with one slice left in the
recursion.  This is an optimization.


\subsection{A few elementary operators}
\label{ssec:ElementaryOperators}

The following operators serve both as an illustration of how to
construct {\tt ManyBodyOperators} as well as building blocks for more
complicate operators.


\subsubsection{Spin operators}

First, we consider some spin operators.  The spin raising and lowering
operators can be implemented as follows:
\begin{code}
-- | Spin raising operator for spin 1/2
sigmaPlus :: ManyBodyOperator
Just sigmaPlus = ketBra d 1 0 1.0
  where
    Just d = toDim 2

-- | Spin lowering operator for spin 1/2
sigmaMinus :: ManyBodyOperator
Just sigmaMinus = ketBra d 0 1 1.0
  where
    Just d = toDim 2
\end{code}

Next, the Pauli spin matrices:
\begin{code}
-- | Pauli sigma^x matrix
sigmaX :: ManyBodyOperator
Just (Just sigmaX) = liftM2 add (ketBra d 0 1 1.0) (ketBra d 1 0 1.0)
  where
    Just d = toDim 2

-- | Pauli sigma^y matrix
sigmaY :: ManyBodyOperator
Just (Just sigmaY) = liftM2 add (ketBra d 0 1 (0.0 :+ 1.0))
                                (ketBra d 1 0 (0.0 :+ (-1.0)))
  where
    Just d = toDim 2

-- | Pauli sigma^y matrix
sigmaZ :: ManyBodyOperator
Just (Just sigmaZ) = liftM2 add (ketBra d 1 1 1.0) (ketBra d 0 0 (-1.0))
  where
    Just d = toDim 2
\end{code}

We can also build angular momentum operators with magnitude greater than
$J=1/2$. \textbf{Todo: Develop a type for $J$ that can only be positive
and write $J_+$, $J_-$, $J_x$, $J_y$, and $J_z$}.


\subsubsection{Harmonic oscillator operators}

The harmonic oscillator operators in a Hilbert space truncated to
a maximum quantum number of $n_{\rm max}$ is:
\begin{code}
-- | Number operator for a harmonic oscillator
numberOperator :: Dim -> ManyBodyOperator
numberOperator d = fromJust $ foldM add (zero d)
  [fromJust $ ketBra d i i (fromIntegral i) | i <- [0..nMax]]
  where
    nMax = fromDim d

-- | Annihilation operator
annihilationOperator :: Dim -> ManyBodyOperator
annihilationOperator d = fromJust $ foldM add (zero d)
  [fromJust $ ketBra d (i - 1) i (sqrt (fromIntegral i)) | i <- [1..nMax]]
  where
    nMax = fromDim d

-- | Creation operator
creationOperator :: Dim -> ManyBodyOperator
creationOperator d = fromJust $ foldM add (zero d)
  [fromJust $ ketBra d (i + 1) i (sqrt (fromIntegral i + 1.0))
  | i <- [0..(nMax - 1)]]
  where
    nMax = fromDim d
\end{code}



\end{code}
