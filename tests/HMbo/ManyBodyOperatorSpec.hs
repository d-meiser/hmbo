module HMbo.ManyBodyOperatorSpec (spec) where

import Test.Hspec
import HMbo
import Data.Maybe
import qualified Data.Vector.Unboxed as VU
import Data.Complex (conjugate, realPart)

matrixElement :: Ket -> ManyBodyOperator -> Ket -> Amplitude
matrixElement psi a phi = VU.foldl1 (+) $ VU.zipWith (*) psi' aPhi
  where
    aPhi = fromJust $ a `apply` phi
    psi' = VU.map conjugate psi

basisState :: Int -> Int -> Ket
basisState d i = VU.fromList [kroneckerDelta i j | j <- [0..(d - 1)]]
  where
    kroneckerDelta m n | m == n = 1.0
                       | otherwise = 0.0

aij :: ManyBodyOperator -> Int -> Int -> Amplitude
aij a i j = matrixElement (basisState d i) a (basisState d j)
  where
    d = fromDim $ getDim a

isClose :: Double -> Amplitude -> Amplitude -> Bool
isClose tol a b = (realPart $ abs (a - b)) < tol

defTol :: Double
defTol = 1.0e-12

isHermitian :: Double -> ManyBodyOperator -> Bool
isHermitian tol a =
  and [isClose tol (aij a i j) (conjugate (aij a j i))
      | i <- [0..(d-1)], j <- [0..i]]
  where
    d = fromDim $ getDim a

isDiagonal :: Double -> ManyBodyOperator -> Bool
isDiagonal tol a =
  and [i == j || isClose tol 0.0 (aij a i j)
      | i <- [0..(d - 1)], j <- [0..(d - 1)]]
  where
    d = fromDim $ getDim a

spec :: Spec
spec = do
  describe "Null Matrix" $ do
    it "Can be built from a dimension." $
      getDim (zero (fromJust $ toDim 2)) `shouldBe` fromJust (toDim 2)
    it "Does not change dimension when scaled." $
      getDim (scale 3.7 $ zero (fromJust (toDim 2))) 
        `shouldBe` fromJust (toDim 2)

  describe "Identity Matrix" $ do
    it "Can be built from a dimension." $
      getDim (eye (fromJust (toDim 2))) `shouldBe` fromJust (toDim 2)
    it "Does not change dimension when scaled." $
      getDim (scale 3.7 $ eye (fromJust (toDim 2)))
        `shouldBe` fromJust (toDim 2)

  describe "ketBra" $ do
    it "Can be constructed from valid indices." $ do
      ketBra (fromJust $ toDim 2) 0 0 1.0 `shouldSatisfy` isJust
      ketBra (fromJust $ toDim 2) 0 1 1.0 `shouldSatisfy` isJust
    it "Cannot be constructed from negative indices." $ do
      ketBra (fromJust $ toDim 2) (-1) 0 1.0 `shouldBe` Nothing
      ketBra (fromJust $ toDim 2) 0 (-1) 1.0 `shouldBe` Nothing
    it "Cannot be constructed from indices that are too large." $
      ketBra (fromJust (toDim 2)) 2 0 1.0 `shouldBe` Nothing

  describe "kron" $ do
    it "Multiplies dimensions." $
      getDim (
          kron (eye (fromJust (toDim 3)))
          (eye (fromJust $ toDim 3))) `shouldBe` fromJust (toDim 9)
    it "Results in ScaledId when both of the operators are identites." $
      kron
        (eye (fromJust (toDim 3)))
        (eye (fromJust (toDim 3)))
          `shouldBe` eye (fromJust (toDim 9))
    it
      "Produces Hermitian operators when applied to Hermitian operators." $
      do
        kron sigmaX sigmaY `shouldSatisfy` isHermitian defTol
        kron sigmaX (eye (fromJust $ toDim 4)) `shouldSatisfy`
          isHermitian defTol
    it
      "Produces a diagonal operator when supplied with diagonal operators." $
      do
        kron sigmaZ sigmaZ `shouldSatisfy` isDiagonal defTol
        kron (eye 3) sigmaZ `shouldSatisfy` isDiagonal defTol
        kron sigmaZ (eye 3) `shouldSatisfy` isDiagonal defTol

  describe "add" $ do
    it "Results in Nothing if dimensions don't match." $
      add
        (eye (fromJust (toDim 2)))
        (eye (fromJust (toDim 3)))
          `shouldBe` Nothing
    it "Adds matrix entries." $ do
      let op1 = sigmaZ `kron` (eye 2)
      let op2 = sigmaX `kron` sigmaY
      let theSum = fromJust $ op1 `add` op2
      aij theSum 0 1 `shouldSatisfy`
        isClose defTol ((aij op1 0 1) + (aij op2 0 1))
      aij theSum 2 1 `shouldSatisfy`
        isClose defTol ((aij op1 2 1) + (aij op2 2 1))
      aij theSum 0 0 `shouldSatisfy`
        isClose defTol ((aij op1 0 0) + (aij op2 0 0))

  describe "apply" $ do
    it "Returns nothing when dimensions don't match." $ do
      let v = VU.fromList [1.3, 3.4] :: Ket
      let d = fromJust $ toDim 3
      apply (eye d) v `shouldBe` Nothing
    it "Returns something when dimensions match." $ do
      let v = VU.fromList [1.3, 3.4] :: Ket
      let d = fromJust $ toDim 2
      apply (eye d) v `shouldSatisfy` isJust

    describe "Identity" $
      it "Returns vectors unchanged." $ do
        let v = VU.fromList [1.3, 3.4] :: Ket
        let d = fromJust $ toDim 2
        apply (eye d) v `shouldBe` Just v

    describe "Zero operator" $
      it "Returns 0." $ do
        let v = VU.fromList [1.3, 3.4] :: Ket
        let w = VU.replicate 2 0 :: Ket
        let d = fromJust $ toDim 2
        fromJust (apply (zero d) v) `shouldBe` w

  describe "transpose" $ do
    it "Returns the same vector when transposed with pivot 1." $ do
      let v = VU.fromList [1.3, 3.4] :: Ket
      transpose 1 v `shouldBe` v
    it "Returns the same vector when transposed with pivot (length v)." $ do
      let v = VU.fromList [1.3, 3.4] :: Ket
      transpose (VU.length v) v `shouldBe` v
    it "Works for a typical case" $ do
      let v = VU.fromList [1, 2, 3, 4] :: Ket
      let w = VU.fromList [1, 3, 2, 4] :: Ket
      transpose 2 v `shouldBe` w

  describe "sigmaZ" $ do
    it "Has no off-diagonal elements." $
      sigmaZ `shouldSatisfy` isDiagonal defTol
    it "Considers the 0 state as spin down." $ do
      aij sigmaZ 0 0 `shouldSatisfy` isClose defTol (-1.0)
    it "Considers the 1 state as spin up." $ do
      aij sigmaZ 1 1 `shouldSatisfy` isClose defTol (1.0)

  describe "sigmaX" $ do
    it "Has no diagonal elements." $ do
      aij sigmaX 0 0 `shouldSatisfy` isClose defTol 0.0
      aij sigmaX 1 1 `shouldSatisfy` isClose defTol 0.0
    it "Is Hermitian." $
      sigmaX `shouldSatisfy` isHermitian defTol

  describe "sigmaY" $ do
    it "Has no diagonal elements." $ do
      aij sigmaY 0 0 `shouldSatisfy` isClose defTol 0.0
      aij sigmaY 1 1 `shouldSatisfy` isClose defTol 0.0
    it "Is Hermitian." $
      sigmaY `shouldSatisfy` isHermitian defTol

--  describe "simplify" $ do
--    it "Produces a single entry when applied to an identity matrix." $
--      length (simplify (eye 2)) `shouldBe` 1
--    it "Produces a single entry when applied to a KetBra." $
--      length (simplify (fromJust $ ketBra 2 0 0 1.0)) `shouldBe` 1
