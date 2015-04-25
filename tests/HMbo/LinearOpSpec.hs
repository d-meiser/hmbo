module HMbo.LinearOpSpec (spec) where

import Test.Hspec
import HMbo
import Data.Maybe
import qualified Data.Vector.Unboxed as VU


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

  describe "add" $
    it "Results in Nothing if dimensions don't match." $
      add
        (eye (fromJust (toDim 2)))
        (eye (fromJust (toDim 3)))
          `shouldBe`
            Nothing

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
    it "Returns the same vector when transposed with pivot 1" $ do
      let v = VU.fromList [1.3, 3.4] :: Ket
      transpose 1 v `shouldBe` v
    it "Returns the same vector when transposed with pivot (length v)" $ do
      let v = VU.fromList [1.3, 3.4] :: Ket
      transpose (VU.length v) v `shouldBe` v
    it "Works for a typical case" $ do
      let v = VU.fromList [1, 2, 3, 4] :: Ket
      let w = VU.fromList [1, 3, 2, 4] :: Ket
      transpose 2 v `shouldBe` w
      
  describe "simplify" $ do
    it "Produces a single entry when applied to an identity matrix" $
      length (simplify (eye 2)) `shouldBe` 1
    it "Produces a single entry when applied to a KetBra" $
      length (simplify (fromJust $ ketBra 2 0 0 1.0)) `shouldBe` 1
      
