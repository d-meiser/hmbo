module HMbo.LinearOpSpec (spec) where

import Test.Hspec
import HMbo
import Data.Maybe

spec :: Spec
spec = do
  describe "Dim" $ do
    it "Can't be built from negative integeter." $ do
      toDim (-1) `shouldBe` Nothing
    it "Can't be built from 0." $ do
      toDim 0 `shouldBe` Nothing
    it "Can be built from positive integer." $ do
      toDim 2 `shouldSatisfy` (Nothing /=)
    it "Has the right dimension." $ do
      toDim 2 `shouldSatisfy` ((Just 2 ==) . (fmap fromDim))

  describe "Null Matrix" $ do
    it "Can be built from a dimension." $ do
      getDim (zero (fromJust $ toDim 2)) `shouldBe` (fromJust $ toDim 2)
    it "Does not change dimension when scaled." $ do
      getDim (scale 3.7 $ zero (fromJust $ toDim 2)) 
        `shouldBe` (fromJust $ toDim 2)

  describe "Identity Matrix" $ do
    it "Can be built from a dimension." $ do
      getDim (identity (fromJust $ toDim 2)) `shouldBe` (fromJust $ toDim 2)
    it "Does not change dimension when scaled." $ do
      getDim (scale 3.7 $ identity (fromJust $ toDim 2))
        `shouldBe` (fromJust $ toDim 2)

  describe "kron" $ do
    it "Multiplies dimensions." $ do
      getDim (
          kron (identity (fromJust $ toDim 3))
          (identity (fromJust $ toDim 3))) `shouldBe`
            (fromJust $ toDim 9)
    it "Results in Null when one of the operators is Null." $ do
      kron
        (identity (fromJust $ toDim 3))
        (zero (fromJust $ toDim 3))
          `shouldBe`
            zero (fromJust $ toDim 9)
    it "Results in ScaledId when both of the operators are identites." $ do
      kron
        (identity (fromJust $ toDim 3))
        (identity (fromJust $ toDim 3))
          `shouldBe`
            identity (fromJust $ toDim 9)
