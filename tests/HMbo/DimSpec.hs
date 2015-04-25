module HMbo.DimSpec (spec) where

import Test.Hspec
import HMbo

spec :: Spec
spec =
  describe "Dim" $ do
    it "Can't be built from negative integeter." $
      toDim (-1) `shouldBe` Nothing
    it "Can't be built from 0." $
      toDim 0 `shouldBe` Nothing
    it "Can be built from positive integer." $
      toDim 2 `shouldSatisfy` (Nothing /=)
    it "Has the right dimension." $
      toDim 2 `shouldSatisfy` ((Just 2 ==) . fmap fromDim)
