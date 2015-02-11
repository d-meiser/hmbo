module HMbo.SparseMatrixSpec (spec) where

import Test.Hspec
import HMbo.SparseMatrix

spec :: Spec
spec = do
  describe "SparseMatrix" $ do
    it "Zero Matrix has no entries" $ do
      numEntries zeroMatrix `shouldBe` 0
