module Hmbo.SparseMatrixSpec (spec) where

import Test.Hspec

spec :: Spec
spec = do
  describe "SparseMatrix" $ do
    it "Zero Matrix has no entries" $ do
      numEntries SparseMatrix `shouldBe` 0
