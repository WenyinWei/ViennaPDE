#pragma once
#include <vector>
#include "gtest/gtest.h"

template<typename NumericT>
void EXPECT_MESH_EQ(
    std::vector<std::vector<std::vector<NumericT>>> & A, 
    std::vector<std::vector<std::vector<NumericT>>> & B, std::string err_info) {
  for (size_t layer_i=0; layer_i< A.size(); layer_i++) {
    for (size_t row_i=0; row_i < A[0].size(); row_i++) {
      for (size_t column_i= 0; column_i < A[0][0].size(); column_i++) {
        if constexpr (std::is_same<NumericT, float>::value) {
          EXPECT_FLOAT_EQ( A[layer_i][row_i][column_i], B[layer_i][row_i][column_i]) << err_info;
        } else if constexpr (std::is_same<NumericT, double>::value) {
          EXPECT_DOUBLE_EQ(A[layer_i][row_i][column_i], B[layer_i][row_i][column_i]) << err_info;
        } else {FAIL() << "Specify your wanted scalar value type to test, float or double?";}
      }
    }
  }
}


template<typename NumericT>
void EXPECT_MESH_NEGATIVE_EQ(
    std::vector<std::vector<std::vector<NumericT>>> & A, 
    std::vector<std::vector<std::vector<NumericT>>> & B, std::string err_info) {
  for (size_t layer_i=0; layer_i< A.size(); layer_i++) {
    for (size_t row_i=0; row_i < A[0].size(); row_i++) {
      for (size_t column_i= 0; column_i < A[0][0].size(); column_i++) {
        if constexpr (std::is_same<NumericT, float>::value) {
          EXPECT_FLOAT_EQ( A[layer_i][row_i][column_i], -B[layer_i][row_i][column_i]) << err_info;
        } else if constexpr (std::is_same<NumericT, double>::value) {
          EXPECT_DOUBLE_EQ(A[layer_i][row_i][column_i], -B[layer_i][row_i][column_i]) << err_info;
        } else {FAIL() << "Specify your wanted scalar value type to test, float or double?";}
      }
    }
  }
}
