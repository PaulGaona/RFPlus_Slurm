#include <Rcpp.h>
using namespace Rcpp;


// Convert Rcpp::String to int
int rcpp_string_to_int(String rcpp_string) {
  std::string std_string = rcpp_string.get_cstring();
  return std::stoi(std_string);
}


// Convert Rcpp::CharacterVector to Rcpp::IntegerVector
IntegerVector char_to_int(CharacterVector char_vector) {
  int n = char_vector.size();
  IntegerVector int_vector(n);
  for (int i = 0; i < n; ++i) {
    std::string str = as<std::string>(char_vector(i));
    int_vector(i) = std::stoi(str);
  }
  return int_vector;
}


// Check if integer is in an Rcpp::IntegerVector
bool isin_int(int element, IntegerVector vector) {
  IntegerVector::iterator it = std::find(vector.begin(), vector.end(), element);
  return it != vector.end();
}


// Split Rcpp::String by comma
CharacterVector split_string_by_comma(String rcpp_string) {
  std::string std_string = rcpp_string.get_cstring();
  std::istringstream ss(std_string);
  std::string token;
  CharacterVector result;
  while (std::getline(ss, token, ',')) {
    result.push_back(token);
  }
  return result;
}


// [[Rcpp::export]]
NumericMatrix extract_psi_cpp(NumericMatrix x,
                              IntegerVector node_preds, 
                              List tree_paths,
                              IntegerVector node_ids, 
                              IntegerVector split_vars,
                              NumericVector split_vals) {
  int n = node_preds.size();
  int n_nodes = node_ids.size();
  
  NumericMatrix psi(n, n_nodes);
  for (int i = 0; i < n; ++i) {
    int terminal_node = node_preds[i];
    IntegerVector tree_path = as<IntegerVector>(
      tree_paths[std::to_string(terminal_node)]
    );
    IntegerVector in_node = match(tree_path, node_ids) - 1;
    for (IntegerVector::iterator j = in_node.begin(); j != in_node.end(); ++j) {
      if (x(i, split_vars(*j)) > split_vals(*j)) {
        psi(i, *j) = 1;
      } else {
        psi(i, *j) = -1;
      }
    }
  }
  
  return psi;
}


// [[Rcpp::export]]
NumericMatrix extract_psi_chr_cpp(NumericMatrix x,
                                  IntegerVector node_preds, 
                                  List tree_paths,
                                  IntegerVector node_ids, 
                                  IntegerVector split_vars,
                                  CharacterVector split_vals,
                                  IntegerVector unordered_factors) {
  int n = node_preds.size();
  int n_nodes = node_ids.size();
  
  NumericMatrix psi(n, n_nodes);
  for (int i = 0; i < n; ++i) {
    int terminal_node = node_preds[i];
    IntegerVector tree_path = as<IntegerVector>(
      tree_paths[std::to_string(terminal_node)]
    );
    IntegerVector in_node = match(tree_path, node_ids) - 1;
    for (IntegerVector::iterator j = in_node.begin(); j != in_node.end(); ++j) {
      String split_val_str = split_vals(*j);
      if (isin_int(split_vars(*j), unordered_factors)) {
        CharacterVector split_val_chr = split_string_by_comma(split_val_str);
        IntegerVector split_val = char_to_int(split_val_chr);
        if (isin_int(x(i, split_vars(*j)), split_val)) {
          psi(i, *j) = 1;
        } else {
          psi(i, *j) = -1;
        }
      } else {
        if (x(i, split_vars(*j)) > rcpp_string_to_int(split_val_str)) {
          psi(i, *j) = 1;
        } else {
          psi(i, *j) = -1;
        }
      }
    }
  }
  
  return psi;
}


