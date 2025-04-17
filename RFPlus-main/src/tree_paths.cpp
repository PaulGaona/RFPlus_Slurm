#include <Rcpp.h>
using namespace Rcpp;


int match_int(int value, IntegerVector vector) {
  IntegerVector::iterator it = std::find(vector.begin(), vector.end(), value);
  if (it != vector.end()) {
    return std::distance(vector.begin(), it);
  } else {
    return -1;
  }
}


// [[Rcpp::export]]
List get_tree_paths_cpp(IntegerVector terminal_node_ids, 
                        IntegerVector left_child_ids,
                        IntegerVector right_child_ids, 
                        IntegerVector node_ids) {
  
  List tree_paths(terminal_node_ids.size());
  for (int i = 0; i < terminal_node_ids.size(); ++i) {
    int terminal_node_id = terminal_node_ids[i];
    int node_id = terminal_node_id;
    IntegerVector tree_path;
    while (node_id != 0) {
      int idx = match_int(node_id, left_child_ids);
      if (idx == -1) {
        idx = match_int(node_id, right_child_ids);
      }
      node_id = node_ids[idx];
      tree_path.push_back(node_id);
    }
    tree_paths[i] = tree_path;
  }
  
  return tree_paths;
}


