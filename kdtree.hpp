/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <vector>
#include <stack>
#include <cmath>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
    if (first[curDim] == second[curDim])
      return first < second;
    return (first[curDim] <= second[curDim]) ? true : false;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
     double currBestDist = 0;
     double potDist = 0;

     // Adds the square of the distances in each dimension
     for (int i = 0; i < Dim; i++) {
       currBestDist += (target[i] - currentBest[i]) * (target[i] - currentBest[i]);
       potDist += (target[i] - potential[i]) * (target[i] - potential[i]);
     }
     if (currBestDist == potDist)
      return potential < currentBest;
     return (currBestDist <= potDist) ? false : true;
}

template <int Dim>
int KDTree<Dim>::partition(vector<Point<Dim>>& vect, int left, int right, int pivot, int dim) {
  // Initializes the variables
  Point<Dim> pivotValue = vect[pivot];
  int storeIndex = left;

  // Swaps the points around
  swap(vect, pivot, right);
  for (int i = left; i < right; i++) {
    if (smallerDimVal(vect[i], pivotValue, dim)) {
      swap (vect, i, storeIndex);
      storeIndex++;
    }
  }
  swap(vect, right, storeIndex);
  return storeIndex;
}

template <int Dim>
void KDTree<Dim>::swap(vector<Point<Dim>>& vect, int one, int two) {
  Point<Dim> p = vect[one];
  vect[one] = vect[two];
  vect[two] = p;
}

template <int Dim>
Point<Dim> KDTree<Dim>::select(vector<Point<Dim>>& vect, int left, int right, int k, int dim) {
  // Initializes variables
  int pivotIndex = left + floor(rand() % (right - left + 1));

  // Base case
  if (left == right)
    return vect[right];

  // Checks if the index is correct
  pivotIndex = partition(vect, left, right, pivotIndex, dim);
  if (k == pivotIndex)
    return vect[pivotIndex];
  else if (k < pivotIndex) {
    return select(vect, left, pivotIndex - 1, k, dim);
  }
  else {
    return select(vect, pivotIndex + 1, right, k, dim);
  }
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::constructorHelper(vector<Point<Dim>>& newPoints,
  int left, int right, int dim) {
  // Base case
  if (left > right)
    return nullptr;

  // Recursive calls and creates new node
  int midpoint = (left + right) / 2;
  KDTreeNode *node = new KDTreeNode(select(newPoints, left, right, midpoint, dim));
  node->left = constructorHelper(newPoints, left, midpoint - 1, (dim + 1) % Dim);
  node->right = constructorHelper(newPoints, midpoint + 1, right, (dim + 1) % Dim);

  return node;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
    // Initializes variables
    vector<Point<Dim>> v = newPoints;

    // Handles empty vector;
    if (newPoints.empty()) {
      root = nullptr;
      size = 0;
      return;
    }

    // std::cout<<std::endl;
    // for (unsigned i = 0; i < newPoints.size(); i++) {
    //   std::cout<<newPoints[i]<<" ";
    // }
    // std::cout<<std::endl;

    root = constructorHelper(v, 0, v.size() - 1, 0);
    size = v.size();
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  // Initializes variables
  std::vector<Point<Dim>> vect;
  std::stack<KDTreeNode*> s;
  KDTreeNode *node;

  // Checks if empty
  if (other.size == 0) {
    root = nullptr;
    size = 0;
    return;
  }

  // Creates vector of points based off of previous tree
  s.push(other.root);
  while (!s.empty()) {
    node = s.top();
    vect.push_back(node);
    s.pop();
    while (node->left != nullptr) {
      s.push(node->left);
      node = node->left;
    }
    if (node->right != nullptr) {
      s.push(node->right);
    }
  }

  // Calls constructor
  return KDTree(vect);

}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  delete this;
  KDTree<Dim>* tree = new KDTree(rhs);
  return tree;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
  // Initializes variables to delete tree
  std::stack<KDTreeNode*> s;
  std::vector<KDTreeNode*> v;
  KDTreeNode* node = root;

  // Checks if KDTree is not empty
  if (size == 0) {
    return;
  }

  // Preorder traversal to go to every node and delete it
  s.push(node);
  while (!s.empty()) {
    node = s.top();
    s.pop();
    if (node->left != nullptr)
      s.push(node->left);
    if (node->right != nullptr)
      s.push(node->right);
    delete node;
   }
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::findRecursive(const Point<Dim>& query,
    KDTreeNode* node, int d) const {
    // Variables
    KDTreeNode *left = nullptr;
    KDTreeNode *right = nullptr;
    KDTreeNode *best = nullptr;

    // Base cases
    if (node->left == nullptr && node->right == nullptr) {
      return node;
    }
    if (node->point == query)
      return node;

    // Checks to see which node is closest
    if (smallerDimVal(query, node->point, d)) {
      if (node->left != nullptr) {
        left = findRecursive(query, node->left, (d + 1) % Dim);
        best = left;
      }
      else {
        right = findRecursive(query, node->right, (d + 1) % Dim);
        best = right;
      }
    }
    else {
      if (node->right != nullptr) {
        right = findRecursive(query, node->right, (d + 1) % Dim);
        best = right;
      }
      else {
        left = findRecursive(query, node->left, (d + 1) % Dim);
        best = left;
      }
    }

    double radius = 0;
    for (int i = 0; i < Dim; i++) {
      radius += (query[i] - best->point[i]) * (query[i] - best->point[i]);
    }
    radius = sqrt(radius);
    // Checks if which node is the closest
    if (shouldReplace(query, best->point, node->point)) {
      if (best == left) {
        if (node->right != nullptr) {
          right = findRecursive(query, node->right, (d + 1) % Dim);
          if (shouldReplace(query, right->point, node->point))
            return node;
          else
            return right;
        }
        else
          return node;
      }
      if (best == right) {
        if (node->left != nullptr) {
          left = findRecursive(query, node->left, (d + 1) % Dim);
          if (shouldReplace(query, left->point, node->point))
            return node;
          else
            return left;
        }
      }
    }
    else if ((((best->point[d] + radius) >= node->point[d]) && (node->point[d] >= best->point[d])) ||
      (((best->point[d] - radius) <= node->point[d]) && (best->point[d] >= node->point[d]))) {
      if (best == left) {
        if (node->right != nullptr) {
          right = findRecursive(query, node->right, (d + 1) % Dim);
          if (shouldReplace(query, best->point, right->point))
            return right;
          else
            return best;
        }
        else
          return best;
      }
      else {
        if (node->left != nullptr) {
          left = findRecursive(query, node->left, (d + 1) % Dim);
          if (shouldReplace(query, best->point, left->point))
            return left;
          else
            return best;
        }
        else
          return best;
      }
    }
    else
      return best;

    return node;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    // Base case
    if (root == nullptr)
      return Point<Dim>();

    // Recursive helper function
    KDTreeNode *node = findRecursive(query, root, 0);

    return node->point;
}
