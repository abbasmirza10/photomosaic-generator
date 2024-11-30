/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <deque>

using namespace std;

template<int Dim>
bool smallerDimVal(const Point<Dim>& first, const Point<Dim>& second, int curDim) {
    if (first[curDim] != second[curDim]) {
        return first[curDim] < second[curDim];
    }
    return first < second;
}


template <int Dim>
bool shouldReplace(const Point<Dim>& target, const Point<Dim>& currentBest, const Point<Dim>& potential) {
    double currDist = 0.0;
    double potentialDist = 0.0;
    for (int i = 0; i < Dim; ++i) {
        currDist += (currentBest[i] - target[i]) * (currentBest[i] - target[i]);
        potentialDist += (potential[i] - target[i]) * (potential[i] - target[i]);
    }
    if (potentialDist < currDist) {
        return true;  
    }
    else if (potentialDist > currDist) {
        return false;  
    }
    else
        for (int i = 0; i < Dim; ++i) {
            if (potential[i] != currentBest[i]) {
                return potential[i] < currentBest[i];
            }
        }
        return false;
}



template<int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints) {
    if (newPoints.empty()) {
        root = nullptr;
    } else {
        vector<Point<Dim>> copyPoints = newPoints;
        root = helperConstructor(copyPoints, 0, newPoints.size() - 1, 0);
    }
}


template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */

  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
}




template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    if (root == nullptr)
        return Point<Dim>();

    Point<Dim> best = root->point;
    helperFindNearest(root, query, 0, best);
    return best;
}

template <int Dim>
void KDTree<Dim>::helperFindNearest(KDTreeNode* current, const Point<Dim>& query, int depth, Point<Dim>& best) const
{
    if (current == nullptr) {
        return;
    }
    Point<Dim> currPoint = current->point;
    int currDim = depth % Dim;
    bool goToLeft = smallerDimVal(query, currPoint, currDim);
    KDTreeNode* nextBranch = (goToLeft) ? current->left : current->right;
    KDTreeNode* otherBranch = (goToLeft) ? current->right : current->left;
    helperFindNearest(nextBranch, query, depth + 1, best);
    if (shouldReplace(query, best, currPoint))
        best = currPoint;
    long radiusSquared = 0;
    for (int i = 0; i < Dim; i++) {
        long diff = query[i] - best[i];
        radiusSquared += diff * diff;
    }
    long dist = query[currDim] - currPoint[currDim];
    if (radiusSquared >= (dist * dist))
        helperFindNearest(otherBranch, query, depth + 1, best);
}



template <typename RandIter, typename Comparator>
void select(RandIter begin, RandIter end, RandIter k, Comparator cmp) {
    if (begin == end) {
        return;
    }
    RandIter pivot = partition(begin, end, begin + rand() % (end - begin), cmp);
    if (k == pivot) {
        return;
    } else if (k < pivot) {
        select(begin, pivot, k, cmp);
    } else {
        select(pivot + 1, end, k, cmp);
    }      
}

template <typename RandIter, typename Comparator>
RandIter partition(RandIter begin, RandIter end, RandIter pivot, Comparator cmp) {
    auto pivVal = *pivot;
    std::iter_swap(pivot, end - 1);
    RandIter tmp = begin;
    for (RandIter iter = begin; iter != end - 1; ++iter) {
        if (cmp(*iter, pivVal)) {
            std::iter_swap(tmp, iter);
            ++tmp;
        }
    }
    std::iter_swap(end - 1, tmp);
    return tmp;
}


template<int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::helperConstructor(vector<Point<Dim>>& points, int left, int right, int dimension) {
    if (left > right) {
        return nullptr;
    }
    int median = (left + right) / 2;
    select(points.begin() + left, points.begin() + right + 1, points.begin() + median,
           [dimension](const Point<Dim>& a, const Point<Dim>& b) {
               return smallerDimVal(a, b, dimension);
           });
    KDTreeNode* val = new KDTreeNode(points[median]);
    int nextDimension = (dimension + 1) % Dim;
    val->left = helperConstructor(points, left, median - 1, nextDimension);
    val->right = helperConstructor(points, median + 1, right, nextDimension);
    return val;
}

