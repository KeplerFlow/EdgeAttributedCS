// The @Treap class implements the treap structure.
#ifndef CORE_GADGET_TREAP_H_
#define CORE_GADGET_TREAP_H_

#include <vector>
#include "../defs.h"

namespace gadget {
class Treap final {
 public:
  explicit Treap(const int n);

  Treap(){};
  Treap(const Treap& other) : n_(other.n_), nd_(other.n_) {
    for (int i = 0; i < n_; ++i) {
      nd_[i].p = other.nd_[i].p; // 拷贝父节点索引
      nd_[i].l = other.nd_[i].l; // 拷贝左子节点索引
      nd_[i].r = other.nd_[i].r; // 拷贝右子节点索引
      nd_[i].s = other.nd_[i].s; // 拷贝子树大小
      nd_[i].w = other.nd_[i].w; // 拷贝优先级
    }
  }

  void Insert(const int x, const bool f, int& r);
  void InsertAfter(const int x, const int y, int& r);
  void Delete(const int x, int& r);
  int Merge(const int r1, const int r2);
  int Rank(const int x) const;
  int Select(const int r, const int rank) const;
  int Root(const int x) const;
  int Minimum(const int x) const;
  int Maximum(const int x) const;
  int Size(const int r) const;
  void Check(const int r) const;

 private:
  struct TreapNode final {
    int p;  // the parent
    int l;  // the left child
    int r;  // the right child
    int s;  // the size of the subtree rooted at this node
    int w;  // the priority
  };

  void LeftRotate(int& x);
  void RightRotate(int& x);
  void SubCheck(const int x) const;

  int n_;
  std::vector<TreapNode> nd_;
};
}  // namespace gadget

#endif
