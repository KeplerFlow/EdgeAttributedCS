#ifndef CORE_GADGET_HEAP_H_
#define CORE_GADGET_HEAP_H_

#include <vector>
#include "../defs.h"

namespace gadget
{
  struct Pair final
  {
    int key;
    int val;
  };
  class MinHeap final
  {
  public:
    explicit MinHeap(const int n);
    // 深拷贝构造函数
    MinHeap(const MinHeap &other)
        : heap_(other.heap_), pos_(other.pos_)
    {
      // 此处，heap_ 和 pos_ 通过它们各自的拷贝构造函数进行复制
      // 拷贝构造函数自动处理了向量的深拷贝
    }

    MinHeap(){};

    void Insert(const int key, const int value);
    void Delete(const int key);
    Pair Top() const;
    bool Contains(const int key) const
    {
      return pos_[key] != -1;
    }
    bool Empty() const
    {
      return heap_.size() == 1;
    }
    int Size() const
    {
      return heap_.size() - 1;
    }

  private:
    void Up(const int h, const int key, const int value);
    void Down(const int h, const int key, const int value);

    std::vector<Pair> heap_;
    std::vector<int> pos_;
  };
} // namespace gadget

#endif
