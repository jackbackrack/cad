#pragma once
#include <tuple>

template<class Fn, class T0> auto tuple_map(Fn&& fn, const std::tuple<T0>& t)
 -> decltype(std::make_tuple(fn(std::get<0>(t))))
    { return std::make_tuple(fn(std::get<0>(t))); }
template<class Fn, class T0, class T1> auto tuple_map(Fn&& fn, const std::tuple<T0, T1>& t) 
 -> decltype(std::make_tuple(fn(std::get<0>(t)),fn(std::get<1>(t))))
    { return std::make_tuple(fn(std::get<0>(t)),fn(std::get<1>(t))); }
template<class Fn, class T0, class T1, class T2> auto tuple_map(Fn&& fn, const std::tuple<T0, T1, T2>& t) 
 -> decltype(std::make_tuple(fn(std::get<0>(t)),fn(std::get<1>(t)),fn(std::get<2>(t))))
    { return std::make_tuple(fn(std::get<0>(t)),fn(std::get<1>(t)),fn(std::get<2>(t))); }


template<class... Types> struct ZippedView {
  std::tuple<Array<Types>&...> arrays;

  using IndexT = decltype(std::get<0>(arrays).size());
  using IndexDiffT = decltype(declval<IndexT>() - declval<IndexT>());

 private:
    struct GetRefHelper {
      IndexT i;
      template<class T> std::reference_wrapper<T> operator()(const Array<T>& array) { return array[i]; }
    };
 public:
  using Val = std::tuple<Types...>;

  struct Ref {
    ZippedView* view;
    IndexT i;
    template<size_t I> auto get() -> decltype(std::get<I>(view->arrays)[i]) { return std::get<I>(view->arrays)[i]; }
    template<size_t I> auto get() const -> decltype(std::get<I>(view->arrays)[i]) { return std::get<I>(view->arrays)[i]; }

    std::tuple<Types&...> as_tuple() const {
      return tuple_map(GetRefHelper{i}, view->arrays);
    }

    Ref& operator=(const Ref& rhs) {
      as_tuple() = rhs.as_tuple();
      return *this;
    }

    Ref& operator=(const std::tuple<Types...>& rhs) {
      as_tuple() = rhs;
      return *this;
    }

    // Implicit conversion to tuple
    operator Val() const
    { return as_tuple(); }

    // Overload that accepts non-lvalue arguments since we are swapping the underlying memory
    inline friend void swap(Ref lhs, Ref rhs) {
      auto lt = lhs.as_tuple();
      auto rt = rhs.as_tuple();
      swap(lt, rt);
    }
  };

  struct Ptr {
    ZippedView* view;
    IndexT i;
    template<size_t I> auto get() -> decltype(&std::get<I>(view->arrays)[i]) { return &std::get<I>(view->arrays)[i]; }
    template<size_t I> auto get() const -> decltype(&std::get<I>(view->arrays)[i]) { return &std::get<I>(view->arrays)[i]; }
  };

  struct Iterator : public std::iterator<std::random_access_iterator_tag, Val, IndexDiffT, Ptr, Ref> {
    ZippedView* view;
    IndexT i;
    Iterator(ZippedView& _view, IndexT _i)
     : view(&_view)
     , i(_i)
    { }
    IndexDiffT operator-(const Iterator& rhs) const { assert(view == rhs.view); return i - rhs.i; }
    Iterator& operator--() { --i; return *this; }
    Iterator& operator++() { ++i; return *this; }
    Iterator operator+(const IndexDiffT rhs) const { return Iterator{*view, i+rhs}; }
    Iterator& operator+=(const IndexDiffT rhs) { i += rhs; return *this; }
    Ref operator*() const { return Ref{view, i}; }

    bool operator==(const Iterator& rhs) const { assert(view == rhs.view); return i == rhs.i; } 
    bool operator!=(const Iterator& rhs) const { assert(view == rhs.view); return i != rhs.i; } 
    bool operator<(const Iterator& rhs) const { assert(view == rhs.view); return i < rhs.i; } 
    bool operator<=(const Iterator& rhs) const { assert(view == rhs.view); return i <= rhs.i; } 
    bool operator>(const Iterator& rhs) const { assert(view == rhs.view); return i > rhs.i; } 
    bool operator>=(const Iterator& rhs) const { assert(view == rhs.view); return i > rhs.i; } 
  };

  IndexT size() const { return std::get<0>(arrays).size(); }
  Iterator begin() { return Iterator{*this, 0}; }
  Iterator end() { return Iterator{*this, size()}; }

  Ref operator[](const IndexT i) { return Ref{this, i}; }
};

template<class... Types> auto zipped_view(Array<Types>&... arrays) -> ZippedView<Types...>
{ return ZippedView<Types...>{std::tuple<Array<Types>&...>{arrays...}}; }
