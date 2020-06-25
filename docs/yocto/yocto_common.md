# Yocto/Common: Common utilities

Yocto/Common is a collection of utilities helpful in implementing other
Yocto/GL libraries. Yocto/Common is implemented in `yocto_common.h`.

**This library is experimental** and will be documented appropriately when
the code reaches stability.

<!--

## Collection helpers

The library contains a set of helpers that make it easier to use the STL.
Here we also define the vocabulary types used in the rest of Yocto/GL.

1. check whether a value is in a container with `contain()`


## Python-like iterators and collection helpers

This library includes a set of functions to help use C++ collections with
more ease, inspired by Python. All functions and operators are defined in
the yocto namespace so they will not affect the code outside. But within
the Yocto/GL collection they are the best way to do this.

1. use `range()` to iterato over an integer sequence
2. use `enumerate()` to iteratare over a vector and number its elements
3. use opeartors + to either concatenate two vectors or a vector and an
   element
4. use operators += to append an element or a vector to a given
vector


## Concurrency utilities

C++ has very basic supprt for concurrency and most of it is still platform
dependent. We provide here very basic support for concurrency utlities
built on top of C++ low-level threading and synchronization.

1. use `concurrent_queue()` for communicationing values between threads
2. use `parallel_for()` for basic parallel for loops

-->
