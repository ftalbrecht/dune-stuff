// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_COMMON_COMPILER_HH
#define DUNE_STUFF_COMMON_COMPILER_HH

#if defined(__GNUC__)
# if defined(__GNUC_PATCHLEVEL__)
#  define DS_GNUC_VERSION (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100 \
                            + __GNUC_PATCHLEVEL__)
# else
#  define DS_GNUC_VERSION (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100)
# endif
#endif

// gtest impl of disabling tests via prefixed name
// seems broken for ICC and gcc < 4.9
#if defined(__INTEL_COMPILER) or (DS_GNUC_VERSION < 40900)
# define DS_COMPILE_DISABLED_TESTS 0
#else
# define DS_COMPILE_DISABLED_TESTS 1
#endif

#endif // DUNE_STUFF_COMMON_COMPILER_HH
