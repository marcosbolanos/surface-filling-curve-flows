# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/marcos/Code/surface-filling-curve-flows/build-libcxx/deps/geometry-central/deps/eigen-src")
  file(MAKE_DIRECTORY "/home/marcos/Code/surface-filling-curve-flows/build-libcxx/deps/geometry-central/deps/eigen-src")
endif()
file(MAKE_DIRECTORY
  "/home/marcos/Code/surface-filling-curve-flows/build-libcxx/deps/geometry-central/deps/eigen-build"
  "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix"
  "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix/tmp"
  "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp"
  "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix/src"
  "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/var/home/marcos/Code/surface-filling-curve-flows/build-libcxx/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
