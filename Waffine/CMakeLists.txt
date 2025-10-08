add_library(waffine WaffineForm.cpp WaffineForm.hpp noise_symbol.cpp noise_symbol.h)
# Note: assume this library was compiled elsewhere.
# We include the required header, but not compiled object files.
# Winterval release: https://github.com/Willmo3/Winterval/tree/main/release/winterval-v0.0.1.tar.gz
target_link_libraries(waffine winterval)
