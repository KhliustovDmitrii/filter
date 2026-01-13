# Utilities library
add_library(utilities STATIC
    mathematics/fourier/fft23.cpp
    mathematics/matrix/matrix.cpp
    mathematics/special_functions/special_functions.cpp
    mathematics/spline/spline.cpp
    parsing/simple_parser/simple_parser.cpp
)

target_include_directories(utilities PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to core if needed
target_link_libraries(utilities PRIVATE core)