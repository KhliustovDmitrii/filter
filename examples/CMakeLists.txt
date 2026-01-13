# Airborne Electromagnetics example
add_executable(airborne_electromagnetics
    airborne_electromagnetics/main.cpp
    airborne_electromagnetics/aggregator/EQUATOR_aggregator.cpp
    airborne_electromagnetics/data_loader/EQUATOR_data_loader.cpp
    airborne_electromagnetics/model/EQUATOR/EQUATOR.cpp
    airborne_electromagnetics/model/EQUATOR_C/EQUATOR_C.cpp
)

target_include_directories(airborne_electromagnetics PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/airborne_electromagnetics
    ${CMAKE_CURRENT_SOURCE_DIR}/airborne_electromagnetics/model
    ${CMAKE_CURRENT_SOURCE_DIR}/airborne_electromagnetics/model/EQUATOR
    ${CMAKE_CURRENT_SOURCE_DIR}/airborne_electromagnetics/model/EQUATOR_C
)

target_link_libraries(airborne_electromagnetics PRIVATE core utilities)

# Polynomial model example (if we keep the Makefile, or we can migrate)
add_executable(polynomial_model
    polynomial_model/main.cpp
    polynomial_model/model/polynomial_model.cpp
)

target_include_directories(polynomial_model PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/polynomial_model
)

target_link_libraries(polynomial_model PRIVATE core utilities)