# Core library containing all implementations
add_library(core STATIC
    # Filter implementations
    filter/kalman_extended/kalman_extended.cpp
    filter/kalman_unscented/kalman_unscented.cpp
    
    # Input filter
    input_filter/averaging_filter/averaging_filter.cpp
    
    # Regularizer
    regularizer/boundary_projector/boundary_projector.cpp
    
    # Updater
    updater/decay_updater/decay_updater.cpp
)

# Interface include directories
target_include_directories(core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link types (headers only)
target_sources(core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/aggregator/aggregator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/data_loader/data_loader.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/filter/probabilistic_filter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/input_filter/input_filter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/model/model.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/reader/reader.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/regularizer/regularizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/updater/updater.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/writer/writer.h
)

# Add type implementations that exist
target_sources(core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/data_loader/data_loader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../types/filter/probabilistic_filter.cpp
)