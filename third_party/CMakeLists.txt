set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(QPOASES_AVOID_LA_NAMING_CONFLICTS ON)
add_subdirectory(qpOASES)

add_subdirectory(matplotlib-cpp)
include_directories(matplotlib-cpp)

set(ENABLE_OSQP TRUE)
set(ENABLE_ECOS TRUE)
add_subdirectory(Epigraph)

# add_subdirectory(qpSWIFT)
