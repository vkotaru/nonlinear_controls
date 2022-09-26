#add_subdirectory(qpOASES)
set(qpOASES_SRC_DIR qpOASES)
set(qpOASES_SRCS
        ${qpOASES_SRC_DIR}/src/BLASReplacement.cpp
        ${qpOASES_SRC_DIR}/src/Bounds.cpp
        ${qpOASES_SRC_DIR}/src/Constraints.cpp
        ${qpOASES_SRC_DIR}/src/Flipper.cpp
        ${qpOASES_SRC_DIR}/src/Indexlist.cpp
        ${qpOASES_SRC_DIR}/src/LAPACKReplacement.cpp
        ${qpOASES_SRC_DIR}/src/Matrices.cpp
        ${qpOASES_SRC_DIR}/src/MessageHandling.cpp
        ${qpOASES_SRC_DIR}/src/Options.cpp
        ${qpOASES_SRC_DIR}/src/OQPinterface.cpp
        ${qpOASES_SRC_DIR}/src/QProblemB.cpp
        ${qpOASES_SRC_DIR}/src/QProblem.cpp
        ${qpOASES_SRC_DIR}/src/SolutionAnalysis.cpp
        ${qpOASES_SRC_DIR}/src/SparseSolver.cpp
        ${qpOASES_SRC_DIR}/src/SQProblem.cpp
        ${qpOASES_SRC_DIR}/src/SQProblemSchur.cpp
        ${qpOASES_SRC_DIR}/src/SubjectTo.cpp
        ${qpOASES_SRC_DIR}/src/Utils.cpp
        )
add_library(qpOASES_cpp ${qpOASES_SRCS})
target_include_directories(qpOASES_cpp PUBLIC . ${qpOASES_SRC_DIR}/include)
# include_directories(matplotlib-cpp)

set(ENABLE_OSQP TRUE)
set(ENABLE_ECOS TRUE)
add_subdirectory(Epigraph)

# add_subdirectory(qpSWIFT)
