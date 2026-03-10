# ------------------------------------------------------------------------------
# Third-party dependencies for nonlinear_controls
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# qpOASES - Quadratic Programming solver
# ------------------------------------------------------------------------------
# Disable qpOASES examples to speed up build
set(QPOASES_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(qpOASES)

# ------------------------------------------------------------------------------
# Epigraph - Convex optimization modeling (optional)
# ------------------------------------------------------------------------------
if(NLC_ENABLE_EPIGRAPH)
    # Configure Epigraph with OSQP and ECOS support
    set(ENABLE_OSQP TRUE CACHE BOOL "" FORCE)
    set(ENABLE_ECOS TRUE CACHE BOOL "" FORCE)
    set(ENABLE_TESTING FALSE CACHE BOOL "" FORCE)

    add_subdirectory(Epigraph)
endif()

# Note: matplotlib-cpp is header-only, included directly where needed
# Note: qpSWIFT integration is not currently used
