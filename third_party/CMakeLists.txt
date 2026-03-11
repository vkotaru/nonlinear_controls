# ------------------------------------------------------------------------------
# Third-party dependencies for nonlinear_controls
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# qpOASES - Quadratic Programming solver
# ------------------------------------------------------------------------------
# Disable qpOASES examples to speed up build
set(QPOASES_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(qpOASES)

# Note: matplotlib-cpp is header-only, included directly where needed
