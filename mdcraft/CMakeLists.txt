configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configuration.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/configuration.h
)

add_subdirectory(data)
if (mdcraft_ENABLE_MPI)
add_subdirectory(decomp)
endif()
add_subdirectory(io)
add_subdirectory(lattice)
add_subdirectory(tools)
add_subdirectory(neibs)
add_subdirectory(solver)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/configuration.h 
        DESTINATION include/${CMAKE_PROJECT_NAME}
)
