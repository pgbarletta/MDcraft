find_package(Python 3 COMPONENTS Interpreter) # Or Python 2 if needed

if(Python_FOUND)
    execute_process(COMMAND "${Python_EXECUTABLE}" -c "import mpi4py; print(mpi4py.__path__[0])"
                    OUTPUT_VARIABLE MPI4PY_PATH
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE _mpi4py_result)

    if(_mpi4py_result EQUAL 0)
        set(MPI4PY_FOUND TRUE CACHE BOOL "Whether mpi4py was found")
        set(MPI4PY_INCLUDE_DIR "${MPI4PY_PATH}/include" CACHE PATH "Path to mpi4py installation")
        # You might need to add logic to find MPI libraries if your project links against them
        # find_package(MPI REQUIRED)
        # set(MPI4PY_LIBRARIES ${MPI_LIBRARIES})
    else()
        set(MPI4PY_FOUND FALSE CACHE BOOL "Whether mpi4py was found")
    endif()
else()
    set(MPI4PY_FOUND FALSE CACHE BOOL "Whether mpi4py was found")
endif()

mark_as_advanced(MPI4PY_INCLUDE_DIR)
