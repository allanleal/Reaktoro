# Configure the setup.py file
configure_file(setup.py.in ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

# Create pybind11-stubgen command string depending if this program has been found or not
if(PYBIND11_STUBGEN)
    set(PYBIND11_STUBGEN_COMMAND ${PYBIND11_STUBGEN} --ignore-invalid=all --no-setup-py --log-level=ERROR reaktoro)
else()
    set(PYBIND11_STUBGEN_COMMAND ${CMAKE_COMMAND} -E true)  # do nothing when it it's time to generate python stubs for reaktoro in the target below
endif()

# Create a custom target to build the python package during build stage
add_custom_target(reaktoro-setuptools ALL
    DEPENDS reaktoro4py
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/reaktoro reaktoro
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:reaktoro4py> reaktoro
    COMMAND ${CMAKE_COMMAND} -E env
        "PATH=${REAKTORO_PATH}"
        "PYTHONPATH=${REAKTORO_PYTHONPATH}"
            ${PYBIND11_STUBGEN_COMMAND}
    COMMAND ${CMAKE_COMMAND} -E copy_directory stubs/reaktoro-stubs/ reaktoro
    COMMAND ${CMAKE_COMMAND} -E rm -rf build  # remove build dir created by previous `python setup.py install` commands (see next) to ensure fresh rebuild since changed python files are not overwritten even with --force option
    COMMAND ${PYTHON_EXECUTABLE} setup.py --quiet install --force --no-compile --prefix=installed/$<CONFIG>
    COMMAND ${CMAKE_COMMAND} -E rm reaktoro/$<TARGET_FILE_NAME:reaktoro4py>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Ensure the path where the python package is installed is not empty
if(NOT DEFINED REAKTORO_PYTHON_INSTALL_PREFIX)
    set(REAKTORO_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# Create an install target for the reaktoro python package
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/installed/$<CONFIG>/
    DESTINATION ${REAKTORO_PYTHON_INSTALL_PREFIX})
