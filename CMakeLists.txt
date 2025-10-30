# --- CMakeLists.txt for BlackHoleProject ---
# Requires CMake version 3.10 or higher for modern features and SFML's FindSFML module.
cmake_minimum_required(VERSION 3.10)

# Project definition
project(BlackHoleProject CXX)

# Set the C++ standard to C++17. This is essential for modern C++ features.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # Require C++17, don't just use it if available
set(CMAKE_CXX_EXTENSIONS OFF)      # Disable compiler-specific extensions for portability

# --- SFML Configuration ---
# Find SFML version 2.6.2 with the required components.
# SFML::* imported targets will be available after this.
# SFML_DIR should be set by the user to the 'lib/cmake/SFML' directory
# of their SFML installation (e.g., C:/SFML_2_6_2_VC17/lib/cmake/SFML).
find_package(SFML 2.6.2 COMPONENTS Graphics Window System REQUIRED)

# --- Source File Definitions ---
# Application source file
set(APP_SOURCES main.cpp)

# ImGui library source files
set(IMGUI_SOURCES
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_tables.cpp
    # imgui/imgui_demo.cpp # Uncomment if you want to include the ImGui demo window
)

# ImGui-SFML integration source file
set(IMGUI_SFML_SOURCES
    imgui-sfml/imgui-SFML.cpp
)

# Combine all source files into one list for the executable
set(ALL_SOURCES
    ${APP_SOURCES}
    ${IMGUI_SOURCES}
    ${IMGUI_SFML_SOURCES}
)

# --- Executable Target ---
# Create the main executable target
add_executable(${PROJECT_NAME} ${ALL_SOURCES})

# --- Target Properties & Linkage ---

# Set include directories for the executable.
# PRIVATE means these include paths are only used for compiling this target.
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/imgui"      # ImGui headers
    "${CMAKE_CURRENT_SOURCE_DIR}/imgui-sfml" # ImGui-SFML headers
    # SFML include directory is usually handled by the SFML::* targets,
    # but explicitly adding SFML_INCLUDE_DIR ensures it.
    "${SFML_INCLUDE_DIR}"
)

# Link necessary libraries to the executable.
# PRIVATE means these libraries are only linked by this target, not exported for others.
target_link_libraries(${PROJECT_NAME} PRIVATE
    SFML::Graphics
    SFML::Window
    SFML::System
)

# --- Output Directory Configuration ---
# Configure where the compiled executable and libraries will be placed.
# This makes the build directory more organized.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# --- Windows Specific: Copy DLLs for Runtime ---
# This block handles copying SFML DLLs (and potentially ImGui-SFML DLL if built separately)
# to the executable's output directory, which is necessary for running on Windows.
if(WIN32)
    # Get the base directory for SFML's libraries.
    # We query the location of SFML::Graphics which should be representative.
    get_target_property(SFML_GRAPHICS_LIB_PATH SFML::Graphics LOCATION)
    if(SFML_GRAPHICS_LIB_PATH)
        get_filename_component(SFML_DLL_DIR "${SFML_GRAPHICS_LIB_PATH}" DIRECTORY)
        message(STATUS "SFML DLL base directory found: ${SFML_DLL_DIR}")

        # List of SFML DLLs to copy. SFML 2.x typically uses a '-2.dll' suffix.
        # We need to find both release (e.g., sfml-graphics-2.dll) and debug (e.g., sfml-graphics-d-2.dll).
        file(GLOB SFML_ALL_DLL_FILES "${SFML_DLL_DIR}/sfml-*-2.dll")

        # If ImGui-SFML was built as a shared library (DLL), add its DLLs.
        # In this project, imgui-SFML.cpp is compiled directly into the executable,
        # so it won't produce a separate DLL unless you change its build type.
        # However, this line is kept for robustness in case you modify your setup.
        file(GLOB IMGUI_SFML_DLL_FILES "${SFML_DLL_DIR}/imgui-sfml-*.dll")
        if(IMGUI_SFML_DLL_FILES)
            list(APPEND SFML_ALL_DLL_FILES ${IMGUI_SFML_DLL_FILES})
        endif()

        if(SFML_ALL_DLL_FILES)
            # Add a custom command to execute after the project builds.
            # This command copies the found DLLs to the directory where the executable is.
            add_custom_command(
                TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${SFML_ALL_DLL_FILES}
                    "$<TARGET_FILE_DIR:${PROJECT_NAME}>/"
                COMMENT "Copying SFML and potentially ImGui-SFML DLLs to executable directory."
            )
        else()
            message(WARNING "No SFML DLLs (sfml-*-2.dll) found in ${SFML_DLL_DIR}. "
                            "Please ensure SFML is correctly installed and its DLLs are present.")
        endif()
    else()
        message(WARNING "Could not determine SFML::Graphics library location. "
                        "Manual DLL copying might be required if SFML DLLs are not found automatically.")
    endif()
endif()