cmake_minimum_required(VERSION 3.28 FATAL_ERROR)

project(SlicerDentalArch)

set(EXTENSION_HOMEPAGE "https://github.com/echoechofu/SlicerDentalArch")
set(EXTENSION_CATEGORY "Dentistry")
set(EXTENSION_CONTRIBUTORS "SlicerDentalArch contributors")
set(EXTENSION_DESCRIPTION "Reviewed dental landmark algorithms and an interactive Slicer editor for width points and a planar reference arch, with live native-unit space analysis and separate revision files.")
set(EXTENSION_ICONURL "")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "NA")
set(EXTENSION_STATUS "Experimental")

find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

add_subdirectory(DentalArchEditor)

include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
