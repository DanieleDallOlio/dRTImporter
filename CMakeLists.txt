cmake_minimum_required(VERSION 3.16.3)

project(dRTImporter)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE
  "https://github.com/DanieleDallOlio/dRTImporter")

set(EXTENSION_CATEGORY
  "DICOM")

set(EXTENSION_CONTRIBUTORS
  "Daniele Dall'Olio (University of Bologna)")

set(EXTENSION_DESCRIPTION
  "DICOM plugin for importing MiM-style dynamic RT Structure Sets as segmentation sequences in 3D Slicer.")

set(EXTENSION_ICONURL
  "https://raw.githubusercontent.com/DanieleDallOlio/dRTImporter/main/dRTImporter.png")

# SlicerRT supplies the RT Structure Set reader used for each temporal frame.
set(EXTENSION_DEPENDS
  "SlicerRT")

#-----------------------------------------------------------------------------
# Find Slicer
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(dRTImporter)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
