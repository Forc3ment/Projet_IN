#############################################################################
#
# This file is part of the ViSP software.
# Copyright (C) 2005 - 2017 by Inria. All rights reserved.
#
# This software is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# ("GPL") version 2 as published by the Free Software Foundation.
# See the file LICENSE.txt at the root directory of this source
# distribution for additional information about the GNU GPL.
#
# For using ViSP with software that can not be combined with the GNU
# GPL, please contact Inria about acquiring a ViSP Professional
# Edition License.
#
# See http://visp.inria.fr for more information.
#
# This software was developed at:
# Inria Rennes - Bretagne Atlantique
# Campus Universitaire de Beaulieu
# 35042 Rennes Cedex
# France
#
# If you have questions regarding the use of this file, please contact
# Inria at visp@inria.fr
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Description:
# Try to find Intel RealSense SDK to work with R200, F200 and SR300 devices.
# Once run this will define: 
#
# REALSENSE_FOUND
# REALSENSE_INCLUDE_DIRS
# REALSENSE_LIBRARIES
#
# Authors:
# Fabien Spindler
#
#############################################################################

set(REALSENSE_INC_SEARCH_PATH /usr/local/include)
set(REALSENSE_LIB_SEARCH_PATH /usr/local/lib)

if(MSVC)
  list(APPEND REALSENSE_INC_SEARCH_PATH "C:/librealsense/include")

  list(APPEND REALSENSE_INC_SEARCH_PATH $ENV{REALSENSE_HOME}/include)

  if(CMAKE_CL_64)
    list(APPEND REALSENSE_LIB_SEARCH_PATH "C:/librealsense/bin/x64")
    list(APPEND REALSENSE_LIB_SEARCH_PATH $ENV{REALSENSE_HOME}/bin/x64)
  else()
    list(APPEND REALSENSE_LIB_SEARCH_PATH "C:/librealsense/bin/Win32")
    list(APPEND REALSENSE_LIB_SEARCH_PATH $ENV{REALSENSE_HOME}/bin/Win32)
  endif()

else()
  list(APPEND REALSENSE_INC_SEARCH_PATH /usr/include)
  list(APPEND REALSENSE_LIB_SEARCH_PATH /usr/lib)

  list(APPEND REALSENSE_INC_SEARCH_PATH $ENV{REALSENSE_HOME}/include)
  list(APPEND REALSENSE_LIB_SEARCH_PATH $ENV{REALSENSE_HOME}/lib)
endif()

find_path(REALSENSE_INCLUDE_DIRS librealsense/rs.hpp
  PATHS
    ${REALSENSE_INC_SEARCH_PATH}
)

if(MSVC)
  find_library(REALSENSE_LIBRARIES_OPT
    NAMES realsense
    PATHS
      ${REALSENSE_LIB_SEARCH_PATH}
  )
  find_library(REALSENSE_LIBRARIES_DBG
    NAMES realsense-d
    PATHS
      ${REALSENSE_LIB_SEARCH_PATH}
  )
  if(REALSENSE_LIBRARIES_OPT)
    set(REALSENSE_LIBRARIES optimized ${REALSENSE_LIBRARIES_OPT})
  endif()
  if(REALSENSE_LIBRARIES_DBG)
    list(APPEND REALSENSE_LIBRARIES debug ${REALSENSE_LIBRARIES_DBG})
  endif()

  mark_as_advanced(REALSENSE_LIBRARIES_OPT REALSENSE_LIBRARIES_DBG)
else()
  find_library(REALSENSE_LIBRARIES
    NAMES realsense
    PATHS
      ${REALSENSE_LIB_SEARCH_PATH}
  )
endif()

if(REALSENSE_LIBRARIES AND REALSENSE_INCLUDE_DIRS)
  set(REALSENSE_FOUND TRUE)
else()
  set(REALSENSE_FOUND FALSE)
endif()
  
mark_as_advanced(
  REALSENSE_INCLUDE_DIRS
  REALSENSE_LIBRARIES
  REALSENSE_INC_SEARCH_PATH
  REALSENSE_LIB_SEARCH_PATH
)