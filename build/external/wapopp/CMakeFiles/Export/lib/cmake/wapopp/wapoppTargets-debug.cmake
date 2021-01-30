#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "wapopp::wapopp" for configuration "Debug"
set_property(TARGET wapopp::wapopp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(wapopp::wapopp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libwapopp.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS wapopp::wapopp )
list(APPEND _IMPORT_CHECK_FILES_FOR_wapopp::wapopp "${_IMPORT_PREFIX}/lib/libwapopp.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
