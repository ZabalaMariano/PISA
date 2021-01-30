include(FindPackageHandleStandardArgs)
set(${CMAKE_FIND_PACKAGE_NAME}_CONFIG ${CMAKE_CURRENT_LIST_FILE})
find_package_handle_standard_args(wapopp CONFIG_MODE)

if(NOT TARGET wapopp::wapopp)
  include("${CMAKE_CURRENT_LIST_DIR}/wapoppTargets.cmake")
endif()
