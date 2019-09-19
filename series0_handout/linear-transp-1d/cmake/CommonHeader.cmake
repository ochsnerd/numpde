
SET(NUMPDE_USE_HUNTER OFF CACHE BOOL "Use Hunter for automatic compiling of additional packages")

SET(HUNTER_ROOT ${CMAKE_BINARY_DIR}/hunter_root)

IF(${NUMPDE_USE_HUNTER})
	include("cmake/HunterGate.cmake")
	HunterGate(
 	   URL "https://github.com/ruslo/hunter/archive/v0.18.3.tar.gz"
    	   SHA1 "3e29b599966d2acc7ddcb10cc50522d794d24d93"
	   )

ELSE()
	MACRO(hunter_add_package args)
	     # do nothing
	ENDMACRO()
      ENDIF()
# set cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include("cmake/PackageFind.cmake")

