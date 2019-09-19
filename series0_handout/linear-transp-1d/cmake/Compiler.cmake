
if(CMAKE_VERSION VERSION_GREATER 3.1.0 OR CMAKE_VERSION VERSION_EQUAL 3.1.0)
    # use c++11
    # only valid for new platforms
    set(CMAKE_CXX_STANDARD 11)
  else()
    if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	# for older cmake versions
	# (note, this CXX flag is only valid for clang and gcc, for MSVC it is not needed)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
endif()





if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # (this is only applicable on clang)			      
    # ignore some mathgl warnings
    add_definitions( -Wno-deprecated-register -Wno-return-type-c-linkage)
endif()


if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    add_definitions(-Wno-deprecated-declarations -Wno-ignored-attributes)
endif()
