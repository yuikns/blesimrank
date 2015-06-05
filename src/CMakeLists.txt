ADD_EXECUTABLE(blesimrank main.cc)

ADD_DEPENDENCIES(blesimrank argcv)


# http://www.gnu.org/software/libc/manual/html_node/Feature-Test-Macros.html#index-_005fPOSIX_005fSOURCE 
#ADD_DEFINITIONS( " -D_POSIX_SOURCE ") 
#INCLUDE(ExternalProject) # IF NEED EXTEND PROJECT
#INCLUDE(ExternalLeveldb)
#MESSAGE(STATUS "$$$$$$$$$$ INCLUDE LEVELDB DIR " ${LEVELDB_INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${LEVELDB_INCLUDE_DIR})

#TARGET_LINK_LIBRARIES(blesimrank argcv_static)
#TARGET_LINK_LIBRARIES(blesimrank argcv_sys_static)
#TARGET_LINK_LIBRARIES(blesimrank argcv_wrapper_leveldb_static)
TARGET_LINK_LIBRARIES(blesimrank m)
TARGET_LINK_LIBRARIES(blesimrank pthread)

#  executable properties
SET_TARGET_PROPERTIES(blesimrank PROPERTIES OUTPUT_NAME "blesimrank")





