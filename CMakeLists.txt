PROJECT(IJKMESHINFO)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${IJK_DIR})
  GET_FILENAME_COMPONENT(IJK_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(IJK_DIR ${IJK_ABSOLUTE_PATH} CACHE PATH "IJK directory")
ENDIF (NOT DEFINED ${IJK_DIR})

SET(CMAKE_INSTALL_PREFIX "${IJK_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${IJK_DIR}/lib CACHE PATH "Library directory")
SET(IJKMESHINFO_DIR "src/ijkmeshinfo")

#---------------------------------------------------------


INCLUDE_DIRECTORIES("${IJK_DIR}/include")
LINK_DIRECTORIES("${IJK_DIR}/lib")

SET(CMAKE_CXX_FLAGS "-std=c++11")

ADD_EXECUTABLE(ijkmeshinfo ijkmeshinfo.cxx ijkmeshinfo_compute.cxx
                           ijkmeshinfoIO.cxx)

SET(CMAKE_INSTALL_PREFIX ${IJK_DIR})
INSTALL(TARGETS ijkmeshinfo DESTINATION "bin/linux")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${IJKMESHINFO_DIR}/ijkmeshinfo.tar ${IJKMESHINFO_DIR}/README ${IJKMESHINFO_DIR}/INSTALL ${IJKMESHINFO_DIR}/*.cxx ${IJKMESHINFO_DIR}/*.h ${IJKMESHINFO_DIR}/*.txx ${IJKMESHINFO_DIR}/CMakeLists.txt )

ADD_CUSTOM_TARGET(doc COMMAND doxygen ijkmeshinfo_doxygen.config)
