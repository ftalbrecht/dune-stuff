MACRO(ADD_CXX_FLAGS)
	FOREACH( ARG ${ARGN} )
# 		ADD_DEFINITIONS(${ARG})
		LIST(APPEND MY_CXX_FLAGS ${ARG} )
	ENDFOREACH( ARG )
ENDMACRO(ADD_CXX_FLAGS)

function(TO_LIST_SPACES _LIST_NAME OUTPUT_VAR)
  set(NEW_LIST_SPACE)
  foreach(ITEM ${${_LIST_NAME}})
    set(NEW_LIST_SPACE ${NEW_LIST_SPACE} ${ITEM})
  endforeach()
#   string(STRIP ${NEW_LIST_SPACE} NEW_LIST_SPACE)
  set(${OUTPUT_VAR} "${NEW_LIST_SPACE}" PARENT_SCOPE)
endfunction()

MACRO(INCLUDE_DIR)
	FOREACH( ARG ${ARGN} )
		IF(IS_DIRECTORY ${ARG} )
			LIST(APPEND MY_CXX_FLAGS "-I${ARG}" )
			LIST(APPEND PLAIN_INCLUDE_DIRS "-I${ARG}" )
			INCLUDE_DIRECTORIES(${ARG})
		ELSE(IS_DIRECTORY ${ARG} )
			MESSAGE( STATUS "Include directory ${ARG} does not exist" )
		ENDIF(IS_DIRECTORY ${ARG} )
    ENDFOREACH( ARG )
ENDMACRO(INCLUDE_DIR)

MACRO(INCLUDE_SYS_DIR)
	FOREACH( ARG ${ARGN} )
		IF(IS_DIRECTORY ${ARG} )
			LIST(APPEND MY_CXX_FLAGS "-isystem${ARG}")
			LIST(APPEND PLAIN_INCLUDE_DIRS "-I${ARG}" )
			INCLUDE_DIRECTORIES(${ARG})
		ELSE(IS_DIRECTORY ${ARG} )
			MESSAGE( STATUS "Include directory ${ARG} does not exist" )
		ENDIF(IS_DIRECTORY ${ARG} )
    ENDFOREACH( ARG )
ENDMACRO(INCLUDE_SYS_DIR)

MACRO( HEADERCHECK )
	ADD_CUSTOM_TARGET( headercheck SOURCES ${ARGN} )
	FOREACH( HEADER ${ARGN} )
		GET_FILENAME_COMPONENT( fn ${HEADER} NAME )
		SET( TEST_NAME "headercheck_${fn}")
		TO_LIST_SPACES( MY_CXX_FLAGS TEST_NAME_FLAGS )
		SET(XARGS ${TEST_NAME_FLAGS} -c ${HEADER} -o ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.o )
		ADD_CUSTOM_TARGET(  ${TEST_NAME} ${CMAKE_CXX_COMPILER} ${XARGS} )
		ADD_TEST( ${TEST_NAME} ${CMAKE_CXX_COMPILER} ${XARGS} )
		add_dependencies( headercheck ${TEST_NAME} )
	ENDFOREACH( HEADER )
ENDMACRO( HEADERCHECK )

MACRO( ADD_CPPCHECK )
	FOREACH( SOURCEFILE ${ARGN} )
		LIST( APPEND CHECKPATHS ${SOURCEFILE} )
	ENDFOREACH( SOURCEFILE )	
	TO_LIST_SPACES( CPPCHECK_FLAGS CPPCHECK_FLAGS_SPLIT )
	ADD_CUSTOM_TARGET(  cppcheck cppcheck "--xml" "--enable=all" "-f" "--quiet" "-j 3" 
						${CPPCHECK_FLAGS_SPLIT} ${CHECKPATHS} "2>cppcheck.xml" )
ENDMACRO( ADD_CPPCHECK )

MACRO( ADD_DUNE_MODULES )
	FOREACH( MODULE ${ARGN} )
		INCLUDE_SYS_DIR( ${CMAKE_CURRENT_SOURCE_DIR}/../dune-${MODULE} )
		LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/.libs )
		FILE( GLOB_RECURSE tmp_header "${CMAKE_CURRENT_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/*.hh" )
		LIST( APPEND DUNE_HEADERS ${tmp_header} )
	ENDFOREACH(MODULE)
ENDMACRO( ADD_DUNE_MODULES )

MACRO( ADD_MY_MODULES )
	FOREACH( MODULE ${ARGN} )
		INCLUDE_DIR( ${CMAKE_CURRENT_SOURCE_DIR}/../dune-${MODULE} )
		LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/.libs )
		FILE( GLOB_RECURSE tmp_header "${CMAKE_CURRENT_SOURCE_DIR}/../dune-${MODULE}/${MODULE}/*.hh" )
		LIST( APPEND DUNE_HEADERS ${tmp_header} )
	ENDFOREACH(MODULE)
ENDMACRO( ADD_MY_MODULES )

MACRO( ADD_DUNE_EXECUTABLE target sources headers libs )
	ADD_EXECUTABLE(${target} ${sources} ${headers} ${DUNE_HEADERS} )
	#for some $@#&& reason these targets DO NOT inherit flags added via TARGET_LINK_LIBRARIES nor INCLUDE_DIRECTORIES
	#so we need some magic to readd them
# 	get_target_property(tmp_flags ${target} COMPILE_FLAGS  )
# 	if( NOT ${tmp_flags} )
# 		set( tmp_flags "" )
# 	endif( NOT ${tmp_flags} )
# 	LIST(APPEND MY_CXX_FLAGS ${tmp_flags} )
# 	foreach(arg ${MY_CXX_FLAGS} )
# 		set(bar "${bar} ${arg}")
# 	endforeach(arg ${MY_CXX_FLAGS} )
# 	set_target_properties( ${target} PROPERTIES COMPILE_FLAGS  ${bar} )
# 	get_target_property(tmp_flags ${target} LINK_FLAGS  )
# 	if( NOT ${tmp_flags} )
# 		set( tmp_flags "" )
# 	endif( NOT ${tmp_flags} )
# # 	message(STATUS ${tmp_flags} )
# 	set( tmp_flags "${bar} ${tmp_flags}" )
# 	set_target_properties( ${target} PROPERTIES LINK_FLAGS  ${tmp_flags} )
	TARGET_LINK_LIBRARIES(${target} ${libs} )
	
	add_dependencies( ${target} config_refresh )
ENDMACRO( ADD_DUNE_EXECUTABLE )

MACRO( SET_CONFIGHEADER_VARS )
	IF( IS_DIRECTORY ${ALUGRID_BASE_PATH} )
		SET( ALUGRID_FOUND "1" )
	ELSE( IS_DIRECTORY ${ALUGRID_BASE_PATH} )
		SET( ALUGRID_FOUND "0" )
	ENDIF( IS_DIRECTORY ${ALUGRID_BASE_PATH} )
ENDMACRO( SET_CONFIGHEADER_VARS )

add_custom_target( config_refresh
				${CMAKE_CURRENT_SOURCE_DIR}/cmake/regen_config_header.sh ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h
				)

ADD_CXX_FLAGS( "-DHAVE_CMAKE_CONFIG -std=c++0x" )

INCLUDE (CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("tr1/array" HAVE_TR1_ARRAY)
CHECK_INCLUDE_FILE_CXX("malloc.h" HAVE_MALLOC_H)

SET( CUSTOM_FLAGS
	"-Wall -Wextra -pedantic -O3 -DDNDEBUG -funroll-loops -g -ggdb -DENABLE_ADAPTIVE=1 -DADAPTIVE_SOLVER -DUSE_BFG_CG_SCHEME -fno-strict-aliasing" CACHE STRING  
	"CUSTOM FLAGS")
FIND_PACKAGE(Boost 1.35.0 REQUIRED)
INCLUDE_DIR(${Boost_INCLUDE_DIR})
FIND_PACKAGE( PkgConfig )
pkg_check_modules( CCGNU REQUIRED libccgnu2 )
ADD_CXX_FLAGS( "${CCGNU_CFLAGS}" )
INCLUDE_SYS_DIR(${Boost_INCLUDE_DIR})

#try to enable link-time-optimisation
if (CMAKE_COMPILER_IS_GNUCC)
	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	if (GCC_VERSION VERSION_GREATER 4.5 OR GCC_VERSION VERSION_EQUAL 4.5)
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
		set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto")
		set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -flto")
	endif()
endif()

ENABLE_TESTING()