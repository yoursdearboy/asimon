find_program(R_COMMAND R DOC "R executable.")

if(R_COMMAND)
    execute_process(WORKING_DIRECTORY .
            COMMAND ${R_COMMAND} RHOME
            OUTPUT_VARIABLE R_ROOT_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    find_path(R_INCLUDE_DIR R.h
            HINTS ${R_ROOT_DIR}
            PATHS /usr/local/lib /usr/local/lib64 /usr/share
            PATH_SUFFIXES include R/include
            DOC "Path to file R.h")
    find_library(R_LIBRARY
            NAMES libR.dylib
            HINTS ${R_ROOT_DIR} ${R_ROOT_DIR}/lib)
    message(${R_LIBRARY})
    add_library(libR UNKNOWN IMPORTED)
    set_target_properties(libR PROPERTIES
            IMPORTED_LOCATION "${R_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE_DIR}"
            IMPORTED_LINK_INTERFACE_LANGUAGES "C")
endif()
