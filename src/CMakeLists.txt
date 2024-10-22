# Compile the library
set(SRC_FILES
    libkd/kdtree.c
    libkd/kdtree_dim.c
    libkd/kdtree_mem.c
    libkd/kdtree_fits_io.c
    libkd/kdint_ddd.c
    libkd/kdint_fff.c
    libkd/kdint_ddu.c
    libkd/kdint_duu.c
    libkd/kdint_dds.c
    libkd/kdint_dss.c
    libkd/kdint_lll.c

    solver/quad-utils.c
    solver/solvedfile.c
    solver/solver.c
    solver/tweak2.c
    solver/verify.c

    util/an-endian.c
    util/bl.c
    util/bl-sort.c
    util/codekd.c
    util/dallpeaks.c
    util/datalog.c
    util/dcen3x3.c
    util/dfind.c
    util/dmedsmooth.c
    util/dobjects.c
    util/dpeaks.c
    util/dselip.c
    util/dsigma.c
    util/dsmooth.c
    util/errors.c
    util/fits_io.c
    util/fit-wcs.cpp
    util/healpix.c
    util/image2xy.c
    util/index.c
    util/ioutils.c
    util/log.c
    util/matchobj.c
    util/mathutil.c
    util/permutedsort.c
    util/quadfile.c
    util/resample.c
    util/simplexy.c
    util/sip.c
    util/sip-utils.cpp
    util/starkd.c
    util/starutil.c
    util/starxy.c
    util/tic.c
)

add_library(astrometry-net-lite ${SRC_FILES})

target_include_directories(astrometry-net-lite
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include/
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include/astrometry
        ${PROJECT_SOURCE_DIR}/include/astrometry/solver
        ${PROJECT_SOURCE_DIR}/include/astrometry/libkd
        ${FETCHCONTENT_BASE_DIR}/eigen3-src/
)

target_link_libraries(astrometry-net-lite
    PUBLIC
        cfitsio
)

target_compile_options(astrometry-net-lite PRIVATE "-w")
