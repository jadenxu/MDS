SET(MDS_CORE_LIBS
    my_node.h
    my_graph.h
    my_graph.cpp
    
    pcoa.h
    pcoa.cpp
    
    smacof.h
    smacof.cpp

    NLDR.h
    NLDR.cpp)
ADD_LIBRARY(mdsCore ${MDS_CORE_LIBS})
TARGET_LINK_LIBRARIES(mdsCore ${CORELIBS})
