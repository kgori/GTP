cmake_minimum_required(VERSION 2.8.4)
project(cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    BipartiteGraph.h
    Bipartition.h
    EdgeAttribute.h
    Geodesic.h
    PhyloTree.h
    PhyloTreeEdge.h
    PolyMain.h
    Ratio.h
    RatioSequence.h
    TreeDistance.h)

add_executable(cpp ${SOURCE_FILES})