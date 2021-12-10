target_sources(${PROJECT_NAME} PRIVATE Domain.cpp Edge.cpp Element.cpp getMesh.cpp Grid.cpp Material.cpp Mesh.cpp MeshAdapt.cpp MeshExtract.cpp MeshUtil.cpp Node.cpp Partition.cpp saveMesh.cpp Side.cpp)

add_subdirectory(bamg)
add_subdirectory(metis)
