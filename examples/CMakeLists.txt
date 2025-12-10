file(GLOB EXAMPLE_SOURCES *.cpp)

foreach(example ${EXAMPLE_SOURCES})
    get_filename_component(example_name ${example} NAME_WE)
    add_executable(${example_name} ${example})
	target_link_libraries(${example_name}
		moebius
		${GLUI_LIBRARIES}
		OpenGL::GL
        OpenGL::GLU
		GLUT::GLUT
	)    

endforeach()
