BUILT_OBJECTS += obj/bedpeSummary.o

obj/bedpeSummary.o: src/bedpeSummary/bedpeSummary.cpp obj/bedpeSummary.d
	$(CXX_COMPILE)
