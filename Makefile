CXX = g++ -O2 -Wall -mfma -Wa,-q -I ./

eights: eight-4th eight-6th eight-8th eight-10th eight-12th eight-14th eight-16th

keplers: kepler-4th kepler-6th kepler-8th kepler-10th kepler-12th kepler-14th kepler-16th

eight-4th: eight.cpp nbodysystem.h hermite4.h
	$(CeightXX) -DHERMITE_FOURTH $< -o $@

eight-6th: eight.cpp nbodysystem.h hermite6.h
	$(CXX) -DHERMITE_SIXTH $< -o $@

eight-8th: eight.cpp nbodysystem.h hermite8.h
	$(CXX) -DHERMITE_EIGTH $< -o $@

eight-10th: eight.cpp nbodysystem.h hermite10.h
	$(CXX) -DHERMITE_TENTH $< -o $@

eight-12th: eight.cpp nbodysystem.h hermite12.h
	$(CXX) -DHERMITE_TWELVETH $< -o $@

eight-14th: eight.cpp nbodysystem.h hermite14.h
	$(CXX) -DHERMITE_FOURTEENTH $< -o $@

eight-16th: eight.cpp nbodysystem.h hermite16.h
	$(CXX) -DHERMITE_SIXTEENTH $< -o $@

kepler-4th: kepler.cpp nbodysystem.h hermite4.h
	$(CXX) -DHERMITE_FOURTH $< -o $@

kepler-6th: kepler.cpp nbodysystem.h hermite6.h
	$(CXX) -DHERMITE_SIXTH $< -o $@

kepler-8th: kepler.cpp nbodysystem.h hermite8.h
	$(CXX) -DHERMITE_EIGTH $< -o $@

kepler-10th: kepler.cpp nbodysystem.h hermite10.h
	$(CXX) -DHERMITE_TENTH $< -o $@

kepler-12th: kepler.cpp nbodysystem.h hermite12.h
	$(CXX) -DHERMITE_TWELVETH $< -o $@

kepler-14th: kepler.cpp nbodysystem.h hermite14.h
	$(CXX) -DHERMITE_FOURTEENTH $< -o $@

kepler-16th: kepler.cpp nbodysystem.h hermite16.h
	$(CXX) -DHERMITE_SIXTEENTH $< -o $@

