CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

profile: profile.cpp
	$(CXX) $(CXXFLAGS) -o profile profile.cpp

test: profile
	./profile --test

plot: profile
	./profile --csv --distance 400 --vmax 300 --amax 1200 --jmax 9000 > move.csv
	python3 plot_profile.py move.csv -o move.png --vmax 300 --amax 1200 --jmax 9000

clean:
	rm -f profile move.csv move.png
