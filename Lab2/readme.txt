smith_waterman.py implements the algorithm in python

bl50.py contains the BLOSUM50 matrix

generator.py generates and runs the test instances described in the report


### How to run the algorithm:


$ python3 smith_waterman.py

computes the optimal alignment of the example strings "HEAGAWGHEE" and "PAWHEAE" with linear gap penalty 8.


$ python3 smith_waterman.py [S1] [S2]

computes the optimal alignment of strings [S1] and [S2] with linear gap penalty 8.


$ python3 smith_waterman.py [S1] [S2] [gap_pen]

computes the optimal alignment of strings [S1] and [S2] with linear gap penalty [gap_pen].


### How to run the tests:

$ python3 generator.py
