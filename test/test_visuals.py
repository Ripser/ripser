import pytest
import numpy as np
import matplotlib.pyplot as plt

from ripser import Rips


"""

    Testing visualization is a little more difficult, but still necessary. An example of how to get started:
    > https://stackoverflow.com/questions/27948126/how-can-i-write-unit-tests-against-code-that-uses-matplotlib

    ```
    def test_plot_square2():
        f, ax = plt.subplots()
        x, y = [0, 1, 2], [0, 1, 2]
        plot_square(x, y)
        x_plot, y_plot = ax.lines[0].get_xydata().T
        np.testing.assert_array_equal(y_plot, np.square(y))

    ```


    Notes
    -----

    ax.get_children() gives all the pieces in the plot, very useful
    Scatter data is of type `PathCollection` and will be a child of ax.

"""


class TestPlotting:
    def test_single(self):
        """ Most just test this doesn't crash
        """
        rips = Rips()

        diagram = np.array([[0,1], [1,1],[2,4], [3,5]])
        
        f, ax = plt.subplots()
        rips.plot(diagram, show=False)

        x_plot, y_plot = ax.lines[0].get_xydata().T

        assert x_plot[0] <= np.min(diagram)
        assert x_plot[1] >= np.max(diagram)

        # get PathCollection
        pathcols = [child for child in ax.get_children() if child.__class__.__name__ == "PathCollection"]
        assert len(pathcols) == 1
        
    def test_multiple(self):

        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, show=False)
        
       
        pathcols = [child for child in ax.get_children() if child.__class__.__name__ == "PathCollection"]

        assert len(pathcols) == 2
        np.testing.assert_array_equal(pathcols[0].get_offsets(), diagrams[0])
        np.testing.assert_array_equal(pathcols[1].get_offsets(), diagrams[1])
        
    def test_legend_true(self):
        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, legend=True, show=False)
        legend = [child for child in ax.get_children() if child.__class__.__name__ == "Legend"]

        assert len(legend) == 1

        
    def test_legend_false(self):
        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, legend=False, show=False)
        legend = [child for child in ax.get_children() if child.__class__.__name__ == "Legend"]

        assert len(legend) == 0



    def test_set_title(self):

        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, title='my title', show=False)
        assert ax.get_title() == 'my title'

        f, ax = plt.subplots()
        rips.plot(diagrams, show=False)
        assert ax.get_title() == ''


    def test_default_square(self):
        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, show=False)
        diagonal = ax.lines[0].get_xydata()

        assert diagonal[0,0] == diagonal[0,1]
        assert diagonal[1,0] == diagonal[1,1]

    def test_default_label(self):
        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, show=False)
        
        assert ax.get_ylabel() == 'Death'
        assert ax.get_xlabel() == 'Birth'

    def test_lifetime(self):
        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, lifetime=True, show=False)
        
        assert ax.get_ylabel() == 'Lifetime'
        assert ax.get_xlabel() == 'Birth'

        line = ax.get_lines()[0]
        np.testing.assert_array_equal(line.get_ydata(), [0,0])

    def test_lifetime_removes_birth(self):
        rips = Rips()

        diagrams = [
            np.array([[0,1], [1,1],[2,4], [3,5]]),
            np.array([[0.5,3], [2,4],[4,5], [10,15]])
        ]

        f, ax = plt.subplots()
        rips.plot(diagrams, lifetime=True, show=False)
        

        pathcols = [child for child in ax.get_children() if child.__class__.__name__ == "PathCollection"]

        modded1 = diagrams[0]
        modded1[:,1] = diagrams[0][:,1] - diagrams[0][:,0]
        modded2 = diagrams[1]
        modded2[:,1] = diagrams[1][:,1] - diagrams[1][:,0]
        assert len(pathcols) == 2
        np.testing.assert_array_equal(pathcols[0].get_offsets(), modded1)
        np.testing.assert_array_equal(pathcols[1].get_offsets(), modded2)
        