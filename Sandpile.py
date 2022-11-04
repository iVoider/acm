import numpy as np

class Sandpile():
    def __init__(self, arr = None ,rows = 3, cols = 3, max_sand = 10):
        """
        arr - 2d array of values
        rows - height of sandpile
        cols - width of sandpile
        max_sand - max count of sandpile grains (must be div by 4)
        """
        self.max_grains = max_sand
        if arr == None:
            self.rows = rows
            self.cols = cols
            self.grid = np.zeros((cols,rows), int)
        else:
            self.rows = len(arr)
            self.cols = len(arr[0])
            self.grid = np.array(arr)


    def topple(self, elem_x, elem_y):
        try:
         p = self.grid[elem_x, elem_y]
         b = p // self.max_grains
         o = p % self.max_grains
         self.grid[elem_x, elem_y] = o

         # increase height of neighbor piles
         self.grid[elem_x-1, elem_y] += b
         self.grid[elem_x+1, elem_y] += b
         self.grid[elem_x, elem_y-1] += b
         self.grid[elem_x, elem_y+1] += b


         self.grid[0] = self.grid[-1] = 0
         self.grid[:, 0] = self.grid[:, -1] = 0
        except:
         return False


    def run(self):

        iterations = 0
        topple = self.topple
        where = np.where

        while np.max(self.grid) >= self.max_grains:
            elem_x, elem_y = where(self.grid >= self.max_grains)
            if topple(elem_x, elem_y) == False:
                return self.grid
            iterations += 1

        return self.grid

    def get_pile(self):
        return self.grid

    def set_sand(self, x, y, number):
        self.grid[x,y] = number

    def __add__(self, other):
        result = Sandpile(rows = self.rows, cols = self.cols)
        try:
            result.grid = self.grid + other.grid
            return result.run()
        except ValueError:
            print("ValueError: sandpile grid sizes must match")