import math


class Vector:
    def __init__(self, lst: list) -> None:
        self.lst = lst

    def __neg__(self):
        return Vector([-i for i in self.lst])

    def __add__(self, other):
        return Vector([a+b for a, b in zip(self.lst, other.lst)])

    def __mul__(self, scalar):
        return Vector([scalar*i for i in self.lst])

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return (1/scalar)*self
    
    def __matmul__(self, other):
        return sum([a*b for a, b in zip(self.lst, other.lst)])

    def __eq__(self, other):
        return self.lst == other.lst
    
    def __getitem__(self, key):
        return self.lst[key]
    
    def __iter__(self):
        for i in self.lst:
            yield i
        
    def __abs__(self):
        return math.sqrt(sum([a**2 for a in self.lst]))
    
    def __round__(self, n):
        return Vector([round(i, n) for i in self.lst])
    
    def __repr__(self):
        return str(tuple(self.lst))
    
    __str__ = __repr__

    def normalized(self):
        return (1/abs(self))*self
