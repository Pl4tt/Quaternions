from __future__ import annotations

import math

from vector import Vector


class Quaternion:
    def __init__(self, scalar: int, vector: Vector) -> None:
        self.scalar = scalar
        self.vector = vector
    
    def __neg__(self):
        return Quaternion(-self.scalar, -self.vector)

    def __add__(self, other):
        return Quaternion(self.scalar+other.scalar, self.vector+other.vector)

    def __sub__(self, other):
        return self.__add__(-other)

    def __mul__(self, other):
        match other:
            case float() | int():
                return Quaternion(other*self.scalar, other*self.vector)
            
            case Vector():
                return Vector([self*i for i in other])

            case Quaternion():
                scalar = self.scalar*other.scalar - self.vector@other.vector
                vector = Vector([
                    self.scalar*other.vector[0] + self.vector[0]*other.scalar + self.vector[1]*other.vector[2] - self.vector[2]*other.vector[1],
                    self.scalar*other.vector[1] - self.vector[0]*other.vector[2] + self.vector[1]*other.scalar + self.vector[2]*other.vector[0],
                    self.scalar*other.vector[2] + self.vector[0]*other.vector[1] - self.vector[1]*other.vector[0] + self.vector[2]*other.scalar,
                ])
                return Quaternion(scalar, vector)
    
    __rmul__ = __mul__

    def __pow__(self, n):
        phi = math.acos(self.scalar/abs(self))
        return round(abs(self)**n*Quaternion(math.cos(n*phi), self.vector.normalized()*math.sin(n*phi)), 10)

    def __truediv__(self, other):
        match other:
            case int() | float():
                return (1/other)*self

            case Quaternion():
                return self*other.inverse()

    def __eq__(self, other):
        return self.scalar == other.scalar and self.vector == other.vector

    def __abs__(self):
        return math.sqrt(self.scalar**2 + self.vector[0]**2 + self.vector[1]**2 + self.vector[2]**2)
    
    def __round__(self, n):
        return Quaternion(round(self.scalar, n), round(self.vector, n))

    def __repr__(self):
        a, b, c, d = self.scalar, *self.vector

        return f"{a if a != 0 else '0'}{'+' if b >= 0 else ''}{b if b != 0 else '0'}i{'+' if c >= 0 else ''}{c if c != 0 else '0'}j{'+' if d >= 0 else ''}{d if d != 0 else '0'}k"

    def conjugate(self):
        return Quaternion(self.scalar, -self.vector)
    
    def inverse(self):
        return round(self.conjugate()*abs(self)**(-2), 10)
    

    def complex_matrix_repr(self):
        """this quaternion written in (complex) matrix notation"""
        return [
            [self.scalar+self.vector[0]*1j, self.vector[1]+self.vector[2]*1j],
            [-self.vector[1]+self.vector[2]*1j, self.scalar-self.vector[0]*1j],
        ]

    def real_matrix_repr(self):
        """this quaternion written in (real) matrix notation"""
        return [
            [self.scalar] + list(self.conjugate().vector),
            [self.vector[0], self.scalar, -self.vector[2], self.vector[1]],
            [self.vector[1], self.vector[2], self.scalar, -self.vector[0]],
            [self.vector[2], -self.vector[1], self.vector[0], self.scalar],
        ]
    
    def rotation_matrix(self):
        """
        returns a matrix M for which:

        let p be a 3-dim vector you want to rotate, let q be this quaternion and q' by q.conjugate():
        q*p*q' = Mp = Vector p rotated by rotation of this quaternion.
        """
        a, b, c, d = self.scalar, *self.vector

        return [
            [1-2*(c**2+d**2), 2*(b*c-a*d), 2*(b*d+a*c)],
            [2*(b*c+a*d), 1-2*(b**2+d**2), 2*(c*d-a*b)],
            [2*(b*d-a*c), 2*(c*d+a*b), 1-2*(b**2+c**2)],
        ]
    
    def rotate_vector(self, point_vector: Vector) -> Vector:
        """rotates a point vector by the rotation represented by the quaternion"""
        hamilton = self*Quaternion(0, point_vector)*self.conjugate()

        return hamilton.vector
    
    def inverse_rotate_vector(self, point_vector: Vector) -> Vector:
        """inverts the process of self.rotate_vector()"""
        return (self.conjugate()*Quaternion(0, point_vector)*self).vector
    
    def rotation_details(self) -> tuple[int, Vector]:
        """returns (angle, vector) of rotation of quaternion (angle in radians)"""
        angle = 2*math.acos(self.scalar)
        vector = (1/math.sin(angle/2))*self.vector

        return round(angle, 5), round(vector, 5)
    
    def normalized(self):
        return self/abs(self)
    
    def polar(self):
        """returns sane quaternion using polar form (only useful for testing)"""
        phi = math.acos(self.scalar/abs(self))
        return round(Quaternion(abs(self)*math.cos(phi), abs(self)*self.vector.normalized()*math.sin(phi)), 5)
    
    @staticmethod
    def rotation_quaternion(vector: Vector, angle: float) -> Quaternion:
        """generates a quaternion that represents a rotation around the given vector of the given angle"""
        # angle in radians
        return Quaternion(math.cos(angle/2), math.sin(angle/2)*vector)
    
    @staticmethod
    def exp(quaternion: Quaternion):
        return round(math.exp(quaternion.scalar)*Quaternion(math.cos(abs(quaternion.vector)), quaternion.vector.normalized()*math.sin(abs(quaternion.vector))), 5)
    
    @staticmethod
    def log(quaternion: Quaternion):
        return round(Quaternion(math.log(abs(quaternion)), quaternion.vector.normalized()*math.acos(quaternion.scalar/abs(quaternion))), 5)


def test():
    q1 = Quaternion(5, Vector([2, 6, 9]))
    q2 = Quaternion(0, Vector([10, -4, -20]))
    q3 = Quaternion(1, Vector([0, 1, 0])).normalized()
    quatlog3 = Quaternion.log(q3)

    assert 2*q1 == q1*2
    assert q1*q2 != q2*q1
    assert str(q1-q2) == "5-8i+10j+29k"
    assert q1+q2 == Quaternion(5, Vector([12, 2, -11]))
    assert q1*q2 == Quaternion(184, Vector([-34, 110, -168]))
    assert str(q1/q2) == "-0.3565891471+0.06589147239999997i-0.2131782948j+0.3255813955k"
    assert q1.conjugate() == Quaternion(5, Vector([-2, -6, -9]))
    assert q1.complex_matrix_repr() == [[(5+2j), (6+9j)], [(-6+9j), (5-2j)]]
    assert q2.complex_matrix_repr() == [[10j, (-4-20j)], [(4-20j), -10j]]
    assert q1.real_matrix_repr() == [[5, -2, -6, -9], [2, 5, -9, 6], [6, 9, 5, -2], [9, -6, 2, 5]]
    assert q2.real_matrix_repr() == [[0, -10, 4, 20], [10, 0, 20, -4], [-4, -20, 0, -10], [-20, 4, 10, 0]]
    assert q2**2 == q2*q2 == Quaternion(-516, Vector([0, 0, 0]))
    assert abs(q2) == 22.715633383201094
    assert str(q2**(-1)) == str(q2.inverse()) == "0-0.019379845i+0.007751938j+0.0387596899k"
    assert q1**(-1) == q1.inverse() == Quaternion(0.0342465753, Vector([-0.0136986301, -0.0410958904, -0.0616438356]))
    assert str(round(q1**(-2), 5)) == str(round((q1**2).inverse(), 5)) == "-0.0045-0.00094i-0.00281j-0.00422k"
    assert abs(q1.vector.normalized()) == 1
    assert Quaternion.log(q3) == quatlog3
    assert Quaternion.exp(quatlog3) == round(q3, 5)
    assert (q1, q2) == (q1.polar(), q2.polar())

    print("ALGEBRAIC TESTS PASSED")

    rotation_q = Quaternion.rotation_quaternion(Vector([0.6, 0.58, 0.55]), math.pi/2)
    rotated_v = rotation_q.rotate_vector(Vector([-3, 5, 1]))

    assert rotation_q == Quaternion(0.7071067811865476, Vector([0.4242640687119285, 0.41012193308819755, 0.3889087296526012]))
    assert rotated_v == Vector([-1.1816500000000003, -1.2902500000000003, 5.6480500000000005])
    assert Vector([-3, 5, 1]) == round(rotation_q.inverse_rotate_vector(Vector([-1.18, -1.3, 5.65])), 1)
    assert rotation_q.rotation_details() == (round(math.pi/2, 5), Vector([0.6, 0.58, 0.55]))
    assert rotation_q.rotation_matrix() == [
        [0.3611, -0.20200000000000018, 0.9100000000000001],
        [0.8980000000000001, 0.3374999999999999, -0.281],
        [-0.25, 0.9190000000000002, 0.3036000000000001]
    ]

    rotation_q2 = Quaternion.rotation_quaternion(Vector([0.13, -0.23, 0.96]), 1)
    rotated_v2 = rotation_q2.rotate_vector(Vector([10, -0.5, 0]))
    
    assert rotation_q2 == Quaternion(0.8775825618903728, Vector([0.06232532001854639, -0.1102678738789667, 0.4602485170600349]))
    assert rotated_v2 == Vector([5.911257523072413, 7.657373336624139, 2.5051409987543676])
    assert rotation_q2.rotation_details() == (1, Vector([0.13, -0.23, 0.96]))
    assert rotation_q2.rotation_matrix() == [
        [0.5520245970685022, -0.8215571064701233, -0.13616805427816003],
        [0.7940671843610381, 0.5685737140572491, -0.21089247888934132],
        [0.2509085987334724, 0.007889977160711797, 0.9679131009495962]
    ]

    print("ROTATION TESTS PASSED")
    
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    test()
