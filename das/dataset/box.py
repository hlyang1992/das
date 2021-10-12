import numpy as np


class LatticeError(Exception):
    """
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    """

    pass


class Box:
    def __init__(self, matrix):
        try:
            m = np.array(matrix, dtype=np.float64).reshape((3, 3))
        except:
            raise LatticeError("new_lattice must be  length 9 sequence or 3x3 matrix!")
        m.setflags(write=False)
        self._matrix = m

    @property
    def matrix(self):
        return self._matrix

    @property
    def volume(self):
        return abs(np.linalg.det(self._matrix))

    @property
    def lengths(self):
        return np.linalg.norm(self._matrix, axis=1)

    @property
    def angles(self):
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        https://en.wikipedia.org/wiki/Lattice_constant
        """
        m = self._matrix
        angles = np.zeros(3)
        for i, ii in enumerate([[1, 2], [2, 0], [0, 1]]):
            angles[i] = self.vector_angle(*self._matrix[ii])
        angles = np.rad2deg(angles)
        return angles

    @property
    def is_orthorhombic(self):
        """Check that lattice matrix is diagonal."""
        return np.all(self._matrix == np.diag(np.diagonal(self._matrix)))

    @property
    def reciprocal_lattice(self):
        """
        Return the reciprocal lattice with a factor of 2 * pi.
        """
        v = np.linalg.inv(self._matrix).T
        return Box(v * 2 * np.pi)

    @property
    def reciprocal_lattice_crystallographic(self):
        """
        Returns the *crystallographic* reciprocal lattice, i.e., no factor of
        2 * pi.
        """
        return Box(self.reciprocal_lattice.matrix / (2 * np.pi))

    def cart_to_frac(self, cart_coords):
        return cart_coords @ np.linalg.inv(self._matrix)

    def frac_to_cart(self, frac_coords):
        return frac_coords @ self._matrix

    def __repr__(self):
        return (
                f"{self.__class__.__name__}("
                + f'[[{",".join(map(repr, self._matrix[0]))}],'
                + f'[{",".join(map(repr, self._matrix[0]))}],'
                + f'[{",".join(map(repr, self._matrix[0]))}]]'
                + ")"
        )

    def __str__(self):
        outs = [
            "Lattice",
            "    abc : " + " ".join(f"{i:.6f}" for i in self.lengths),
            " angles : " + " ".join(f"{i:.6f}" for i in self.angles),
            f" volume : {self.volume:.6f}",
            "      A : " + " ".join(f"{i:.6f}" for i in self._matrix[0]),
            "      B : " + " ".join(f"{i:.6f}" for i in self._matrix[1]),
            "      C : " + " ".join(f"{i:.6f}" for i in self._matrix[2]),
        ]
        return "\n".join(outs)

    def __eq__(self, other):
        return False if other is None else np.allclose(self._matrix, other.matrix)

    def __ne__(self, other):
        return not self.__eq__(other)

    def copy(self):
        return Box(self._matrix.copy())

    @classmethod
    def cubic(cls, a):
        return cls(np.eye(3) * a)

    @classmethod
    def tetragonal(cls, a, c):
        return cls(np.diag([a, a, c]))

    @classmethod
    def orthorhombic(cls, a, b, c):
        return cls(np.diag([a, b, c]))

    @staticmethod
    def vector_angle(v1, v2):
        return np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2))

    def get_lmp_box(self):
        new_matrix = np.zeros((3, 3))
        a, b, c = self.lengths
        alpha, beta, gamma = np.deg2rad(self.angles)
        a_x = a
        b_x = b * np.cos(gamma)
        b_y = b * np.sin(gamma)
        c_x = c * np.cos(beta)
        c_y = (np.dot(self.matrix[1, :], self.matrix[2, :]) - b_x * c_x) / b_y
        c_z = np.sqrt(c ** 2 - c_x ** 2 - c_y ** 2)
        new_matrix[0] = [a_x, 0, 0]
        new_matrix[1] = [b_x, b_y, 0]
        new_matrix[2] = [c_x, c_y, c_z]
        return new_matrix

    @classmethod
    def from_lmp_box(cls, lx, ly, lz, xy, xz, yz):
        matrix = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])
        return cls(matrix)
