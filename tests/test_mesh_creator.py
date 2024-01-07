import unittest
from mesh_creator import mesh_creation


class TestMeshCreation(unittest.TestCase):

    def setUp(self) -> None:
        None

    def tearDown(self) -> None:
        None

    def test_output(self):
        self.assertEqual(mesh_creation(1, 1), [[0]])
        self.assertEqual(mesh_creation(2, 2), [[0, 1], [0, 1]])
        self.assertEqual(mesh_creation(0, 0), [])

    def test_input(self):
        self.assertRaises(ValueError, mesh_creation, -2, 2)
        self.assertRaises(ValueError, mesh_creation, 2, -2)

    def test_type(self):
        self.assertRaises(TypeError, mesh_creation, [1, 2], -2)
        self.assertRaises(TypeError, mesh_creation, '1', 2)
        self.assertRaises(TypeError, mesh_creation, '112', {2, 3, 4})
        self.assertRaises(TypeError, mesh_creation, True, 2)
        self.assertRaises(TypeError, mesh_creation, 5 + 2j, 2)
