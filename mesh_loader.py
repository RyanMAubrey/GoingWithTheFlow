from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np


@dataclass
class TriMesh:
    vertices: np.ndarray
    faces: np.ndarray
    normals: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    texcoords: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))


def load_obj(filepath, scale=1.0, center=False):
    # Parse a wavefront obj file into a triangulated TriMesh
    verts = []
    vnorms = []
    vtexs = []
    faces = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            key = parts[0]

            if key == "v" and len(parts) >= 4:
                verts.append([float(parts[1]), float(parts[2]), float(parts[3])])

            elif key == "vn" and len(parts) >= 4:
                vnorms.append([float(parts[1]), float(parts[2]), float(parts[3])])

            elif key == "vt" and len(parts) >= 3:
                vtexs.append([float(parts[1]), float(parts[2])])

            elif key == "f":
                # Grab vertex indices, handles v, v/vt, v/vt/vn, v//vn
                face_vi = []
                for token in parts[1:]:
                    vi = int(token.split("/")[0]) - 1
                    face_vi.append(vi)

                # Fan triangulate for quads and ngons
                for k in range(1, len(face_vi) - 1):
                    faces.append([face_vi[0], face_vi[k], face_vi[k + 1]])

    vertices = np.array(verts, dtype=np.float64) * scale
    if center:
        vertices -= vertices.mean(axis=0)

    normals = np.array(vnorms, dtype=np.float64) if vnorms else np.empty((0, 3))
    texcoords = np.array(vtexs, dtype=np.float64) if vtexs else np.empty((0, 2))
    faces_arr = np.array(faces, dtype=np.int32)

    return TriMesh(vertices=vertices, faces=faces_arr, normals=normals, texcoords=texcoords)


def load_pose_sequence(filepaths, scale=1.0, center=False):
    # Load multiple obj files that share the same topology
    meshes = [load_obj(fp, scale=scale, center=center) for fp in filepaths]

    # Make sure all frames match
    ref = meshes[0]
    for i, m in enumerate(meshes):
        assert m.vertices.shape == ref.vertices.shape, f"Frame {i} vertex count mismatch"
        assert m.faces.shape == ref.faces.shape, f"Frame {i} face count mismatch"

    return meshes