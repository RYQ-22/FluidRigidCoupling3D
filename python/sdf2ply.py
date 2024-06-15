import numpy as np
from skimage import measure
import argparse
import open3d as o3d

def marching_cubes(sdf, voxel_size, threshold=0.0):
    verts, faces, normals, values = measure.marching_cubes(sdf, level=threshold, spacing=(voxel_size, voxel_size, voxel_size))
    return verts, faces, normals

def main():
    parser = argparse.ArgumentParser(description='Process some arguments.')
    parser.add_argument('-frame', type=int, help='Frame number')
    parser.add_argument('-project_path', type=str, help='project path')
    parser.add_argument('-l', type=float, help = 'voxel size')
    args = parser.parse_args()
    frame_number = args.frame
    project_path = args.project_path
    l = args.l
    # read from .npy file
    sdf = np.load(project_path + f"/python/sdf/sdf_{frame_number}.npy")

    # sdf to mesh
    vertices, faces, normals = marching_cubes(sdf, l)
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(vertices)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    #mesh.vertex_normals = o3d.utility.Vector3dVector(normals)
    #mesh.filter_smooth_taubin(number_of_iterations=100, lambda_filter=0.9, mu=-0.8)
    mesh = mesh.filter_smooth_simple(number_of_iterations=6)
    mesh = mesh.filter_smooth_laplacian(number_of_iterations=2)
    mesh = mesh.subdivide_loop(number_of_iterations=2)
    o3d.io.write_triangle_mesh(project_path + f"/python/ply/liquid_{frame_number}.ply", mesh)

if __name__ == "__main__":
    main()
    