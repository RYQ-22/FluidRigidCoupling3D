import os
import argparse
import subprocess

def ply_to_png(ply_path, rigidBody_path, png_path, scene_path, pbrt_folder):
    pbrt_text = f"""
    Transform [ -0.9999538 0.0036316987 -0.008376651 0 0.004065331 0.99861884 -0.05233363 0 0.008175063 -0.052366752 -0.9985882 0 31.001501 -37.19425 235.93471 1 ]
    Camera "perspective"
        "float fov" [ 50 ]
    Film "rgb"
        "string filename" [ "dambreak0.exr" ]
        "float maxcomponentvalue" [ 10 ]
        "integer yresolution" [ 1080 ]
        "integer xresolution" [ 1920 ]
        "string sensor" "canon_eos_100d"
        "float iso" 150
    Sampler "halton"
        "integer pixelsamples" [ 128 ]
    Integrator "volpath"
        "integer maxdepth" [ 64 ]

    WorldBegin

    MakeNamedMaterial "WaterAir"
        "string type" [ "dielectric" ]
        "float eta" [ 1.33 ]
        "bool remaproughness" [ false ] 

    MakeNamedMaterial "Solid"
        "string type" [ "diffuse" ]    	    	
    	"rgb reflectance" [ 0.527 0.805 0.918 ]

    AttributeBegin
        Rotate 60 0 1 0
        LightSource "infinite"
            "float scale" [3]
            "string filename" "textures/sky.exr"
    AttributeEnd

    AttributeBegin
        # boundary 1:1:2        
        Material "dielectric"
        "spectrum eta" "glass-BAF10"        
        Translate -40 1 100
        Scale 78 78 78 
        Rotate 90 0 1 0 
        AttributeBegin	        
            Shape "plymesh"
                "string filename" [ "geometry/boundary11.ply" ]
        AttributeEnd    
        AttributeBegin	        
            Shape "plymesh"
                "string filename" [ "geometry/boundary22.ply" ]
        AttributeEnd    
        AttributeBegin	        
            Shape "plymesh"
                "string filename" [ "geometry/boundary33.ply" ]
        AttributeEnd    
        AttributeBegin	        
            Shape "plymesh"
                "string filename" [ "geometry/boundary44.ply" ]
        AttributeEnd    
        AttributeBegin	        
            Shape "plymesh"
                "string filename" [ "geometry/boundary55.ply" ]
        AttributeEnd    
        AttributeBegin	        
            Shape "plymesh"
                "string filename" [ "geometry/boundary66.ply" ]
        AttributeEnd    
    AttributeEnd

    AttributeBegin
        Material "dielectric"
        "spectrum eta" "glass-F11"

        Translate -42.5 -4.5 103
        Scale 161 161.37 167.14 
        Rotate 90 0 1 0

        AttributeBegin
            NamedMaterial "WaterAir"
            Shape "plymesh"
                "string filename" [ "{ply_path}" ]
        AttributeEnd
    AttributeEnd

    AttributeBegin    

        Translate -42.5 -4.5 103
        Scale 161 161.37 167.14 
        Rotate 90 0 1 0

        AttributeBegin
            NamedMaterial "Solid"
            Shape "plymesh"
                "string filename" [ "{rigidBody_path}" ]
        AttributeEnd
    AttributeEnd

    AttributeBegin      
        Translate 0 -0.04 0        
        
        AttributeBegin
            Shape "plymesh"
                "string filename" [ "geometry/mesh_00008.ply" ]
        AttributeEnd
    AttributeEnd
    """
    with open(scene_path, "w") as scene_file:
        scene_file.write(pbrt_text)
    command_render = pbrt_folder + "/pbrt --gpu --outfile " + png_path + " " + scene_path 
    command_denoise = pbrt_folder + "/imgtool denoise-optix --outfile " + png_path + " " + png_path
    try:
        subprocess.run(command_render, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        subprocess.run(command_denoise, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        print("Error Output:")
        print(e.stderr)        

def main():
    parser = argparse.ArgumentParser(description='Process some arguments.')
    parser.add_argument('-project_path', type=str, help = 'project folder')
    parser.add_argument('-pbrt_folder', type=str, help='pbrt folder')
    parser.add_argument('-frame', type=int, help='total frame number')
    args = parser.parse_args()
    pbrt_folder = args.pbrt_folder
    frame = args.frame
    project_path = args.project_path

    input_folder = project_path + "/python/ply"
    output_folder = project_path + "/python/png"
    scene_folder = project_path + "/python/scene"
    for i in range(frame):
        ply_file_path = input_folder + f"/liquid_{i}.ply"
        png_file_path = output_folder + f"/liquid_{i}.png"
        scene_file_path = scene_folder + f"/scene_{i}.pbrt"
        rigidBody_path = input_folder + f"/rigidbody_{i}.ply"
        ply_to_png(ply_file_path, rigidBody_path, png_file_path, scene_file_path, pbrt_folder)
        

if __name__ == "__main__":
    main()