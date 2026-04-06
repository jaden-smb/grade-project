"""
Blender script: import an OBJ sequence and set up per-frame visibility keyframes.

Run this from inside Blender's scripting workspace, or via:
  blender --background --python src/lbm/export/blender_import.py

Edit SEQ_DIR below to point at your obj_sequence output folder.
"""

import bpy
import os

SEQ_DIR = "C:/path/to/output/scenario_a_equilibrium/obj_sequence"

files = sorted(f for f in os.listdir(SEQ_DIR) if f.endswith(".obj"))

for frame_idx, filename in enumerate(files):
    filepath = os.path.join(SEQ_DIR, filename)

    bpy.ops.wm.obj_import(filepath=filepath)
    obj = bpy.context.selected_objects[0]
    obj.name = f"frame_{frame_idx:05d}"

    obj.hide_render = True
    obj.hide_viewport = True
    obj.keyframe_insert("hide_render", frame=0)
    obj.keyframe_insert("hide_viewport", frame=0)

    obj.hide_render = False
    obj.hide_viewport = False
    obj.keyframe_insert("hide_render", frame=frame_idx + 1)
    obj.keyframe_insert("hide_viewport", frame=frame_idx + 1)

    obj.hide_render = True
    obj.hide_viewport = True
    obj.keyframe_insert("hide_render", frame=frame_idx + 2)
    obj.keyframe_insert("hide_viewport", frame=frame_idx + 2)

bpy.context.scene.frame_end = len(files)
