#!/usr/bin/env python3
"""
PyMOL을 이용한 WT vs Humanized nanobody 구조 시각화

논문 스타일로 nanobody 구조를 렌더링:
  - Cartoon representation with FR/CDR region coloring
  - Mutation sites highlighted as sticks
  - WT와 Humanized를 superimpose하여 비교

Input:
    - ../01_structure_prediction/wt_nanobody.pdb
    - ../01_structure_prediction/humanized_nanobody.pdb
Output:
    - figures/04_structure_wt.png/pdf           : WT 단독
    - figures/04_structure_humanized.png/pdf     : Humanized 단독
    - figures/04_structure_overlay.png/pdf       : 두 구조 overlay

Usage:
    python 04_pymol_structure.py
"""

import os
import pymol
from pymol import cmd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures")
STRUCT_DIR = os.path.join(SCRIPT_DIR, "..", "01_structure_prediction")
os.makedirs(FIGURES_DIR, exist_ok=True)

WT_PDB = os.path.join(STRUCT_DIR, "wt_nanobody.pdb")
HUM_PDB = os.path.join(STRUCT_DIR, "humanized_nanobody.pdb")

# Mutation positions (PDB sequential numbering)
MUT_POSITIONS = [5, 11, 35, 50, 74, 86, 92, 120]

# Region definitions (PDB sequential)
REGIONS = {
    'FR1':  (1, 25),
    'CDR1': (26, 33),
    'FR2':  (34, 50),
    'CDR2': (51, 57),
    'FR3':  (58, 95),
    'CDR3': (96, 114),
    'FR4':  (115, 125),
}

REGION_COLORS = {
    'FR1':  [0.75, 0.78, 0.85],   # light steel blue
    'CDR1': [0.90, 0.30, 0.30],   # red
    'FR2':  [0.70, 0.85, 0.82],   # light teal
    'CDR2': [0.95, 0.75, 0.20],   # amber
    'FR3':  [0.82, 0.75, 0.88],   # light purple
    'CDR3': [0.30, 0.60, 0.90],   # blue
    'FR4':  [0.72, 0.88, 0.72],   # light green
}


def setup_scene():
    """Common scene settings for publication quality."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_shadow_decay_factor", 0.1)
    cmd.set("ray_shadow_decay_range", 2.0)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("ray_trace_mode", 0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_oval_length", 1.2)
    cmd.set("cartoon_oval_width", 0.25)
    cmd.set("stick_radius", 0.15)
    cmd.set("sphere_scale", 0.2)
    cmd.set("spec_reflect", 0.4)
    cmd.set("spec_power", 200)
    cmd.set("ambient", 0.3)
    cmd.set("direct", 0.6)
    cmd.set("reflect", 0.5)
    cmd.set("light_count", 2)
    cmd.set("ray_trace_gain", 0.1)


def color_regions(obj_name):
    """Color FR/CDR regions."""
    for region, (start, end) in REGIONS.items():
        sel_name = f"{obj_name}_{region}"
        cmd.select(sel_name, f"{obj_name} and resi {start}-{end}")
        cmd.set_color(f"col_{region}", REGION_COLORS[region])
        cmd.color(f"col_{region}", sel_name)
        cmd.delete(sel_name)


def show_mutations(obj_name, label_color=[0.8, 0.0, 0.0]):
    """Show mutation sites as sticks (disabled by default)."""
    pass


def set_cdr_up_view(obj_name):
    """Orient using default orient (no rotation)."""
    cmd.orient(obj_name)
    cmd.zoom(obj_name, 3)


def render_single(obj_name, out_prefix):
    """Render a single structure."""
    cmd.hide("everything")
    cmd.show("cartoon", obj_name)
    color_regions(obj_name)
    show_mutations(obj_name)
    set_cdr_up_view(obj_name)

    # Render
    cmd.ray(2400, 2400)
    cmd.png(os.path.join(FIGURES_DIR, f"{out_prefix}.png"), dpi=300)
    print(f"Saved: figures/{out_prefix}.png")


def main():
    pymol.finish_launching(['pymol', '-cq'])
    setup_scene()

    # Load structures
    cmd.load(WT_PDB, "wt")
    cmd.load(HUM_PDB, "humanized")

    # Align humanized to WT
    cmd.align("humanized", "wt")

    # --- WT single ---
    cmd.disable("humanized")
    cmd.enable("wt")
    render_single("wt", "04_structure_wt")

    # --- Humanized single ---
    cmd.disable("wt")
    cmd.enable("humanized")
    render_single("humanized", "04_structure_humanized")

    # --- Overlay ---
    cmd.enable("wt")
    cmd.enable("humanized")
    cmd.hide("everything")

    # WT: gray cartoon
    cmd.show("cartoon", "wt")
    cmd.color("gray70", "wt")
    cmd.set("cartoon_transparency", 0.5, "wt")

    # Humanized: colored cartoon
    cmd.show("cartoon", "humanized")
    color_regions("humanized")
    cmd.set("cartoon_transparency", 0.0, "humanized")

    # Show mutations on humanized
    show_mutations("humanized")

    set_cdr_up_view("humanized")

    cmd.ray(2400, 2400)
    cmd.png(os.path.join(FIGURES_DIR, "04_structure_overlay.png"), dpi=300)
    print("Saved: figures/04_structure_overlay.png")

    cmd.quit()


if __name__ == "__main__":
    main()
