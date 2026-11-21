"""Expected logical hierarchies, shared by layout and VTK readback checks."""
main_tree = {
    "Regions": {"Fluid_Name_2": {}, "Fluid_Name": 0, "Fluid_Name_1": 1,
                "Part": {"Common": {}, "Fluid_Name": 0}},
    "Regions_1": {"Fluid_Name": 0, "Common": {}},
    "Fluid_Name": 0,
    "Fluid_Name_1": {},  # A remote block name is available to this root node.
    "Other_1": {},
    "Other_2": {},
    "Other": 2,  # Unplaced blocks follow all explicit children.
    "Tail": 3,
}
trees = {
    "explicit": main_tree,
    "implicit": main_tree,
    "flat": {"A": 0, "B": 1},
    "empty": {"Empty": {"Unnamed_Node": {}}, "Empty_1": {}},
    "rootrefs": {"B": 1, "A": 0},
    "early": {"A": {"A": 0}, "B_1": {}, "B": 1},
    "late": {"A": {"A": 0}, "B_1": {}, "B": 1},
    "growth": {**{f"Node_{i:03}": {"Inside": ({"A": 0} if i == 1 else
                                                {"B": 1} if i == 96 else {})}
                  for i in range(1, 97)}, "A": 0},
    "rootgrowth": {"B_1": {"A": 0}, "C_1": {}, "D_1": {}, "E_1": {},
                   "B": 1, "C": 2, "D": 3, "E": 4},
}
