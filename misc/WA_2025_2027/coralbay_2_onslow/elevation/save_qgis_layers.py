# Script to save file paths of all layers in a group in a QGIS session.
#
# Run from QGIS by
#   going to Plugins -> Python Console.
#   show editor, open this file and click run.
#
# Used to take elevation layers visible in QGIS and save them in order for SWALS.

group_name = "Raw Elevation"
out_file = "./swals_elevation_files_in_preference_order.txt"

root = QgsProject.instance().layerTreeRoot()
group = root.findGroup(group_name)
if group is None:
    raise ValueError(f"Group '{group_name}' not found in project.")

def write_layers_in_order(node, file_obj):
    for child in node.children():
        if isinstance(child, QgsLayerTreeLayer):
            layer = child.layer()
            file_obj.write(f"{layer.source()}\n")
        elif isinstance(child, QgsLayerTreeGroup):
            write_layers_in_order(child, file_obj)

with open(out_file, "w", encoding="utf-8") as f:
    write_layers_in_order(group, f)

print(f"Saved layer sources in exact panel order to {out_file}")
