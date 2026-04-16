import argparse
import xml.etree.ElementTree as ET
import pandas as pd

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Extract coordinates from an XML file and save to CSV.")
parser.add_argument("-i", "--input", required=True, help="Input XML file.")
parser.add_argument("-o", "--output", required=True, help="Output CSV file.")

# Parse the arguments
args = parser.parse_args()

# Load and parse the XML file
tree = ET.parse(args.input)
root = tree.getroot()

# Extract the coordinates
x_vals = []
y_vals = []
for coordinate in root.findall("coordinate"):
    x_vals.append(float(coordinate.find("x").text))
    y_vals.append(float(coordinate.find("y").text))

# Save to the specified CSV file
data = pd.DataFrame({"x": x_vals, "y": y_vals})
data.to_csv(args.output, index=False)
print(f"Data saved to {args.output}")
