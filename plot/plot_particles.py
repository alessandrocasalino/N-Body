import argparse
import pandas as pd
import matplotlib.pyplot as plt

# Create the parser
parser = argparse.ArgumentParser(description='Plot data from a CSV file.')
parser.add_argument('csv_file', type=str, help='The path to the CSV file')

# Parse the command line arguments
args = parser.parse_args()

# Load the data
data = pd.read_csv(args.csv_file, comment='#')

print(data)

# Create a scatter plot
plt.figure(figsize=(10, 10))
if 'z' in data.columns:
    # 3D plot
    ax = plt.axes(projection='3d')
    ax.scatter(data['x'], data['y'], data['z'])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
else:
    # 2D plot
    plt.scatter(data['x'], data['y'])
    plt.xlabel('X')
    plt.ylabel('Y')

# Save the plot as a PNG file
plt.savefig('plot.png')
