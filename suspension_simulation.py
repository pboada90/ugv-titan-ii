import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation
import ezdxf

# Lengths (adjustable according to design)
l_a = 30.0      # length of the upper arm
l_b = 51.8      # length of the lower arm
l_c = 63.8      # length of the load-bearing arm
t_w = 15.0      # wheel width (rectangle)
t_h = 30.0      # wheel height (rectangle)
offset = 8.0    # Additional offset to move the wheel forward

# Anchor points on the chassis
x1, y1 = 4.0, 50.0
x2, y2 = 0.0, 0.0

# Define initial angles (in radians)
theta_b_0 = np.radians(-20)  # The lower arm also starts 20° below
theta_b = theta_b_0
angle_range = np.radians(20)

# Function to get the coordinates of the points
def get_positions(theta_b):    
    # Coordinates of the end of the lower arm
    x_b = x2 + l_b * np.cos(theta_b)
    y_b = y2 + l_b * np.sin(theta_b)

    m = -(x1 - x_b) / (y1 - y_b)
    b = (l_c**2 - l_a**2 - x_b**2 + x1**2 - y_b**2 + y1**2) / (2 * (y1 - y_b))

    a = 1 + m**2
    b_2 = 2 * (m * (b - y_b) - x_b)
    c = x_b**2 + (b - y_b)**2 - l_c**2

    root_1 = (-b_2 + np.sqrt(b_2**2 - 4 * a * c)) / (2 * a)
    root_2 = (-b_2 - np.sqrt(b_2**2 - 4 * a * c)) / (2 * a)

    root = max([root_1, root_2])

    x_a = root
    y_a = m * root + b

    theta_a = np.arctan2(y_a - y1, x_a - x1)
    
    return (x_a, y_a), (x_b, y_b), theta_a

from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge

# Updated function to load DXF points and organize them
def load_dxf_points(dxf_path):
    doc = ezdxf.readfile(dxf_path)
    msp = doc.modelspace()
    line_segments = []
    for entity in msp:
        if entity.dxftype() == 'LINE':
            start_point = (entity.dxf.start.x, entity.dxf.start.y)
            end_point = (entity.dxf.end.x, entity.dxf.end.y)
            line_segments.append([start_point, end_point])
        elif entity.dxftype() in ('LWPOLYLINE', 'POLYLINE'):
            polyline_points = [(point[0], point[1]) for point in entity.get_points()]
            line_segments.extend([[polyline_points[i], polyline_points[i+1]] for i in range(len(polyline_points)-1)])

    merged = linemerge([LineString(segment) for segment in line_segments])

    if isinstance(merged, LineString):
        merged_lines = [merged]
    elif isinstance(merged, MultiLineString):
        merged_lines = [line for line in merged.geoms]
    else:
        merged_lines = []

    # Convert to a flat list of coordinates
    coords = []
    for geom in merged_lines:
        coords.extend(list(geom.coords))  # Use extend to flatten the list

    return np.array(coords)  # Create a numpy array from the flat list

# Load the DXF file points
steering_fist_points = load_dxf_points("Fist assembly.dxf")

# Create the figure
fig, ax = plt.subplots()
ax.set_xlim(-10, 70)
ax.set_ylim(-40, 80)
ax.set_aspect('equal')  # Keep equal proportions on the axes

line, = ax.plot([], [], 'o-', lw=2)  # Line of the arms
theta_a_text = ax.text(5, 60, '', fontsize=12, color='red')  # Theta_a angle text
theta_b_text = ax.text(5, -35, '', fontsize=12, color='red')  # Theta_b angle text

# Create a polygon for the DXF steering fist and add it to the plot
steering_fist_polygon = Polygon(np.zeros((0, 2)), closed=True, edgecolor='black', fill=False)
#steering_fist_polygon = Polygon([], closed=True, edgecolor='black', fill=False)
ax.add_patch(steering_fist_polygon)

# Coordinates of the pivot points in the DXF file
dxf_lower_joint = np.array([0.0, 0.0])  # Lower pivot point in the DXF
dxf_upper_joint = np.array([-13.5339, 62.3721])  # Upper pivot point in the DXF

def update():
    global theta_b

    # Get arm positions
    (x_a, y_a), (x_b, y_b), theta_a = get_positions(theta_b)

    line.set_data([x1, x_a, x_b, x2], [y1, y_a, y_b, y2])

    # Calculate the rotation and translation vector for the steering fist
    dxf_vector = dxf_upper_joint - dxf_lower_joint
    actual_vector = np.array([x_a - x_b, y_a - y_b])
    scale_factor = np.linalg.norm(actual_vector) / np.linalg.norm(dxf_vector)
    rotation_angle = np.arctan2(actual_vector[1], actual_vector[0]) - np.arctan2(dxf_vector[1], dxf_vector[0])

    # Rotate and scale the steering fist
    rotation_matrix = np.array([
        [np.cos(rotation_angle), -np.sin(rotation_angle)],
        [np.sin(rotation_angle), np.cos(rotation_angle)]
    ])
    fist_transformed = (steering_fist_points - dxf_lower_joint) * scale_factor @ rotation_matrix.T + [x_b, y_b]
    steering_fist_polygon.set_xy(fist_transformed)

    # Update angle text
    theta_a_text.set_text(f'θ_a: {np.degrees(theta_a):.2f}°')
    theta_b_text.set_text(f'θ_b: {np.degrees(theta_b):.2f}°')

    fig.canvas.draw()

# Function to handle key presses
def on_key(event):
    global theta_b
    if event.key == "right":
        if theta_b < theta_b_0 + angle_range:
            theta_b += np.radians(0.5)
    elif event.key == "left":
        if theta_b > theta_b_0 - angle_range:
            theta_b -= np.radians(0.5)
    update()

# Connect the key function to the canvas
fig.canvas.mpl_connect("key_press_event", on_key)

# Show the figure
plt.show()