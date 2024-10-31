import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Adjusted lengths for suspension arms
l_a = 40.0      # Length of the upper arm
l_b = 61.8      # Length of the lower arm
l_c = 63.8      # Length of the load-bearing arm

# Width and height of the wheel (for visual reference only)
t_w = 15.0      
t_h = 30.0      
offset = 8.0    

# Anchor points on the chassis
x1, y1 = 5.0, 50.0
x2, y2 = 0.0, 0.0

# Initial angles (in radians)
theta_b_0 = np.radians(-15)  
theta_b = theta_b_0
angle_range = np.radians(25)

# Reference attachment points in the original design
dxf_lower_joint = np.array([0.0, 0.0])
dxf_upper_joint = np.array([-13.5339, 62.3721])

# Knuckle shape points according to the reference image
knuckle_shape = np.array([
    [0.0, 0.0],
    [2.35, 0.0],
    [2.35, 3.0],
    [12.75, 9.25],
    [15.5, 9.25],
    [15.5, 17.25],
    [15.5, 25.25],
    [12.75, 25.25],
    [2.35, 31.5],
    [-11.1839, 59.3721],
    [-11.1839, 62.3721],
    [-13.5339, 62.3721],
    [-15.8839, 62.3721],
    [-15.8839, 50.0],
    [-2.35, 20.0],
    [-2.35, 0.0]
])

# Additional knuckle point
cvd_additional_point = np.array([-9.0, 17.2642])
cvd_initial_point = np.array([-9.0, 16.0])

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

# Create the figure
fig, ax = plt.subplots()
ax.set_xlim(-55, 80)
ax.set_ylim(-60, 80)
ax.set_aspect('equal')

bottom_arm, = ax.plot([], [], 'o-', lw=2, color="blue")
upper_arm, = ax.plot([], [], 'o-', lw=2, color="blue")
cvd_line, = ax.plot([], [], 'o-', lw=1, color="green")
theta_a_text = ax.text(5, 65, '', fontsize=12, color='blue')
theta_b_text = ax.text(5, -45, '', fontsize=12, color='blue')
cvd_text = ax.text(-50, 20, '', fontsize=12, color='green')

# Draw the knuckle contour with rotation and scaling
def draw_knuckle(x_b, y_b, x_a, y_a):
    # Rotation and scaling vector
    dxf_vector = dxf_upper_joint - dxf_lower_joint
    actual_vector = np.array([x_a - x_b, y_a - y_b])
    scale_factor = np.linalg.norm(actual_vector) / np.linalg.norm(dxf_vector)
    rotation_angle = np.arctan2(actual_vector[1], actual_vector[0]) - np.arctan2(dxf_vector[1], dxf_vector[0])

    rotation_matrix = np.array([
        [np.cos(rotation_angle), -np.sin(rotation_angle)],
        [np.sin(rotation_angle), np.cos(rotation_angle)]
    ])
    
    # Transform knuckle shape points
    transformed_shape = (knuckle_shape - dxf_lower_joint) * scale_factor @ rotation_matrix.T + [x_b, y_b]
    cvd_transformed_point = (cvd_additional_point - dxf_lower_joint) * scale_factor @ rotation_matrix.T + [x_b, y_b]
    
    # Draw the knuckle as a polygon
    polygon = Polygon(transformed_shape, closed=True, fill=None, edgecolor='black', linewidth=2)
    ax.add_patch(polygon)
    cvd_line.set_data([cvd_initial_point[0], cvd_transformed_point[0]], [cvd_initial_point[1], cvd_transformed_point[1]])

    theta_cvd = np.arctan2(cvd_transformed_point[1] - cvd_initial_point[1], cvd_transformed_point[0] - cvd_initial_point[0])
    length_cvd = np.linalg.norm(cvd_transformed_point - cvd_initial_point)
    cvd_text.set_text(f'θ_cvd: {np.degrees(theta_cvd):.2f}°\nl_cvd: {length_cvd:.2f}mm')

# Update function for the plot
def update():
    global theta_b

    (x_a, y_a), (x_b, y_b), theta_a = get_positions(theta_b)
    upper_arm.set_data([x1, x_a], [y1, y_a])
    bottom_arm.set_data([x_b, x2], [y_b, y2])
    
    # Remove existing patches for knuckle contour, if any
    for patch in ax.patches:
        patch.remove()

    # Draw the new knuckle shape
    draw_knuckle(x_b, y_b, x_a, y_a)

    theta_a_text.set_text(f'θ_a: {np.degrees(theta_a):.2f}°')
    theta_b_text.set_text(f'θ_b: {np.degrees(theta_b):.2f}°')
    fig.canvas.draw()

# Handle keys for rotation
def on_key(event):
    global theta_b
    if event.key == "right":
        if theta_b < theta_b_0 + angle_range:
            theta_b += np.radians(0.5)
    elif event.key == "left":
        if theta_b > theta_b_0 - angle_range:
            theta_b -= np.radians(0.5)
    update()

fig.canvas.mpl_connect("key_press_event", on_key)

# Show the figure
plt.show()