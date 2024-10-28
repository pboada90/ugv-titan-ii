import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.animation import FuncAnimation

# Longitudes (ajustables según el diseño)
l_a = 30.0      # longitud del brazo superior
l_b = 51.8      # longitud del brazo inferior
l_c = 63.8      # longitud del brazo portamasa
t_w = 15.0      # ancho de la rueda (rectángulo)
t_h = 30.0      # altura de la rueda (rectángulo)
offset = 8.0    # Desplazamiento adicional para mover la rueda hacia adelante

# Puntos de anclaje en el chasis
x1, y1 = 4.0, 50.0
x2, y2 = 0.0, 0.0

# Definir ángulos iniciales (en radianes)
theta_b_0 = np.radians(-20)  # El brazo inferior también empieza 20° abajo
theta_b = theta_b_0
angle_range = np.radians(20)

# Función para obtener las coordenadas de los puntos
def get_positions(theta_b):    
    # Coordenadas del extremo del brazo inferior
    x_b = x2 + l_b * np.cos(theta_b)
    y_b = y2 + l_b * np.sin(theta_b)

    m = -(x1 - x_b) / (y1 - y_b)
    b = (l_c**2 - l_a**2 - x_b**2 + x1**2 - y_b**2 + y1**2) / (2 * (y1 - y_b))

    a = 1 + m**2
    b_2 = 2*(m*(b - y_b) - x_b)
    c = x_b**2 + (b - y_b)**2 - l_c**2

    root_1 = (-b_2 + np.sqrt(b_2**2 - 4*a*c))/(2*a)
    root_2 = (-b_2 - np.sqrt(b_2**2 - 4*a*c))/(2*a)

    root = max([root_1, root_2])

    x_a = root
    y_a = m*root + b

    theta_a = np.arctan2(y_a - y1, x_a - x1)
    
    return (x_a, y_a), (x_b, y_b), theta_a

# Función para obtener las esquinas del rectángulo (rueda)
def get_rectangle(x_c, y_c, angle, width, height):
    # Coordenadas de las esquinas del rectángulo sin rotar
    rect = np.array([
        [height / 2, -width / 2],
        [height / 2, width / 2],
        [-height / 2, width / 2],
        [-height / 2, -width / 2],
        [height / 2, -width / 2]
    ])
    
    # Matriz de rotación
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    
    # Rotar y trasladar el rectángulo a su posición
    rect_rotated = rect @ rotation_matrix.T
    rect_rotated[:, 0] += x_c
    rect_rotated[:, 1] += y_c
    
    return rect_rotated[:, 0], rect_rotated[:, 1]

# Crear la figura
fig, ax = plt.subplots()
ax.set_xlim(-10, 70)
ax.set_ylim(-40, 80)
ax.set_aspect('equal')  # Mantener proporciones iguales en los ejes

line, = ax.plot([], [], 'o-', lw=2)  # Línea de los brazos
wheel, = ax.plot([], [], 'b-')  # Dibujo de la rueda (rectángulo)
theta_a_text = ax.text(5, 60, '', fontsize=12, color='red')  # Texto del ángulo theta_a
theta_b_text = ax.text(5, -35, '', fontsize=12, color='red')  # Texto del ángulo theta_b

# Función para actualizar la animación
def update():
    global theta_b

    # Obtener posiciones de los brazos
    (x_a, y_a), (x_b, y_b), theta_a = get_positions(theta_b)
    
    # Coordenadas del centro de la rueda (punto medio entre los extremos de los brazos)
    x_c = (x_a + x_b) / 2
    y_c = (y_a + y_b) / 2

    # Calcular el ángulo de la línea que conecta los extremos de los brazos
    angle_line = np.arctan2(y_b - y_a, x_b - x_a)

    # Desplazar la rueda hacia adelante en la dirección de la línea que conecta los brazos
    x_c += offset * np.cos(angle_line)
    y_c += offset * np.sin(angle_line)

    # Obtener las esquinas del rectángulo (rueda)
    rect_x, rect_y = get_rectangle(x_c, y_c, angle_line, t_w, t_h)

    # Actualizar la línea de los brazos
    line.set_data([x1, x_a, x_b, x2], [y1, y_a, y_b, y2])

    # Actualizar el rectángulo (rueda)
    wheel.set_data(rect_x, rect_y)

    theta_a_deg = np.degrees(theta_a)
    theta_a_text.set_text(f'θ_a: {theta_a_deg:.2f}°')

    theta_b_deg = np.degrees(theta_b)
    theta_b_text.set_text(f'θ_b: {theta_b_deg:.2f}°')

    fig.canvas.draw()

# Función para manejar las teclas
def on_key(event):
    global theta_b
    if event.key == "right":
        if theta_b < theta_b_0 + angle_range:
            theta_b += np.radians(0.5)
    elif event.key == "left":
        if theta_b > theta_b_0 - angle_range:
            theta_b -= np.radians(0.5)
    update()

# Conectar la función de teclas al canvas
fig.canvas.mpl_connect("key_press_event", on_key)

# Mostrar la figura
plt.show()