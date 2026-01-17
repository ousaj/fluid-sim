import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import defaultdict

# Función para leer los valores de Average FPS y Total particles de un archivo
def extract_data_from_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        # Usamos expresiones regulares para encontrar los valores
        fps_match = re.search(r'Average FPS:\s*(\d+)', content)
        particles_match = re.search(r'Total particles:\s*(\d+)', content)
        
        if fps_match and particles_match:
            fps = int(fps_match.group(1))
            particles = int(particles_match.group(1))
            return fps, particles
        else:
            return None, None

directory = './benchmarks'

files = [f for f in os.listdir(directory) if f.endswith('.txt')]

fps_by_particles = defaultdict(lambda: {'sum_fps': 0, 'count': 0})

for file in files:
    file_path = os.path.join(directory, file)
    fps, particles = extract_data_from_file(file_path)
    if fps is not None and particles is not None:
        print(f"Archivo: {file}, FPS: {fps}, Partículas: {particles}")  # Línea de depuración
        fps_by_particles[particles]['sum_fps'] += fps
        fps_by_particles[particles]['count'] += 1

avg_fps_values = []
particles_values = []

for particles, data in fps_by_particles.items():
    avg_fps = data['sum_fps'] / data['count']
    avg_fps_values.append(avg_fps)
    particles_values.append(particles)

print(f"Valores de partículas: {particles_values}")
print(f"Promedios de FPS: {avg_fps_values}")

if particles_values and avg_fps_values:
    df = pd.DataFrame({'particles_values': particles_values, 'avg_fps_values': avg_fps_values})
    sns.barplot(x='particles_values', y='avg_fps_values', data=df)
    plt.title("Distribución de FPS por total de partículas")
    plt.xlabel("Total de partículas")
    plt.ylabel("Promedio de FPS")
    plt.show()
else:
    print("No se encontraron datos para graficar.")
plt.show()
