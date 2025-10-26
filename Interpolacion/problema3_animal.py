import numpy as np
import matplotlib.pyplot as plt

# ========== IMPLEMENTACIÓN DE SPLINES NATURALES ==========

def spline_natural_burden_faires(x, y, xx):
    """
    Interpolación por splines cúbicos naturales (método Burden & Faires)
    
    Parámetros:
    x, y: puntos de datos
    xx: puntos donde evaluar el spline
    
    Retorna:
    yy: valores interpolados en xx
    a, b, c, d: coeficientes del spline
    h: diferencias entre puntos x
    """
    n = len(x) - 1
    
    # Paso 1: Calcular h_i = x_{i+1} - x_i
    h = np.zeros(n)
    for i in range(n):
        h[i] = x[i+1] - x[i]
    
    # Paso 2: Calcular alpha_i
    alpha = np.zeros(n+1)
    for i in range(1, n):
        alpha[i] = (3/h[i]) * (y[i+1] - y[i]) - (3/h[i-1]) * (y[i] - y[i-1])
    
    # Paso 3: Resolver sistema tridiagonal para c
    l = np.zeros(n+1)
    mu = np.zeros(n+1)
    z = np.zeros(n+1)
    c = np.zeros(n+1)
    
    # Condiciones naturales
    l[0] = 1
    mu[0] = 0
    z[0] = 0
    
    for i in range(1, n):
        l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]
    
    l[n] = 1
    z[n] = 0
    c[n] = 0
    
    # Sustitución hacia atrás
    for j in range(n-1, -1, -1):
        c[j] = z[j] - mu[j] * c[j+1]
    
    # Paso 4: Calcular b y d
    b = np.zeros(n)
    d = np.zeros(n)
    for i in range(n):
        b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (c[i+1] + 2 * c[i]) / 3
        d[i] = (c[i+1] - c[i]) / (3 * h[i])
    
    # Coeficientes a son simplemente los valores y
    a = y[:-1].copy()
    
    # Evaluar el spline en los puntos xx
    yy = np.zeros(len(xx))
    
    for idx, x_val in enumerate(xx):
        # Encontrar el intervalo correcto
        i = 0
        for j in range(n):
            if x[j] <= x_val <= x[j+1]:
                i = j
                break
            elif x_val < x[0]:
                i = 0
            elif x_val > x[n]:
                i = n-1
        
        # Evaluar spline cúbico
        dx = x_val - x[i]
        yy[idx] = a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3
    
    return yy, a, b, c, d, h

# ========== PROBLEMA 3: CURVA PLANA NO FUNCIONAL ==========

# Puntos proporcionados (sustituyen a los anteriores)
puntos_animal = np.array([
    [6.848, 4.44964871],
    [6.416, 4.07494145],
    [6.16, 3.62997658],
    [5.936, 3.1381733],
    [5.632, 2.64637002],
    [5.344, 2.22482436],
    [5.12, 1.59250585],
    [4.864, 0.91334895],
    [5.456, 1.8501171],
    [5.936, 2.4824356],
    [6.48, 3.06791569],
    [6.976, 3.34894614],
    [7.344, 3.62997658],
    [7.6, 3.93442623],
    [7.84, 4.33255269],
    [8.32, 4.44964871],
    [8.608, 4.44964871],
    [8.752, 4.73067916],
    [8.48, 5.19906323],
    [8.208, 5.50351288],
    [7.824, 5.38641686],
    [7.536, 5.22248244],
    [6.96, 5.03512881],
    [6.432, 4.84777518],
    [6.016, 4.75409836],
    [5.552, 4.56674473],
    [5.216, 4.4028103],
    [4.864, 4.23887588],
    [4.336, 4.02810304],
    [3.712, 3.7236534],
    [3.104, 3.60655738],
    [2.56, 3.39578454],
    [2.128, 3.32552693],
    [1.76, 3.37236534],
    [1.392, 3.44262295],
    [0.976, 3.7704918],
    [1.152, 4.07494145],
    [1.536, 4.02810304],
    [2.048, 4.00468384],
    [2.528, 3.95784543],
    [3.04, 3.95784543],
    [3.488, 4.09836066],
    [3.856, 4.3793911],
    [4.352, 4.87119438],
    [4.768, 5.36299766],
    [5.28, 5.6206089],
    [5.76, 5.6206089],
    [6.128, 5.05854801],
    [6.032, 4.49648712],
    [5.776, 3.91100703],
    [5.76, 3.46604215],
    [6.016, 3.20843091],
    [6.368, 3.02107728],
    [6.576, 2.76346604],
    [6.752, 2.34192037],
    [6.304, 2.55269321],
    [5.984, 2.66978923],
    [5.504, 2.85714286],
    [5.024, 2.92740047],
    [4.944, 3.30210773],
    [5.056, 3.84074941]
])

print(f"Número de puntos: {puntos_animal.shape[0]}")

# Separar coordenadas
x_puntos = puntos_animal[:, 0]
y_puntos = puntos_animal[:, 1]

# 2. CALCULAR PARÁMETRO t (longitud acumulada del polígono)
t_param = [0.0]
for i in range(1, len(puntos_animal)):
    dx = x_puntos[i] - x_puntos[i-1]
    dy = y_puntos[i] - y_puntos[i-1]
    distancia = np.sqrt(dx**2 + dy**2)
    t_param.append(t_param[-1] + distancia)

t_param = np.array(t_param)

# 3. CREAR PUNTOS t PARA INTERPOLACIÓN
t_fino = np.linspace(t_param[0], t_param[-1], 500)

# 4. CONSTRUIR SPLINES PARA x(t) y y(t)
print("Calculando spline para x(t)...")
x_spline, a_x, b_x, c_x, d_x, h_x = spline_natural_burden_faires(t_param, x_puntos, t_fino)

print("Calculando spline para y(t)...")
y_spline, a_y, b_y, c_y, d_y, h_y = spline_natural_burden_faires(t_param, y_puntos, t_fino)

# ========== GRÁFICAS ==========

# Gráfica 1: Puntos originales
plt.figure(figsize=(12, 10))
plt.plot(x_puntos, y_puntos, 'ro', markersize=4, linewidth=1, label='Puntos originales')
plt.title('PROBLEMA 3: Puntos Originales del Animal\n(60 puntos proporcionados)', fontsize=14, fontweight='bold')
plt.xlabel('Coordenada X')
plt.ylabel('Coordenada Y')
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.legend()
plt.savefig('problema3_puntos_originales.png', dpi=300, bbox_inches='tight')

# Gráfica 2: Curva interpolada
plt.figure(figsize=(12, 10))
plt.plot(x_spline, y_spline, 'b-', linewidth=2.5, label='Curva spline interpolada')
plt.plot(x_puntos, y_puntos, 'ro', markersize=4, alpha=0.6, label='Puntos de control')
plt.title('PROBLEMA 3: Curva del Animal - Spline Interpolante\n(Interpolación paramétrica x(t), y(t))', fontsize=14, fontweight='bold')
plt.xlabel('Coordenada X')
plt.ylabel('Coordenada Y')
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.legend()
plt.savefig('problema3_curva_interpolada.png', dpi=300, bbox_inches='tight')

# Gráfica 3: Comparación lado a lado
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Subgráfica 1: Puntos originales
ax1.plot(x_puntos, y_puntos, 'ro', markersize=4, linewidth=1)
ax1.set_title('Puntos Originales (60 puntos)')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.grid(True, alpha=0.3)
ax1.axis('equal')

# Subgráfica 2: Curva interpolada
ax2.plot(x_spline, y_spline, 'b-', linewidth=2)
ax2.plot(x_puntos, y_puntos, 'ro', markersize=3, alpha=0.5)
ax2.set_title('Curva Spline Interpolada')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.grid(True, alpha=0.3)
ax2.axis('equal')

plt.tight_layout()
plt.savefig('problema3_comparacion.png', dpi=300, bbox_inches='tight')

plt.show()

print("¡Interpolación completada! Se han generado las gráficas:")
print("- problema3_puntos_originales.png")
print("- problema3_curva_interpolada.png") 
print("- problema3_comparacion.png")