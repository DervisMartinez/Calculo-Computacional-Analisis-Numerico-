# =============================================================================
# PROBLEMA 4: ANÁLISIS TEÓRICO - SPLINE CON FRONTERAS FIJAS
# Universidad de Carabobo - Facultad de Ciencia y Tecnología
# Estudiante: Dervis Martínez - C.I: 31.456.326
# =============================================================================

import numpy as np

def f(x, d=6):
    """
    Función a interpolar: f(x) = 1 / (1 + (5x/(d+1))^2)
    """
    return 1 / (1 + (5*x/(d+1))**2)

def f_derivative(x, d=6):
    """Derivada exacta de f(x)"""
    return - (50*x/(d+1)**2) / (1 + (5*x/(d+1))**2)**2

def clamped_cubic_spline(x, y, k0, kn):
    """
    Spline cúbico con condiciones de frontera fijas
    S'(x0) = k0, S'(xn) = kn
    """
    n = len(x) - 1
    h = np.diff(x)
    
    # Construir sistema tridiagonal
    A = np.zeros((n+1, n+1))
    b = np.zeros(n+1)
    
    # Condiciones de frontera fijas
    A[0,0] = 2 * h[0]
    A[0,1] = h[0]
    b[0] = 3 * ((y[1] - y[0]) / h[0] - k0)
    
    A[n,n-1] = h[n-1]
    A[n,n] = 2 * h[n-1]
    b[n] = 3 * (kn - (y[n] - y[n-1]) / h[n-1])
    
    # Ecuaciones internas
    for i in range(1, n):
        A[i,i-1] = h[i-1]
        A[i,i] = 2 * (h[i-1] + h[i])
        A[i,i+1] = h[i]
        b[i] = 3 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
    
    # Resolver para c (segundas derivadas)
    c = np.linalg.solve(A, b)
    
    # Calcular coeficientes b y d
    b_coef = np.zeros(n)
    d_coef = np.zeros(n)
    
    for i in range(n):
        b_coef[i] = (y[i+1] - y[i]) / h[i] - h[i] * (2*c[i] + c[i+1]) / 3
        d_coef[i] = (c[i+1] - c[i]) / (3 * h[i])
    
    return y[:-1], b_coef, c[:-1], d_coef

def main():
    """Análisis teórico del Problema 4"""
    print("=" * 70)
    print("PROBLEMA 4: ANÁLISIS TEÓRICO - SPLINE CON FRONTERAS FIJAS")
    print("=" * 70)
    print("UNIVERSIDAD DE CARABOBO - FACULTAD DE CIENCIA Y TECNOLOGÍA")
    print("ESTUDIANTE: Dervis Martínez - C.I: 31.456.326")
    print("=" * 70)
    
    d = 6
    a, b = -(d+1), d+1
    
    # Calcular derivadas en los extremos
    k0 = f_derivative(a, d)  # f'(-7)
    kn = f_derivative(b, d)  # f'(7)
    
    print(f"Derivadas exactas en extremos:")
    print(f"f'({a}) = {k0:.6f}")
    print(f"f'({b}) = {kn:.6f}")
    
    # Ejemplo de aplicación con n=6
    x_nodes = np.linspace(a, b, 7)  # n=6
    y_nodes = f(x_nodes, d)
    
    # Aplicar spline con fronteras fijas
    a_coef, b_coef, c_coef, d_coef = clamped_cubic_spline(x_nodes, y_nodes, k0, kn)
    
    print(f"\nEjemplo con n=6 - Spline con fronteras fijas:")
    for i in range(len(a_coef)):
        print(f"Intervalo [{x_nodes[i]:.2f}, {x_nodes[i+1]:.2f}]:")
        print(f"  a_{i} = {a_coef[i]:.6f}, b_{i} = {b_coef[i]:.6f}, "
              f"c_{i} = {c_coef[i]:.6f}, d_{i} = {d_coef[i]:.6f}")
    
    print("\nSISTEMA DE ECUACIONES PARA SPLINE CON FRONTERAS FIJAS")
    print("=" * 60)
    print("Para i = 1, 2, ..., n-1:")
    print("h_{i-1}c_{i-1} + 2(h_{i-1} + h_i)c_i + h_i c_{i+1} = 3[(a_{i+1} - a_i)/h_i - (a_i - a_{i-1})/h_{i-1}]")
    print("\nCondiciones de frontera fijas:")
    print("2h_0 c_0 + h_0 c_1 = 3[(a_1 - a_0)/h_0 - k_0]")
    print("h_{n-1} c_{n-1} + 2h_{n-1} c_n = 3[k_n - (a_n - a_{n-1})/h_{n-1}]")
    
    print("\nALGORITMO MODIFICADO:")
    print("1. Calcular h_i = x_{i+1} - x_i para i = 0,1,...,n-1")
    print("2. Calcular α_0 = 3[(a_1 - a_0)/h_0 - k_0]")
    print("3. Calcular α_n = 3[k_n - (a_n - a_{n-1})/h_{n-1}]")
    print("4. Para i = 1,2,...,n-1: α_i = 3[(a_{i+1} - a_i)/h_i - (a_i - a_{i-1})/h_{i-1}]")
    print("5. Resolver el sistema tridiagonal")
    print("6. Calcular b_i y d_i como en el spline natural")
    
    return {
        'k0': k0,
        'kn': kn,
        'coeficientes_clamped': (a_coef, b_coef, c_coef, d_coef)
    }

if __name__ == "__main__":
    resultados = main()