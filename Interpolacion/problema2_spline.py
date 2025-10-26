# =============================================================================
# PROBLEMA 2: SPLINE CÚBICO NATURAL
# Universidad de Carabobo - Facultad de Ciencia y Tecnología
# Estudiante: Dervis Martínez - C.I: 31.456.326
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def f(x, d=6):
    """
    Función a interpolar: f(x) = 1 / (1 + (5x/(d+1))^2)
    """
    return 1 / (1 + (5*x/(d+1))**2)

def natural_cubic_spline(x, y):
    """
    Implementación del spline cúbico natural
    """
    n = len(x) - 1
    h = np.diff(x)
    
    # Construir sistema tridiagonal
    A = np.zeros((n+1, n+1))
    b = np.zeros(n+1)
    
    # Condiciones naturales: S''(x0) = S''(xn) = 0
    A[0,0] = 1
    A[n,n] = 1
    
    # Ecuaciones internas
    for i in range(1, n):
        A[i,i-1] = h[i-1]
        A[i,i] = 2*(h[i-1] + h[i])
        A[i,i+1] = h[i]
        b[i] = 3*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
    
    # Resolver sistema para c (segundas derivadas)
    c = np.linalg.solve(A, b)
    
    # Calcular coeficientes b y d
    b_coef = np.zeros(n)
    d_coef = np.zeros(n)
    
    for i in range(n):
        b_coef[i] = (y[i+1]-y[i])/h[i] - h[i]*(2*c[i] + c[i+1])/3
        d_coef[i] = (c[i+1]-c[i])/(3*h[i])
    
    return y[:-1], b_coef, c[:-1], d_coef

def main():
    """Implementación completa del Problema 2"""
    print("=" * 70)
    print("PROBLEMA 2: SPLINE CÚBICO NATURAL")
    print("=" * 70)
    print("UNIVERSIDAD DE CARABOBO - FACULTAD DE CIENCIA Y TECNOLOGÍA")
    print("ESTUDIANTE: Dervis Martínez - C.I: 31.456.326")
    print("=" * 70)
    
    d = 6
    a, b = -(d+1), d+1
    n_values = [4, 6, 8]
    
    x_plot = np.linspace(a, b, 1000)
    y_true = f(x_plot, d)
    
    # Configuración de gráficas
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Splines interpolantes
    ax1.plot(x_plot, y_true, 'k-', linewidth=2, label='f(x)')
    
    coeficientes_spline = {}
    
    for n in n_values:
        x_nodes = np.linspace(a, b, n+1)
        y_nodes = f(x_nodes, d)
        
        # Usar scipy para spline cúbico natural
        spline = CubicSpline(x_nodes, y_nodes, bc_type='natural')
        y_spline = spline(x_plot)
        
        ax1.plot(x_plot, y_spline, '--', linewidth=1.5, label=f'Spline n={n}')
        ax1.plot(x_nodes, y_nodes, 'o', markersize=6, alpha=0.7)
        
        # Calcular coeficientes manualmente
        a_coef, b_coef, c_coef, d_coef = natural_cubic_spline(x_nodes, y_nodes)
        coeficientes_spline[n] = {
            'a': a_coef, 'b': b_coef, 'c': c_coef, 'd': d_coef,
            'x_nodes': x_nodes
        }
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Interpolación Spline Cúbico Natural')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Errores del spline
    for n in n_values:
        x_nodes = np.linspace(a, b, n+1)
        y_nodes = f(x_nodes, d)
        
        spline = CubicSpline(x_nodes, y_nodes, bc_type='natural')
        y_spline = spline(x_plot)
        
        error = y_true - y_spline
        ax2.plot(x_plot, error, label=f'Error n={n}')
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('f(x) - S(x)')
    ax2.set_title('Error de Interpolación (Spline Cúbico Natural)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('spline_interpolacion.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Coeficientes detallados
    print("\nCOEFICIENTES SPLINES CÚBICOS NATURALES")
    print("=" * 60)
    for n in n_values:
        coefs = coeficientes_spline[n]
        x_nodes = coefs['x_nodes']
        
        print(f"\nSpline n={n}:")
        for i in range(len(coefs['a'])):
            print(f"Intervalo [{x_nodes[i]:.2f}, {x_nodes[i+1]:.2f}]:")
            print(f"  a_{i} = {coefs['a'][i]:.6f}, b_{i} = {coefs['b'][i]:.6f}, "
                  f"c_{i} = {coefs['c'][i]:.6f}, d_{i} = {coefs['d'][i]:.6f}")
    
    print(f"\nArchivo generado: spline_interpolacion.png")
    
    return coeficientes_spline

if __name__ == "__main__":
    resultados = main()