# =============================================================================
# PROBLEMA 1: INTERPOLACIÓN POLINOMIAL DE NEWTON
# Universidad de Carabobo - Facultad de Ciencia y Tecnología
# Estudiante: Dervis Martínez - C.I: 31.456.326
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

def f(x, d=6):
    """
    Función a interpolar: f(x) = 1 / (1 + (5x/(d+1))^2)
    Parámetros:
        x: punto(s) de evaluación
        d: último dígito de la cédula (6)
    """
    return 1 / (1 + (5*x/(d+1))**2)

def newton_divided_differences(x, y):
    """
    Calcula las diferencias divididas para el polinomio de Newton
    """
    n = len(x)
    F = np.zeros((n, n))
    F[:,0] = y
    
    for j in range(1, n):
        for i in range(n - j):
            F[i,j] = (F[i+1,j-1] - F[i,j-1]) / (x[i+j] - x[i])
    
    return F[0,:]  # Retorna los coeficientes de la primera fila

def newton_polynomial(x, nodes, coefs):
    """
    Evalúa el polinomio de Newton en el punto x
    """
    n = len(coefs)
    result = coefs[0]
    product = 1.0
    
    for i in range(1, n):
        product *= (x - nodes[i-1])
        result += coefs[i] * product
    
    return result

def main():
    """Implementación completa del Problema 1"""
    print("=" * 70)
    print("PROBLEMA 1: INTERPOLACIÓN POLINOMIAL DE NEWTON")
    print("=" * 70)
    print("UNIVERSIDAD DE CARABOBO - FACULTAD DE CIENCIA Y TECNOLOGÍA")
    print("ESTUDIANTE: Dervis Martínez - C.I: 31.456.326")
    print("=" * 70)
    
    # Parámetros del problema
    d = 6
    a, b = -(d+1), d+1  # Intervalo [-7, 7]
    n_values = [4, 6, 8]
    
    # Generar puntos para graficar
    x_plot = np.linspace(a, b, 1000)
    y_true = f(x_plot, d)
    
    # Configuración de gráficas
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Gráfica de los polinomios interpolantes
    ax1.plot(x_plot, y_true, 'k-', linewidth=2, label='f(x)')
    
    coeficientes_newton = {}
    
    for n in n_values:
        # Calcular nodos equiespaciados
        x_nodes = np.linspace(a, b, n+1)
        y_nodes = f(x_nodes, d)
        
        # Calcular coeficientes del polinomio de Newton
        coefs = newton_divided_differences(x_nodes, y_nodes)
        coeficientes_newton[n] = coefs
        
        # Evaluar polinomio en puntos de graficación
        y_interp = np.array([newton_polynomial(xi, x_nodes, coefs) for xi in x_plot])
        
        ax1.plot(x_plot, y_interp, '--', linewidth=1.5, label=f'P_{n}(x)')
        ax1.plot(x_nodes, y_nodes, 'o', markersize=6, alpha=0.7)
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Interpolación de Newton - Polinomios Interpolantes')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Gráfica de errores
    for n in n_values:
        x_nodes = np.linspace(a, b, n+1)
        y_nodes = f(x_nodes, d)
        
        coefs = newton_divided_differences(x_nodes, y_nodes)
        y_interp = np.array([newton_polynomial(xi, x_nodes, coefs) for xi in x_plot])
        
        error = y_true - y_interp
        ax2.plot(x_plot, error, label=f'Error n={n}')
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('f(x) - P(x)')
    ax2.set_title('Error de Interpolación (Método de Newton)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('newton_interpolacion.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Impresión de coeficientes
    print("\nCOEFICIENTES POLINOMIOS NEWTON")
    print("=" * 50)
    for n in n_values:
        x_nodes = np.linspace(a, b, n+1)
        print(f"\nPolinomio n={n}:")
        print(f"Nodos x: {x_nodes}")
        print(f"Coeficientes diferencias divididas: {coeficientes_newton[n]}")
    
    print(f"\nArchivo generado: newton_interpolacion.png")
    
    return coeficientes_newton

if __name__ == "__main__":
    resultados = main()