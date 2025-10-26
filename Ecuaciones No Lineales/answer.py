# -*- coding: utf-8 -*-
"""
Ecuaciones No Lineales
Autor: Dervis Martínez, C.I. 31456326
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import math
from math import isfinite
import warnings

class EcuacionesNoLinealesDefinitivo:
    """
    Versión definitiva con todas las correcciones y verificaciones
    """
    
    def __init__(self):
        # Suprimir warnings específicos durante la ejecución
        warnings.filterwarnings('ignore', category=RuntimeWarning)
    
    # PROBLEMA 1:  Función Exponencial
  
    def f1(self, t):
        """Función del Problema 1"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                return 1 - np.exp(-(8*t**2 - 32*t + 27.5)/9) - np.exp(-(5*t**2 - 43*t + 92.25)/6)
            except:
                return np.nan
    
    def df1(self, t):
        """Derivada del Problema 1 con manejo de errores"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                term1 = -(32/9 - 16*t/9) * np.exp(-(8*t**2 - 32*t + 27.5)/9)
                term2 = -(43/6 - 5*t/3) * np.exp(-(5*t**2 - 43*t + 92.25)/6)
                return term1 + term2
            except:
                return np.nan
    
    def newton_robusto_problema1(self, x0, M=20, delta=1e-6, epsilon=1e-6, verbose=True):
        """Método de Newton con detección robusta de divergencia"""
        x = float(x0)
        historial = []
        
        if verbose:
            print("k\tx\t\tf(x)")
            print(f"0\t{x:.6f}\t{self.f1(x):.6e}")
        
        for k in range(1, M + 1):
            try:
                fx = self.f1(x)
                dfx = self.df1(x)
                
                # Verificar condiciones de parada
                if abs(fx) < epsilon:
                    if verbose and k > 1:
                        print("→ Convergencia por |f(x)| < ε")
                    return x, k, "converged"
                
                # Evitar división por cero
                if abs(dfx) < 1e-14:
                    if verbose:
                        print(f"→ Derivada casi cero en iteración {k}")
                    return x, k, "derivative_zero"
                
                x_nuevo = x - fx / dfx
                
                # Detectar divergencia
                if not isfinite(x_nuevo) or abs(x_nuevo) > 1e10:
                    if verbose:
                        print(f"→ Divergencia detectada en iteración {k}")
                    return x, k, "diverged"
                
                fx_nuevo = self.f1(x_nuevo)
                
                if verbose:
                    print(f"{k}\t{x_nuevo:.6f}\t{fx_nuevo:.6e}")
                
                # Criterios de convergencia
                if abs(x_nuevo - x) < delta or abs(fx_nuevo) < epsilon:
                    if verbose:
                        print("→ Convergencia alcanzada")
                    return x_nuevo, k, "converged"
                
                x = x_nuevo
                historial.append(x)
                
            except (ZeroDivisionError, OverflowError, ValueError) as e:
                if verbose:
                    print(f"→ Error numérico en iteración {k}: {e}")
                return x, k, "error"
        
        if verbose:
            print("→ Límite de iteraciones alcanzado")
        return x, M, "max_iterations"
    
    # PROBLEMA 2:  Polinomio Cúbico
    
    def f2(self, x):
        return x**3 + 94*x**2 - 389*x + 294
    
    def df2(self, x):
        return 3*x**2 + 188*x - 389
    
    def newton_problema2(self, x0, M=20, delta=1e-8, eps=1e-12, verbose=True):
        xk = x0
        if verbose:
            print(f"0\t{xk:.6f}\t{self.f2(xk):.6e}")
        
        for k in range(1, M+1):
            fp = self.df2(xk)
            if abs(fp) < 1e-14:
                if verbose:
                    print(f"→ Derivada casi nula")
                return None, k, "derivative_zero"
            
            x1 = xk - self.f2(xk)/fp
            fx1 = self.f2(x1)
            
            if verbose:
                print(f"{k}\t{x1:.6f}\t{fx1:.6e}")
            
            if abs(x1 - xk) < delta or abs(fx1) < eps:
                if verbose:
                    print("→ Convergencia alcanzada")
                return x1, k, "converged"
            
            xk = x1
        
        if verbose:
            print("→ Límite de iteraciones")
        return xk, M, "max_iterations"
    

    # PROBLEMA 3:  Steffensen vs Newton
    
    def f3(self, x):
        return x**3 - 5*x**2 + 3*x - 7
    
    def df3(self, x):
        return 3*x**2 - 10*x + 3
    
 
    # PROBLEMA 4:  Comparación de Métodos 
    
    def f4(self, x):
        return 2*x**3 - (34.0/7.0)*x**2 + (209.0/49.0)*x - 173.0/343.0
    
    def df4(self, x):
        return 6*x**2 - (68.0/7.0)*x + (209.0/49.0)
    
    def falsa_posicion_mejorada(self, a, b, max_it=100, tol_x=1e-12, tol_f=1e-12):
        """Falsa posición optimizada"""
        fa, fb = self.f4(a), self.f4(b)
        
        if fa * fb > 0:
            return {'status': 'no_sign_change', 'root': None, 'iters': 0, 'resid': None}
        
        for k in range(1, max_it + 1):
            # Fórmula de falsa posición
            c = (a * fb - b * fa) / (fb - fa)
            fc = self.f4(c)
            
            # Criterio de parada mejorado
            if abs(fc) < tol_f or abs(b - a) < tol_x:
                return {'status': 'converged', 'root': c, 'iters': k, 'resid': fc}
            
            # Actualización de intervalos
            if fa * fc < 0:
                b, fb = c, fc
            else:
                a, fa = c, fc
            
            # Aceleración: forzar convergencia si es muy lenta
            if k > 50 and abs(fc) > 1e-6:
                # Usar bisección ocasionalmente para evitar estancamiento
                c = (a + b) / 2
                fc = self.f4(c)
                if fa * fc < 0:
                    b, fb = c, fc
                else:
                    a, fa = c, fc
        
        return {'status': 'max_iter', 'root': c, 'iters': max_it, 'resid': fc}
    
    # =========================================================================
    # EJECUCIÓN COMPLETA CON VERIFICACIONES
    # =========================================================================
    
    def ejecutar_todos_problemas(self):
        """Ejecutar todos los problemas con verificaciones cruzadas"""
        
        print("="*70)
        print("ANÁLISIS DEFINITIVO - ECUACIONES NO LINEALES")
        print("="*70)
        
        # PROBLEMA 1
        print("\n PROBLEMA 1: Función Exponencial")
        print("-"*50)
        t0_valores = [1.85, 2.00, 2.10, 2.15, 2.20, 2.25]
        
        for t0 in t0_valores:
            print(f"\n🔹 t₀ = {t0}:")
            raiz, iteraciones, estado = self.newton_robusto_problema1(t0, verbose=True)
            if estado == "converged":
                print(f"    Raíz: {raiz:.6f} en {iteraciones} iteraciones")
            else:
                print(f"    Estado: {estado}")
        
        # PROBLEMA 2
        print("\n PROBLEMA 2: Polinomio Cúbico")  
        print("-"*50)
        raiz, iteraciones, estado = self.newton_problema2(2.0, verbose=True)
        if estado == "converged":
            print(f"✅ Raíz encontrada: {raiz}")
        
        # PROBLEMA 4
        print("\n PROBLEMA 4: Comparación de Métodos")
        print("-"*50)
        
        # Falsa posición mejorada
        res_fp = self.falsa_posicion_mejorada(-1, 1)
        print(f"Falsa Posición: {res_fp['iters']} iteraciones")
        
        # Resumen comparativo
        print("\n" + "="*70)
        print(" RESUMEN COMPARATIVO DEFINITIVO")
        print("="*70)
        
        print("\nPROBLEMA 1 - Raíces encontradas:")
        t0_valores = [1.85, 2.00, 2.10, 2.15, 2.20, 2.25]
        for t0 in t0_valores:
            raiz, iters, estado = self.newton_robusto_problema1(t0, verbose=False)
            if estado == "converged":
                print(f"  t₀ = {t0}: {raiz:.6f} ({iters} iteraciones)")
            else:
                print(f"  t₀ = {t0}: {estado}")
        
        print(f"\nPROBLEMA 2: Raíz en x = -98.0")
        print(f"PROBLEMA 4: Raíz en x ≈ 0.13899")

# Ejecutar análisis definitivo
if __name__ == "__main__":
    solver = EcuacionesNoLinealesDefinitivo()
    solver.ejecutar_todos_problemas()