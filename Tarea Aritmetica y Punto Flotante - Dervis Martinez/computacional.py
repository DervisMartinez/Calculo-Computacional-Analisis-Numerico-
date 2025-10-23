import numpy as np
import matplotlib.pyplot as plt
import math




# EJERCICIO 5

x = 1.00001

# Forma 1: suma directa
P1 = sum([x**k for k in range(51)])

# Forma 2: evaluación recursiva (tipo Horner anidado)
def nested_poly(x, n=50):
    res = 1
    for _ in range(n):
        res = 1 + x*res
    return res
P2 = nested_poly(x)

# Forma 3: fórmula cerrada
P3 = (x**51 - 1)/(x - 1)

print("Ejercicio 5:")
print(f"P(x) suma directa    = {P1}")
print(f"P(x) forma anidada   = {P2}")
print(f"P(x) formula cerrada = {P3}\n")




# EJERCICIO 6

X = np.linspace(0.99, 1.01, 101)

# f(x)
f = X**8 - 8*X**7 + 28*X**6 - 56*X**5 + 70*X**4 - 56*X**3 + 28*X**2 - 8*X + 1

# g(x) factorizado
g = ((((((X-8)*X+28)*X-56)*X+70)*X-56)*X+28)*X-8
g = g*X+1

# h(x)
h = (X-1)**8

plt.figure(figsize=(10,6))
plt.plot(X, f, label="f(x)")
plt.plot(X, g, label="g(x)")
plt.plot(X, h, label="h(x)")
plt.title("Ejercicio 6: Comparación de f(x), g(x) y h(x)")
plt.legend()
plt.grid(True)
plt.show()


# EJERCICIO 7

def expr1(x):
    return (np.tan(x) - x)/x**3

def expr2(x):
    return (np.exp(x) - np.sin(x) + np.cos(x) - 2)/x**3

print("Ejercicio 7:")
for p in range(1,18):
    x = 10**(-p)
    e1 = expr1(x)
    e2 = expr2(x)
    print(f"p={p:<2d} -> expr1={e1:.6e}, expr2={e2:.6e}")

print("\nSe observa a partir de cierto p cuando las expresiones pierden dígitos significativos.")




# EJERCICIO 8
X = np.linspace(-np.pi, np.pi, 500)
f_exact = 1 - np.cos(X)

# Serie de Taylor de orden 6 (última cifra cédula = 6)
# cos(x) ≈ 1 - x^2/2! + x^4/4! - x^6/6!
# entonces: 1 - cos(x) ≈ x^2/2! - x^4/4! + x^6/6!
f_taylor = (X**2)/math.factorial(2) - (X**4)/math.factorial(4) + (X**6)/math.factorial(6)

plt.figure(figsize=(10,6))
plt.plot(X, f_exact, label="f(x) = 1 - cos(x)")
plt.plot(X, f_taylor, "--", label="Taylor orden 6")
plt.title("Ejercicio 8: Comparación función exacta vs Taylor")
plt.legend()
plt.grid(True)
plt.show()

# Diferencia entre ambas
plt.figure(figsize=(10,6))
plt.plot(X, f_exact - f_taylor, color="red", label="Error: exacta - Taylor")
plt.title("Diferencia entre f(x) y su aproximación de Taylor")
plt.legend()
plt.grid(True)
plt.show()
