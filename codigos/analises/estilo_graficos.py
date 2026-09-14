"""Cores e interpolacao compartilhadas pelas analises."""
from __future__ import annotations
import numpy as np
from matplotlib import colormaps

ORDEM_N = (1, 3, 5, 10, 25, 50)
CORES_N = {n: colormaps["Blues"](v) for n, v in zip(ORDEM_N, np.linspace(0.92, 0.38, len(ORDEM_N)))}
CORES_H = {
    0.01: "#2166ac",  # azul
    0.50: "#d6604d",  # laranja-avermelhado
    0.99: "#1b9e77",  # verde
}
MAPAS_H = {
    0.01: colormaps["Blues"],
    0.50: colormaps["Oranges"],
    0.99: colormaps["Greens"],
}

def cor_n(valor):
    return CORES_N.get(int(float(valor)), colormaps["Blues"](0.55))

def cor_h(valor):
    valor = float(valor)
    chave = min(CORES_H, key=lambda x: abs(x - valor))
    return CORES_H[chave]

def tons_h(valor, quantidade):
    """Retorna tons da familia de H, do escuro ao claro."""
    valor = float(valor)
    chave = min(MAPAS_H, key=lambda x: abs(x - valor))
    return MAPAS_H[chave](np.linspace(0.92, 0.38, int(quantidade)))

def curva_suave(x, y, pontos=400, log_x=False, log_y=False):
    """PCHIP preservador de forma; mantém os originais se houver menos de 3 pontos."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    mascara = np.isfinite(x) & np.isfinite(y)
    if log_x: mascara &= x > 0
    if log_y: mascara &= y > 0
    x, y = x[mascara], y[mascara]
    ordem = np.argsort(x); x, y = x[ordem], y[ordem]
    x, indices = np.unique(x, return_index=True); y = y[indices]
    if len(x) < 3: return x, y
    xi = np.log10(x) if log_x else x
    yi = np.log10(y) if log_y else y
    xd = np.linspace(xi.min(), xi.max(), pontos)
    # Derivadas PCHIP de Fritsch-Carlson: preservam monotonicidade e forma.
    h = np.diff(xi)
    delta = np.diff(yi) / h
    m = np.zeros_like(yi)
    if len(yi) == 3:
        interior = (delta[:-1] * delta[1:]) > 0
    else:
        interior = (delta[:-1] * delta[1:]) > 0
    for k in range(1, len(yi) - 1):
        if interior[k - 1]:
            w1, w2 = 2 * h[k] + h[k - 1], h[k] + 2 * h[k - 1]
            m[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])
    m[0] = ((2*h[0] + h[1])*delta[0] - h[0]*delta[1]) / (h[0] + h[1])
    m[-1] = ((2*h[-1] + h[-2])*delta[-1] - h[-1]*delta[-2]) / (h[-1] + h[-2])
    if np.sign(m[0]) != np.sign(delta[0]): m[0] = 0
    elif abs(m[0]) > 3*abs(delta[0]): m[0] = 3*delta[0]
    if np.sign(m[-1]) != np.sign(delta[-1]): m[-1] = 0
    elif abs(m[-1]) > 3*abs(delta[-1]): m[-1] = 3*delta[-1]
    j = np.clip(np.searchsorted(xi, xd) - 1, 0, len(xi) - 2)
    s = (xd - xi[j]) / h[j]
    yd = ((2*s**3 - 3*s**2 + 1)*yi[j] + (s**3 - 2*s**2 + s)*h[j]*m[j]
          + (-2*s**3 + 3*s**2)*yi[j+1] + (s**3 - s**2)*h[j]*m[j+1])
    return (10**xd if log_x else xd), (10**yd if log_y else yd)
