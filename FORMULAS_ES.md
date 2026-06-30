# Simulador de turbina eólica BEM — referencia de fórmulas

Este documento enumera todas las ecuaciones físicas utilizadas en la
simulación, en el orden en que el solver las aplica. Corresponde
directamente al código en `physics/bem.py` en el repositorio de GitHub —
cada sección a continuación indica la función donde vive esa ecuación,
para poder cotejarla con el código.

---

## 1. Magnitudes iniciales

**Relación de velocidad en punta (tip-speed ratio)**

```
lambda = omega * R / V
```
- `omega` — velocidad angular del rotor, rad/s (`omega = RPM * 2*pi / 60`)
- `R` — radio del rotor, m
- `V` — velocidad del viento, m/s

**Área barrida por el rotor**

```
A = pi * R^2
```

---

## 2. Geometría de la pala (en cada estación radial r)

La pala se discretiza en N estaciones radiales entre el radio del buje
`R_hub` y la punta `R`. En cada estación r:

**Cuerda (estrechamiento lineal)**
```
c(r) = c0 + (0.25*c0 - c0) * (r - R_hub)/(R - R_hub)
```
`c0` es la cuerda en la raíz; la cuerda en la punta se fija en el 25% de
la cuerda de raíz.

**Torsión (alabeo lineal)**
```
beta(r) = beta0 * (1 - (r - R_hub)/(R - R_hub))
```
`beta0` es el ángulo de torsión en la raíz (grados); la torsión disminuye
linealmente hasta 0° en la punta.

**Solidez local**
```
sigma(r) = B * c(r) / (2 * pi * r)
```
`B` es el número de palas. La solidez es la fracción del área anular
barrida en el radio r que está físicamente ocupada por material de pala.

---

## 3. Aerodinámica del perfil (coeficientes de sustentación y resistencia)

Este es el **modelo analítico simplificado** utilizado en lugar de datos
medidos en túnel de viento (ver la nota de "Precisión" en el README
principal — esta es la mayor fuente de desviación respecto a los números
de una turbina real).

**Antes de la pérdida (|alpha| <= 14°): teoría de perfil delgado**
```
Cl(alpha) = 2*pi * sin(alpha - alpha0)
Cd(alpha) = Cd_min + k * Cl^2
```
- `alpha0 = -2°` (ángulo de ataque de sustentación nula)
- `Cd_min = 0.008`, `k = 0.02` (polar parabólica de resistencia inducida)

**Después de la pérdida (|alpha| > 14°): transición hacia comportamiento
tipo placa plana**
```
Cl(alpha) = signo(alpha) * 1.1 * cos(min(|alpha|-14°, 76°))
Cd(alpha) = Cd_min + 4k + 1.4 * sin(min(|alpha|-14°, 90°))
```
Esto reproduce el comportamiento característico tras la pérdida
aerodinámica (caída de la sustentación y aumento brusco de la
resistencia), sin modelar las características de pérdida medidas de un
perfil real específico.

---

## 4. Iteración de la Teoría del Elemento de Pala y Cantidad de Movimiento (BEM)

Este es el núcleo de la simulación, resuelto de forma independiente en
cada estación radial y luego integrado a lo largo de toda la pala.

**Ángulo de entrada de flujo** (ángulo entre el plano del rotor y el
viento relativo local)
```
phi = atan2( (1-a)*V , (1+a')*omega*r )
```

**Ángulo de ataque**
```
alpha = phi - (beta(r) + theta)
```
`theta` es el ángulo de paso de pala (pitch, en grados), aplicado de
forma uniforme a lo largo de la envergadura.

**Coeficientes de fuerza proyectados sobre el plano del rotor**
```
Cn = Cl*cos(phi) + Cd*sin(phi)     (normal al plano del rotor -> empuje)
Ct = Cl*sin(phi) - Cd*cos(phi)     (tangencial -> par motor)
```

**Factor de pérdida en punta de pala (Prandtl)**
```
f_tip = (B/2) * (R - r) / (r * sin(phi))
F_tip = (2/pi) * arccos( exp(-f_tip) )
```

**Factor de pérdida en el buje (Prandtl)** (misma forma, reflejada en el
buje)
```
f_hub = (B/2) * (r - R_hub) / (r * sin(phi))
F_hub = (2/pi) * arccos( exp(-f_hub) )
```

**Factor de corrección combinado**
```
F = F_tip * F_hub
```
Esto tiene en cuenta que el rotor tiene un número finito de palas
discretas, en lugar de comportarse como un disco actuador continuo
idealizado — las velocidades inducidas disminuyen cerca de la punta y del
buje, y F captura ese efecto.

**Factor de inducción axial** (a partir del equilibrio entre cantidad de
movimiento y elemento de pala)
```
a = 1 / ( 4*F*sin^2(phi) / (sigma*Cn) + 1 )
```

**Corrección de Glauert para inducción alta** (aplicada cuando a > 0.4,
donde la teoría simple de cantidad de movimiento deja de ser válida —
el "estado de estela turbulenta")
```
K = 4*F*sin^2(phi) / (sigma*Cn)
a_c = 0.2
a = 0.5 * ( 2 + K*(1-2*a_c) - sqrt( (K*(1-2*a_c)+2)^2 + 4*(K*a_c^2 - 1) ) )
```

**Factor de inducción tangencial**
```
a' = 1 / ( 4*F*sin(phi)*cos(phi) / (sigma*Ct) - 1 )
```

Estas ecuaciones se resuelven de forma **iterativa**: se parte de un
valor inicial estimado para `a, a'` (0.2 y 0.02), se recalculan `phi`,
`alpha`, `Cl`, `Cd`, `Cn`, `Ct`, `F`, luego nuevos valores de `a, a'`, y
se repite hasta que dejan de cambiar (convergencia), típicamente en
20–60 iteraciones por estación.

---

## 5. Integración de fuerzas y potencia

Una vez que `a, a'` han convergido en todas las estaciones:

**Velocidad relativa al cuadrado en la estación**
```
Vrel^2 = ((1-a)*V)^2 + ((1+a')*omega*r)^2
```

**Contribución al empuje de esta estación** (anillo de ancho dr)
```
dT = 0.5 * rho * Vrel^2 * B * c(r) * Cn * dr
```

**Contribución al par motor de esta estación**
```
dQ = 0.5 * rho * Vrel^2 * B * c(r) * Ct * r * dr
```

Sumadas sobre todas las estaciones para obtener el empuje total `T` y el
par motor total `Q`.

**Potencia mecánica**
```
P = Q * omega
```

**Coeficiente de potencia**
```
Cp = P / (0.5 * rho * A * V^3)
```

**Coeficiente de empuje**
```
Ct_total = T / (0.5 * rho * A * V^2)
```

---

## 6. El límite de Betz

```
Cp_max = 16/27 ≈ 0.5926
```

Esta es la fracción máxima teórica de la energía cinética del viento que
*cualquier* rotor idealizado (disco actuador, sin pérdidas) puede
extraer — se deriva directamente de la teoría de cantidad de movimiento
1D y se aplica independientemente del diseño de la pala. No se impone
como una regla en el código; es una propiedad estructural de las
ecuaciones anteriores. Si un resultado lo superara, eso indicaría un
error en el código, no una turbina mejor.

Las turbinas reales suelen alcanzar un Cp entre 0.35 y 0.48 en su
relación de velocidad de punta de diseño, claramente por debajo del
límite de Betz, debido a las pérdidas en punta, la resistencia
aerodinámica y los efectos de número finito de palas que la corrección
de Glauert tiene en cuenta.

---

## 7. Objetivo del optimizador

El optimizador (ver `optimizer/optimize.py`) mantiene fija la geometría
de la pala y la velocidad del viento, y busca sobre el ángulo de paso
`theta` y la relación de velocidad de punta `lambda` para maximizar:

```
maximizar  P(theta, lambda)  sujeto a:
    theta  en [-10°, 20°]
    lambda en [2, 13]
```

resuelto mediante `scipy.optimize.differential_evolution` (un optimizador
global que no requiere gradiente), verificado mediante una búsqueda en
cuadrícula exhaustiva sobre los mismos límites. Ver los comentarios en el
archivo del optimizador para entender por qué no se utilizó
deliberadamente un método basado en gradiente (por ejemplo,
`scipy.optimize.minimize`): la función objetivo incorpora la resolución
iterativa descrita arriba, que no es suficientemente suave para una
convergencia fiable basada en gradientes.

---

## Dónde encontrar datos de turbinas reales

Todo lo anterior es teoría BEM internamente coherente y dimensionalmente
correcta. Lo que **no** es, es una validación frente a una turbina real
específica, porque la polar del perfil aerodinámico (sección 3) es
analítica y no medida. Para obtener números que coincidan con una máquina
real (por ejemplo, la turbina de referencia NREL de 5MW), las ecuaciones
de las secciones 1, 2 y 4–7 se mantienen igual — solo sería necesario
sustituir la polar del perfil aerodinámico (sección 3) y las funciones de
cuerda/torsión de la sección 2 por los datos publicados reales de esa
turbina.
