from gurobipy import * 
import numpy as np
import matplotlib.pyplot as plt

# CREACIÓN DE CASOS
np.random.seed(5)

# Horizonte temporal:
T = 60
t = [i for i in range(1,T+1)]
g = 120/T

# Camiones:
K  = 5
k  = [i for i in range(1,K+1)]
Q  = 25 # Misma capacidad para todos los camiones
Q0 = [np.random.randint(10,16) for i in k] # Carga inicial del camión k

# Estaciones
N  = 50
n  = [i for i in range(1,N+1)]
C  = [np.random.randint(20,31) for i in n]
C0 = [np.random.randint(10,16) for i in n]

# Demanda
dem = {}
for i in n:
    # Clasificamos las estaciones
    tipo = np.random.randint(1,3)
    if tipo == 1:
        for j in t:
            dem[(i,j)] = np.random.randint(1,2*g)
    if tipo == 2:
        for j in t:
            dem[(i,j)] = -np.random.randint(1,2*g)

# Construcción de la red espacio temporal
    # Se incluyen 'u0' como posición inicial de cada camión y 'dummy'
    # como nodo auxiliar que representa el final de los trayectos
u0 = [-i for i in k]
dummy = ('dum','dum')
    # Conjunto de estados
S1 = [(i,0) for i in u0]
S2 = [(i,j) for i in n for j in t]
S = S1 + S2 + [dummy]
    # Conjunto de predecesores y sucesores
Pred = {(i,j):(i,j-1) for i in n for j in t if j > 1 }
Succ = {(i,j):(i,j+1) for i in n for j in t if j < T}
    # Localización de los nodos:
U = np.random.rand(N+K)*20
V = np.random.rand(N+K)*20
    # Distancia (== Tiempo) discretizada entre estaciones
dist ={(i,j): round(np.hypot(((abs(U[i-1]-U[j-1]))//g)+1,((abs(V[i-1]-V[j-1]))//g)+1)) 
     for i in n for j in n if i!=j}
    # Se añade distancia vehículo - estación
for i in u0:
    for j in n:
        dist[(i,j)]= round(np.hypot(((abs(U[i]-U[j-1]))//g)+1,((abs(V[i]-V[j-1]))//g)+1))
    # Conjunto de arcos
        # Todos los pares de viajes directos -> ((i,j),(i',j')) tq j'-j > d((i,j),(i',j')) > j'-j-1
A11 = [(i,j) for i in S2 for j in S2 if j[0]!=i[0] 
       and j[1]-i[1]>=dist[(i[0],j[0])] and dist[(i[0],j[0])]>j[1]-i[1]-1]
A12 = [(i,j) for i in S1 for j in S2 if j[1]-i[1]>=dist[(i[0],j[0])] and dist[(i[0],j[0])]>j[1]-i[1]-1]
A1 = A11 + A12
        # Esperar un tiempo en la misma estación
A2 = [(i,Succ[i]) for i in S2 if i[1] < T]
        # Fin de ruta desde cualquier estado
A3 = [(i,dummy) for i in S if i!= dummy]

A = A1 + A2 + A3
AA = A1 + A2
# Algunos conjuntos adicionales
A_V = [(a[0],a[1],v) for a in A for v in k] # (ArcoSalida,ArcoLlegada,k)

IncidPred = {s:[] for s in S if s != dummy} # {Arco:[ArcoPred1,ArcoPred2...]}
IncidSucc = {s:[] for s in S if s != dummy}
for a in AA:
    if a[0] not in IncidPred[a[1]]:
            IncidPred[a[1]].append(a[0])
    if a[1] not in IncidSucc[a[0]]:
            IncidSucc[a[0]].append(a[1])
for a in A3:
    IncidSucc[a[0]].append(a[1])

#CONSTRUCCIÓN DEL MODELO ARCO-FLUJO
# Declaración del modelo
mod = Model('RSBCD')

# Variables de decisión
x = mod.addVars(S2, vtype=GRB.CONTINUOUS, name = 'x')     # bicicletas que quedan en cada estado
ypos = mod.addVars(S2, vtype=GRB.CONTINUOUS, name = 'y+') # deficit de bicicletas en cada estado
yneg = mod.addVars(S2, vtype=GRB.CONTINUOUS, name = 'y-') # exceso de bicicletas en cada estado
delta = mod.addVars(A_V, vtype=GRB.BINARY, name = 'd')    # 1 si el vehículo k utiliza el arco a
q = mod.addVars(A_V, vtype=GRB.CONTINUOUS, name = 'q')    # carga del vehículo k al pasar por el arco a

# Función objetivo
mod.setObjective(quicksum(ypos[s]+yneg[s] for s in S2),GRB.MINIMIZE)

# sujeto a
    # Conservación de flujo (en cada estación, en cada período)
mod.addConstrs(quicksum(q[a,s,j] for a in IncidPred[s] for j in k) -
               quicksum(q[s,a,j] for a in IncidSucc[s] for j in k) -
               x[s] + ypos[s] - yneg[s] == dem[s] - C[s[0]-1] 
               for s in S2 if s[1] == 1)
mod.addConstrs(quicksum(q[a,s,j] for a in IncidPred[s] for j in k) -
               quicksum(q[s,a,j] for a in IncidSucc[s] for j in k) +
               x[Pred[s]] - x[s] + ypos[s] + yneg[s] == dem[s]
               for s in S2 if s[1] > 1)
    # Cada nodo es visitado, como máximo, una vez por período
mod.addConstrs(quicksum(delta[a,s,j] for a in IncidPred[s] for j in k) <= 1 
               for s in S2)
    # Conservación de flujo de los vehículos
mod.addConstrs(quicksum(delta[a,s,j] for a in IncidPred[s]) == 
               quicksum(delta[s,a,j] for a in IncidSucc[s]) 
               for j in k for s in S2)
    # Cada vehículo se utiliza una única vez
mod.addConstrs(quicksum(delta[s,a,j]  for s in S1 for a in IncidSucc[s]) == 1 for j in k)
    # Los vehículos salen con su carga inicial
mod.addConstrs(quicksum(q[s,a,j]  for s in S1 for a in IncidSucc[s]) == Q0[j-1] for j in k)
    # No negatividad y máximos de las variables continuas
mod.addConstrs(x[s] >= 0 for s in S2)
mod.addConstrs(x[s] <= C[s[0]-1] for s in S2)

mod.addConstrs(ypos[s] >= 0 for s in S2)
mod.addConstrs(yneg[s] >= 0 for s in S2)

mod.addConstrs(q[s] >= 0 for s in A_V)
mod.addConstrs(q[s] <= Q*delta[s] for s in A_V)

# Configuración de Gurobipy
mod.Params.TimeLimit = 3600
mod.Params.MIPGap = 0.1
# mod.feasRelaxS(0,False,False,True)
mod.Params.Heuristics = 0.1
mod.optimize()