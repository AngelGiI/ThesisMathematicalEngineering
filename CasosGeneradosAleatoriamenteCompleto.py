import time
import numpy as np
import copy

# Semilla del caso
np.random.seed(600)

# Registro del tiempo
start = time.time()

# Horizonte temporal
T = 60
g = 120/T
t = [i for i in range(1,T+1)]

# Vehículos
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
tipo={}
for i in n:
    # Clasificamos las estaciones
    tip = np.random.randint(1,3)
    tipo[i] = tip
    if tip == 1:
        for j in t:
            dem[(i,j)] = np.random.randint(1,2*g)
    if tip == 2:
        for j in t:
            dem[(i,j)] = -np.random.randint(1,2*g)

# Estados      
S1 = [(-i,0) for i in k]
S2 = [(i,j) for i in n for j in t]
SS = S1+S2

end = time.time()
print('Estados',end - start)
start = time.time()

# Conjunto de predecesores y sucesores
Pred = {(i,j):(i,j-1) for i in n for j in t if j > 1 }
Suc = {(i,j):(i,j+1) for i in n for j in t if j < T}
for i in k:
    Suc[(-i,0)]=[]
    
# Localización de los nodos:
U = np.random.rand(N+K)*20
V = np.random.rand(N+K)*20

# Distancia (== unidades de tiempo) discretizada entre estaciones
dist ={(i,j): round(np.hypot(((abs(U[i-1]-U[j-1]))//g)+1,((abs(V[i-1]-V[j-1]))//g)+1)) 
     for i in n for j in n if i!=j}
for i in n:
    dist[i,i]=1
# Se añade distancia vehículo - estación
for i in k:
    for j in n:
        dist[(-i,j)]= round(np.hypot(((abs(U[i]-U[j-1]))//3)+1,((abs(V[i]-V[j-1]))//3)+1))

end = time.time()
print('Distancias',end - start)
start = time.time()     

# Cálculo de los arcos
A1 = [(i,j) for i in SS for j in S2 if j[0]!=i[0] 
       and j[1]-i[1]>=dist[(i[0],j[0])] and dist[(i[0],j[0])]>j[1]-i[1]-1]
        # Esperar un tiempo en la misma estación
A2 = [(i,Suc[i]) for i in S2 if i[1] < T]
AA=A1+A2

end = time.time()
print('arcos',end - start)
start = time.time()

# Cálculo de Incidencias (deltas).
IncidPred = {s:[] for s in SS} # {Arco:[ArcoPred1,ArcoPred2...]}
IncidSucc = {s:[] for s in SS}
for a in AA:
    if a[0] not in IncidPred[a[1]]:
        IncidPred[a[1]].append(a[0])
    if a[1] not in IncidSucc[a[0]]:
        IncidSucc[a[0]].append(a[1])

end = time.time()
print('Incidencias',end - start)
start = time.time()

# Inicialización con idénticos niveles de llenado, capacidades y pesos
# para todos los procedimientos 
C = [np.random.randint(20,26) for i in n]
w1=0.01
w2=0.1
# AVET
Ct = [np.random.randint(10,16) for i in n]
Qt = [np.random.randint(10,16) for i in k]
# PILOT
QtM = Qt.copy()
CtM = Ct.copy()
# NO
Nivel = Ct.copy()

# Solución obtenida si no se opera (NO)
start = time.time() 
Noop = 0
tt = 1
while tt <= T:
    for s in n:
            if tipo[s]==1: # est de dejada, dem > 0
                Noop = Noop + max(0,dem[(s,tt)]-Nivel[s-1])
                Nivel[s-1] = Nivel[s-1] - min(dem[(s,tt)],Nivel[s-1])
            elif tipo[s]==2: # est de recogida, dem < 0
                Noop = Noop + max(0,Nivel[s-1] - dem[(s,tt)] - C[s-1])
                Nivel[s-1] = min(C[s-1],Nivel[s-1] - dem[(s,tt)])
    tt+=1
end = time.time()
print('NO',Noop,end - start)
start = time.time()


# Procedimiento heurístico AVET.
# Inicialización
rutas = {i:[(-i,0)] for i in k}
cantidad = {i:0 for i in k}
Fobj = 0
temp = 0
ocupados = []
ytemp ={i:0 for i in n}
# Programa principal
start = time.time()
while temp <= T:
    if temp != 0:
        for i in k: # Operaciones de carga y descarga
            if temp != 0 and rutas[i][-1][1] == temp:
                ocupados.remove(rutas[i][-1][0])
                ytemp[rutas[i][-1][0]]=0
                Ct[rutas[i][-1][0]-1] += cantidad[i]
                Qt[i-1] -= cantidad[i]
                
        for s in n: # Actualización de las demandas
            if tipo[s]==1: # est de dejada, dem > 0
                Fobj = Fobj + max(0,dem[(s,temp)]-Ct[s-1])
                Ct[s-1] = Ct[s-1] - min(dem[(s,temp)],Ct[s-1])
                ytemp[s]+= max(0,dem[(s,temp)]-Ct[s-1])
            elif tipo[s]==2: # est de recogida, dem < 0
                Fobj = Fobj + max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                Ct[s-1] = min(C[s-1],Ct[s-1] - dem[(s,temp)])
                ytemp[s]+= max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                
    for i in k: # Elección del siguiente estado a visitar
        if rutas[i][-1][1] == temp:
                pos = rutas[i][-1] # Se renombra por comodidad
                valor = 0
                elegido = ['dum']
                nbicis=0
                for s in IncidSucc[pos]:
                    if s[0] not in ocupados or s == Suc[pos]:
                        Nvalor = 0
                        demanda = 0
                        mrg = 0
                        invol = 0
                        
                        # Demanda acumulada del candidato en el momento de llegada
                        for j in range(1,dist[(pos[0],s[0])]+1):
                            if temp+j <= T:
                                demanda+= dem[s[0],temp+j]
                                
                        # Cálculo del márgen       
                        aux = demanda
                        tp = temp+dist[(pos[0],s[0])]+1
                        if tipo[s[0]]==2: # est de recogida, dem < 0
                            while Ct[s[0]-1] - aux < C[s[0]-1] and tp <= T:
                                mrg += 1
                                aux+= dem[s[0],tp]
                                tp+=1
                            # Si ya no hace falta balancear la estación:
                            if tp > T and Ct[s[0]-1] - aux <= C[s[0]-1]:
                                Criterio = 0
                            else:
                                invol = (ytemp[s[0]]+ Ct[s[0]-1] - demanda - C[s[0]-1])/dist[(pos[0],s[0])]
                                Nvalor = min(Ct[s[0]-1]-demanda,Q-Qt[i-1])/dist[(pos[0],s[0])]
                                Criterio = Nvalor - w1*mrg + w2*invol
                        
                        # Procedimiento análogo para estaciones de dejada
                        elif tipo[s[0]]==1:
                            while Ct[s[0]-1] - aux > 0 and tp <= T:
                                mrg += 1
                                aux+= dem[s[0],tp]
                                tp+=1
                            if tp > T and Ct[s[0]-1] - aux >= 0:
                                Criterio = 0
                            else:
                                invol = (ytemp[s[0]]+ demanda - Ct[s[0]-1])/dist[(pos[0],s[0])]
                                Nvalor = min(C[s[0]-1]-Ct[s[0]-1]+demanda,Qt[i-1])/dist[(pos[0],s[0])]
                                Criterio = Nvalor - w1*mrg + w2*invol
                                
                        if Criterio > valor:
                            elegido = [s]
                            valor = Criterio
                            nbicis = int(Nvalor*dist[(pos[0],s[0])])
                s = elegido[0]
                if s == 'dum':
                    pass
                elif tipo[s[0]]==2: # est de recogida, dem < 0
                    cantidad[i] = -nbicis
                elif tipo[s[0]]==1:
                    cantidad[i] = nbicis
                rutas[i].append(s)
                ocupados.append(s[0])
    temp+=1
end = time.time()
print('AVET',Fobj,end - start)
start = time.time()


# PILOT
# Inicialización
rutasM = {i:[(-i,0)] for i in k}
llegadaM = {i:0 for i in k}
FobjM = 0
tempM = 0
ocupadosM = []
ytempM ={i:0 for i in n}

# Programa principal
while tempM <= T:
    if tempM != 0:
        # Bloques 1 y 2 igual que en AVET, actualizando la solución 'M' de PILOT.
        for i in k:
            if rutasM[i][-1][1] == tempM:
                ocupadosM.remove(rutasM[i][-1][0])
                CtM[rutasM[i][-1][0]-1] += llegadaM[i]
                QtM[i-1] -= llegadaM[i]
                ytempM[rutasM[i][-1][0]]=0
        for s in n:
            if tipo[s]==1: # est de dejada, dem > 0
                FobjM = FobjM + max(0,dem[(s,tempM)]-CtM[s-1])
                CtM[s-1] = CtM[s-1] - min(dem[(s,tempM)],CtM[s-1])
                ytempM[s]+= max(0,dem[(s,tempM)]-CtM[s-1])
            elif tipo[s]==2: # est de recogida, dem < 0
                FobjM = FobjM + max(0,CtM[s-1] - dem[(s,tempM)] - C[s-1])
                CtM[s-1] = min(C[s-1],CtM[s-1] - dem[(s,tempM)])
                ytempM[s]+= max(0,dem[(s,tempM)]-CtM[s-1])
    for i in k:
        if rutasM[i][-1][1] == tempM and tempM < T:
            FobjTemp = 99999
            elegidoT = ['dum']
            llegadaM[i] = 0
            for suc in IncidSucc[rutasM[i][-1]]:
                if suc[0] not in ocupadosM or suc == Suc[rutasM[i][-1]]:
                    # Se inicializa AVET con los parámetros actuales de PILOT
                    rutas = copy.deepcopy(rutasM)
                    llegada = copy.deepcopy(llegadaM)
                    Qt = QtM.copy()
                    Ct = CtM.copy()
                    temp = tempM
                    Fobj=0
                    ocupados = ocupadosM.copy()
                    ytemp = ytempM.copy()
                    
                    # Se añade manualmente el candidato a evaluar
                    rutas[i].append(suc)
                    ocupados.append(suc[0])
                    demanda = 0
                    for j in range(1,dist[(rutas[i][-2][0],suc[0])]+1):
                        demanda+= dem[suc[0],tempM+j]
                    if tipo[suc[0]]==2: # est de recogida, dem < 0
                        nbicis = min(Ct[suc[0]-1]-demanda,Q-Qt[i-1])
                        llegada[i] = -nbicis
                    elif tipo[suc[0]]==1:
                        nbicis = min(C[suc[0]-1]-Ct[suc[0]-1]+demanda,Qt[i-1])
                        llegada[i] = nbicis    
                        
                    while temp <= T: 
                        if temp != tempM:
                            for ii in k:
                                if rutas[ii][-1][1] == temp:
                                    ocupados.remove(rutas[ii][-1][0])
                                    Ct[rutas[ii][-1][0]-1] += llegada[ii]
                                    Qt[ii-1] -= llegada[ii]
                                    ytemp[rutas[ii][-1][0]]=0
                            for s in n:
                                if tipo[s]==1: # est de dejada, dem > 0
                                    Fobj = Fobj + max(0,dem[(s,temp)]-Ct[s-1])
                                    Ct[s-1] = Ct[s-1] - min(dem[(s,temp)],Ct[s-1])
                                    ytemp[s]+= max(0,dem[(s,temp)]-Ct[s-1])
                                elif tipo[s]==2: # est de recogida, dem < 0
                                    Fobj = Fobj + max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                                    Ct[s-1] = min(C[s-1],Ct[s-1] - dem[(s,temp)])
                                    ytemp[s]+= max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                        for ii in k:
                            if rutas[ii][-1][1] == temp:
                                    valor = -0.1
                                    elegido = ['dum']
                                    nbicis=0
                                    for s in IncidSucc[rutas[ii][-1]]:
                                        if s[0] not in ocupados or s == Suc[rutas[ii][-1]]:
                                            Nvalor = 0
                                            demanda = 0
                                            mrg = 0
                                            invol = 0
                                            
                                            # Demanda acumulada del candidato en el momento de llegada
                                            for j in range(1,dist[(rutas[ii][-1][0],s[0])]+1):
                                                demanda+= dem[s[0],temp+j]
                                            
                                            # Margen
                                            aux = demanda
                                            tp = temp+dist[(rutas[ii][-1][0],s[0])]+1
                                            if tipo[s[0]]==2: # est de recogida, dem < 0
                                                while Ct[s[0]-1] - aux < C[s[0]-1] and tp <= T:
                                                    mrg += 1
                                                    aux+= dem[s[0],tp]
                                                    tp+=1
                                                if tp > T and Ct[s[0]-1] - aux <= C[s[0]-1]:
                                                    Criterio = 0
                                                else:
                                                    invol = (Ct[s[0]-1] - demanda - C[s[0]-1])/dist[(rutas[ii][-1][0],s[0])]
                                                    Nvalor = min(Ct[s[0]-1]-demanda,Q-Qt[ii-1])/dist[(rutas[ii][-1][0],s[0])]
                                                    Criterio = Nvalor - w1*mrg + w2*invol
                                            elif tipo[s[0]]==1:
                                                while Ct[s[0]-1] - aux > 0 and tp <= T:
                                                    mrg += 1
                                                    aux+= dem[s[0],tp]
                                                    tp+=1
                                                if tp > T and Ct[s[0]-1] - aux >= 0:
                                                    Criterio = 0
                                                else:    
                                                    invol = (demanda - Ct[s[0]-1])/dist[(rutas[ii][-1][0],s[0])]
                                                    Nvalor = min(C[s[0]-1]-Ct[s[0]-1]+demanda,Qt[ii-1])/dist[(rutas[ii][-1][0],s[0])]
                                                    Criterio = Nvalor - w1*mrg + w2*invol
                                            if Criterio > valor:
                                                elegido = [s]
                                                valor = Criterio
                                                nbicis = int(Nvalor*dist[(rutas[ii][-1][0],s[0])])
                                    s = elegido[0]
                                    if s == 'dum':
                                        pass
                                    elif tipo[s[0]]==2: # est de recogida, dem < 0
                                        llegada[ii] = -nbicis
                                    elif tipo[s[0]]==1:
                                        llegada[ii] = nbicis
                                    rutas[ii].append(s)
                                    ocupados.append(s[0])
                        temp +=1
                    # Una vez acaba AVET, se escoge el candidato con mejor función objetivo
                    if Fobj < FobjTemp:
                        elegidoT = suc
                        FobjTemp = Fobj
            # Se actualiza la solución con el estado elegido
            if tempM == T:
                pass
            else: # Se rehacen los cálculos perdidos del número de bicicletas
                demanda = 0
                for j in range(1,dist[(rutasM[i][-1][0],elegidoT[0])]+1):
                    demanda+= dem[elegidoT[0],tempM+j]
                if tipo[elegidoT[0]]==2: # est de recogida, dem < 0
                    nbicis = min(CtM[elegidoT[0]-1]-demanda,Q-QtM[i-1])
                    llegadaM[i] = -nbicis
                elif tipo[elegidoT[0]]==1:
                    nbicis = min(C[elegidoT[0]-1]-CtM[elegidoT[0]-1]+demanda,QtM[i-1])
                    llegadaM[i] = nbicis
                rutasM[i].append(elegidoT)
                ocupadosM.append(elegidoT[0])
    tempM+=1

end = time.time()
print('PILOT',FobjM,end - start)