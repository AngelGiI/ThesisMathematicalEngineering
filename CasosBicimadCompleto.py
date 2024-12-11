import time
import json
import numpy as np
import copy

f1 = open('C:\\Users\Angel\Desktop\Estaciones.json')
data = json.load(f1)
f1.close()

MinLong = 180
MaxLong = -180
MinLat = 180
MaxLat = -180

Estaciones = {}
for i in data["stations"]:
    Estaciones[i["id"]]=[i["name"],i["total_bases"],i["dock_bikes"],float(i["longitude"]),float(i["latitude"])]
    if Estaciones[i["id"]][3] < MinLong:
        MinLong = Estaciones[i["id"]][3]
    elif Estaciones[i["id"]][3] > MaxLong:
        MaxLong = Estaciones[i["id"]][3]
    elif Estaciones[i["id"]][4] < MinLat:
        MinLat = Estaciones[i["id"]][4]
    elif Estaciones[i["id"]][4] > MaxLat:
        MaxLat = Estaciones[i["id"]][4]

NoF = []
for i in range(1,269):
    if i not in Estaciones.keys():
        NoF.append(i)

Movimientos = {}
for i in Estaciones.keys():
    Estaciones[i][3] = round((Estaciones[i][3]-MinLong)*100/(MaxLong-MinLong),2)
    Estaciones[i][4] = round((Estaciones[i][4]-MinLat)*100/(MaxLat-MinLat),2)
    Movimientos[i]= np.zeros((24,), dtype=int)


f2 = open('C:\\Users\Angel\Desktop\\202105_movements.json')
mov = [json.loads(line) for line in f2]
f2.close()


Dist = {}
for i in mov:
    Dist[(i["idunplug_station"],i["idplug_station"])] = i["travel_time"]//60
    if i["idplug_station"] != 2009 and i["idunplug_station"] != 2009:
        if  'T0' in i["unplug_hourTime"]:
            if 'T00' in i["unplug_hourTime"]:
                Movimientos[i["idplug_station"]][0]+=1
                Movimientos[i["idunplug_station"]][0]-=1
            elif 'T01' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][1]+=1
                Movimientos[i["idunplug_station"]][1]-=1  
            elif 'T02' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][2]+=1
                Movimientos[i["idunplug_station"]][2]-=1
            elif 'T03' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][3]+=1
                Movimientos[i["idunplug_station"]][3]-=1
            elif 'T04' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][4]+=1
                Movimientos[i["idunplug_station"]][4]-=1
            elif 'T05' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][5]+=1
                Movimientos[i["idunplug_station"]][5]-=1
            elif 'T06' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][6]+=1
                Movimientos[i["idunplug_station"]][6]-=1
            elif 'T07' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][7]+=1
                Movimientos[i["idunplug_station"]][7]-=1
            elif 'T08' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][8]+=1
                Movimientos[i["idunplug_station"]][8]-=1
            elif 'T09' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][9]+=1
                Movimientos[i["idunplug_station"]][9]-=1
        if  'T1' in i["unplug_hourTime"]:
            if 'T10' in i["unplug_hourTime"]:
                Movimientos[i["idplug_station"]][10]+=1
                Movimientos[i["idunplug_station"]][10]-=1
            elif 'T11' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][11]+=1
                Movimientos[i["idunplug_station"]][11]-=1  
            elif 'T12' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][12]+=1
                Movimientos[i["idunplug_station"]][12]-=1
            elif 'T13' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][13]+=1
                Movimientos[i["idunplug_station"]][13]-=1
            elif 'T14' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][14]+=1
                Movimientos[i["idunplug_station"]][14]-=1
            elif 'T15' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][15]+=1
                Movimientos[i["idunplug_station"]][15]-=1
            elif 'T16' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][16]+=1
                Movimientos[i["idunplug_station"]][16]-=1
            elif 'T17' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][17]+=1
                Movimientos[i["idunplug_station"]][17]-=1
            elif 'T18' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][18]+=1
                Movimientos[i["idunplug_station"]][18]-=1
            elif 'T19' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][19]+=1
                Movimientos[i["idunplug_station"]][19]-=1
        if  'T2' in i["unplug_hourTime"]:
            if 'T20' in i["unplug_hourTime"]:
                Movimientos[i["idplug_station"]][20]+=1
                Movimientos[i["idunplug_station"]][20]-=1
            elif 'T21' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][21]+=1
                Movimientos[i["idunplug_station"]][21]-=1  
            elif 'T22' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][22]+=1
                Movimientos[i["idunplug_station"]][22]-=1
            elif 'T23' in i['unplug_hourTime']:
                Movimientos[i["idplug_station"]][23]+=1
                Movimientos[i["idunplug_station"]][23]-=1


tipo = {}
Def = 0
Sob = 0
Nul = 0
for i in Movimientos.keys():
    tipo[i] = sum(Movimientos[i][j] for j in range(0,20))
    if tipo[i] > 0 :
        Sob+=1
    elif tipo[i] < 0:
        Def+=1
    else:
        Nul+=1
    for j in range(0,24):
        if tipo[i] > 0 and Movimientos[i][j] < 0:
            Movimientos[i][j] = 0
        elif tipo[i] < 0 and Movimientos[i][j] > 0:
            Movimientos[i][j] = 0
        elif tipo[i] == 0 and Movimientos[i][j] > 0:
            Movimientos[i][j] = 0
 
mov = {}
for i in Movimientos.keys():
    mov[i] = Movimientos[i][0:20]


np.random.seed(13)
# Horizonte temporal
T = 60
g = 120/T
t = [i for i in range(1,T+1)]

for i in Dist:
    Dist[i]=round(Dist[i]//g)
    if Dist[i]==0:
        Dist[i]+=1

# Vehículos
K  = 5
k  = [i for i in range(1,K+1)]
Q  = 25 # Misma capacidad para todos los camiones
Q0 = [np.random.randint(10,16) for i in k] # Carga inicial del camión k

# Estaciones
N  = 263 # Con 263 va bien, pero con 264 nunca termina el bucle
         # porque se omitió una estación con datos extremadamente atípicos
n  = [i for i in range(1,N+1)]

# INSTANCIAS 
instancia = []
while len(instancia) < N:
    i = np.random.randint(1,270)
    if i not in Estaciones.keys():
        pass
    elif i not in instancia:
        instancia.append(i)
Estaciones[i][3] = round((Estaciones[i][3]-MinLong)*100/(MaxLong-MinLong),2)
Estaciones[i][4] = round((Estaciones[i][4]-MinLat)*100/(MaxLat-MinLat),2)
C=[]
C0=[]
dem={}
U = []
V = []
ttipo={}

for i in instancia:
    C.append(Estaciones[i][1])
    C0.append(Estaciones[i][2])
    U.append(Estaciones[i][3])
    V.append(Estaciones[i][4])

for i in range(1,len(instancia)+1):
    for j in t:
        dem[(i,j)]=int((mov[instancia[i-1]][(j-1)//(T//20)]*g)//60)
        ttipo[i] = (tipo[instancia[i-1]])//abs(tipo[instancia[i-1]])
        if ttipo[i]== -1 or ttipo[i]==0:
            ttipo[i] = 2

U.append(np.random.rand(K)*20)
V.append(np.random.rand(K)*20)

dist ={(i,j): round(np.hypot(((abs(U[i-1]-U[j-1]))//g)+1,((abs(V[i-1]-V[j-1]))//g)+1)) 
     for i in n for j in n if i!=j}
for i in n:
    dist[i,i]=1
    # Se añade distancia vehículo - estación
for i in k:
    for j in n:
        dist[(-i,j)]= round(np.hypot(((abs(U[i]-U[j-1]))//3)+1,((abs(V[i]-V[j-1]))//3)+1))
        
S1 = [(-i,0) for i in k]
S2 = [(i,j) for i in n for j in t]
SS = S1+S2

Pred = {(i,j):(i,j-1) for i in n for j in t if j > 1 }
Succ = {(i,j):(i,j+1) for i in n for j in t if j < T}
for i in k:
    Succ[(-i,0)]=[]
    
A1 = [(i,j) for i in SS for j in S2 if j[0]!=i[0] 
       and j[1]-i[1]>=dist[(i[0],j[0])] and dist[(i[0],j[0])]>j[1]-i[1]-1]
        # Esperar un tiempo en la misma estación
A2 = [(i,Succ[i]) for i in S2 if i[1] < T]
AA=A1+A2

IncidPred = {s:[] for s in SS} # {Arco:[ArcoPred1,ArcoPred2...]}
IncidSucc = {s:[] for s in SS}
for a in AA:
    if a[0] not in IncidPred[a[1]]:
        IncidPred[a[1]].append(a[0])
    if a[1] not in IncidSucc[a[0]]:
        IncidSucc[a[0]].append(a[1])

rutas = {i:[(-i,0)] for i in k}
llegada = {i:0 for i in k}
Q = 25
Qt = [np.random.randint(10,16) for i in k]
Ct = C0.copy()
Fobj = 0
temp = 0
ocupados = []
ytemp ={i:0 for i in n}
QtN = Qt.copy()
CN = C.copy()
CtN = Ct.copy()
CtM = Ct.copy()
QtM = Qt.copy()
w1=0.01
w2=0.1
# Programa principal
start = time.time()
while temp <= T:
    if temp != 0:
        for i in k:
            if temp != 0 and rutas[i][-1][1] == temp:
                ocupados.remove(rutas[i][-1][0])
                ytemp[rutas[i][-1][0]]=0
                Ct[rutas[i][-1][0]-1] += llegada[i]
                Qt[i-1] -= llegada[i]
        for s in n:
            if ttipo[s]==1: # est de dejada, dem > 0
                Fobj = Fobj + max(0,dem[(s,temp)]-Ct[s-1])
                Ct[s-1] = Ct[s-1] - min(dem[(s,temp)],Ct[s-1])
                ytemp[s]+= max(0,dem[(s,temp)]-Ct[s-1])
            elif ttipo[s]==2: # est de recogida, dem < 0
                Fobj = Fobj + max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                Ct[s-1] = min(C[s-1],Ct[s-1] - dem[(s,temp)])
                ytemp[s]+= max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
    for i in k:
        if rutas[i][-1][1] == temp and temp < T:
                pos = rutas[i][-1] # Se renombra por comodidad
                valor = 0
                elegido = ['dum']
                nbicis=0
                for s in IncidSucc[pos]:
                    if s[0] not in ocupados or s == Succ[pos]:
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
                        if ttipo[s[0]]==2: # est de recogida, dem < 0
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
                        elif ttipo[s[0]]==1:
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
                elif ttipo[s[0]]==2: # est de recogida, dem < 0
                    llegada[i] = -nbicis
                elif ttipo[s[0]]==1:
                    llegada[i] = nbicis
                rutas[i].append(s)
                ocupados.append(s[0])
    temp+=1
end = time.time()
print('AVET',Fobj,end-start)
FobjN = 0
temp = 0
while temp <= T:
    if temp != 0:
        for s in n:
            if ttipo[s]==1: # est de dejada, dem > 0
                FobjN = FobjN + max(0,dem[(s,temp)]-CtN[s-1])
                CtN[s-1] = CtN[s-1] - min(dem[(s,temp)],CtN[s-1])
            elif ttipo[s]==2: # est de recogida, dem < 0
                FobjN = FobjN + max(0,CtN[s-1] - dem[(s,temp)] - CN[s-1])
                CtN[s-1] = min(CN[s-1],CtN[s-1] - dem[(s,temp)])
    temp+=1
print('NO',FobjN)


start = time.time()

# PILOT

rutasM = {i:[(-i,0)] for i in k}
llegadaM = {i:0 for i in k}
FobjM = 0
tempM = 0
ocupadosM = []
ytempM ={i:0 for i in n}

while tempM <= T:
    if tempM != 0:
        for i in k:
            if rutasM[i][-1][1] == tempM:
                ocupadosM.remove(rutasM[i][-1][0])
                CtM[rutasM[i][-1][0]-1] += llegadaM[i]
                QtM[i-1] -= llegadaM[i]
                ytempM[rutasM[i][-1][0]]=0
        for s in n:
            if ttipo[s]==1: # est de dejada, dem > 0
                FobjM = FobjM + max(0,dem[(s,tempM)]-CtM[s-1])
                CtM[s-1] = CtM[s-1] - min(dem[(s,tempM)],CtM[s-1])
                ytempM[s]+= max(0,dem[(s,tempM)]-CtM[s-1])
            elif ttipo[s]==2: # est de recogida, dem < 0
                FobjM = FobjM + max(0,CtM[s-1] - dem[(s,tempM)] - C[s-1])
                CtM[s-1] = min(C[s-1],CtM[s-1] - dem[(s,tempM)])
                ytempM[s]+= max(0,dem[(s,tempM)]-CtM[s-1])
    for i in k:
        if rutasM[i][-1][1] == tempM and tempM < T:
            FobjTemp = 99999
            elegidoT = ['dum']
            llegadaM[i] = 0
            for suc in IncidSucc[rutasM[i][-1]]:
                if suc[0] not in ocupadosM or suc == Succ[rutasM[i][-1]]:
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
                    if ttipo[suc[0]]==2: # est de recogida, dem < 0
                        nbicis = min(Ct[suc[0]-1]-demanda,Q-Qt[i-1])
                        llegada[i] = -nbicis
                    elif ttipo[suc[0]]==1:
                        nbicis = min(C[suc[0]-1]-Ct[suc[0]-1]+demanda,Qt[i-1])
                        llegada[i] = nbicis    
                    # T1 = temp + T//10 <-- limitación de profundidad
                    
                    # Programa principal
                    while temp <= T: #and temp <= T1:
                        if temp != tempM:
                            for ii in k:
                                if rutas[ii][-1][1] == temp:
                                    ocupados.remove(rutas[ii][-1][0])
                                    Ct[rutas[ii][-1][0]-1] += llegada[ii]
                                    Qt[ii-1] -= llegada[ii]
                                    ytemp[rutas[ii][-1][0]]=0
                            for s in n:
                                if ttipo[s]==1: # est de dejada, dem > 0
                                    Fobj = Fobj + max(0,dem[(s,temp)]-Ct[s-1])
                                    Ct[s-1] = Ct[s-1] - min(dem[(s,temp)],Ct[s-1])
                                    ytemp[s]+= max(0,dem[(s,temp)]-Ct[s-1])
                                elif ttipo[s]==2: # est de recogida, dem < 0
                                    Fobj = Fobj + max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                                    Ct[s-1] = min(C[s-1],Ct[s-1] - dem[(s,temp)])
                                    ytemp[s]+= max(0,Ct[s-1] - dem[(s,temp)] - C[s-1])
                        for ii in k:
                            if rutas[ii][-1][1] == temp:
                                    valor = -0.1
                                    elegido = ['dum']
                                    nbicis=0
                                    for s in IncidSucc[rutas[ii][-1]]:
                                        if s[0] not in ocupados or s == Succ[rutas[ii][-1]]:
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
                                            if ttipo[s[0]]==2: # est de recogida, dem < 0
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
                                            elif ttipo[s[0]]==1:
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
                                    elif ttipo[s[0]]==2: # est de recogida, dem < 0
                                        llegada[ii] = -nbicis
                                    elif ttipo[s[0]]==1:
                                        llegada[ii] = nbicis
                                    rutas[ii].append(s)
                                    ocupados.append(s[0])
                        temp+=1
                    if Fobj < FobjTemp:
                        elegidoT = suc
                        FobjTemp = Fobj
            if tempM == T:
                pass
            else:
                demanda = 0
                for j in range(1,dist[(rutasM[i][-1][0],elegidoT[0])]+1):
                    demanda+= dem[elegidoT[0],tempM+j]
                if ttipo[elegidoT[0]]==2: # est de recogida, dem < 0
                    nbicis = min(CtM[elegidoT[0]-1]-demanda,Q-QtM[i-1])
                    llegadaM[i] = -nbicis
                elif ttipo[elegidoT[0]]==1:
                    nbicis = min(C[elegidoT[0]-1]-CtM[elegidoT[0]-1]+demanda,QtM[i-1])
                    llegadaM[i] = nbicis
                rutasM[i].append(elegidoT)
                ocupadosM.append(elegidoT[0])
    tempM+=1

end = time.time()
print('PILOT',FobjM,end - start)